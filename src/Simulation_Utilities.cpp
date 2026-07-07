#include "Simulation_Utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Halo_Models.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr double KEPLER_DOMAIN_TOLERANCE = 1.0e-10;

double Clamp_Cosine(double value)
{
	return std::max(-1.0, std::min(1.0, value));
}

void Abort_Hyperbolic_Kepler_Shift(const std::string& message)
{
	std::cerr << "Error in Hyperbolic_Kepler_Shift(): " << message << std::endl;
	std::exit(EXIT_FAILURE);
}

bool Finite_Positive(double value)
{
	return std::isfinite(value) && value > 0.0;
}

}	// namespace

// 1. Event Class
Event::Event()
: time(0.0), position({0, 0, 0}), velocity({0, 0, 0})
{
}

Event::Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel)
: time(t), position(pos), velocity(vel)
{
}

double Event::Radius() const
{
	return position.Norm();
}

double Event::Speed() const
{
	return velocity.Norm();
}

double Event::Angular_Momentum() const
{
	return position.Cross(velocity).Norm();
}

double Event::Asymptotic_Speed_Sqr(Solar_Model& solar_model) const
{
	double r	 = Radius();
	double v	 = Speed();
	double vesc	 = solar_model.Local_Escape_Speed(r);
	double u_sqr = v * v - vesc * vesc;
	return u_sqr;
}

double Event::Isoreflection_Angle(const libphysica::Vector& vel_sun) const
{
	return acos(position.Normalized().Dot(vel_sun.Normalized()));
}

int Event::Isoreflection_Ring(const libphysica::Vector& vel_sun, unsigned int number_of_rings) const
{
	double theta					= Isoreflection_Angle(vel_sun);
	std::vector<double> ring_angles = Isoreflection_Ring_Angles(number_of_rings);
	for(unsigned int ring = 0; ring < ring_angles.size(); ring++)
		if(theta <= ring_angles[ring])
			return ring;
	std::cerr << "Error in Event::Isoreflection_Ring(): Angle = " << theta << " out of bound." << std::endl;
	std::exit(EXIT_FAILURE);
}

Event Event::In_Units(double unit_distance, double unit_time) const
{
	return Event(libphysica::natural_units::In_Units(time, unit_time), libphysica::natural_units::In_Units(position, unit_distance), libphysica::natural_units::In_Units(velocity, unit_distance / unit_time));
}

// Overload <<
std::ostream& operator<<(std::ostream& output, const Event& event)
{
	return output << "{"
				  << event.time
				  << ","
				  << event.position
				  << ","
				  << event.velocity
				  << "}";
}

// 2. Generator of initial conditions
double PDF_Initial_Speed(double v, obscura::DM_Distribution& halo_model, Solar_Model& solar_model)
{
	double v_esc			 = solar_model.Local_Escape_Speed(rSun);
	double v_average		 = halo_model.Average_Speed();
	double v_inverse_average = halo_model.Eta_Function(0.0);
	return halo_model.PDF_Speed(v) * (v + v_esc * v_esc / v) / (v_average + v_esc * v_esc * v_inverse_average);
}

// Conditional pdf for cos_theta given a speed value v
double PDF_Cos_Theta(double cos_theta, double v, obscura::DM_Distribution& halo_model)
{
	double normalization = 2.0 * M_PI * v * v / halo_model.PDF_Speed(v);

	// Construct velocity vector with angle theta to the Sun's velocity
	libphysica::Vector vel_sun = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Get_Observer_Velocity();
	libphysica::Vector vel	   = libphysica::Spherical_Coordinates(v, acos(cos_theta), 0.0, vel_sun);

	return normalization * halo_model.PDF_Velocity(vel);
}

Event Initial_Conditions(obscura::DM_Distribution& halo_model, Solar_Model& solar_model, std::mt19937& PRNG)
{
	// 1. Initial velocity
	// 1.1. Sample initial speed u asymptotically far from the Sun.
	std::function<double(double)> pdf_v = [&halo_model, &solar_model](double v) {
		return PDF_Initial_Speed(v, halo_model, solar_model);
	};
	double u = libphysica::Rejection_Sampling(pdf_v, halo_model.Minimum_DM_Speed(), halo_model.Maximum_DM_Speed(), 1200.0, PRNG);

	// 1.2. Sample cos(theta) where theta is the angle between v and v_sun.
	std::function<double(double)> pdf_cos_theta = [u, &halo_model](double cos_theta) {
		return PDF_Cos_Theta(cos_theta, u, halo_model);
	};
	// Determine domain of cos theta
	libphysica::Vector vel_sun = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Get_Observer_Velocity();
	double v_sun			   = vel_sun.Norm();
	double v_gal			   = halo_model.Maximum_DM_Speed() - v_sun;
	double cos_theta_max	   = std::min(1.0, (v_gal * v_gal - v_sun * v_sun - u * u) / (2.0 * u * v_sun));

	double y_max	 = PDF_Cos_Theta(-1.0, u, halo_model);
	double cos_theta = libphysica::Rejection_Sampling(pdf_cos_theta, -1.0, cos_theta_max, y_max, PRNG);
	// double cos_theta = libphysica::Sample_Uniform(PRNG,-1.0,1.0);// to test isotropic initial conditions

	// 1.3. Construct velocity vector
	double phi							= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector initial_velocity = libphysica::Spherical_Coordinates(u, acos(cos_theta), phi, vel_sun);

	// 1.4. Blue-shift the speed
	double asymptotic_distance = 1000.0 * AU;
	double vesc_asymptotic	   = solar_model.Local_Escape_Speed(asymptotic_distance);
	double v				   = sqrt(u * u + vesc_asymptotic * vesc_asymptotic);
	initial_velocity		   = v * initial_velocity.Normalized();

	// 2. Initial position
	// 2.1 Find the maximum impact parameter such that the particle still hits the Sun.
	double v_esc				= solar_model.Local_Escape_Speed(rSun);
	double impact_parameter_max = sqrt(u * u + v_esc * v_esc) / v * rSun;
	libphysica::Vector e_z		= (-1.0) * initial_velocity.Normalized();
	libphysica::Vector e_x({0, e_z[2], -e_z[1]});
	e_x.Normalize();
	libphysica::Vector e_y = e_z.Cross(e_x);

	// 2.2 Find a random point in the plane.
	double phi_disk						= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double xi							= libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	double impact_parameter				= sqrt(xi) * impact_parameter_max;
	libphysica::Vector initial_position = asymptotic_distance * e_z + impact_parameter * (cos(phi_disk) * e_x + sin(phi_disk) * e_y);

	return Event(0.0, initial_position, initial_velocity);
}

// 3. Analytically propagate a particle at event on a hyperbolic Kepler orbit to a radius R (without passing the periapsis)
void Hyperbolic_Kepler_Shift(Event& event, double R_final)
{
	const double mu = G_Newton * mSun;
	const double R_initial = event.Radius();
	const double speed = event.Speed();

	if(R_final < rSun || R_initial < rSun)
		Abort_Hyperbolic_Kepler_Shift("orbits inside the Sun cannot be described analytically.");
	if(!Finite_Positive(R_final) || !Finite_Positive(R_initial) || !Finite_Positive(speed))
		Abort_Hyperbolic_Kepler_Shift("initial or final radius/speed is non-finite or non-positive.");
	if(std::fabs(R_final - R_initial) <= KEPLER_DOMAIN_TOLERANCE * std::max(R_initial, R_final))
		return;

	const double radial_velocity = event.position.Dot(event.velocity) / R_initial;
	if(!std::isfinite(radial_velocity) || radial_velocity == 0.0)
		Abort_Hyperbolic_Kepler_Shift("cannot choose a forward Kepler branch at zero or non-finite radial velocity.");

	const double radial_direction = (radial_velocity < 0.0) ? -1.0 : 1.0;
	if(radial_direction < 0.0 && R_final > R_initial)
		Abort_Hyperbolic_Kepler_Shift("requested outward shift while particle is on the inbound branch.");
	if(radial_direction > 0.0 && R_final < R_initial)
		Abort_Hyperbolic_Kepler_Shift("requested inward shift while particle is on the outbound branch.");

	const double energy = 0.5 * speed * speed - mu / R_initial;
	if(!Finite_Positive(energy))
		Abort_Hyperbolic_Kepler_Shift("orbit is not hyperbolic.");

	const libphysica::Vector r_hat = event.position / R_initial;
	const libphysica::Vector h_vec = event.position.Cross(event.velocity);
	const double h = h_vec.Norm();
	const double h_scale = R_initial * speed;

	if(!std::isfinite(h) || h <= KEPLER_DOMAIN_TOLERANCE * h_scale)
	{
		const double final_speed_sqr = 2.0 * (energy + mu / R_final);
		if(!Finite_Positive(final_speed_sqr))
			Abort_Hyperbolic_Kepler_Shift("radial branch has non-finite or non-positive final speed.");
		event.position = R_final * r_hat;
		event.velocity = radial_direction * sqrt(final_speed_sqr) * r_hat;
		return;
	}

	const double p = h * h / mu;
	const libphysica::Vector eccentricity_vector = event.velocity.Cross(h_vec) / mu - r_hat;
	const double eccentricity = eccentricity_vector.Norm();
	if(!Finite_Positive(p) || !std::isfinite(eccentricity) || eccentricity <= 1.0)
		Abort_Hyperbolic_Kepler_Shift("invalid hyperbolic orbital elements.");

	const double cos_theta_final_raw = (p / R_final - 1.0) / eccentricity;
	if(!std::isfinite(cos_theta_final_raw)
	   || cos_theta_final_raw > 1.0 + KEPLER_DOMAIN_TOLERANCE
	   || cos_theta_final_raw < -1.0 - KEPLER_DOMAIN_TOLERANCE)
		Abort_Hyperbolic_Kepler_Shift("target radius is not reachable on this Kepler orbit.");

	libphysica::Vector axis_x = eccentricity_vector / eccentricity;
	libphysica::Vector axis_z = h_vec / h;
	libphysica::Vector axis_y = axis_z.Cross(axis_x);
	const double axis_y_norm = axis_y.Norm();
	if(!Finite_Positive(axis_y_norm))
		Abort_Hyperbolic_Kepler_Shift("failed to construct orbital basis.");
	axis_y = axis_y / axis_y_norm;
	axis_x = axis_y.Cross(axis_z).Normalized();

	const double theta_final_abs = acos(Clamp_Cosine(cos_theta_final_raw));
	const double theta_final = radial_direction * theta_final_abs;

	event.position = R_final * (cos(theta_final) * axis_x + sin(theta_final) * axis_y);
	event.velocity = sqrt(mu / p) * (-sin(theta_final) * axis_x + (eccentricity + cos(theta_final)) * axis_y);
}

// 4. Equiareal isodetection rings
std::vector<double> Isoreflection_Ring_Angles(unsigned int number_of_rings)
{
	std::vector<double> thetas;
	double theta = 0;
	for(unsigned int i = 0; i < number_of_rings - 1; i++)
	{
		theta = acos(cos(theta) - 2.0 / number_of_rings);
		thetas.push_back(theta);
	}
	thetas.push_back(180. * deg);
	return thetas;
}

}	// namespace DaMaSCUS_SUN

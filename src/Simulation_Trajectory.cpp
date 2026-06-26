#include "Simulation_Trajectory.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/Astronomy.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr double RK45_MIN_STEP_FACTOR = 0.1;
constexpr double RK45_MAX_STEP_FACTOR = 4.0;
// 单次 Runge_Kutta_45_Step 内层 while(!accepted) 的最大重试次数。
// 在正常物理场景下 <100 次即可收敛；设为 2000 是为了在极端数值情形下
// 保证外层循环一定能在有限时间内拿回控制权，从而触发 snapshot 回调。
constexpr int RK45_MAX_INNER_RETRIES = 2000;
constexpr double FREE_ENERGY_DRIFT_REL_E_SCALE_EV = 1.0e-30;
constexpr double MAX_OPTICAL_DEPTH_STEP = 0.05;
constexpr double CAPTURE_MODE_MAX_OPTICAL_DEPTH_STEP = 0.10;
constexpr double OPTICAL_DEPTH_RELATIVE_TOLERANCE = 1.0e-2;
constexpr double OPTICAL_DEPTH_ABSOLUTE_TOLERANCE = 1.0e-12 * MAX_OPTICAL_DEPTH_STEP;
constexpr unsigned int MAX_OPTICAL_DEPTH_RETRIES = 100;
constexpr std::size_t MAX_OPTICAL_DEPTH_PIECES = 4;
constexpr std::size_t CAPTURE_MODE_OPTICAL_DEPTH_PIECES = 2;
constexpr unsigned long int SNAPSHOT_PROGRESS_STEP_INTERVAL = 16;

struct OpticalDepthPiece
{
	double start;
	double end;
	double rate_start;
	double rate_end;
};

double RK45_Absolute_Max_Time_Step()
{
	return 1.0e6 * sec;  // 约 11.6 天；防止太阳外弱加速度区域步长增长到 inf。
}

double RK45_Sanitized_Time_Step(double step)
{
	if(!std::isfinite(step) || step <= 0.0)
		return 0.1 * sec;
	return std::min(step, RK45_Absolute_Max_Time_Step());
}

double Free_Propagation_Time_Step_Cap(double radius, double speed, double maximum_distance)
{
	double cap = RK45_Absolute_Max_Time_Step();
	const double safe_speed = std::max(std::fabs(speed), 1.0e-12 * km / sec);

	// 避免单步跨越过大的径向距离，否则可能直接跳过 2R_sun 边界并把 r 推到非物理值。
	const double crossing_scale = std::max(0.25 * maximum_distance, 10.0 * km);
	cap = std::min(cap, crossing_scale / safe_speed);

	// 太阳外近似 Kepler 运动，步长不应大于局域动力学时间的一个固定比例。
	const double safe_radius = std::max(radius, 1.0 * km);
	const double dynamical_time = sqrt(safe_radius * safe_radius * safe_radius / (G_Newton * mSun));
	if(std::isfinite(dynamical_time) && dynamical_time > 0.0)
		cap = std::min(cap, 0.1 * dynamical_time);

	return std::max(cap, 1.0e-8 * sec);
}

bool Outward_Escaping_At_Boundary(const Event& event, Solar_Model& solar_model, double boundary_radius)
{
	const double radius = event.Radius();
	if(radius < boundary_radius)
		return false;

	const double radial_velocity = (radius > 0.0) ? event.position.Dot(event.velocity) / radius : 0.0;
	return radial_velocity > 0.0 && event.Speed() > solar_model.Local_Escape_Speed(radius);
}

double Clamp_Unit_Interval(double value)
{
	return std::max(0.0, std::min(1.0, value));
}

Event Interpolate_Event(const Event& before, const Event& after, double fraction)
{
	fraction = Clamp_Unit_Interval(fraction);
	return Event(before.time + fraction * (after.time - before.time),
	             before.position + fraction * (after.position - before.position),
	             before.velocity + fraction * (after.velocity - before.velocity));
}

double Radius_At_Fraction(const Event& before, const Event& after, double fraction)
{
	return Interpolate_Event(before, after, fraction).Radius();
}

bool Surface_Crossing_Fractions(const Event& before, const Event& after, double& first, double& second)
{
	const libphysica::Vector displacement = after.position - before.position;
	const double a = displacement.Dot(displacement);
	if(a <= 0.0)
		return false;

	const double b = 2.0 * before.position.Dot(displacement);
	const double c = before.position.Dot(before.position) - rSun * rSun;
	const double discriminant = b * b - 4.0 * a * c;
	if(!std::isfinite(discriminant) || discriminant < 0.0)
		return false;

	const double sqrt_discriminant = std::sqrt(std::max(0.0, discriminant));
	first = (-b - sqrt_discriminant) / (2.0 * a);
	second = (-b + sqrt_discriminant) / (2.0 * a);
	if(first > second)
		std::swap(first, second);
	return true;
}

bool Fraction_In_Unit_Interval(double fraction)
{
	constexpr double tolerance = 1.0e-12;
	return fraction >= -tolerance && fraction <= 1.0 + tolerance;
}

double Find_Surface_Crossing_Fraction(const Event& before, const Event& after, double lower, double upper)
{
	double r_lower = Radius_At_Fraction(before, after, lower) - rSun;
	for(int iteration = 0; iteration < 60; iteration++)
	{
		const double mid = 0.5 * (lower + upper);
		const double r_mid = Radius_At_Fraction(before, after, mid) - rSun;
		if((r_lower <= 0.0 && r_mid <= 0.0) || (r_lower >= 0.0 && r_mid >= 0.0))
		{
			lower = mid;
			r_lower = r_mid;
		}
		else
			upper = mid;
	}
	return Clamp_Unit_Interval(0.5 * (lower + upper));
}

bool Solar_Interior_Fraction_Interval(const Event& before, const Event& after, double r_before, double r_after, double& start, double& end)
{
	const bool before_inside = r_before < rSun;
	const bool after_inside = r_after < rSun;
	double first_crossing = 0.0;
	double second_crossing = 0.0;
	const bool has_crossings = Surface_Crossing_Fractions(before, after, first_crossing, second_crossing);
	if(before_inside && after_inside)
	{
		start = 0.0;
		end = 1.0;
		return true;
	}
	if(!before_inside && after_inside)
	{
		if(has_crossings && Fraction_In_Unit_Interval(first_crossing))
			start = Clamp_Unit_Interval(first_crossing);
		else if(has_crossings && Fraction_In_Unit_Interval(second_crossing))
			start = Clamp_Unit_Interval(second_crossing);
		else
			start = Find_Surface_Crossing_Fraction(before, after, 0.0, 1.0);
		end = 1.0;
		return end > start;
	}
	if(before_inside && !after_inside)
	{
		start = 0.0;
		if(has_crossings && Fraction_In_Unit_Interval(second_crossing))
			end = Clamp_Unit_Interval(second_crossing);
		else if(has_crossings && Fraction_In_Unit_Interval(first_crossing))
			end = Clamp_Unit_Interval(first_crossing);
		else
			end = Find_Surface_Crossing_Fraction(before, after, 0.0, 1.0);
		return end > start;
	}

	if(has_crossings && Fraction_In_Unit_Interval(first_crossing) && Fraction_In_Unit_Interval(second_crossing))
	{
		start = Clamp_Unit_Interval(first_crossing);
		end = Clamp_Unit_Interval(second_crossing);
		return end > start;
	}
	return false;
}

double Scattering_Rate_At_Fraction(const Event& before, const Event& after, double fraction, Solar_Model& solar_model, obscura::DM_Particle& DM)
{
	const Event event = Interpolate_Event(before, after, fraction);
	const double radius = event.Radius();
	if(radius >= rSun)
		return 0.0;
	return solar_model.Total_DM_Scattering_Rate(DM, radius, event.Speed());
}

bool Build_Optical_Depth_Pieces(const Event& before,
                                const Event& after,
                                double interior_start,
                                double interior_end,
                                double actual_dt,
                                Solar_Model& solar_model,
                                obscura::DM_Particle& DM,
                                bool use_cached_start_rate,
                                double cached_start_rate,
                                std::size_t requested_piece_count,
                                std::array<OpticalDepthPiece, MAX_OPTICAL_DEPTH_PIECES>& pieces,
                                std::size_t& piece_count,
                                double& tau_two_piece,
                                double& tau_four_piece,
                                bool& invalid_rate)
{
	piece_count = 0;
	tau_two_piece = 0.0;
	tau_four_piece = 0.0;
	invalid_rate = false;
	if(!(interior_end > interior_start) || actual_dt <= 0.0)
		return false;

	requested_piece_count = std::max<std::size_t>(1, std::min(requested_piece_count, MAX_OPTICAL_DEPTH_PIECES));
	const double span = interior_end - interior_start;
	std::array<double, MAX_OPTICAL_DEPTH_PIECES + 1> fractions;
	std::array<double, MAX_OPTICAL_DEPTH_PIECES + 1> rates;
	for(std::size_t i = 0; i <= requested_piece_count; i++)
		fractions[i] = interior_start + static_cast<double>(i) / static_cast<double>(requested_piece_count) * span;

	for(std::size_t i = 0; i <= requested_piece_count; i++)
	{
		if(i == 0 && use_cached_start_rate)
			rates[i] = cached_start_rate;
		else
			rates[i] = Scattering_Rate_At_Fraction(before, after, fractions[i], solar_model, DM);

		if(!std::isfinite(rates[i]) || rates[i] < 0.0)
		{
			invalid_rate = true;
			return true;
		}
	}

	const auto trapezoid_tau = [&](std::size_t start_index, std::size_t end_index)
	{
		const double piece_dt = actual_dt * (fractions[end_index] - fractions[start_index]);
		return 0.5 * (rates[start_index] + rates[end_index]) * piece_dt;
	};

	if(requested_piece_count >= 4)
		tau_two_piece = trapezoid_tau(0, 2) + trapezoid_tau(2, 4);
	else
	{
		for(std::size_t i = 0; i < requested_piece_count; i++)
			tau_two_piece += trapezoid_tau(i, i + 1);
	}
	for(std::size_t i = 0; i < requested_piece_count; i++)
	{
		if(fractions[i + 1] > fractions[i])
		{
			pieces[piece_count++] = OpticalDepthPiece{fractions[i], fractions[i + 1], rates[i], rates[i + 1]};
			tau_four_piece += trapezoid_tau(i, i + 1);
		}
	}

	return piece_count > 0;
}

double Solve_Linear_Rate_Optical_Depth_Fraction(double rate_start, double rate_end, double interval_dt, double target_tau)
{
	if(target_tau <= 0.0)
		return 0.0;
	if(!std::isfinite(target_tau) || interval_dt <= 0.0)
		return 1.0;
	if(!std::isfinite(rate_start) || !std::isfinite(rate_end) || rate_start < 0.0 || rate_end < 0.0)
		return 1.0;

	const double total_tau = 0.5 * (rate_start + rate_end) * interval_dt;
	if(!std::isfinite(total_tau) || total_tau <= 0.0)
		return 1.0;
	if(target_tau >= total_tau)
		return 1.0;

	const double target = target_tau / interval_dt;
	const double a = 0.5 * (rate_end - rate_start);
	const double b = rate_start;
	if(std::fabs(a) < 1.0e-300)
		return (b > 0.0) ? Clamp_Unit_Interval(target / b) : 1.0;

	const double discriminant = b * b + 4.0 * a * target;
	if(std::isfinite(discriminant) && discriminant >= 0.0)
	{
		const double sqrt_discriminant = sqrt(discriminant);
		const double roots[2] = {
		    (-b + sqrt_discriminant) / (2.0 * a),
		    (-b - sqrt_discriminant) / (2.0 * a)
		};
		for(double root : roots)
		{
			if(std::isfinite(root) && root >= 0.0 && root <= 1.0)
				return Clamp_Unit_Interval(root);
		}
	}

	double low = 0.0;
	double high = 1.0;
	for(int iteration = 0; iteration < 30; iteration++)
	{
		const double mid = 0.5 * (low + high);
		const double tau_mid = interval_dt * (rate_start * mid + 0.5 * (rate_end - rate_start) * mid * mid);
		if(tau_mid < target_tau)
			low = mid;
		else
			high = mid;
	}
	return Clamp_Unit_Interval(0.5 * (low + high));
}

bool RK45_Errors_Within_Tolerance(const double errors[3], const double tolerances[3])
{
	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]) || errors[i] > tolerances[i])
			return false;
	}

	return true;
}

double RK45_Next_Step_Size(double current_step, const double errors[3], const double tolerances[3])
{
	// 若任何一个 error 为 NaN/Inf，说明 RK 中间阶段出现数值崩溃。
	// 此时必须让步长"缩小"而不是"放大"，否则接下来会继续崩，外层
	// while(!accepted) 会陷入无限循环（见 DaMaSCUS-SUN-EVAP#rank1-stuck）。
	for(int i = 0; i < 3; i++)
	{
		if(!std::isfinite(errors[i]))
			return RK45_MIN_STEP_FACTOR * current_step;
	}

	double factor = RK45_MAX_STEP_FACTOR;

	for(int i = 0; i < 3; i++)
	{
		if(errors[i] <= 0.0)
			continue;

		const double candidate = 0.84 * pow(tolerances[i] / errors[i], 0.25);
		if(std::isfinite(candidate))
			factor = std::min(factor, candidate);
	}

	factor = std::max(RK45_MIN_STEP_FACTOR, std::min(factor, RK45_MAX_STEP_FACTOR));
	return RK45_Sanitized_Time_Step(factor * current_step);
}
}

// 1. Result of one trajectory
bool TrajectoryTerminationInvalidatesSurvival(TrajectoryTerminationReason reason)
{
	switch(reason)
	{
		case TrajectoryTerminationReason::EnergyDriftEscape:
		case TrajectoryTerminationReason::NumericalFailure:
		case TrajectoryTerminationReason::NonFiniteState:
		case TrajectoryTerminationReason::SpeedLimit:
		case TrajectoryTerminationReason::WallTimeLimit:
		case TrajectoryTerminationReason::MaxFreeSteps:
		case TrajectoryTerminationReason::MaxScatterings:
		case TrajectoryTerminationReason::Unknown:
			return true;
		default:
			return false;
	}
}

Trajectory_Result::Trajectory_Result(const Event& event_ini, const Event& event_final, unsigned long int nScat, TrajectoryBincount bc)
: initial_event(event_ini), final_event(event_final), number_of_scatterings(nScat), bincount(std::move(bc))
{
}

bool Trajectory_Result::Particle_Reflected() const
{
	if(bincount.termination_reason != TrajectoryTerminationReason::OutwardEscape)
		return false;
	double r	= final_event.Radius();
	double vesc = sqrt(2 * G_Newton * mSun / r);
	return r > rSun && final_event.Speed() > vesc && number_of_scatterings > 0;
}

bool Trajectory_Result::Particle_Free() const
{
	return bincount.termination_reason == TrajectoryTerminationReason::OutwardEscape
	    && number_of_scatterings == 0;
}

bool Trajectory_Result::Particle_Captured(Solar_Model& solar_model) const
{
	double r	= final_event.Radius();
	double vesc = solar_model.Local_Escape_Speed(r);
	return final_event.Speed() < vesc;
}

void Trajectory_Result::Print_Summary(Solar_Model& solar_model, unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Trajectory result summary" << std::endl
				  << std::endl
				  << "Number of scatterings:\t" << number_of_scatterings << std::endl
				  << "Simulation time [days]:\t" << libphysica::Round(In_Units(final_event.time, day)) << std::endl
				  << "Final radius [rSun]:\t" << libphysica::Round(In_Units(final_event.Radius(), rSun)) << std::endl
				  << "Final speed [km/sec]:\t" << libphysica::Round(In_Units(final_event.Speed(), km / sec)) << std::endl
				  << "Free particle:\t\t[" << (Particle_Free() ? "x" : " ") << "]" << std::endl
				  << "Captured:\t\t[" << (Particle_Captured(solar_model) ? "x" : " ") << "]" << std::endl
				  << "Reflection:\t\t[" << (Particle_Reflected() ? "x" : " ") << "]";

		if(Particle_Reflected())
		{
			double u_i_sqr = initial_event.Asymptotic_Speed_Sqr(solar_model);
			double u_f_sqr = final_event.Asymptotic_Speed_Sqr(solar_model);
			// 检查渐近速度平方是否为正
			if(u_i_sqr > 0.0 && u_f_sqr > 0.0)
			{
				double u_i = sqrt(u_i_sqr);
				double u_f = sqrt(u_f_sqr);
				std::cout << "\t(ratio u_f/u_i = " << libphysica::Round(u_f / u_i) << ")" << std::endl;
			}
			else
				std::cout << "\t(asymptotic speeds unavailable)" << std::endl;
		}
		else
			std::cout << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

// 2. Simulator
Trajectory_Simulator::Trajectory_Simulator(const Solar_Model& model, unsigned long int max_time_steps, unsigned long int max_scatterings, double max_distance)
: solar_model(model), free_flight_reference_energy_eV(std::numeric_limits<double>::quiet_NaN()), current_physical_bound_state(false), terminate_on_capture(false), trajectory_in_progress(false), accumulated_snapshot_overhead_sec(0.0), maximum_time_steps(max_time_steps), maximum_scatterings(max_scatterings), maximum_distance(max_distance), current_mpi_rank(0), current_trajectory_id(0)
{
	std::random_device rd;
	PRNG.seed(rd());
	rate_nuclei_cache.resize(solar_model.target_isotopes.size());
}

void Trajectory_Simulator::Publish_Snapshot_Progress() const
{
	if(snapshot_progress_callback)
	{
		const auto callback_start = std::chrono::steady_clock::now();
		snapshot_progress_callback(*this);
		const auto callback_end = std::chrono::steady_clock::now();
		if(trajectory_in_progress)
			accumulated_snapshot_overhead_sec += 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(callback_end - callback_start).count();
	}
}

void Trajectory_Simulator::Set_Snapshot_Progress_Callback(std::function<void(const Trajectory_Simulator&)> callback)
{
	snapshot_progress_callback = std::move(callback);
}

void Trajectory_Simulator::Enable_Capture_Mode(bool enabled)
{
	terminate_on_capture = enabled;
}

bool Trajectory_Simulator::Trajectory_In_Progress() const
{
	return trajectory_in_progress;
}

unsigned long int Trajectory_Simulator::Current_Trajectory_ID() const
{
	return current_trajectory_id;
}

double Trajectory_Simulator::Current_Trajectory_Wall_Time_Seconds() const
{
	if(!trajectory_in_progress)
		return 0.0;

	const double total_wall_time_sec = 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - current_trajectory_wall_start).count();
	return std::max(0.0, total_wall_time_sec - accumulated_snapshot_overhead_sec);
}

const TrajectoryBincount& Trajectory_Simulator::Current_Trajectory_Bincount() const
{
	return current_bincount;
}

// Accumulate one step into the current bincount
void Trajectory_Simulator::Accumulate_Bincount_Step(double r_km, double v2_km2s2, double dt_sec)
{
	if(dt_sec <= 0.0 || r_km < 0.0 || r_km >= BIN_MAX_KM)
		return;
	int bin_idx = static_cast<int>(r_km / BIN_WIDTH_KM);
	if(bin_idx < 0) bin_idx = 0;
	if(bin_idx >= NUM_BINS) return;
	current_bincount.dt_hist[bin_idx] += dt_sec;
	current_bincount.v2dt_hist[bin_idx] += v2_km2s2 * dt_sec;
}

void Trajectory_Simulator::Reset_Bincount_Anchor(const Event& event)
{
	if(terminate_on_capture)
		return;

	prev_time_sec = In_Units(event.time, sec);
	prev_r_km = In_Units(event.Radius(), km);
	const double v_kms = In_Units(event.Speed(), km / sec);
	prev_v2_km2s2 = v_kms * v_kms;
}

double Trajectory_Simulator::Capture_Energy_eV(double radius, double speed, obscura::DM_Particle& DM)
{
	double vesc = solar_model.Local_Escape_Speed(radius);
	double E = 0.5 * DM.mass * (speed * speed - vesc * vesc);
	return In_Units(E, eV);
}

bool Trajectory_Simulator::Update_Capture_State(double radius, double speed, double time, obscura::DM_Particle& DM, bool allow_new_capture)
{
	if(terminate_on_capture && !allow_new_capture && !current_bincount.is_captured)
		return false;

	double E_eV = Capture_Energy_eV(radius, speed, DM);
	double dE_from_prev_eV = std::numeric_limits<double>::quiet_NaN();
	if(std::isfinite(previous_capture_energy_eV))
		dE_from_prev_eV = E_eV - previous_capture_energy_eV;
	previous_capture_energy_eV = E_eV;
	bool negative_energy = E_eV < 0.0;
	double t_now_sec = In_Units(time, sec);

	if(!allow_new_capture)
	{
		if(current_bincount.is_captured && std::isfinite(free_flight_reference_energy_eV))
		{
			double drift_eV = std::fabs(E_eV - free_flight_reference_energy_eV);
			if(std::isfinite(drift_eV) && drift_eV > current_bincount.max_free_energy_drift_eV)
				current_bincount.max_free_energy_drift_eV = drift_eV;
			double scale_eV = std::max(std::fabs(free_flight_reference_energy_eV), FREE_ENERGY_DRIFT_REL_E_SCALE_EV);
			double drift_rel = drift_eV / scale_eV;
			if(std::isfinite(drift_rel) && drift_rel > current_bincount.max_free_energy_drift_rel)
				current_bincount.max_free_energy_drift_rel = drift_rel;
		}

		if(current_physical_bound_state && current_bincount.is_captured)
		{
			current_bincount.t_last_bound = t_now_sec;
			current_bincount.t_last_negative = t_now_sec;
		}
		return false;
	}

	bool was_bound = current_physical_bound_state;
	current_physical_bound_state = negative_energy;
	free_flight_reference_energy_eV = E_eV;

	if(negative_energy)
	{
		if(!current_bincount.is_captured)
		{
			current_bincount.is_captured = true;
			current_bincount.t_capture = t_now_sec;
			current_bincount.t_first_negative = t_now_sec;
			current_bincount.t_last_bound = t_now_sec;
			current_bincount.r_first_negative_km = In_Units(radius, km);
			current_bincount.E_first_negative_eV = E_eV;
			current_bincount.dE_first_negative_from_prev_eV = dE_from_prev_eV;
		}
		else if(!was_bound)
		{
			current_bincount.t_final_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();
		}
		current_bincount.t_last_bound = t_now_sec;
		current_bincount.t_last_negative = t_now_sec;
		return true;
	}

	if(was_bound && current_bincount.is_captured)
	{
		current_bincount.t_final_unbinding_scatter = t_now_sec;
	}

	return false;
}

TrajectoryTerminationReason Trajectory_Simulator::Propagate_Freely(Event& current_event, obscura::DM_Particle& DM)
{
	if(Outward_Escaping_At_Boundary(current_event, solar_model, maximum_distance))
		return TrajectoryTerminationReason::OutwardEscape;

	Free_Particle_Propagator particle_propagator(current_event);
	double minus_log_xi = -log(libphysica::Sample_Uniform(PRNG));
	TrajectoryTerminationReason outcome = TrajectoryTerminationReason::MaxFreeSteps;
	unsigned long int time_steps = 0;
	unsigned long int step_attempts = 0;
	unsigned int optical_depth_retries = 0;
	const bool fast_capture_mode = terminate_on_capture;
	const double optical_depth_step_limit = fast_capture_mode ? CAPTURE_MODE_MAX_OPTICAL_DEPTH_STEP : MAX_OPTICAL_DEPTH_STEP;
	const std::size_t optical_depth_piece_target = fast_capture_mode ? CAPTURE_MODE_OPTICAL_DEPTH_PIECES : MAX_OPTICAL_DEPTH_PIECES;

	auto commit_accepted_event = [&](const Event& accepted_event)
	{
		const double t_now_sec = In_Units(accepted_event.time, sec);
		if(!terminate_on_capture)
		{
			const double r_now_km  = In_Units(accepted_event.Radius(), km);
			const double v_now_kms = In_Units(accepted_event.Speed(), km / sec);
			const double v2_now    = v_now_kms * v_now_kms;

			if(prev_time_sec >= 0.0)
			{
				const double dt_sec = t_now_sec - prev_time_sec;
				Accumulate_Bincount_Step(prev_r_km, prev_v2_km2s2, dt_sec);
				prev_dt_sec = dt_sec;
			}

			prev_time_sec = t_now_sec;
			prev_r_km     = r_now_km;
			prev_v2_km2s2 = v2_now;
		}

		return Update_Capture_State(accepted_event.Radius(), accepted_event.Speed(), accepted_event.time, DM, false);
	};

	while(time_steps < maximum_time_steps && outcome == TrajectoryTerminationReason::MaxFreeSteps)
	{
		step_attempts++;
		Event event_before = particle_propagator.Event_In_3D();
		double r_before = particle_propagator.Current_Radius();
		double v_before = particle_propagator.Current_Speed();
		double rate_before = 0.0;
		if(r_before >= rSun)
		{
			const double step_cap = Free_Propagation_Time_Step_Cap(r_before, v_before, maximum_distance);
			particle_propagator.time_step = std::min(RK45_Sanitized_Time_Step(particle_propagator.time_step), step_cap);
			const double radial_velocity_before =
			    (r_before > 0.0) ? event_before.position.Dot(event_before.velocity) / r_before : 0.0;
			if(r_before > rSun && radial_velocity_before < 0.0)
			{
				const double surface_time_estimate = (r_before - rSun) / (-radial_velocity_before);
				if(std::isfinite(surface_time_estimate) && surface_time_estimate > 0.0)
					particle_propagator.time_step =
					    std::min(particle_propagator.time_step, std::max(surface_time_estimate, 1.0e-8 * sec));
			}
		}
		else
		{
			particle_propagator.time_step = RK45_Sanitized_Time_Step(particle_propagator.time_step);
			rate_before = solar_model.Total_DM_Scattering_Rate(DM, r_before, v_before);
			if(!std::isfinite(rate_before) || rate_before < 0.0)
			{
				std::cerr << "\nWarning in Propagate_Freely(): invalid pre-step scattering rate (rank "
				          << current_mpi_rank << ", traj " << current_trajectory_id
				          << ", rate=" << rate_before << "). Marking trajectory as numerically failed." << std::endl;
				Publish_Snapshot_Progress();
				current_event = particle_propagator.Event_In_3D();
				return TrajectoryTerminationReason::NumericalFailure;
			}
			if(rate_before > 0.0)
			{
				particle_propagator.time_step = std::min(particle_propagator.time_step, optical_depth_step_limit / rate_before);
			}
		}

		const Free_Particle_Propagator propagator_before = particle_propagator;
		double t_before = particle_propagator.Current_Time();
		bool rk_step_ok = particle_propagator.Runge_Kutta_45_Step(solar_model);
		double actual_dt = particle_propagator.Current_Time() - t_before;
		double r_after = particle_propagator.Current_Radius();
		double v_after = particle_propagator.Current_Speed();
		Event event_after = particle_propagator.Event_In_3D();

		if(!rk_step_ok)
		{
			std::cerr << "\nWarning in Propagate_Freely(): RK45 step failed tolerance after retry limits (rank "
			          << current_mpi_rank << ", traj " << current_trajectory_id
			          << "). Marking trajectory as numerically failed." << std::endl;
			Publish_Snapshot_Progress();
			current_event = particle_propagator.Event_In_3D();
			return TrajectoryTerminationReason::NumericalFailure;
		}

		if(!std::isfinite(r_after) || !std::isfinite(v_after) || !std::isfinite(actual_dt) || actual_dt <= 0.0)
		{
			// 出现 NaN/Inf，轨迹已无法继续。向外传一次进度，再中止这条轨迹，
			// 避免 rank 被单条坏轨迹卡住导致 snapshot / MPI_Barrier 死锁。
			std::cerr << "\nWarning in Propagate_Freely(): non-finite state (rank " << current_mpi_rank
			          << ", traj " << current_trajectory_id << ", r=" << r_after
			          << ", v=" << v_after << ", dt=" << actual_dt << "). Aborting trajectory." << std::endl;
			Publish_Snapshot_Progress();
			current_event = particle_propagator.Event_In_3D();
			return TrajectoryTerminationReason::NonFiniteState;
		}

		if(v_after > v_max)
		{
			std::cerr << "\nWarning in Propagate_Freely(): DM speed exceeds the maximum of v_max = " << v_max << std::endl
					  << "\tAbort simulation." << std::endl;
			Publish_Snapshot_Progress();
			current_event = particle_propagator.Event_In_3D();
			return TrajectoryTerminationReason::SpeedLimit;
		}

		// 单条轨迹 wall-clock 预算：防止任何一条轨迹把整个 rank 永远拖住。
		// 正式统计默认不限制；若配置为正值，每 256 步检查一次。
		if(max_trajectory_wall_time_sec > 0.0 && (step_attempts & 0xFFu) == 0u)
		{
			double traj_wall = Current_Trajectory_Wall_Time_Seconds();
			if(traj_wall > max_trajectory_wall_time_sec)
			{
				std::cerr << "\nWarning in Propagate_Freely(): trajectory wall-time budget exceeded (rank "
				          << current_mpi_rank << ", traj " << current_trajectory_id << ", traj_wall="
				          << traj_wall << "s > " << max_trajectory_wall_time_sec
				          << "s). Aborting trajectory." << std::endl;
				Publish_Snapshot_Progress();
				current_event = particle_propagator.Event_In_3D();
				return TrajectoryTerminationReason::WallTimeLimit;
			}
		}

		// Check for scatterings and reflection
		bool scattering = false;
		bool reflection = false;
		Event accepted_event = event_after;
		double interior_start = 0.0;
		double interior_end = 0.0;
		if(Solar_Interior_Fraction_Interval(event_before, event_after, r_before, r_after, interior_start, interior_end))
		{
			if(r_after < rSun && v_after < 0.0)
			{
				std::cerr << "Warning: Negative velocity detected (v = " << v_after << ") at r = " << r_after << ", skipping scattering calculation." << std::endl;
				break;
			}
			std::array<OpticalDepthPiece, MAX_OPTICAL_DEPTH_PIECES> optical_pieces;
			std::size_t optical_piece_count = 0;
			double tau_two_piece = 0.0;
			double delta_tau = 0.0;
			bool invalid_rate = false;
			const bool use_cached_start_rate = r_before < rSun && interior_start <= 1.0e-12;
				Build_Optical_Depth_Pieces(event_before,
				                           event_after,
				                           interior_start,
				                           interior_end,
				                           actual_dt,
				                           solar_model,
				                           DM,
				                           use_cached_start_rate,
				                           rate_before,
				                           optical_depth_piece_target,
				                           optical_pieces,
				                           optical_piece_count,
				                           tau_two_piece,
				                           delta_tau,
				                           invalid_rate);
			if(invalid_rate || !std::isfinite(delta_tau) || delta_tau < 0.0)
			{
				std::cerr << "\nWarning in Propagate_Freely(): invalid scattering rate (rank "
				          << current_mpi_rank << ", traj " << current_trajectory_id
				          << "). Marking trajectory as numerically failed." << std::endl;
				Publish_Snapshot_Progress();
				current_event = event_after;
				return TrajectoryTerminationReason::NumericalFailure;
			}

				const double tau_error_scale = std::max(delta_tau, OPTICAL_DEPTH_ABSOLUTE_TOLERANCE);
				const double tau_relative_error = std::fabs(delta_tau - tau_two_piece) / tau_error_scale;
				const bool reject_for_optical_depth = delta_tau > optical_depth_step_limit;
				const bool reject_for_tau_accuracy = !fast_capture_mode && tau_relative_error > OPTICAL_DEPTH_RELATIVE_TOLERANCE;
				if(reject_for_optical_depth || reject_for_tau_accuracy)
				{
					optical_depth_retries++;
					if(optical_depth_retries > MAX_OPTICAL_DEPTH_RETRIES)
					{
						std::cerr << "\nWarning in Propagate_Freely(): optical-depth step rejected more than "
						          << MAX_OPTICAL_DEPTH_RETRIES << " consecutive times (rank "
						          << current_mpi_rank << ", traj " << current_trajectory_id
						          << ", delta_tau=" << delta_tau
						          << ", relative_error=" << tau_relative_error
						          << "). Marking trajectory as numerically failed." << std::endl;
						Publish_Snapshot_Progress();
						current_event = event_before;
						return TrajectoryTerminationReason::NumericalFailure;
					}

					particle_propagator = propagator_before;
					double factor = RK45_MAX_STEP_FACTOR;
					if(delta_tau > optical_depth_step_limit && delta_tau > 0.0)
						factor = std::min(factor, 0.8 * optical_depth_step_limit / delta_tau);
					if(reject_for_tau_accuracy && tau_relative_error > 0.0)
						factor = std::min(factor, 0.8 * std::sqrt(OPTICAL_DEPTH_RELATIVE_TOLERANCE / tau_relative_error));
					factor = std::max(RK45_MIN_STEP_FACTOR, std::min(factor, 1.0));
					particle_propagator.time_step = RK45_Sanitized_Time_Step(factor * actual_dt);
					continue;
				}

			if(delta_tau >= minus_log_xi)
			{
				double remaining_tau = minus_log_xi;
				for(std::size_t piece_index = 0; piece_index < optical_piece_count; piece_index++)
				{
					const OpticalDepthPiece& piece = optical_pieces[piece_index];
					const double piece_dt = actual_dt * (piece.end - piece.start);
					const double piece_tau = 0.5 * (piece.rate_start + piece.rate_end) * piece_dt;
					if(piece_tau >= remaining_tau)
					{
						const double local_fraction =
						    Solve_Linear_Rate_Optical_Depth_Fraction(piece.rate_start, piece.rate_end, piece_dt, remaining_tau);
						const double fraction = piece.start + local_fraction * (piece.end - piece.start);
						accepted_event = Interpolate_Event(event_before, event_after, fraction);
						scattering = true;
						outcome = TrajectoryTerminationReason::Scatter;
						break;
					}
					remaining_tau -= piece_tau;
				}
				if(!scattering)
				{
					accepted_event = Interpolate_Event(event_before, event_after, interior_end);
					scattering = true;
					outcome = TrajectoryTerminationReason::Scatter;
				}
			}
			else
				minus_log_xi -= delta_tau;
		}

		if(!scattering)
		{
			if(r_after >= maximum_distance && r_after >= r_before && v_after > solar_model.Local_Escape_Speed(r_after))
			{
				reflection = true;
				outcome = TrajectoryTerminationReason::OutwardEscape;
			}
		}

		time_steps++;
		optical_depth_retries = 0;
		const bool captured_now = commit_accepted_event(accepted_event);
		current_event = accepted_event;
		if((time_steps % SNAPSHOT_PROGRESS_STEP_INTERVAL) == 0u)
			Publish_Snapshot_Progress();
		if(terminate_on_capture && captured_now)
			return TrajectoryTerminationReason::CaptureMode;

		if(reflection || scattering)
			break;

	}

	return outcome;
}


int Trajectory_Simulator::Sample_Target(obscura::DM_Particle& DM, double r, double DM_speed)
{
	if(r > rSun)
	{
		throw std::runtime_error("Sample_Target(): r > rSun.");
	}
	else
	{
		// C: 复用预分配的 rate_nuclei_cache，避免每次散射事件堆分配
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			rate_nuclei_cache[i] = solar_model.DM_Scattering_Rate_Nucleus(DM, r, DM_speed, i);
			if(!std::isfinite(rate_nuclei_cache[i]) || rate_nuclei_cache[i] < 0.0)
				throw std::runtime_error("Sample_Target(): nucleus scattering rate is negative or non-finite.");
		}
		double rate_electron = solar_model.DM_Scattering_Rate_Electron(DM, r, DM_speed);
		if(!std::isfinite(rate_electron) || rate_electron < 0.0)
			throw std::runtime_error("Sample_Target(): electron scattering rate is negative or non-finite.");
		double total_rate	 = std::accumulate(rate_nuclei_cache.begin(), rate_nuclei_cache.end(), rate_electron);
		if(!std::isfinite(total_rate) || total_rate <= 0.0)
			throw std::runtime_error("Sample_Target(): total scattering rate is non-positive or non-finite.");

		double xi = libphysica::Sample_Uniform(PRNG);
		// Electron
		double sum = rate_electron / total_rate;
		if(sum > xi)
			return -1;
		// Nuclei
		for(unsigned int i = 0; i < solar_model.target_isotopes.size(); i++)
		{
			sum += rate_nuclei_cache[i] / total_rate;
			if(sum > xi || i == solar_model.target_isotopes.size() - 1)
				return i;
		}
		throw std::runtime_error("Sample_Target(): no target could be sampled.");
	}
}

libphysica::Vector Trajectory_Simulator::Sample_Target_Velocity(double temperature, double target_mass, const libphysica::Vector& vel_DM)
{
	// Sampling algorithm taken from Romano & Walsh, "An improved target velocity sampling algorithm for free gas elastic scattering"
	double kappa = sqrt(target_mass / 2.0 / temperature);
	double vDM	 = vel_DM.Norm();
	if(!std::isfinite(vDM) || vDM <= 0.0)
		throw std::runtime_error("Sample_Target_Velocity(): DM speed is non-positive or non-finite.");
	// 1. Sample target speed vT and mu = cos alpha
	double y  = kappa * vDM;
	double x  = y;
	double mu = 1.0;
	while(sqrt(x * x + y * y - 2.0 * x * y * mu) / (x + y) < libphysica::Sample_Uniform(PRNG, 0.0, 1.0))
	{
		mu			= libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		double xi_1 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
		if(xi_1 < 2.0 / (sqrt(M_PI) * y + 2.0))
		{
			double xi_2 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_3 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double z	= -log(xi_2 * xi_3);
			x			= sqrt(z);
		}
		else
		{
			double xi_2 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_3 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double xi_4 = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
			double z	= -log(xi_2) - pow(cos(M_PI / 2.0 * xi_3), 2.0) * log(xi_4);
			x			= sqrt(z);
		}
	}
	// 2. Construct the target velocity vel_T
	double vT		 = x / kappa;
	double cos_theta = mu;
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	double phi		 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double cos_phi	 = cos(phi);
	double sin_phi	 = sin(phi);

	libphysica::Vector u = vel_DM.Normalized();
	libphysica::Vector reference = (std::fabs(u[2]) < 0.9)
	                             ? libphysica::Vector({0.0, 0.0, 1.0})
	                             : libphysica::Vector({1.0, 0.0, 0.0});
	libphysica::Vector e1 = reference.Cross(u).Normalized();
	libphysica::Vector e2 = u.Cross(e1).Normalized();
	libphysica::Vector unit_vector_T = cos_theta * u
	                                 + sin_theta * cos_phi * e1
	                                 + sin_theta * sin_phi * e2;

	return vT * unit_vector_T;
}

libphysica::Vector Trajectory_Simulator::New_DM_Velocity(double cos_scattering_angle, double DM_mass, double target_mass, libphysica::Vector& vel_DM, libphysica::Vector& vel_target)
{
	// Construction of n, the unit vector pointing into the direction of vfinal.
	double phi			 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	libphysica::Vector n = libphysica::Spherical_Coordinates(1.0, acos(cos_scattering_angle), phi, vel_DM);

	double relative_speed = (vel_target - vel_DM).Norm();

	return target_mass * relative_speed / (target_mass + DM_mass) * n + (DM_mass * vel_DM + target_mass * vel_target) / (target_mass + DM_mass);
}

void Trajectory_Simulator::Scatter(Event& current_event, obscura::DM_Particle& DM)
{
	double r = current_event.Radius();
	double v = current_event.Speed();
	// 1. Find target properties.
	int target_index = Sample_Target(DM, r, v);

	double target_mass;
	if(target_index == -1)
		target_mass = mElectron;
	else
		target_mass = solar_model.target_isotopes[target_index].mass;

	libphysica::Vector vel_target = Sample_Target_Velocity(solar_model.Temperature(r), target_mass, current_event.velocity);

	// 2. Sample the scattering angle
	double cos_alpha = (target_index == -1) ? DM.Sample_Scattering_Angle_Electron(PRNG, v, r) : DM.Sample_Scattering_Angle_Nucleus(PRNG, solar_model.target_isotopes[target_index], v, r);

	// 3. Construct the final DM velocity
	current_event.velocity = New_DM_Velocity(cos_alpha, DM.mass, target_mass, current_event.velocity, vel_target);
}

void Trajectory_Simulator::Fix_PRNG_Seed(unsigned int fixed_seed)
{
	PRNG.seed(fixed_seed);
}

Trajectory_Result Trajectory_Simulator::Simulate(const Event& initial_condition, obscura::DM_Particle& DM, unsigned int mpi_rank)
{
	Event current_event = initial_condition;
	long unsigned int number_of_scatterings = 0;

	// Initialize per-trajectory bincount
	current_bincount = TrajectoryBincount();
	prev_time_sec = -1.0;
	prev_r_km = 0.0;
	prev_v2_km2s2 = 0.0;
	prev_dt_sec = 0.0;
	previous_capture_energy_eV = Capture_Energy_eV(current_event.Radius(), current_event.Speed(), DM);
	free_flight_reference_energy_eV = std::numeric_limits<double>::quiet_NaN();
	current_physical_bound_state = false;
	current_mpi_rank = mpi_rank;
	current_trajectory_id++;
	trajectory_in_progress = true;
	current_trajectory_wall_start = std::chrono::steady_clock::now();
	accumulated_snapshot_overhead_sec = 0.0;
	Publish_Snapshot_Progress();

	TrajectoryTerminationReason termination_reason = TrajectoryTerminationReason::Unknown;
	while(number_of_scatterings < maximum_scatterings)
	{
		TrajectoryTerminationReason propagation_reason = Propagate_Freely(current_event, DM);
		if(propagation_reason == TrajectoryTerminationReason::Scatter && current_event.Radius() < rSun)
		{
			try
			{
				Scatter(current_event, DM);
			}
			catch(const std::exception& error)
			{
				std::cerr << "\nWarning in Trajectory_Simulator::Simulate(): " << error.what()
				          << " Marking trajectory as numerically failed." << std::endl;
				termination_reason = TrajectoryTerminationReason::NumericalFailure;
				break;
			}
			number_of_scatterings++;
			bool captured_after_scatter = Update_Capture_State(current_event.Radius(), current_event.Speed(), current_event.time, DM, true);
			Reset_Bincount_Anchor(current_event);
			if(terminate_on_capture && captured_after_scatter)
			{
				termination_reason = TrajectoryTerminationReason::CaptureMode;
				break;
			}
		}
		else
		{
			termination_reason = propagation_reason;
			break;
		}
	}

	if(termination_reason == TrajectoryTerminationReason::Unknown && number_of_scatterings >= maximum_scatterings)
		termination_reason = TrajectoryTerminationReason::MaxScatterings;

	trajectory_in_progress = false;
	current_bincount.termination_reason = termination_reason;
	current_bincount.number_of_scatterings = number_of_scatterings;
	current_bincount.t_termination = In_Units(current_event.time, sec);

	if(current_bincount.is_captured)
	{
		double r_final = current_event.Radius();
		double v_final = current_event.Speed();
		double vesc_final = solar_model.Local_Escape_Speed(r_final);
		double E_final = 0.5 * DM.mass * (v_final * v_final - vesc_final * vesc_final);
		double E_final_eV = In_Units(E_final, eV);
		if(termination_reason == TrajectoryTerminationReason::OutwardEscape && E_final_eV >= 0.0)
		{
			if(std::isfinite(current_bincount.t_final_unbinding_scatter))
			{
				current_bincount.boundary_escape_observed = true;
				current_bincount.t_boundary_escape = current_bincount.t_termination;
				current_bincount.event_observed = true;
			}
			else
			{
				termination_reason = TrajectoryTerminationReason::EnergyDriftEscape;
				current_bincount.termination_reason = termination_reason;
				current_bincount.survival_valid = false;
				current_bincount.numerically_invalid_escape = true;
				current_bincount.boundary_escape_observed = false;
				current_bincount.event_observed = false;
			}
		}
		if(!current_bincount.event_observed)
			current_bincount.t_final_unbinding_scatter = std::numeric_limits<double>::quiet_NaN();
		current_bincount.truncated = !current_bincount.event_observed;
	}

	if(TrajectoryTerminationInvalidatesSurvival(termination_reason))
		current_bincount.survival_valid = false;

	return Trajectory_Result(initial_condition, current_event, number_of_scatterings, current_bincount);
}  

// 3. Equation of motion solution with Runge-Kutta-Fehlberg
Free_Particle_Propagator::Free_Particle_Propagator(const Event& event)
{
	// 1. Coordinate system
	axis_x = event.position.Normalized();
	axis_z = event.position.Cross(event.velocity).Normalized();
	axis_y = axis_z.Cross(axis_x);

	// 2. Coordinates
	time			 = event.time;
	radius			 = event.Radius();
	phi				 = 0.0;
	v_radial		 = (radius == 0) ? event.Speed() : event.position.Dot(event.velocity) / radius;
	angular_momentum = (event.position.Cross(event.velocity)).Dot(axis_z);

	// 3. Error tolerances (fixed-size array, no heap allocation)
	error_tolerances[0] = 1.0 * km;
	error_tolerances[1] = 1.0e-3 * km / sec;
	error_tolerances[2] = 1.0e-7;
}

double Free_Particle_Propagator::dr_dt(double v)
{
	return v;
}

double Free_Particle_Propagator::dv_dt(double r, double mass)
{
	// 防止 RK45 中间阶段 r_i ≤ 0（向内飞越球心时的数值越界）
	// 导致离心项 L²/r³ 符号翻转从而使误差发散、步长无限缩小。
	// 1 km 远小于任何真实轨迹的转折半径（~10⁴ km），不影响物理结果。
	if(r < 1.0 * km) r = 1.0 * km;
	return angular_momentum * angular_momentum / r / r / r - G_Newton * mass / r / r;
}

double Free_Particle_Propagator::dphi_dt(double r)
{
	return angular_momentum / r / r;
}

bool Free_Particle_Propagator::Runge_Kutta_45_Step(Solar_Model& solar_model)
{
	time_step = RK45_Sanitized_Time_Step(time_step);
	bool accepted = false;
	int inner_iter = 0;
	while(!accepted)
	{
	inner_iter++;
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];
	double r_i;  // intermediate radius for each RK stage

	// Stage 0
	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, solar_model.Mass(radius));
	k_p[0] = time_step * dphi_dt(radius);

	// Stage 1
	r_i = radius + k_r[0] / 4.0;
	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	// Stage 2
	r_i = radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1];
	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[2] = time_step * dphi_dt(r_i);

	// Stage 3
	r_i = radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2];
	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[3] = time_step * dphi_dt(r_i);

	// Stage 4
	r_i = radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3];
	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[4] = time_step * dphi_dt(r_i);

	// Stage 5
	r_i = radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4];
	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(r_i, solar_model.Mass(r_i));
	k_p[5] = time_step * dphi_dt(r_i);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	// Error and adapting the time step (B: 栈数组替代堆分配vector)
	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	// Check if errors fall below the tolerance
	// 自适应步长方法：根据误差调整 time_step
	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		// 边界检查：防止半径变为负数（数值误差）
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
		return true;
	}
	else
	{
		// 绝对最小步长保护：防止步长无限缩小导致死循环
		// 当粒子轨迹中间阶段 r_i 为负值时（角动量 L≠0 且接近 r≈0），
		// dv_dt = L²/r³ 会发散，误差无法满足容差，步长每次乘 0.1 直到下溢。
		// 强制接受：确保内层循环始终能返回，snapshot callback 才能被触发。
		static const double abs_min_step = 1.0e-8 * sec;  // 10 纳秒物理时间下限
		const bool reached_min_step  = time_step_new < abs_min_step;
		const bool reached_iter_cap  = inner_iter >= RK45_MAX_INNER_RETRIES;
		if(reached_min_step || reached_iter_cap)
		{
			// NaN 守卫：若本轮 radius_4/v_radial_4/phi_4 已经是非有限数，
			// 绝不能把 NaN 写回传播器状态，否则后续的每一步都会继续崩。
			// 这时保留原状态，只把时间推进 abs_min_step，由外层 Propagate_Freely
			// 的 is_finite 检查来中止本条轨迹。
		const bool state_ok =
		    std::isfinite(radius_4) && std::isfinite(v_radial_4) && std::isfinite(phi_4);
		time = time + abs_min_step;
		if(state_ok)
		{
			radius = radius_4;
			if(radius < 0.0)
				radius = 0.0;
			v_radial = v_radial_4;
			phi      = phi_4;
		}
		time_step = abs_min_step;
		accepted  = true;
		return false;
	}
	else
	{
		time_step = time_step_new;
		// Loop continues with smaller time_step (was recursive call before Q1 fix)
	}
	}
	}   // end while(!accepted)
	return true;
}

// Overload with constant mass (for Kepler orbits outside the sun / testing)
bool Free_Particle_Propagator::Runge_Kutta_45_Step(double constant_mass)
{
	time_step = RK45_Sanitized_Time_Step(time_step);
	bool accepted = false;
	int inner_iter = 0;
	while(!accepted)
	{
	inner_iter++;
	double k_r[6];
	double k_v[6];
	double k_p[6];

	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, constant_mass);
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, constant_mass);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], constant_mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], constant_mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], constant_mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], constant_mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4104.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4104.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4104.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	double errors[3] = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	double time_step_new = RK45_Next_Step_Size(time_step, errors, error_tolerances);

	if(RK45_Errors_Within_Tolerance(errors, error_tolerances))
	{
		time	  = time + time_step;
		radius	  = radius_4;
		if(radius < 0.0)
			radius = 0.0;
		v_radial  = v_radial_4;
		phi		  = phi_4;
		time_step = time_step_new;
		accepted  = true;
		return true;
	}
	else
	{
		static const double abs_min_step = 1.0e-8 * sec;  // 10 纳秒物理时间下限（同 Solar Model 版本）
		const bool reached_min_step  = time_step_new < abs_min_step;
		const bool reached_iter_cap  = inner_iter >= RK45_MAX_INNER_RETRIES;
		if(reached_min_step || reached_iter_cap)
		{
		const bool state_ok =
		    std::isfinite(radius_4) && std::isfinite(v_radial_4) && std::isfinite(phi_4);
		time = time + abs_min_step;
		if(state_ok)
		{
			radius = radius_4;
			if(radius < 0.0)
				radius = 0.0;
			v_radial = v_radial_4;
			phi      = phi_4;
		}
		time_step = abs_min_step;
		accepted  = true;
		return false;
	}
	else
	{
		time_step = time_step_new;
	}
	}
	}
	return true;
}

double Free_Particle_Propagator::Current_Time()
{
	return time;
}

double Free_Particle_Propagator::Current_Radius()
{
	// 边界检查：确保返回非负半径
	return (radius < 0.0) ? 0.0 : radius;
}

double Free_Particle_Propagator::Current_Speed()
{
	// 边界检查：确保 radius 非负
	double r = (radius < 0.0) ? 0.0 : radius;
	if(r == 0 || angular_momentum == 0)
		return fabs(v_radial);  // 返回速度绝对值
	else
		return sqrt(v_radial * v_radial + angular_momentum * angular_momentum / r / r);
}

Event Free_Particle_Propagator::Event_In_3D()
{
	double v_phi			= angular_momentum / pow(radius, 2);
	libphysica::Vector xNew = radius * (cos(phi) * axis_x + sin(phi) * axis_y);
	libphysica::Vector vNew = (v_radial * cos(phi) - v_phi * radius * sin(phi)) * axis_x + (v_radial * sin(phi) + radius * v_phi * cos(phi)) * axis_y;

	return Event(time, xNew, vNew);
}

}	// namespace DaMaSCUS_SUN

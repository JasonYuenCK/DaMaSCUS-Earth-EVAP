#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <vector>
#include <string>
#include <array>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/DM_Particle.hpp"

#include "Simulation_Trajectory.hpp"

namespace DaMaSCUS_SUN
{

// Survival-analysis record for every captured trajectory.
struct EvaporationRecord
{
	int rank = -1;
	unsigned long int trajectory_id = 0;
	double t_evap = 0.0;    // finite only for observed unbinding events [seconds]
	double t_capture = -1.0;
	double t_final_unbinding_scatter = -1.0;
	double t_boundary_escape = -1.0;
	double t_termination = -1.0;
	double observed_lifetime = 0.0;
	double lifetime_unbinding = -1.0;
	double lifetime_boundary = -1.0;
	double r_first_negative_km = -1.0;
	double E_first_negative_eV = 0.0;
	double dE_first_negative_from_prev_eV = 0.0;
	bool event_observed = false;
	bool boundary_escape_observed = false;
	bool survival_valid = true;
	bool numerically_invalid_escape = false;
	bool censored = true;
	bool truncated = true; // compatibility alias for censored
	TrajectoryTerminationReason termination_reason = TrajectoryTerminationReason::Unknown;
	double max_free_energy_drift_eV = 0.0;
	double max_free_energy_drift_rel = 0.0;
	unsigned long int number_of_scatterings = 0;
};

struct EvaporationModeBincount
{
	unsigned long int count = 0;
	unsigned long int truncated_count = 0;
	std::array<double, NUM_BINS> dt_hist{};
	std::array<double, NUM_BINS> v2dt_hist{};
	std::array<double, NUM_BINS> dt_sq_hist{};
	std::array<double, NUM_BINS> v2dt_sq_hist{};
};

class Simulation_Data
{
  private:
	// Configuration
	unsigned int target_captured_per_rank;   // ceil(sample_size / N_ranks)
	unsigned long int max_trajectories_per_rank;
	double initial_and_final_radius = 2.0 * libphysica::natural_units::rSun;
	unsigned int minimum_number_of_scatterings = 1;
	unsigned long int maximum_number_of_scatterings = DEFAULT_MAXIMUM_SCATTERINGS;
	unsigned long int maximum_free_time_steps = DEFAULT_MAXIMUM_FREE_TIME_STEPS;

	// Results
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_reflected_particles;
	unsigned long int number_of_captured_particles;
	unsigned long int number_of_complete_evaporation_particles;
	unsigned long int number_of_censored_captured_particles;
	unsigned long int number_of_invalid_survival_captured_particles;
	std::array<unsigned long int, TRAJECTORY_TERMINATION_REASON_COUNT> termination_reason_counts;
	double average_number_of_scatterings;
	double computing_time;
	bool early_stopped;

	// Aggregated bincount histograms
	std::array<double, NUM_BINS> captured_dt_hist;
	std::array<double, NUM_BINS> captured_v2dt_hist;
	std::array<double, NUM_BINS> not_captured_dt_hist;
	std::array<double, NUM_BINS> not_captured_v2dt_hist;

	// Per-bin sum of squares for error estimation
	std::array<double, NUM_BINS> captured_dt_sq_hist;      // Σ (per-traj dt)²
	std::array<double, NUM_BINS> captured_v2dt_sq_hist;    // Σ (per-traj v²dt)²
	std::array<double, NUM_BINS> not_captured_dt_sq_hist;
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist;

	// Evaporation records
	std::vector<EvaporationRecord> evaporation_records;
	bool evaporation_mode_bincount_enabled = false;
	std::vector<double> evaporation_mode_boundaries_log10_s;
	std::vector<std::string> evaporation_mode_labels;
	bool evaporation_mode_include_truncated = false;
	std::vector<EvaporationModeBincount> evaporation_mode_bincounts;

	// --- Per-trajectory computation time statistics ---
	double total_wall_time_captured;
	double total_wall_time_not_captured;

	// Wall-clock time log-histogram: 10^{-4} ~ 10^{8} s, 10 bins per decade, 120 bins
	static constexpr int WALL_TIME_BINS = 120;
	static constexpr double WALL_TIME_LOG_MIN = -4.0;
	static constexpr double WALL_TIME_LOG_MAX = 8.0;
	std::array<unsigned long int, WALL_TIME_BINS> wall_time_hist_captured;
	std::array<unsigned long int, WALL_TIME_BINS> wall_time_hist_not_captured;
	unsigned long int wall_time_overflow_captured;
	unsigned long int wall_time_overflow_not_captured;

	// --- Per-trajectory RK45 step count statistics ---
	unsigned long int total_rk45_steps_captured;
	unsigned long int total_rk45_steps_not_captured;

	// RK45 step count log-histogram: 10^{0} ~ 10^{14}, 10 bins per decade, 140 bins
	static constexpr int STEP_COUNT_BINS = 140;
	static constexpr double STEP_COUNT_LOG_MIN = 0.0;
	static constexpr double STEP_COUNT_LOG_MAX = 14.0;
	std::array<unsigned long int, STEP_COUNT_BINS> step_count_hist_captured;
	std::array<unsigned long int, STEP_COUNT_BINS> step_count_hist_not_captured;
	unsigned long int step_count_overflow_captured;
	unsigned long int step_count_overflow_not_captured;

	// MPI
	int mpi_rank, mpi_processes;
	void Perform_MPI_Reductions();

	// For reflection spectrum compatibility (kept but not actively used in new logic)
	unsigned int isoreflection_rings;
	double minimum_speed_threshold;
	std::vector<unsigned long int> number_of_data_points;
	double KDE_boundary_correction_factor = 0.75;

  public:
	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min = 0.0, unsigned int iso_rings = 1);

	void Configure(double initial_radius, unsigned int min_scattering, unsigned long int max_scattering, unsigned long int max_free_steps = DEFAULT_MAXIMUM_FREE_TIME_STEPS);
	void Configure_Evaporation_Mode_Bincount(bool enabled, const std::vector<double>& boundaries_log10_s, const std::vector<std::string>& labels, bool include_truncated);

	void Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg = SnapshotConfig(), unsigned int fixed_seed = 0, bool capture_mode = false);

	// Output files
	void Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM);

	double Free_Ratio() const;
	double Capture_Ratio() const;
	double Reflection_Ratio(int isoreflection_ring = -1) const;

	double Minimum_Speed() const;
	double Lowest_Speed(unsigned int iso_ring = 0) const;
	double Highest_Speed(unsigned int iso_ring = 0) const;

	void Print_Summary(unsigned int mpi_rank = 0);
	void Print_Capture_Mode_Summary(unsigned int mpi_rank = 0);
};
}	// namespace DaMaSCUS_SUN
#endif

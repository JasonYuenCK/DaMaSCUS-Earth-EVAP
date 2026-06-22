#include "Data_Generation.hpp"

#include <dirent.h>
#include <errno.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"

namespace DaMaSCUS_SUN
{

using namespace libphysica::natural_units;

namespace
{
constexpr int SNAPSHOT_WALL_TIME_BINS = 120;
constexpr int SNAPSHOT_STEP_COUNT_BINS = 140;

bool Has_Positive_Evaporation_Time(double t_evap)
{
	return std::isfinite(t_evap) && t_evap > 0.0;
}

bool Has_Positive_Evaporation_Time(const EvaporationRecord& rec)
{
	return Has_Positive_Evaporation_Time(rec.observed_lifetime);
}

bool Is_Completed_Evaporation_Record(const EvaporationRecord& rec)
{
	return rec.survival_valid
	    && rec.event_observed
	    && std::isfinite(rec.lifetime_unbinding)
	    && rec.lifetime_unbinding >= 0.0;
}

const char* TerminationReason_Name(TrajectoryTerminationReason reason)
{
	switch(reason)
	{
		case TrajectoryTerminationReason::OutwardEscape: return "outward_escape";
		case TrajectoryTerminationReason::Scatter: return "scatter";
		case TrajectoryTerminationReason::WallTimeLimit: return "wall_time_limit";
		case TrajectoryTerminationReason::MaxFreeSteps: return "max_free_steps";
		case TrajectoryTerminationReason::MaxScatterings: return "max_scatterings";
		case TrajectoryTerminationReason::NonFiniteState: return "non_finite_state";
		case TrajectoryTerminationReason::SpeedLimit: return "speed_limit";
		case TrajectoryTerminationReason::NumericalFailure: return "numerical_failure";
		case TrajectoryTerminationReason::CaptureMode: return "capture_mode";
		case TrajectoryTerminationReason::EnergyDriftEscape: return "energy_drift_escape";
		case TrajectoryTerminationReason::Unknown:
		default: return "unknown";
	}
}

int TerminationReason_Index(TrajectoryTerminationReason reason)
{
	int index = static_cast<int>(reason);
	if(index < 0 || index >= TRAJECTORY_TERMINATION_REASON_COUNT)
		return static_cast<int>(TrajectoryTerminationReason::Unknown);
	return index;
}

bool Completed_Outward_Escape(TrajectoryTerminationReason reason)
{
	return reason == TrajectoryTerminationReason::OutwardEscape;
}

bool Incomplete_Physical_Termination(TrajectoryTerminationReason reason)
{
	return reason != TrajectoryTerminationReason::OutwardEscape
	    && reason != TrajectoryTerminationReason::CaptureMode;
}

unsigned long int Count_Incomplete_Physical_Terminations(const std::array<unsigned long int, TRAJECTORY_TERMINATION_REASON_COUNT>& counts)
{
	unsigned long int total = 0;
	for(int reason_index = 0; reason_index < TRAJECTORY_TERMINATION_REASON_COUNT; reason_index++)
	{
		TrajectoryTerminationReason reason = static_cast<TrajectoryTerminationReason>(reason_index);
		if(Incomplete_Physical_Termination(reason))
			total += counts[reason_index];
	}
	return total;
}

int Classify_Evaporation_Mode(double log10_t_evap, const std::vector<double>& boundaries_log10_s)
{
	if(!std::isfinite(log10_t_evap))
		return -1;
	int mode = 0;
	while(mode < static_cast<int>(boundaries_log10_s.size()) && log10_t_evap >= boundaries_log10_s[mode])
		mode++;
	return mode;
}

void Accumulate_Mode_Bincount(EvaporationModeBincount& mode_bincount, const TrajectoryBincount& bincount)
{
	mode_bincount.count++;
	if(bincount.truncated)
		mode_bincount.truncated_count++;
	for(int bin = 0; bin < NUM_BINS; bin++)
	{
		mode_bincount.dt_hist[bin] += bincount.dt_hist[bin];
		mode_bincount.v2dt_hist[bin] += bincount.v2dt_hist[bin];
		mode_bincount.dt_sq_hist[bin] += bincount.dt_hist[bin] * bincount.dt_hist[bin];
		mode_bincount.v2dt_sq_hist[bin] += bincount.v2dt_hist[bin] * bincount.v2dt_hist[bin];
	}
}

bool Build_Evaporation_Record(const TrajectoryBincount& bincount, int mpi_rank, unsigned long int trajectory_id, EvaporationRecord& rec)
{
	if(!bincount.is_captured || !std::isfinite(bincount.t_capture))
		return false;

	const bool numerically_invalid_escape =
	    bincount.numerically_invalid_escape || bincount.termination_reason == TrajectoryTerminationReason::EnergyDriftEscape;
	const bool survival_valid = bincount.survival_valid && !numerically_invalid_escape;
	const bool event_observed = survival_valid && bincount.event_observed;
	const double t_termination = std::isfinite(bincount.t_termination) ? bincount.t_termination : bincount.t_capture;
	const double lifetime_unbinding = event_observed ? (bincount.t_final_unbinding_scatter - bincount.t_capture) : -1.0;
	const double lifetime_boundary = bincount.boundary_escape_observed ? (bincount.t_boundary_escape - bincount.t_capture) : -1.0;
	double observed_lifetime = event_observed ? lifetime_unbinding : (t_termination - bincount.t_capture);
	if(!survival_valid)
		observed_lifetime = std::numeric_limits<double>::quiet_NaN();
	else if(!std::isfinite(observed_lifetime) || observed_lifetime < 0.0)
		observed_lifetime = 0.0;

	rec.rank = mpi_rank;
	rec.trajectory_id = trajectory_id;
	rec.t_evap = event_observed ? lifetime_unbinding : std::numeric_limits<double>::quiet_NaN();
	rec.t_capture = bincount.t_capture;
	rec.t_final_unbinding_scatter = std::isfinite(bincount.t_final_unbinding_scatter) ? bincount.t_final_unbinding_scatter : -1.0;
	rec.t_boundary_escape = std::isfinite(bincount.t_boundary_escape) ? bincount.t_boundary_escape : -1.0;
	rec.t_termination = t_termination;
	rec.observed_lifetime = observed_lifetime;
	rec.lifetime_unbinding = (std::isfinite(lifetime_unbinding) && lifetime_unbinding >= 0.0) ? lifetime_unbinding : -1.0;
	rec.lifetime_boundary = (std::isfinite(lifetime_boundary) && lifetime_boundary >= 0.0) ? lifetime_boundary : -1.0;
	rec.r_first_negative_km = bincount.r_first_negative_km;
	rec.E_first_negative_eV = bincount.E_first_negative_eV;
	rec.dE_first_negative_from_prev_eV = bincount.dE_first_negative_from_prev_eV;
	rec.event_observed = event_observed;
	rec.boundary_escape_observed = survival_valid && bincount.boundary_escape_observed;
	rec.survival_valid = survival_valid;
	rec.numerically_invalid_escape = numerically_invalid_escape;
	rec.censored = survival_valid && !event_observed;
	rec.truncated = rec.censored;
	rec.termination_reason = bincount.termination_reason;
	rec.max_free_energy_drift_eV = bincount.max_free_energy_drift_eV;
	rec.max_free_energy_drift_rel = bincount.max_free_energy_drift_rel;
	rec.number_of_scatterings = bincount.number_of_scatterings;
	return true;
}

struct Rank_Snapshot_State
{
	uint64_t run_id = 0;
	int32_t snapshot_index = 0;
	int32_t rank = 0;
	int32_t done = 0;
	int32_t trajectory_in_progress = 0;
	uint64_t local_captured = 0;
	uint64_t local_total = 0;
	uint64_t current_trajectory_id = 0;
	double rank_elapsed_wall_sec = 0.0;
	double current_trajectory_wall_sec = 0.0;
	double current_trajectory_physical_sec = 0.0;
	int32_t current_trajectory_captured = 0;
	std::array<double, NUM_BINS> current_trajectory_dt_hist{};
	std::array<double, NUM_BINS> current_trajectory_v2dt_hist{};
	std::array<double, NUM_BINS> captured_dt_hist{};
	std::array<double, NUM_BINS> captured_v2dt_hist{};
	std::array<double, NUM_BINS> captured_dt_sq_hist{};
	std::array<double, NUM_BINS> captured_v2dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_dt_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_hist{};
	std::array<double, NUM_BINS> not_captured_dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist{};
	double total_wall_time_captured = 0.0;
	double total_wall_time_not_captured = 0.0;
	uint64_t total_rk45_steps_captured = 0;
	uint64_t total_rk45_steps_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_captured{};
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_not_captured{};
	uint64_t wall_time_overflow_captured = 0;
	uint64_t wall_time_overflow_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_captured{};
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_not_captured{};
	uint64_t step_count_overflow_captured = 0;
	uint64_t step_count_overflow_not_captured = 0;
	std::vector<EvaporationRecord> completed_evaporation_records;
};

struct Snapshot_Report_State
{
	long long snapshot_time_label = 0;
	double snapshot_interval_seconds = 0.0;
	uint64_t total_trajectories = 0;
	uint64_t captured_particles = 0;
	uint64_t snapshot_bincount_captured_samples = 0;
	uint64_t snapshot_bincount_not_captured_samples = 0;
	std::array<double, NUM_BINS> captured_dt_hist{};
	std::array<double, NUM_BINS> captured_v2dt_hist{};
	std::array<double, NUM_BINS> captured_dt_sq_hist{};
	std::array<double, NUM_BINS> captured_v2dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_dt_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_hist{};
	std::array<double, NUM_BINS> not_captured_dt_sq_hist{};
	std::array<double, NUM_BINS> not_captured_v2dt_sq_hist{};
	double total_wall_time_captured = 0.0;
	double total_wall_time_not_captured = 0.0;
	uint64_t total_rk45_steps_captured = 0;
	uint64_t total_rk45_steps_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_captured{};
	std::array<uint64_t, SNAPSHOT_WALL_TIME_BINS> wall_time_hist_not_captured{};
	uint64_t wall_time_overflow_captured = 0;
	uint64_t wall_time_overflow_not_captured = 0;
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_captured{};
	std::array<uint64_t, SNAPSHOT_STEP_COUNT_BINS> step_count_hist_not_captured{};
	uint64_t step_count_overflow_captured = 0;
	uint64_t step_count_overflow_not_captured = 0;
	std::vector<Rank_Snapshot_State> rank_states;
	std::vector<EvaporationRecord> completed_evaporation_records;
};

struct Snapshot_Load_Diagnostics
{
	int ready_ranks = 0;
	std::vector<std::string> per_rank_status;
};

std::string Format_Rank_Status(const Rank_Snapshot_State& state);

template<typename T, size_t N>
bool Find_Nonzero_Bin_Range(const std::array<T, N>& captured_hist, const std::array<T, N>& not_captured_hist, int& first_bin, int& last_bin)
{
	first_bin = 0;
	while(first_bin < static_cast<int>(N) && captured_hist[first_bin] == 0 && not_captured_hist[first_bin] == 0)
		first_bin++;

	last_bin = static_cast<int>(N) - 1;
	while(last_bin >= first_bin && captured_hist[last_bin] == 0 && not_captured_hist[last_bin] == 0)
		last_bin--;

	return first_bin <= last_bin;
}

template<typename T>
void Write_Binary_Value(std::ofstream& file, const T& value)
{
	file.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void Write_Binary_Array(std::ofstream& file, const std::array<T, N>& values)
{
	file.write(reinterpret_cast<const char*>(values.data()), N * sizeof(T));
}

template<typename T>
void Read_Binary_Value(std::ifstream& file, T& value)
{
	file.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void Read_Binary_Array(std::ifstream& file, std::array<T, N>& values)
{
	file.read(reinterpret_cast<char*>(values.data()), N * sizeof(T));
}

void Write_Evaporation_Record_Binary(std::ofstream& file, const EvaporationRecord& rec)
{
	const int32_t rank = static_cast<int32_t>(rec.rank);
	const uint64_t trajectory_id = static_cast<uint64_t>(rec.trajectory_id);
	const int32_t event_observed = rec.event_observed ? 1 : 0;
	const int32_t boundary_escape_observed = rec.boundary_escape_observed ? 1 : 0;
	const int32_t survival_valid = rec.survival_valid ? 1 : 0;
	const int32_t numerically_invalid_escape = rec.numerically_invalid_escape ? 1 : 0;
	const int32_t censored = rec.censored ? 1 : 0;
	const int32_t truncated = rec.truncated ? 1 : 0;
	const int32_t termination_reason = static_cast<int32_t>(TerminationReason_Index(rec.termination_reason));
	const uint64_t number_of_scatterings = static_cast<uint64_t>(rec.number_of_scatterings);

	Write_Binary_Value(file, rank);
	Write_Binary_Value(file, trajectory_id);
	Write_Binary_Value(file, rec.t_evap);
	Write_Binary_Value(file, rec.t_capture);
	Write_Binary_Value(file, rec.t_final_unbinding_scatter);
	Write_Binary_Value(file, rec.t_boundary_escape);
	Write_Binary_Value(file, rec.t_termination);
	Write_Binary_Value(file, rec.observed_lifetime);
	Write_Binary_Value(file, rec.lifetime_unbinding);
	Write_Binary_Value(file, rec.lifetime_boundary);
	Write_Binary_Value(file, rec.r_first_negative_km);
	Write_Binary_Value(file, rec.E_first_negative_eV);
	Write_Binary_Value(file, rec.dE_first_negative_from_prev_eV);
	Write_Binary_Value(file, event_observed);
	Write_Binary_Value(file, boundary_escape_observed);
	Write_Binary_Value(file, survival_valid);
	Write_Binary_Value(file, numerically_invalid_escape);
	Write_Binary_Value(file, censored);
	Write_Binary_Value(file, truncated);
	Write_Binary_Value(file, termination_reason);
	Write_Binary_Value(file, rec.max_free_energy_drift_eV);
	Write_Binary_Value(file, rec.max_free_energy_drift_rel);
	Write_Binary_Value(file, number_of_scatterings);
}

void Read_Evaporation_Record_Binary(std::ifstream& file, EvaporationRecord& rec)
{
	int32_t rank = 0;
	uint64_t trajectory_id = 0;
	int32_t event_observed = 0;
	int32_t boundary_escape_observed = 0;
	int32_t survival_valid = 0;
	int32_t numerically_invalid_escape = 0;
	int32_t censored = 0;
	int32_t truncated = 0;
	int32_t termination_reason = 0;
	uint64_t number_of_scatterings = 0;

	Read_Binary_Value(file, rank);
	Read_Binary_Value(file, trajectory_id);
	Read_Binary_Value(file, rec.t_evap);
	Read_Binary_Value(file, rec.t_capture);
	Read_Binary_Value(file, rec.t_final_unbinding_scatter);
	Read_Binary_Value(file, rec.t_boundary_escape);
	Read_Binary_Value(file, rec.t_termination);
	Read_Binary_Value(file, rec.observed_lifetime);
	Read_Binary_Value(file, rec.lifetime_unbinding);
	Read_Binary_Value(file, rec.lifetime_boundary);
	Read_Binary_Value(file, rec.r_first_negative_km);
	Read_Binary_Value(file, rec.E_first_negative_eV);
	Read_Binary_Value(file, rec.dE_first_negative_from_prev_eV);
	Read_Binary_Value(file, event_observed);
	Read_Binary_Value(file, boundary_escape_observed);
	Read_Binary_Value(file, survival_valid);
	Read_Binary_Value(file, numerically_invalid_escape);
	Read_Binary_Value(file, censored);
	Read_Binary_Value(file, truncated);
	Read_Binary_Value(file, termination_reason);
	Read_Binary_Value(file, rec.max_free_energy_drift_eV);
	Read_Binary_Value(file, rec.max_free_energy_drift_rel);
	Read_Binary_Value(file, number_of_scatterings);

	rec.rank = static_cast<int>(rank);
	rec.trajectory_id = static_cast<unsigned long int>(trajectory_id);
	rec.event_observed = (event_observed != 0);
	rec.boundary_escape_observed = (boundary_escape_observed != 0);
	rec.survival_valid = (survival_valid != 0);
	rec.numerically_invalid_escape = (numerically_invalid_escape != 0);
	rec.censored = (censored != 0);
	rec.truncated = (truncated != 0);
	rec.termination_reason = static_cast<TrajectoryTerminationReason>(TerminationReason_Index(static_cast<TrajectoryTerminationReason>(termination_reason)));
	rec.number_of_scatterings = static_cast<unsigned long int>(number_of_scatterings);
}

double Capture_Ratio_From_Counts(uint64_t total_trajectories, uint64_t captured_particles)
{
	return (total_trajectories > 0) ? static_cast<double>(captured_particles) / static_cast<double>(total_trajectories) : 0.0;
}

double Snapshot_Bin_Error(double sum, double sum_sq, double count)
{
	if(count <= 1.0)
		return 0.0;

	double mean = sum / count;
	double variance = sum_sq / count - mean * mean;
	if(variance < 0.0)
		variance = 0.0;

	return sqrt(count * variance);
}

bool Has_Bincount_Contribution(const std::array<double, NUM_BINS>& dt_hist, const std::array<double, NUM_BINS>& v2dt_hist)
{
	return std::any_of(dt_hist.begin(), dt_hist.end(), [](double value) { return value != 0.0; })
	    || std::any_of(v2dt_hist.begin(), v2dt_hist.end(), [](double value) { return value != 0.0; });
}

void Write_Report_Header(std::ofstream& file, double mass_gev, double sigma_cm2, uint64_t total_trajectories, uint64_t captured_particles, bool early_stopped, long long snapshot_time_label = -1, double snapshot_interval_seconds = 0.0)
{
	if(snapshot_time_label >= 0)
	{
		file << "# snapshot_target_wall_time_s = " << snapshot_time_label << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << snapshot_interval_seconds << "\n";
	}

	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# total_trajectories = " << total_trajectories << "\n";
	file << "# captured_particles = " << captured_particles << "\n";

	double p_cap = Capture_Ratio_From_Counts(total_trajectories, captured_particles);
	file << "# capture_rate = " << std::fixed << std::setprecision(8) << p_cap << "\n";

	{
		double N = static_cast<double>(total_trajectories);
		double p = p_cap;
		double z = 1.96;
		double sigma_p = (N > 0.0) ? sqrt(p * (1.0 - p) / N) : 0.0;
		file << "# capture_rate_err = " << std::fixed << std::setprecision(8) << sigma_p << "\n";
		if(N > 0.0)
		{
			double denom = 1.0 + z * z / N;
			double center = p + z * z / (2.0 * N);
			double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
			double ci_lower = (center - spread) / denom;
			double ci_upper = (center + spread) / denom;
			if(ci_lower < 0.0) ci_lower = 0.0;
			if(ci_upper > 1.0) ci_upper = 1.0;
			file << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << ci_lower << "\n";
			file << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << ci_upper << "\n";
		}
	}

	if(early_stopped)
		file << "# EARLY_STOP: max_trajectories reached\n";
}

void Write_Evaporation_Record_List(std::ofstream& file, const std::vector<EvaporationRecord>& records)
{
	size_t event_count = 0;
	size_t censored_count = 0;
	size_t boundary_escape_count = 0;
	size_t survival_valid_count = 0;
	size_t numerically_invalid_escape_count = 0;
	size_t positive_observed_count = 0;
	for(const auto& rec : records)
	{
		if(rec.event_observed)
			event_count++;
		else if(rec.censored)
			censored_count++;
		if(rec.boundary_escape_observed)
			boundary_escape_count++;
		if(rec.survival_valid)
			survival_valid_count++;
		if(rec.numerically_invalid_escape)
			numerically_invalid_escape_count++;
		if(rec.survival_valid && Has_Positive_Evaporation_Time(rec))
			positive_observed_count++;
	}

	file << "# evaporation_record_count = " << records.size() << "\n";
	file << "# evaporation_event_observed_count = " << event_count << "\n";
	file << "# evaporation_censored_count = " << censored_count << "\n";
	file << "# evaporation_boundary_escape_count = " << boundary_escape_count << "\n";
	file << "# evaporation_survival_valid_count = " << survival_valid_count << "\n";
	file << "# evaporation_numerically_invalid_escape_count = " << numerically_invalid_escape_count << "\n";
	file << "# evaporation_positive_observed_lifetime_count = " << positive_observed_count << "\n";
	file << "# WARNING: t_evap is a physical event time only when event_observed=1; use observed_lifetime plus event_observed for survival analysis, and require survival_valid=1.\n";
	file << "# rank  trajectory_id  t_evap[s]  observed_lifetime[s]  lifetime_unbinding[s]  lifetime_boundary[s]  t_capture[s]  t_final_unbinding_scatter[s]  t_boundary_escape[s]  t_termination[s]  event_observed(0/1)  boundary_escape_observed(0/1)  survival_valid(0/1)  numerically_invalid_escape(0/1)  censored(0/1)  truncated(0/1)  r_capture[km]  E_capture[eV]  dE_capture_from_prev[eV]  termination_reason  max_free_energy_drift[eV]  max_free_energy_drift_rel  number_of_scatterings\n";
	for(const auto& rec : records)
	{
		file << rec.rank << "\t" << rec.trajectory_id << "\t" << std::scientific << std::setprecision(10)
		     << rec.t_evap << "\t" << rec.observed_lifetime << "\t" << rec.lifetime_unbinding
		     << "\t" << rec.lifetime_boundary << "\t" << rec.t_capture << "\t" << rec.t_final_unbinding_scatter
		     << "\t" << rec.t_boundary_escape << "\t" << rec.t_termination
		     << "\t" << (rec.event_observed ? 1 : 0) << "\t" << (rec.boundary_escape_observed ? 1 : 0)
		     << "\t" << (rec.survival_valid ? 1 : 0) << "\t" << (rec.numerically_invalid_escape ? 1 : 0)
		     << "\t" << (rec.censored ? 1 : 0) << "\t" << (rec.truncated ? 1 : 0)
		     << "\t" << rec.r_first_negative_km << "\t" << rec.E_first_negative_eV
		     << "\t" << rec.dE_first_negative_from_prev_eV << "\t" << TerminationReason_Name(rec.termination_reason)
		     << "\t" << rec.max_free_energy_drift_eV << "\t" << rec.max_free_energy_drift_rel
		     << "\t" << rec.number_of_scatterings << "\n";
	}
}

std::string Format_Physical_Time_Scientific(double physical_time_sec)
{
	std::ostringstream stream;

	if(!std::isfinite(physical_time_sec))
		stream << "nan";
	else
		stream << std::uppercase << std::scientific << std::setprecision(3) << physical_time_sec;

	return stream.str();
}

bool Path_Exists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0;
}

std::string Join_Path(const std::string& directory, const std::string& name)
{
	if(directory.empty() || directory.back() == '/')
		return directory + name;
	return directory + "/" + name;
}

bool Remove_Path_Recursive(const std::string& path)
{
	struct stat info;
	if(lstat(path.c_str(), &info) != 0)
		return errno == ENOENT;

	if(S_ISDIR(info.st_mode))
	{
		DIR* dir = opendir(path.c_str());
		if(dir == NULL)
			return false;

		bool success = true;
		struct dirent* entry = NULL;
		while((entry = readdir(dir)) != NULL)
		{
			const std::string name = entry->d_name;
			if(name == "." || name == "..")
				continue;
			if(!Remove_Path_Recursive(Join_Path(path, name)))
				success = false;
		}
		closedir(dir);

		if(rmdir(path.c_str()) != 0 && errno != ENOENT)
			success = false;

		return success;
	}

	return std::remove(path.c_str()) == 0 || errno == ENOENT;
}

bool Clear_Directory_Contents(const std::string& directory)
{
	DIR* dir = opendir(directory.c_str());
	if(dir == NULL)
		return errno == ENOENT;

	bool success = true;
	struct dirent* entry = NULL;
	while((entry = readdir(dir)) != NULL)
	{
		const std::string name = entry->d_name;
		if(name == "." || name == "..")
			continue;
		if(!Remove_Path_Recursive(Join_Path(directory, name)))
			success = false;
	}
	closedir(dir);
	return success;
}

void Ensure_Directory_Exists(const std::string& directory)
{
	struct stat info;
	if(stat(directory.c_str(), &info) != 0)
		mkdir(directory.c_str(), 0755);
}

long long Snapshot_Time_Label_Seconds(int snapshot_index, double interval_seconds)
{
	return static_cast<long long>(std::llround(snapshot_index * interval_seconds));
}

std::string Snapshot_Text_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s.txt";
}

std::string Snapshot_Completed_Evaporation_Diagnostic_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_completed_evaporation_diagnostic.txt";
}

void Remove_Stale_Snapshot_Evaporation_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	std::string path = snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_evaporation.txt";
	std::remove(path.c_str());
}

void Write_Snapshot_Diagnostics(std::ofstream& file, int snapshot_index, double interval_seconds, int caller_rank, const Snapshot_Load_Diagnostics& diagnostics, bool merged)
{
	file << "# snapshot_status = " << (merged ? "merged" : "waiting") << "\n";
	file << "# snapshot_target_wall_time_s = " << Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds) << "\n";
	file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
	file << "# attempted_by_rank = " << caller_rank << "\n";
	file << "# ready_ranks = " << diagnostics.ready_ranks << " / " << diagnostics.per_rank_status.size() << "\n";
	file << "#\n";
	file << "# rank_index  source_or_wait_reason_and_status\n";
	for(size_t rank = 0; rank < diagnostics.per_rank_status.size(); rank++)
		file << rank << "\t" << diagnostics.per_rank_status[rank] << "\n";
}

bool Snapshot_Text_File_Is_Merged(const std::string& path)
{
	std::ifstream file(path);
	if(!file.is_open())
		return false;

	std::string line;
	while(std::getline(file, line))
	{
		if(line == "# snapshot_status = merged")
			return true;
		if(!line.empty() && line[0] != '#')
			break;
	}

	return false;
}

std::string Rank_Snapshot_Checkpoint_Path(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds)
{
	return rank_snapshot_dir + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_rank" + std::to_string(rank) + ".bin";
}

std::string Rank_Snapshot_Final_Path(const std::string& rank_snapshot_dir, int rank)
{
	return rank_snapshot_dir + "rank" + std::to_string(rank) + "_final.bin";
}

void Cleanup_Snapshot_Checkpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds).c_str());
}

void Cleanup_Final_Snapshot_States(const std::string& rank_snapshot_dir, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(Rank_Snapshot_Final_Path(rank_snapshot_dir, rank).c_str());

	rmdir(rank_snapshot_dir.c_str());
}

bool Write_Text_File_Atomically(const std::string& path, int unique_tag, const std::function<void(std::ofstream&)>& writer)
{
	std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(unique_tag);
	std::ofstream file(tmp_path, std::ios::trunc);
	if(!file.is_open())
		return false;

	writer(file);
	file.close();

	if(!file.good())
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	if(std::rename(tmp_path.c_str(), path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	return true;
}

bool Write_Rank_Snapshot_State(const std::string& path, const Rank_Snapshot_State& state)
{
	std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(state.rank);
	std::ofstream file(tmp_path, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		return false;

	Write_Binary_Value(file, state.run_id);
	Write_Binary_Value(file, state.snapshot_index);
	Write_Binary_Value(file, state.rank);
	Write_Binary_Value(file, state.done);
	Write_Binary_Value(file, state.trajectory_in_progress);
	Write_Binary_Value(file, state.local_captured);
	Write_Binary_Value(file, state.local_total);
	Write_Binary_Value(file, state.current_trajectory_id);
	Write_Binary_Value(file, state.rank_elapsed_wall_sec);
	Write_Binary_Value(file, state.current_trajectory_wall_sec);
	Write_Binary_Value(file, state.current_trajectory_physical_sec);
	Write_Binary_Value(file, state.current_trajectory_captured);
	Write_Binary_Array(file, state.current_trajectory_dt_hist);
	Write_Binary_Array(file, state.current_trajectory_v2dt_hist);
	Write_Binary_Array(file, state.captured_dt_hist);
	Write_Binary_Array(file, state.captured_v2dt_hist);
	Write_Binary_Array(file, state.captured_dt_sq_hist);
	Write_Binary_Array(file, state.captured_v2dt_sq_hist);
	Write_Binary_Array(file, state.not_captured_dt_hist);
	Write_Binary_Array(file, state.not_captured_v2dt_hist);
	Write_Binary_Array(file, state.not_captured_dt_sq_hist);
	Write_Binary_Array(file, state.not_captured_v2dt_sq_hist);
	Write_Binary_Value(file, state.total_wall_time_captured);
	Write_Binary_Value(file, state.total_wall_time_not_captured);
	Write_Binary_Value(file, state.total_rk45_steps_captured);
	Write_Binary_Value(file, state.total_rk45_steps_not_captured);
	Write_Binary_Array(file, state.wall_time_hist_captured);
	Write_Binary_Array(file, state.wall_time_hist_not_captured);
	Write_Binary_Value(file, state.wall_time_overflow_captured);
	Write_Binary_Value(file, state.wall_time_overflow_not_captured);
	Write_Binary_Array(file, state.step_count_hist_captured);
	Write_Binary_Array(file, state.step_count_hist_not_captured);
	Write_Binary_Value(file, state.step_count_overflow_captured);
	Write_Binary_Value(file, state.step_count_overflow_not_captured);
	const uint64_t completed_evaporation_record_count = static_cast<uint64_t>(state.completed_evaporation_records.size());
	Write_Binary_Value(file, completed_evaporation_record_count);
	for(const auto& rec : state.completed_evaporation_records)
		Write_Evaporation_Record_Binary(file, rec);
	file.close();

	if(!file.good())
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	if(std::rename(tmp_path.c_str(), path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	return true;
}

bool Read_Rank_Snapshot_State(const std::string& path, uint64_t expected_run_id, Rank_Snapshot_State& state)
{
	std::ifstream file(path, std::ios::binary);
	if(!file.is_open())
		return false;

	Read_Binary_Value(file, state.run_id);
	Read_Binary_Value(file, state.snapshot_index);
	Read_Binary_Value(file, state.rank);
	Read_Binary_Value(file, state.done);
	Read_Binary_Value(file, state.trajectory_in_progress);
	Read_Binary_Value(file, state.local_captured);
	Read_Binary_Value(file, state.local_total);
	Read_Binary_Value(file, state.current_trajectory_id);
	Read_Binary_Value(file, state.rank_elapsed_wall_sec);
	Read_Binary_Value(file, state.current_trajectory_wall_sec);
	Read_Binary_Value(file, state.current_trajectory_physical_sec);
	Read_Binary_Value(file, state.current_trajectory_captured);
	Read_Binary_Array(file, state.current_trajectory_dt_hist);
	Read_Binary_Array(file, state.current_trajectory_v2dt_hist);
	Read_Binary_Array(file, state.captured_dt_hist);
	Read_Binary_Array(file, state.captured_v2dt_hist);
	Read_Binary_Array(file, state.captured_dt_sq_hist);
	Read_Binary_Array(file, state.captured_v2dt_sq_hist);
	Read_Binary_Array(file, state.not_captured_dt_hist);
	Read_Binary_Array(file, state.not_captured_v2dt_hist);
	Read_Binary_Array(file, state.not_captured_dt_sq_hist);
	Read_Binary_Array(file, state.not_captured_v2dt_sq_hist);
	Read_Binary_Value(file, state.total_wall_time_captured);
	Read_Binary_Value(file, state.total_wall_time_not_captured);
	Read_Binary_Value(file, state.total_rk45_steps_captured);
	Read_Binary_Value(file, state.total_rk45_steps_not_captured);
	Read_Binary_Array(file, state.wall_time_hist_captured);
	Read_Binary_Array(file, state.wall_time_hist_not_captured);
	Read_Binary_Value(file, state.wall_time_overflow_captured);
	Read_Binary_Value(file, state.wall_time_overflow_not_captured);
	Read_Binary_Array(file, state.step_count_hist_captured);
	Read_Binary_Array(file, state.step_count_hist_not_captured);
	Read_Binary_Value(file, state.step_count_overflow_captured);
	Read_Binary_Value(file, state.step_count_overflow_not_captured);
	uint64_t completed_evaporation_record_count = 0;
	Read_Binary_Value(file, completed_evaporation_record_count);
	state.completed_evaporation_records.clear();
	state.completed_evaporation_records.resize(static_cast<size_t>(completed_evaporation_record_count));
	for(size_t i = 0; i < state.completed_evaporation_records.size(); i++)
		Read_Evaporation_Record_Binary(file, state.completed_evaporation_records[i]);

	if(!file)
		return false;

	return state.run_id == expected_run_id;
}

bool Write_Snapshot_Text_Status(const std::string& snapshot_root, int snapshot_index, double interval_seconds, int caller_rank, const Snapshot_Load_Diagnostics& diagnostics, bool merged)
{
	return Write_Text_File_Atomically(Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index * 1000 + caller_rank, [&](std::ofstream& file)
	{
		Write_Snapshot_Diagnostics(file, snapshot_index, interval_seconds, caller_rank, diagnostics, merged);
	});
}

bool Write_Snapshot_Completed_Evaporation_Diagnostic_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds, double mass_gev, double sigma_cm2, const Snapshot_Report_State& report, int caller_rank, const Snapshot_Load_Diagnostics& diagnostics)
{
	return Write_Text_File_Atomically(Snapshot_Completed_Evaporation_Diagnostic_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index * 1000 + caller_rank + 500000, [&](std::ofstream& file)
	{
		Write_Snapshot_Diagnostics(file, snapshot_index, interval_seconds, caller_rank, diagnostics, true);
		file << "#\n";
		file << "# DIAGNOSTIC_ONLY = 1\n";
		file << "# NOT_FOR_SURVIVAL_ANALYSIS = 1\n";
		file << "# record_selection = survival_valid && event_observed && lifetime_unbinding >= 0\n";
		file << "# selection_bias = completed_by_snapshot_wall_time; unfinished long-lived trajectories are absent\n";
		file << "# final_survival_input = evaporation_summary.txt\n";
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);
		file << "# snapshot_completed_evaporation_record_count = " << report.completed_evaporation_records.size() << "\n";
		file << "# Use this file only to inspect already completed evaporation times during a running simulation.\n";
		Write_Evaporation_Record_List(file, report.completed_evaporation_records);
	});
}

void Accumulate_Snapshot_Report_State(Snapshot_Report_State& report, const Rank_Snapshot_State& state)
{
	report.total_trajectories += state.local_total;
	report.captured_particles += state.local_captured;
	report.snapshot_bincount_captured_samples += state.local_captured;
	report.snapshot_bincount_not_captured_samples += (state.local_total - state.local_captured);
	report.total_wall_time_captured += state.total_wall_time_captured;
	report.total_wall_time_not_captured += state.total_wall_time_not_captured;
	report.total_rk45_steps_captured += state.total_rk45_steps_captured;
	report.total_rk45_steps_not_captured += state.total_rk45_steps_not_captured;
	report.wall_time_overflow_captured += state.wall_time_overflow_captured;
	report.wall_time_overflow_not_captured += state.wall_time_overflow_not_captured;
	report.step_count_overflow_captured += state.step_count_overflow_captured;
	report.step_count_overflow_not_captured += state.step_count_overflow_not_captured;
	report.completed_evaporation_records.insert(report.completed_evaporation_records.end(),
	                                           state.completed_evaporation_records.begin(),
	                                           state.completed_evaporation_records.end());

	for(int bin = 0; bin < NUM_BINS; bin++)
	{
		report.captured_dt_hist[bin] += state.captured_dt_hist[bin];
		report.captured_v2dt_hist[bin] += state.captured_v2dt_hist[bin];
		report.captured_dt_sq_hist[bin] += state.captured_dt_sq_hist[bin];
		report.captured_v2dt_sq_hist[bin] += state.captured_v2dt_sq_hist[bin];
		report.not_captured_dt_hist[bin] += state.not_captured_dt_hist[bin];
		report.not_captured_v2dt_hist[bin] += state.not_captured_v2dt_hist[bin];
		report.not_captured_dt_sq_hist[bin] += state.not_captured_dt_sq_hist[bin];
		report.not_captured_v2dt_sq_hist[bin] += state.not_captured_v2dt_sq_hist[bin];
	}

	if(state.trajectory_in_progress && Has_Bincount_Contribution(state.current_trajectory_dt_hist, state.current_trajectory_v2dt_hist))
	{
		if(state.current_trajectory_captured)
		{
			report.snapshot_bincount_captured_samples++;
			for(int bin = 0; bin < NUM_BINS; bin++)
			{
				report.captured_dt_hist[bin] += state.current_trajectory_dt_hist[bin];
				report.captured_v2dt_hist[bin] += state.current_trajectory_v2dt_hist[bin];
				report.captured_dt_sq_hist[bin] += state.current_trajectory_dt_hist[bin] * state.current_trajectory_dt_hist[bin];
				report.captured_v2dt_sq_hist[bin] += state.current_trajectory_v2dt_hist[bin] * state.current_trajectory_v2dt_hist[bin];
			}
		}
		else
		{
			report.snapshot_bincount_not_captured_samples++;
			for(int bin = 0; bin < NUM_BINS; bin++)
			{
				report.not_captured_dt_hist[bin] += state.current_trajectory_dt_hist[bin];
				report.not_captured_v2dt_hist[bin] += state.current_trajectory_v2dt_hist[bin];
				report.not_captured_dt_sq_hist[bin] += state.current_trajectory_dt_hist[bin] * state.current_trajectory_dt_hist[bin];
				report.not_captured_v2dt_sq_hist[bin] += state.current_trajectory_v2dt_hist[bin] * state.current_trajectory_v2dt_hist[bin];
			}
		}
	}

	for(int bin = 0; bin < SNAPSHOT_WALL_TIME_BINS; bin++)
	{
		report.wall_time_hist_captured[bin] += state.wall_time_hist_captured[bin];
		report.wall_time_hist_not_captured[bin] += state.wall_time_hist_not_captured[bin];
	}

	for(int bin = 0; bin < SNAPSHOT_STEP_COUNT_BINS; bin++)
	{
		report.step_count_hist_captured[bin] += state.step_count_hist_captured[bin];
		report.step_count_hist_not_captured[bin] += state.step_count_hist_not_captured[bin];
	}

}

bool Load_Snapshot_Report_State(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, Snapshot_Report_State& report, Snapshot_Load_Diagnostics* diagnostics = NULL)
{
	report = Snapshot_Report_State();
	report.snapshot_time_label = Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds);
	report.snapshot_interval_seconds = interval_seconds;
	report.rank_states.reserve(mpi_processes);
	if(diagnostics != NULL)
	{
		diagnostics->ready_ranks = 0;
		diagnostics->per_rank_status.assign(mpi_processes, "pending");
	}

	bool all_ranks_ready = true;

	for(int rank = 0; rank < mpi_processes; rank++)
	{
		Rank_Snapshot_State state;
		const std::string checkpoint_path = Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		const bool checkpoint_exists = Path_Exists(checkpoint_path);
		if(Read_Rank_Snapshot_State(checkpoint_path, run_id, state))
		{
			Accumulate_Snapshot_Report_State(report, state);
			report.rank_states.push_back(state);
			if(diagnostics != NULL)
			{
				diagnostics->ready_ranks++;
				diagnostics->per_rank_status[rank] = "checkpoint, " + Format_Rank_Status(state);
			}
			continue;
		}

		const std::string final_path = Rank_Snapshot_Final_Path(rank_snapshot_dir, rank);
		const bool final_exists = Path_Exists(final_path);
		if(Read_Rank_Snapshot_State(final_path, run_id, state) && state.done && state.rank_elapsed_wall_sec <= report.snapshot_time_label)
		{
			Accumulate_Snapshot_Report_State(report, state);
			report.rank_states.push_back(state);
			if(diagnostics != NULL)
			{
				diagnostics->ready_ranks++;
				std::ostringstream status;
				status << "final(done_at=" << std::fixed << std::setprecision(3) << state.rank_elapsed_wall_sec << "s), " << Format_Rank_Status(state);
				diagnostics->per_rank_status[rank] = status.str();
			}
			continue;
		}

		all_ranks_ready = false;
		if(diagnostics != NULL)
		{
			std::ostringstream status;
			if(checkpoint_exists)
				status << "checkpoint_unreadable_or_stale";
			else
				status << "no_checkpoint";

			status << "; ";
			if(!final_exists)
				status << "no_final";
			else if(!Read_Rank_Snapshot_State(final_path, run_id, state))
				status << "final_unreadable_or_stale";
			else if(!state.done)
				status << "final_not_done";
			else if(state.rank_elapsed_wall_sec > report.snapshot_time_label)
				status << "final_after_target(done_at=" << std::fixed << std::setprecision(3) << state.rank_elapsed_wall_sec << "s)";
			else
				status << "final_unknown";

			diagnostics->per_rank_status[rank] = status.str();
		}
	}

	return all_ranks_ready;
}

bool Write_Snapshot_Report_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds, double mass_gev, double sigma_cm2, const Snapshot_Report_State& report, int caller_rank, const Snapshot_Load_Diagnostics& diagnostics)
{
	Remove_Stale_Snapshot_Evaporation_File(snapshot_root, snapshot_index, interval_seconds);

	if(!Write_Snapshot_Completed_Evaporation_Diagnostic_File(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report, caller_rank, diagnostics))
		return false;

	if(!Write_Text_File_Atomically(Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index, [&](std::ofstream& file)
	{
		Write_Snapshot_Diagnostics(file, snapshot_index, interval_seconds, caller_rank, diagnostics, true);
		file << "#\n";
		file << "# Cumulative snapshot report\n";
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);

		double snapshot_captured_samples = static_cast<double>(report.snapshot_bincount_captured_samples);
		double snapshot_not_captured_samples = static_cast<double>(report.snapshot_bincount_not_captured_samples);
		double mean_wall_cap = (report.captured_particles > 0) ? report.total_wall_time_captured / static_cast<double>(report.captured_particles) : 0.0;
		double not_captured_particles = static_cast<double>(report.total_trajectories - report.captured_particles);
		double captured_particles = static_cast<double>(report.captured_particles);
		double mean_wall_nc = (not_captured_particles > 0.0) ? report.total_wall_time_not_captured / not_captured_particles : 0.0;
		double mean_steps_cap = (report.captured_particles > 0) ? static_cast<double>(report.total_rk45_steps_captured) / captured_particles : 0.0;
		double mean_steps_nc = (not_captured_particles > 0.0) ? static_cast<double>(report.total_rk45_steps_not_captured) / not_captured_particles : 0.0;

		file << "# total_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << report.total_wall_time_captured << "\n";
		file << "# total_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << report.total_wall_time_not_captured << "\n";
		file << "# mean_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_cap << "\n";
		file << "# mean_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_nc << "\n";
		file << "# total_rk45_steps_captured = " << report.total_rk45_steps_captured << "\n";
		file << "# total_rk45_steps_not_captured = " << report.total_rk45_steps_not_captured << "\n";
		file << "# mean_rk45_steps_captured = " << std::scientific << std::setprecision(6) << mean_steps_cap << "\n";
		file << "# mean_rk45_steps_not_captured = " << std::scientific << std::setprecision(6) << mean_steps_nc << "\n";
		file << "# snapshot_bincount_captured_samples = " << report.snapshot_bincount_captured_samples << "\n";
		file << "# snapshot_bincount_not_captured_samples = " << report.snapshot_bincount_not_captured_samples << "\n";
		file << "# ranks_ready = " << diagnostics.ready_ranks << " / " << diagnostics.per_rank_status.size() << "\n";
		file << "#\n";
		file << "# [Bincount histogram]\n";
		file << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			double cap_err_dt = Snapshot_Bin_Error(report.captured_dt_hist[bin], report.captured_dt_sq_hist[bin], snapshot_captured_samples);
			double cap_err_v2dt = Snapshot_Bin_Error(report.captured_v2dt_hist[bin], report.captured_v2dt_sq_hist[bin], snapshot_captured_samples);
			double not_cap_err_dt = Snapshot_Bin_Error(report.not_captured_dt_hist[bin], report.not_captured_dt_sq_hist[bin], snapshot_not_captured_samples);
			double not_cap_err_v2dt = Snapshot_Bin_Error(report.not_captured_v2dt_hist[bin], report.not_captured_v2dt_sq_hist[bin], snapshot_not_captured_samples);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << report.captured_dt_hist[bin] << "\t" << report.captured_v2dt_hist[bin]
			     << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			     << "\t" << report.not_captured_dt_hist[bin] << "\t" << report.not_captured_v2dt_hist[bin]
			     << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
	}))
		return false;

	return true;
}

std::string Format_Rank_Status(const Rank_Snapshot_State& state)
{
	std::ostringstream stream;
	stream << std::fixed << std::setprecision(3);
	if(state.done)
	{
		stream << "done";
	}
	else
	{
		stream << "running";
	}

	stream << " (captured=" << state.local_captured
	       << ", completed=" << state.local_total
	       << ", completed_evap=" << state.completed_evaporation_records.size()
	       << ", rank_wall=" << state.rank_elapsed_wall_sec << "s";

	if(state.trajectory_in_progress)
	{
		stream << ", current_traj=" << state.current_trajectory_id
		       << ", traj_wall=" << state.current_trajectory_wall_sec << "s"
		       << ", traj_physical_s=" << Format_Physical_Time_Scientific(state.current_trajectory_physical_sec)
		       << ", partial_bucket=" << (state.current_trajectory_captured ? "captured" : "not_captured");
	}
	else if(!state.done)
	{
		stream << ", next_traj=" << (state.local_total + 1);
	}

	stream << ")";
	return stream.str();
}

bool Try_Write_Merged_Snapshot(const std::string& snapshot_root, const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, double mass_gev, double sigma_cm2, int caller_rank)
{
	const std::string snapshot_text_path = Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds);
	const std::string completed_evaporation_path = Snapshot_Completed_Evaporation_Diagnostic_File_Path(snapshot_root, snapshot_index, interval_seconds);
	if(Snapshot_Text_File_Is_Merged(snapshot_text_path) && Path_Exists(completed_evaporation_path))
	{
		Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
		return true;
	}

	Snapshot_Report_State report;
	Snapshot_Load_Diagnostics diagnostics;
	if(!Load_Snapshot_Report_State(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes, run_id, report, &diagnostics))
	{
		Write_Snapshot_Text_Status(snapshot_root, snapshot_index, interval_seconds, caller_rank, diagnostics, false);
		return false;
	}

	if(!Write_Snapshot_Report_File(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report, caller_rank, diagnostics))
		return false;

	Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
	return true;
}
}

Simulation_Data::Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min, unsigned int iso_rings)
: number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0),
  number_of_complete_evaporation_particles(0), number_of_censored_captured_particles(0),
  number_of_invalid_survival_captured_particles(0),
  average_number_of_scatterings(0.0), computing_time(0.0), early_stopped(false),
  total_wall_time_captured(0.0), total_wall_time_not_captured(0.0),
  wall_time_overflow_captured(0), wall_time_overflow_not_captured(0),
  total_rk45_steps_captured(0), total_rk45_steps_not_captured(0),
  step_count_overflow_captured(0), step_count_overflow_not_captured(0),
  mpi_rank(0), mpi_processes(1), isoreflection_rings(iso_rings), minimum_speed_threshold(u_min),
  number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)),
  data(iso_rings, std::vector<libphysica::DataPoint>())
{
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    target_captured_per_rank = (sample_size + mpi_processes - 1) / mpi_processes;  // ceil division
    if(max_trajectories == 0)
        max_trajectories_per_rank = std::numeric_limits<unsigned long int>::max();  // no limit
    else
        max_trajectories_per_rank = (max_trajectories + mpi_processes - 1) / mpi_processes;

    captured_dt_hist.fill(0.0);
    captured_v2dt_hist.fill(0.0);
    not_captured_dt_hist.fill(0.0);
    not_captured_v2dt_hist.fill(0.0);
	captured_dt_sq_hist.fill(0.0);
	captured_v2dt_sq_hist.fill(0.0);
	not_captured_dt_sq_hist.fill(0.0);
	not_captured_v2dt_sq_hist.fill(0.0);
	termination_reason_counts.fill(0);

	// Initialize computation time statistics
    total_wall_time_captured = 0.0;
    total_wall_time_not_captured = 0.0;
    wall_time_hist_captured.fill(0);
    wall_time_hist_not_captured.fill(0);
    wall_time_overflow_captured = 0;
    wall_time_overflow_not_captured = 0;
    total_rk45_steps_captured = 0;
    total_rk45_steps_not_captured = 0;
    step_count_hist_captured.fill(0);
    step_count_hist_not_captured.fill(0);
    step_count_overflow_captured = 0;
    step_count_overflow_not_captured = 0;
}

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, unsigned long int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius      = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps       = max_free_steps;
}

void Simulation_Data::Configure_Evaporation_Mode_Bincount(bool enabled, const std::vector<double>& boundaries_log10_s, const std::vector<std::string>& labels, bool include_truncated)
{
	evaporation_mode_bincount_enabled = enabled;
	evaporation_mode_boundaries_log10_s = boundaries_log10_s;
	evaporation_mode_labels = labels;
	evaporation_mode_include_truncated = include_truncated;
	evaporation_mode_bincounts.clear();
	if(evaporation_mode_bincount_enabled)
		evaporation_mode_bincounts.resize(evaporation_mode_boundaries_log10_s.size() + 1);
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg, unsigned int fixed_seed, bool capture_mode)
{
	if(capture_mode)
		snapshot_cfg.enabled = false;

	auto time_start = std::chrono::system_clock::now();
	unsigned long int local_captured = 0;
	unsigned long int local_total = 0;

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	simulator.max_trajectory_wall_time_sec = snapshot_cfg.max_trajectory_wall_time_sec;
	simulator.Enable_Capture_Mode(capture_mode);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	// Snapshot configuration
	const double snapshot_interval = (snapshot_cfg.interval_seconds > 0.0) ? snapshot_cfg.interval_seconds : 60.0;
	const double snapshot_mass_gev = In_Units(DM.mass, GeV);
	const double snapshot_sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
	std::string snapshot_root;
	std::string rank_snapshot_dir;
	uint64_t snapshot_run_id = 0;
	int next_snapshot_index = 1;
	if(snapshot_cfg.enabled && mpi_rank == 0 && snapshot_cfg.max_trajectory_wall_time_sec <= 0.0)
		std::cerr << "Warning in Generate_Data(): snapshot_enabled=true but max_trajectory_wall_time_sec="
		          << snapshot_cfg.max_trajectory_wall_time_sec
		          << " disables the per-trajectory wall-time guard. A single slow trajectory can keep snapshot files in waiting status. "
		          << "Consider setting max_trajectory_wall_time_sec to a finite value below snapshot_interval="
		          << snapshot_interval << " s." << std::endl;
	if(snapshot_cfg.enabled)
	{
		snapshot_root = g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		Ensure_Directory_Exists(g_top_level_dir + "results_" + std::to_string(log10(In_Units(DM.mass, GeV))) + "_" + std::to_string(log10(In_Units(DM.Sigma_Proton(), cm * cm))) + "/");
		if(mpi_rank == 0)
		{
			Ensure_Directory_Exists(snapshot_root);
			Clear_Directory_Contents(snapshot_root);
			Ensure_Directory_Exists(snapshot_root);
			Ensure_Directory_Exists(rank_snapshot_dir);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		Ensure_Directory_Exists(snapshot_root);
		Ensure_Directory_Exists(rank_snapshot_dir);
		if(mpi_rank == 0)
			snapshot_run_id = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		MPI_Bcast(&snapshot_run_id, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	}

	auto elapsed_since_start = [&]()
	{
		return 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
	};

	auto build_rank_snapshot_state = [&](bool done)
	{
		Rank_Snapshot_State state;
		state.run_id = snapshot_run_id;
		state.rank = mpi_rank;
		state.done = done ? 1 : 0;
		state.trajectory_in_progress = (!done && simulator.Trajectory_In_Progress()) ? 1 : 0;
		state.local_captured = static_cast<uint64_t>(local_captured);
		state.local_total = static_cast<uint64_t>(local_total);
		state.current_trajectory_id = state.trajectory_in_progress ? static_cast<uint64_t>(simulator.Current_Trajectory_ID()) : 0;
		state.rank_elapsed_wall_sec = done ? computing_time : elapsed_since_start();
		state.current_trajectory_wall_sec = state.trajectory_in_progress ? simulator.Current_Trajectory_Wall_Time_Seconds() : 0.0;
		state.current_trajectory_physical_sec = state.trajectory_in_progress ? simulator.Current_Trajectory_Physical_Time_Seconds() : 0.0;
		if(state.trajectory_in_progress)
		{
			const TrajectoryBincount& current_bincount = simulator.Current_Trajectory_Bincount();
			state.current_trajectory_captured = current_bincount.is_captured ? 1 : 0;
			state.current_trajectory_dt_hist = current_bincount.dt_hist;
			state.current_trajectory_v2dt_hist = current_bincount.v2dt_hist;
		}
		state.captured_dt_hist = captured_dt_hist;
		state.captured_v2dt_hist = captured_v2dt_hist;
		state.captured_dt_sq_hist = captured_dt_sq_hist;
		state.captured_v2dt_sq_hist = captured_v2dt_sq_hist;
		state.not_captured_dt_hist = not_captured_dt_hist;
		state.not_captured_v2dt_hist = not_captured_v2dt_hist;
		state.not_captured_dt_sq_hist = not_captured_dt_sq_hist;
		state.not_captured_v2dt_sq_hist = not_captured_v2dt_sq_hist;
		state.total_wall_time_captured = total_wall_time_captured;
		state.total_wall_time_not_captured = total_wall_time_not_captured;
		state.total_rk45_steps_captured = static_cast<uint64_t>(total_rk45_steps_captured);
		state.total_rk45_steps_not_captured = static_cast<uint64_t>(total_rk45_steps_not_captured);
		std::copy(wall_time_hist_captured.begin(), wall_time_hist_captured.end(), state.wall_time_hist_captured.begin());
		std::copy(wall_time_hist_not_captured.begin(), wall_time_hist_not_captured.end(), state.wall_time_hist_not_captured.begin());
		state.wall_time_overflow_captured = static_cast<uint64_t>(wall_time_overflow_captured);
		state.wall_time_overflow_not_captured = static_cast<uint64_t>(wall_time_overflow_not_captured);
		std::copy(step_count_hist_captured.begin(), step_count_hist_captured.end(), state.step_count_hist_captured.begin());
		std::copy(step_count_hist_not_captured.begin(), step_count_hist_not_captured.end(), state.step_count_hist_not_captured.begin());
		state.step_count_overflow_captured = static_cast<uint64_t>(step_count_overflow_captured);
		state.step_count_overflow_not_captured = static_cast<uint64_t>(step_count_overflow_not_captured);
		state.completed_evaporation_records.reserve(evaporation_records.size());
		for(const auto& rec : evaporation_records)
		{
			if(Is_Completed_Evaporation_Record(rec))
				state.completed_evaporation_records.push_back(rec);
		}
		return state;
	};

	auto try_write_ready_snapshots = [&](int max_snapshot_index)
	{
		bool all_snapshots_merged = true;
		for(int snapshot_index = 1; snapshot_index <= max_snapshot_index; snapshot_index++)
		{
			if(!Try_Write_Merged_Snapshot(snapshot_root, rank_snapshot_dir, snapshot_index, snapshot_interval, mpi_processes, snapshot_run_id, snapshot_mass_gev, snapshot_sigma_cm2, mpi_rank))
				all_snapshots_merged = false;
		}
		return all_snapshots_merged;
	};

	auto publish_checkpoint_snapshots = [&]()
	{
		if(!snapshot_cfg.enabled)
			return;

		double elapsed = elapsed_since_start();
		while(elapsed >= next_snapshot_index * snapshot_interval)
		{
			Rank_Snapshot_State state = build_rank_snapshot_state(false);
			state.snapshot_index = next_snapshot_index;
			Write_Rank_Snapshot_State(Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, mpi_rank, next_snapshot_index, snapshot_interval), state);
			try_write_ready_snapshots(next_snapshot_index);
			next_snapshot_index++;
		}
	};
	early_stopped = false;
	int last_progress_milestone = -1;
	if(snapshot_cfg.enabled)
		simulator.Set_Snapshot_Progress_Callback([&](const Trajectory_Simulator&)
		{
			publish_checkpoint_snapshots();
		});

	while(local_captured < target_captured_per_rank && local_total < max_trajectories_per_rank)
	{
		Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);

		auto traj_t0 = std::chrono::high_resolution_clock::now();
		Trajectory_Result trajectory = simulator.Simulate(IC, DM, mpi_rank);
		auto traj_t1 = std::chrono::high_resolution_clock::now();
		double traj_wall_sec = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(traj_t1 - traj_t0).count();

		local_total++;
		number_of_trajectories++;
		termination_reason_counts[TerminationReason_Index(trajectory.bincount.termination_reason)]++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);
		const bool completed_outward_escape = Completed_Outward_Escape(trajectory.bincount.termination_reason);

		// Helper lambdas for log-histogram binning
		auto wall_time_bin = [](double t) -> int {
			if(t <= 0.0) return -1;
			return static_cast<int>((log10(t) - WALL_TIME_LOG_MIN) * 10.0);
		};
		auto step_count_bin = [](unsigned long int s) -> int {
			if(s == 0) return -1;
			return static_cast<int>((log10(static_cast<double>(s)) - STEP_COUNT_LOG_MIN) * 10.0);
		};

		if(trajectory.bincount.is_captured)
		{
			number_of_captured_particles++;
			local_captured++;
			if(trajectory.bincount.event_observed)
				number_of_complete_evaporation_particles++;
			else if(!trajectory.bincount.survival_valid)
				number_of_invalid_survival_captured_particles++;
			else
				number_of_censored_captured_particles++;

			if(!capture_mode)
			{
				// Accumulate captured bincount
				for(int b = 0; b < NUM_BINS; b++)
				{
					captured_dt_hist[b]   += trajectory.bincount.dt_hist[b];
					captured_v2dt_hist[b] += trajectory.bincount.v2dt_hist[b];
					captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
					captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
				}

				// Record every captured trajectory, including right-censored ones.
				EvaporationRecord rec;
				if(Build_Evaporation_Record(trajectory.bincount, mpi_rank, number_of_trajectories, rec))
				{
					evaporation_records.push_back(rec);
					if(evaporation_mode_bincount_enabled
						&& rec.survival_valid
						&& rec.event_observed
						&& Has_Positive_Evaporation_Time(rec.lifetime_unbinding))
					{
						int mode = Classify_Evaporation_Mode(log10(rec.lifetime_unbinding), evaporation_mode_boundaries_log10_s);
						if(mode >= 0 && mode < static_cast<int>(evaporation_mode_bincounts.size()))
							Accumulate_Mode_Bincount(evaporation_mode_bincounts[mode], trajectory.bincount);
					}
				}
			}

			// Computation time & step count statistics (captured)
			total_wall_time_captured += traj_wall_sec;
			total_rk45_steps_captured += trajectory.total_rk45_steps;
			int wb = wall_time_bin(traj_wall_sec);
			if(wb >= 0 && wb < WALL_TIME_BINS) wall_time_hist_captured[wb]++;
			else if(wb >= WALL_TIME_BINS) wall_time_overflow_captured++;
			int sb = step_count_bin(trajectory.total_rk45_steps);
			if(sb >= 0 && sb < STEP_COUNT_BINS) step_count_hist_captured[sb]++;
			else if(sb >= STEP_COUNT_BINS) step_count_overflow_captured++;
		}
		else
		{
			if(!capture_mode && completed_outward_escape)
			{
				// Accumulate non-captured bincount
				for(int b = 0; b < NUM_BINS; b++)
				{
					not_captured_dt_hist[b]   += trajectory.bincount.dt_hist[b];
					not_captured_v2dt_hist[b] += trajectory.bincount.v2dt_hist[b];
					not_captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
					not_captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
				}
			}

			// Computation time & step count statistics (not captured)
			total_wall_time_not_captured += traj_wall_sec;
			total_rk45_steps_not_captured += trajectory.total_rk45_steps;
			int wb = wall_time_bin(traj_wall_sec);
			if(wb >= 0 && wb < WALL_TIME_BINS) wall_time_hist_not_captured[wb]++;
			else if(wb >= WALL_TIME_BINS) wall_time_overflow_not_captured++;
			int sb = step_count_bin(trajectory.total_rk45_steps);
			if(sb >= 0 && sb < STEP_COUNT_BINS) step_count_hist_not_captured[sb]++;
			else if(sb >= STEP_COUNT_BINS) step_count_overflow_not_captured++;

			if(!capture_mode && completed_outward_escape)
			{
				if(trajectory.Particle_Free())
					number_of_free_particles++;
				else if(trajectory.Particle_Reflected())
				{
					number_of_reflected_particles++;
					// Keep reflection data for compatibility (if needed)
					Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
					double v_final = trajectory.final_event.Speed();
					if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > KDE_boundary_correction_factor * minimum_speed_threshold)
					{
						unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
						data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
					}
				}
			}
		}

		// Progress bar (every 20%)
		if(mpi_rank == 0)
		{
			double progress = std::min(1.0, 1.0 * local_captured / target_captured_per_rank);
			int milestone = static_cast<int>(progress * 5);  // 0=0%,1=20%,...,5=100%
			if(milestone > last_progress_milestone)
			{
				last_progress_milestone = milestone;
				double time_elapsed = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
				libphysica::Print_Progress_Bar(progress, 0, 44, time_elapsed);
			}
		}

		publish_checkpoint_snapshots();
	}

	if(local_total >= max_trajectories_per_rank && local_captured < target_captured_per_rank)
		early_stopped = true;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

	publish_checkpoint_snapshots();

	// Final snapshot state: keep per-rank binary under snapshot/rank_snapshot and merge on interval boundaries.
	if(snapshot_cfg.enabled)
	{
		Rank_Snapshot_State final_state = build_rank_snapshot_state(true);
		final_state.snapshot_index = next_snapshot_index - 1;
		Write_Rank_Snapshot_State(Rank_Snapshot_Final_Path(rank_snapshot_dir, mpi_rank), final_state);

		int final_snapshot_index = std::max(1, static_cast<int>(std::ceil(computing_time / snapshot_interval)));
		if(try_write_ready_snapshots(final_snapshot_index))
			Cleanup_Final_Snapshot_States(rank_snapshot_dir, mpi_processes);
	}

	if(mpi_rank == 0)
	{
		libphysica::Print_Progress_Bar(1.0, 0, 44, computing_time);
		std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Perform_MPI_Reductions();
}

void Simulation_Data::Perform_MPI_Reductions()
{
	average_number_of_scatterings *= number_of_trajectories;
	MPI_Allreduce(MPI_IN_PLACE, &number_of_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_free_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_reflected_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_complete_evaporation_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_censored_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_invalid_survival_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, termination_reason_counts.data(), TRAJECTORY_TERMINATION_REASON_COUNT, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &average_number_of_scatterings, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings /= number_of_trajectories;

	// Reduce bincount histograms across all ranks
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for(auto& mode_bincount : evaporation_mode_bincounts)
	{
		MPI_Allreduce(MPI_IN_PLACE, &mode_bincount.count, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &mode_bincount.truncated_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, mode_bincount.dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, mode_bincount.v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, mode_bincount.dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, mode_bincount.v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	// Reduce early_stopped flag (any rank early stopped => global flag)
	int local_es = early_stopped ? 1 : 0;
	int global_es = 0;
	MPI_Allreduce(&local_es, &global_es, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	early_stopped = (global_es > 0);

	// Gather evaporation records from all ranks
	// First gather counts
	int local_evap_count = evaporation_records.size();
	std::vector<int> evap_counts(mpi_processes);
	MPI_Allgather(&local_evap_count, 1, MPI_INT, evap_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

	int total_evap = std::accumulate(evap_counts.begin(), evap_counts.end(), 0);

	// Pack local evaporation data. Use doubles for MPI transfer; the textual output documents the field meanings.
	constexpr int EVAPORATION_MPI_FIELDS = 23;
	std::vector<double> local_evap_data(local_evap_count * EVAPORATION_MPI_FIELDS);
	for(int i = 0; i < local_evap_count; i++)
	{
		local_evap_data[EVAPORATION_MPI_FIELDS*i]     = static_cast<double>(evaporation_records[i].rank);
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 1] = static_cast<double>(evaporation_records[i].trajectory_id);
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 2] = evaporation_records[i].t_evap;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 3] = evaporation_records[i].t_capture;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 4] = evaporation_records[i].t_final_unbinding_scatter;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 5] = evaporation_records[i].t_boundary_escape;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 6] = evaporation_records[i].t_termination;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 7] = evaporation_records[i].observed_lifetime;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 8] = evaporation_records[i].lifetime_unbinding;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 9] = evaporation_records[i].lifetime_boundary;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 10] = evaporation_records[i].r_first_negative_km;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 11] = evaporation_records[i].E_first_negative_eV;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 12] = evaporation_records[i].dE_first_negative_from_prev_eV;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 13] = evaporation_records[i].event_observed ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 14] = evaporation_records[i].boundary_escape_observed ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 15] = evaporation_records[i].survival_valid ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 16] = evaporation_records[i].numerically_invalid_escape ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 17] = evaporation_records[i].censored ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 18] = evaporation_records[i].truncated ? 1.0 : 0.0;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 19] = static_cast<double>(static_cast<int>(evaporation_records[i].termination_reason));
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 20] = evaporation_records[i].max_free_energy_drift_eV;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 21] = evaporation_records[i].max_free_energy_drift_rel;
		local_evap_data[EVAPORATION_MPI_FIELDS*i + 22] = static_cast<double>(evaporation_records[i].number_of_scatterings);
	}

	std::vector<int> recv_counts(mpi_processes), displacements(mpi_processes);
	for(int j = 0; j < mpi_processes; j++)
	{
		recv_counts[j] = evap_counts[j] * EVAPORATION_MPI_FIELDS;
		displacements[j] = (j == 0) ? 0 : displacements[j-1] + recv_counts[j-1];
	}

	std::vector<double> global_evap_data(total_evap * EVAPORATION_MPI_FIELDS);
	MPI_Allgatherv(local_evap_data.data(), local_evap_count * EVAPORATION_MPI_FIELDS, MPI_DOUBLE,
	               global_evap_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	// Unpack into evaporation_records
	evaporation_records.clear();
	evaporation_records.resize(total_evap);
	for(int i = 0; i < total_evap; i++)
	{
		evaporation_records[i].rank          = static_cast<int>(global_evap_data[EVAPORATION_MPI_FIELDS*i]);
		evaporation_records[i].trajectory_id = static_cast<unsigned long int>(global_evap_data[EVAPORATION_MPI_FIELDS*i + 1]);
		evaporation_records[i].t_evap        = global_evap_data[EVAPORATION_MPI_FIELDS*i + 2];
		evaporation_records[i].t_capture = global_evap_data[EVAPORATION_MPI_FIELDS*i + 3];
		evaporation_records[i].t_final_unbinding_scatter = global_evap_data[EVAPORATION_MPI_FIELDS*i + 4];
		evaporation_records[i].t_boundary_escape = global_evap_data[EVAPORATION_MPI_FIELDS*i + 5];
		evaporation_records[i].t_termination = global_evap_data[EVAPORATION_MPI_FIELDS*i + 6];
		evaporation_records[i].observed_lifetime = global_evap_data[EVAPORATION_MPI_FIELDS*i + 7];
		evaporation_records[i].lifetime_unbinding = global_evap_data[EVAPORATION_MPI_FIELDS*i + 8];
		evaporation_records[i].lifetime_boundary = global_evap_data[EVAPORATION_MPI_FIELDS*i + 9];
		evaporation_records[i].r_first_negative_km = global_evap_data[EVAPORATION_MPI_FIELDS*i + 10];
		evaporation_records[i].E_first_negative_eV = global_evap_data[EVAPORATION_MPI_FIELDS*i + 11];
		evaporation_records[i].dE_first_negative_from_prev_eV = global_evap_data[EVAPORATION_MPI_FIELDS*i + 12];
		evaporation_records[i].event_observed = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 13] > 0.5);
		evaporation_records[i].boundary_escape_observed = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 14] > 0.5);
		evaporation_records[i].survival_valid = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 15] > 0.5);
		evaporation_records[i].numerically_invalid_escape = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 16] > 0.5);
		evaporation_records[i].censored = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 17] > 0.5);
		evaporation_records[i].truncated = (global_evap_data[EVAPORATION_MPI_FIELDS*i + 18] > 0.5);
		evaporation_records[i].termination_reason = static_cast<TrajectoryTerminationReason>(TerminationReason_Index(static_cast<TrajectoryTerminationReason>(static_cast<int>(global_evap_data[EVAPORATION_MPI_FIELDS*i + 19]))));
		evaporation_records[i].max_free_energy_drift_eV = global_evap_data[EVAPORATION_MPI_FIELDS*i + 20];
		evaporation_records[i].max_free_energy_drift_rel = global_evap_data[EVAPORATION_MPI_FIELDS*i + 21];
		evaporation_records[i].number_of_scatterings = static_cast<unsigned long int>(global_evap_data[EVAPORATION_MPI_FIELDS*i + 22]);
	}

	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	// Reduce computation time & step count statistics
	MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_not_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_not_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, step_count_hist_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, step_count_hist_not_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
}

void Simulation_Data::Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM)
{
	if(mpi_rank != 0)
		return;

	// Create output directory
	struct stat info;
	if(stat(output_dir.c_str(), &info) != 0)
		mkdir(output_dir.c_str(), 0755);

	double mass_gev = In_Units(DM.mass, GeV);
	double sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);

	auto write_header = [&](std::ofstream& f) {
		Write_Report_Header(f, mass_gev, sigma_cm2, number_of_trajectories, number_of_captured_particles, early_stopped);
	};

	std::remove((output_dir + "/captured_bincount.txt").c_str());
	std::remove((output_dir + "/not_captured_bincount.txt").c_str());
	std::remove((output_dir + "/evaporation_mode_summary.txt").c_str());
	std::remove((output_dir + "/evaporation_mode_bincount.txt").c_str());

	// 1. Merged bincount
	{
		std::ofstream f(output_dir + "/bincount.txt");
		write_header(f);
		const unsigned long int physical_not_captured_particles = number_of_free_particles + number_of_reflected_particles;
		const unsigned long int raw_not_captured_particles = number_of_trajectories - number_of_captured_particles;
		const unsigned long int excluded_not_captured_particles =
		    (raw_not_captured_particles > physical_not_captured_particles) ? (raw_not_captured_particles - physical_not_captured_particles) : 0;
		f << "# not_captured_bincount_samples = " << physical_not_captured_particles << "\n";
		if(excluded_not_captured_particles > 0)
			f << "# excluded_incomplete_not_captured = " << excluded_not_captured_particles << "\n";
		f << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		double N_cap = static_cast<double>(number_of_captured_particles);
		double N_nc = static_cast<double>(physical_not_captured_particles);
		for(int b = 0; b < NUM_BINS; b++)
		{
			double cap_err_dt = Snapshot_Bin_Error(captured_dt_hist[b], captured_dt_sq_hist[b], N_cap);
			double cap_err_v2dt = Snapshot_Bin_Error(captured_v2dt_hist[b], captured_v2dt_sq_hist[b], N_cap);
			double not_cap_err_dt = Snapshot_Bin_Error(not_captured_dt_hist[b], not_captured_dt_sq_hist[b], N_nc);
			double not_cap_err_v2dt = Snapshot_Bin_Error(not_captured_v2dt_hist[b], not_captured_v2dt_sq_hist[b], N_nc);

			f << b << "\t" << std::scientific << std::setprecision(10)
			  << captured_dt_hist[b] << "\t" << captured_v2dt_hist[b]
			  << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			  << "\t" << not_captured_dt_hist[b] << "\t" << not_captured_v2dt_hist[b]
			  << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
		f.close();
	}

	// 2. Evaporation summary
	{
		std::ofstream f(output_dir + "/evaporation_summary.txt");
		write_header(f);
		Write_Evaporation_Record_List(f, evaporation_records);
		f.close();
	}

	if(evaporation_mode_bincount_enabled)
	{
		std::ofstream f(output_dir + "/evaporation_mode_summary.txt");
		write_header(f);
		f << "# evaporation_mode_include_truncated = " << (evaporation_mode_include_truncated ? 1 : 0) << "\n";
		f << "# evaporation_mode_include_truncated_status = deprecated_ignored\n";
		f << "# mode_assignment = event_observed && survival_valid && lifetime_unbinding > 0\n";
		f << "# mode_count = " << evaporation_mode_bincounts.size() << "\n";
		f << "# boundary_count = " << evaporation_mode_boundaries_log10_s.size() << "\n";
		f << "# boundaries_log10_s";
			for(double boundary : evaporation_mode_boundaries_log10_s)
				f << " " << std::scientific << std::setprecision(10) << boundary;
			f << "\n";
			f << "# mode_index  label  log10_t_lower  log10_t_upper  count  censored_count\n";
		for(size_t mode = 0; mode < evaporation_mode_bincounts.size(); mode++)
		{
			double lower = (mode == 0) ? -std::numeric_limits<double>::infinity() : evaporation_mode_boundaries_log10_s[mode - 1];
			double upper = (mode < evaporation_mode_boundaries_log10_s.size()) ? evaporation_mode_boundaries_log10_s[mode] : std::numeric_limits<double>::infinity();
			const std::string label = (mode < evaporation_mode_labels.size()) ? evaporation_mode_labels[mode] : ("mode_" + std::to_string(mode));
			f << mode << "\t" << label << "\t" << std::scientific << std::setprecision(10)
			  << lower << "\t" << upper << "\t" << evaporation_mode_bincounts[mode].count
			  << "\t" << evaporation_mode_bincounts[mode].truncated_count << "\n";
		}
		f.close();
	}

	if(evaporation_mode_bincount_enabled)
	{
		std::ofstream f(output_dir + "/evaporation_mode_bincount.txt");
		write_header(f);
		f << "# bin_index";
		for(size_t mode = 0; mode < evaporation_mode_bincounts.size(); mode++)
		{
			const std::string label = (mode < evaporation_mode_labels.size()) ? evaporation_mode_labels[mode] : ("mode_" + std::to_string(mode));
			f << "  " << label << "_dt[s]"
			  << "  " << label << "_v2dt[km2/s]"
			  << "  " << label << "_err_dt[s]"
			  << "  " << label << "_err_v2dt[km2/s]";
		}
		f << "\n";
		for(int b = 0; b < NUM_BINS; b++)
		{
			f << b;
			for(const auto& mode_bincount : evaporation_mode_bincounts)
			{
				double N_mode = static_cast<double>(mode_bincount.count);
				double err_dt = Snapshot_Bin_Error(mode_bincount.dt_hist[b], mode_bincount.dt_sq_hist[b], N_mode);
				double err_v2dt = Snapshot_Bin_Error(mode_bincount.v2dt_hist[b], mode_bincount.v2dt_sq_hist[b], N_mode);
				f << "\t" << std::scientific << std::setprecision(10)
				  << mode_bincount.dt_hist[b] << "\t" << mode_bincount.v2dt_hist[b]
				  << "\t" << err_dt << "\t" << err_v2dt;
			}
			f << "\n";
		}
		f.close();
	}

	// 3. Computation time & RK45 step count summary
	{
		std::ofstream f(output_dir + "/computation_time_summary.txt");
		write_header(f);

		double mean_wall_cap = (number_of_captured_particles > 0) ? total_wall_time_captured / number_of_captured_particles : 0.0;
		double mean_wall_nc  = (number_of_trajectories > number_of_captured_particles) ? total_wall_time_not_captured / (number_of_trajectories - number_of_captured_particles) : 0.0;
		double mean_steps_cap = (number_of_captured_particles > 0) ? static_cast<double>(total_rk45_steps_captured) / number_of_captured_particles : 0.0;
		double mean_steps_nc  = (number_of_trajectories > number_of_captured_particles) ? static_cast<double>(total_rk45_steps_not_captured) / (number_of_trajectories - number_of_captured_particles) : 0.0;

		f << "# ever_captured_count = " << number_of_captured_particles << "\n";
		f << "# complete_evaporation_count = " << number_of_complete_evaporation_particles << "\n";
		f << "# censored_captured_count = " << number_of_censored_captured_particles << "\n";
		f << "# invalid_survival_captured_count = " << number_of_invalid_survival_captured_particles << "\n";
		f << "# evaporation_record_count = " << evaporation_records.size() << "\n";
		f << "#\n";
		f << "# total_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << total_wall_time_captured << "\n";
		f << "# total_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << total_wall_time_not_captured << "\n";
		f << "# mean_wall_time_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_cap << "\n";
		f << "# mean_wall_time_not_captured_sec = " << std::scientific << std::setprecision(6) << mean_wall_nc << "\n";
		f << "# total_rk45_steps_captured = " << total_rk45_steps_captured << "\n";
		f << "# total_rk45_steps_not_captured = " << total_rk45_steps_not_captured << "\n";
		f << "# mean_rk45_steps_captured = " << std::scientific << std::setprecision(6) << mean_steps_cap << "\n";
		f << "# mean_rk45_steps_not_captured = " << std::scientific << std::setprecision(6) << mean_steps_nc << "\n";
		f << "#\n";
		f << "# incomplete_or_guarded_terminations = " << Count_Incomplete_Physical_Terminations(termination_reason_counts) << "\n";
		f << "# [Termination reason counts]\n";
		f << "# reason  count\n";
		for(int reason_index = 0; reason_index < TRAJECTORY_TERMINATION_REASON_COUNT; reason_index++)
		{
			TrajectoryTerminationReason reason = static_cast<TrajectoryTerminationReason>(reason_index);
			f << TerminationReason_Name(reason) << "\t" << termination_reason_counts[reason_index] << "\n";
		}
		f << "#\n";

		// Wall-clock time histogram
		int wall_first_bin = 0;
		int wall_last_bin = -1;
		bool has_wall_bins = Find_Nonzero_Bin_Range(wall_time_hist_captured, wall_time_hist_not_captured, wall_first_bin, wall_last_bin);
		int wall_output_first_bin = has_wall_bins ? wall_first_bin - 1 : 0;
		f << "# [Wall-clock time histogram]\n";
		if(has_wall_bins)
		{
			double observed_log10_t_min = WALL_TIME_LOG_MIN + 0.1 * wall_output_first_bin;
			double observed_log10_t_max = WALL_TIME_LOG_MIN + 0.1 * (wall_last_bin + 1);
			f << "# occupied_log10_t_min = " << std::fixed << std::setprecision(1) << observed_log10_t_min << "\n";
			f << "# occupied_log10_t_max = " << std::fixed << std::setprecision(1) << observed_log10_t_max << "\n";
		}
		else
		{
			f << "# occupied_log10_t_min = nan\n";
			f << "# occupied_log10_t_max = nan\n";
		}
		f << "# bin_index  log10_t_lower  log10_t_upper  count_captured  count_not_captured\n";
		for(int b = wall_output_first_bin, compact_index = 0; has_wall_bins && b <= wall_last_bin; b++, compact_index++)
		{
			double lo = WALL_TIME_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << compact_index << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
			  << "\t" << ((b >= 0) ? wall_time_hist_captured[b] : 0)
			  << "\t" << ((b >= 0) ? wall_time_hist_not_captured[b] : 0) << "\n";
		}
		f << "# overflow_captured = " << wall_time_overflow_captured << "\n";
		f << "# overflow_not_captured = " << wall_time_overflow_not_captured << "\n";
		f << "#\n";

		// RK45 step count histogram
		int step_first_bin = 0;
		int step_last_bin = -1;
		bool has_step_bins = Find_Nonzero_Bin_Range(step_count_hist_captured, step_count_hist_not_captured, step_first_bin, step_last_bin);
		f << "# [RK45 step count histogram]\n";
		if(has_step_bins)
		{
			double observed_log10_steps_min = STEP_COUNT_LOG_MIN + 0.1 * step_first_bin;
			double observed_log10_steps_max = STEP_COUNT_LOG_MIN + 0.1 * (step_last_bin + 1);
			f << "# occupied_log10_steps_min = " << std::fixed << std::setprecision(1) << observed_log10_steps_min << "\n";
			f << "# occupied_log10_steps_max = " << std::fixed << std::setprecision(1) << observed_log10_steps_max << "\n";
		}
		else
		{
			f << "# occupied_log10_steps_min = nan\n";
			f << "# occupied_log10_steps_max = nan\n";
		}
		f << "# bin_index  log10_steps_lower  log10_steps_upper  count_captured  count_not_captured\n";
		for(int b = step_first_bin, compact_index = 0; has_step_bins && b <= step_last_bin; b++, compact_index++)
		{
			double lo = STEP_COUNT_LOG_MIN + b * 0.1;
			double hi = lo + 0.1;
			f << compact_index << "\t" << std::fixed << std::setprecision(1) << lo << "\t" << hi
			  << "\t" << step_count_hist_captured[b] << "\t" << step_count_hist_not_captured[b] << "\n";
		}
		f << "# overflow_captured = " << step_count_overflow_captured << "\n";
		f << "# overflow_not_captured = " << step_count_overflow_not_captured << "\n";
		f.close();
	}
}

double Simulation_Data::Free_Ratio() const
{
	return 1.0 * number_of_free_particles / number_of_trajectories;
}
double Simulation_Data::Capture_Ratio() const
{
	return 1.0 * number_of_captured_particles / number_of_trajectories;
}
double Simulation_Data::Reflection_Ratio(int isoreflection_ring) const
{
	if(isoreflection_ring < 0)
		return 1.0 * number_of_reflected_particles / number_of_trajectories;
	else
		return 1.0 * data[isoreflection_ring].size() / number_of_trajectories;
}

double Simulation_Data::Minimum_Speed() const
{
	return KDE_boundary_correction_factor * minimum_speed_threshold;
}

double Simulation_Data::Lowest_Speed(unsigned int iso_ring) const
{
	return (*std::min_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

double Simulation_Data::Highest_Speed(unsigned int iso_ring) const
{
	return (*std::max_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

void Simulation_Data::Print_Capture_Mode_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		double N = static_cast<double>(number_of_trajectories);
		double p = (N > 0.0) ? static_cast<double>(number_of_captured_particles) / N : 0.0;
		double z = 1.96;
		double ci_lower = p;
		double ci_upper = p;
		if(N > 0.0)
		{
			double denom = 1.0 + z * z / N;
			double center = p + z * z / (2.0 * N);
			double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
			ci_lower = (center - spread) / denom;
			ci_upper = (center + spread) / denom;
			if(ci_lower < 0.0) ci_lower = 0.0;
			if(ci_upper > 1.0) ci_upper = 1.0;
		}

		std::cout << SEPARATOR
		          << "CAPTURE MODE summary" << std::endl
		          << std::endl
		          << "Termination condition:\t\tpost-scatter E < 0" << std::endl
		          << "File output:\t\t\tdisabled" << std::endl
		          << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
		          << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
		          << "Capture rate:\t\t\t" << std::fixed << std::setprecision(8) << p << std::endl
		          << "Capture rate -error (95%):\t" << std::fixed << std::setprecision(8) << (p - ci_lower) << std::endl
		          << "Capture rate +error (95%):\t" << std::fixed << std::setprecision(8) << (ci_upper - p) << std::endl
		          << "Capture rate 95% CI:\t\t[" << std::fixed << std::setprecision(8) << ci_lower << ", " << ci_upper << "]" << std::endl;

		if(early_stopped)
			std::cout << "*** EARLY STOP: max_trajectories reached ***" << std::endl;

		std::cout << "Simulation time:\t\t" << libphysica::Time_Display(computing_time) << std::endl
		          << SEPARATOR << std::endl;
	}
}

void Simulation_Data::Print_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Simulation data summary" << std::endl
				  << std::endl
				  << "Results:" << std::endl
				  << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles [%]:\t\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl
				  << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
				  << "Complete evaporation count:\t" << number_of_complete_evaporation_particles << std::endl
				  << "Censored captured count:\t" << number_of_censored_captured_particles << std::endl
				  << "Invalid survival count:\t\t" << number_of_invalid_survival_captured_particles << std::endl;
		std::cout << "Termination reasons:" << std::endl;
		for(int reason_index = 0; reason_index < TRAJECTORY_TERMINATION_REASON_COUNT; reason_index++)
		{
			if(termination_reason_counts[reason_index] == 0)
				continue;
			std::cout << "  " << TerminationReason_Name(static_cast<TrajectoryTerminationReason>(reason_index))
			          << ":\t" << termination_reason_counts[reason_index] << std::endl;
		}
		unsigned long int incomplete_terminations = Count_Incomplete_Physical_Terminations(termination_reason_counts);
		if(incomplete_terminations > 0)
		{
			std::cout << "WARNING: incomplete/guarded terminations:\t" << incomplete_terminations << std::endl
			          << "         Incomplete non-captured trajectories are excluded from free/reflected and not_captured bincount statistics." << std::endl;
		}

		// Capture rate error (Wilson interval)
		{
			double N = static_cast<double>(number_of_trajectories);
			double p = Capture_Ratio();
			double z = 1.96;
			double sigma_p = (N > 0) ? sqrt(p * (1.0 - p) / N) : 0.0;
			std::cout << "Capture rate error (1σ):\t" << std::fixed << std::setprecision(6) << sigma_p << std::endl;
			if(N > 0)
			{
				double denom = 1.0 + z * z / N;
				double center = p + z * z / (2.0 * N);
				double spread = z * sqrt(p * (1.0 - p) / N + z * z / (4.0 * N * N));
				double ci_lower = (center - spread) / denom;
				double ci_upper = (center + spread) / denom;
				if(ci_lower < 0.0) ci_lower = 0.0;
				if(ci_upper > 1.0) ci_upper = 1.0;
				std::cout << "Capture rate 95% CI:\t\t[" << std::fixed << std::setprecision(6) << ci_lower << ", " << ci_upper << "]" << std::endl;
			}
		}

		if(early_stopped)
			std::cout << "*** EARLY STOP: max_trajectories reached ***" << std::endl;

		// Median for observed unbinding events only; censored records belong in survival analysis.
		std::vector<double> observed_unbinding_lifetimes;
		for(const auto& rec : evaporation_records)
		{
			if(rec.survival_valid && rec.event_observed && Has_Positive_Evaporation_Time(rec.lifetime_unbinding))
				observed_unbinding_lifetimes.push_back(rec.lifetime_unbinding);
		}
		if(!observed_unbinding_lifetimes.empty())
		{
			std::sort(observed_unbinding_lifetimes.begin(), observed_unbinding_lifetimes.end());
			double median;
			size_t n = observed_unbinding_lifetimes.size();
			if(n % 2 == 0)
				median = 0.5 * (observed_unbinding_lifetimes[n/2 - 1] + observed_unbinding_lifetimes[n/2]);
			else
				median = observed_unbinding_lifetimes[n/2];
			std::cout << "Observed unbinding lifetime median [s]:\t" << std::scientific << std::setprecision(4) << median << " (" << observed_unbinding_lifetimes.size() << " events)" << std::endl;
		}
		else
		{
			std::cout << "Observed unbinding lifetime median:\tN/A (no positive observed unbinding events)" << std::endl;
		}

		std::cout << std::endl
				  << "Trajectory rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_trajectories / computing_time) << std::endl
				  << "Capture rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_captured_particles / computing_time) << std::endl
				  << "Simulation time:\t\t" << libphysica::Time_Display(computing_time) << std::endl;

		// Per-trajectory computation time & RK45 step count summary
		std::cout << std::endl
				  << "Captured wall-clock time:\t" << libphysica::Time_Display(total_wall_time_captured) << std::endl
				  << "Not-captured wall-clock time:\t" << libphysica::Time_Display(total_wall_time_not_captured) << std::endl;
		if(number_of_captured_particles > 0)
			std::cout << "Mean captured traj time:\t\t" << std::scientific << std::setprecision(4) << total_wall_time_captured / number_of_captured_particles << " s" << std::endl;
		if(number_of_trajectories > number_of_captured_particles)
			std::cout << "Mean not-captured traj time:\t" << std::scientific << std::setprecision(4) << total_wall_time_not_captured / (number_of_trajectories - number_of_captured_particles) << " s" << std::endl;
		std::cout << "Total RK45 steps (captured):\t" << total_rk45_steps_captured << std::endl
				  << "Total RK45 steps (not captured):\t" << total_rk45_steps_not_captured << std::endl;
		if(number_of_captured_particles > 0)
			std::cout << "Mean RK45 steps/traj (captured):\t" << std::scientific << std::setprecision(4) << static_cast<double>(total_rk45_steps_captured) / number_of_captured_particles << std::endl;
		if(number_of_trajectories > number_of_captured_particles)
			std::cout << "Mean RK45 steps/traj (not cap):\t" << std::scientific << std::setprecision(4) << static_cast<double>(total_rk45_steps_not_captured) / (number_of_trajectories - number_of_captured_particles) << std::endl;

		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN

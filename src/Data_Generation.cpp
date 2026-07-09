#include "Data_Generation.hpp"

#include <dirent.h>
#include <errno.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <unordered_set>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
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
bool Has_Positive_Evaporation_Time(double t_evap)
{
	return std::isfinite(t_evap) && t_evap > 0.0;
}

void MPI_Trace_Point(int mpi_rank, const std::string& label)
{
	if(std::getenv("DAMASCUS_MPI_TRACE") == nullptr)
		return;
	std::cerr << "[mpi-trace rank " << mpi_rank << "] " << label << std::endl;
}

constexpr unsigned long int CAPTURE_MODE_MPI_SYNC_INTERVAL = 32UL;
constexpr unsigned long int NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK = 128UL;
constexpr double NUMERICAL_FAILURE_WARNING_FRACTION = 1.0e-4;
constexpr double NUMERICAL_FAILURE_ABORT_FRACTION = 1.0e-2;

unsigned long int Normal_Mode_MPI_Sync_Interval(double sigma_cm2)
{
	if(!std::isfinite(sigma_cm2) || sigma_cm2 <= 0.0)
		return NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK;
	if(sigma_cm2 >= 1.0e-35)
		return 64UL;
	if(sigma_cm2 >= 1.0e-36)
		return 128UL;
	if(sigma_cm2 >= 1.0e-37)
		return 1024UL;
	if(sigma_cm2 >= 1.0e-38)
		return 8192UL;
	if(sigma_cm2 >= 1.0e-39)
		return 65536UL;
	return 1048576UL;
}

bool Is_Completed_Evaporation_Record(const EvaporationRecord& rec)
{
	return rec.survival_valid
	    && rec.event_observed
	    && std::isfinite(rec.lifetime_unbinding)
	    && rec.lifetime_unbinding >= 0.0;
}

bool Is_Completed_Evaporation_Event(const CompactEvaporationEvent& event)
{
	return std::isfinite(event.lifetime_unbinding) && event.lifetime_unbinding >= 0.0;
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

bool Build_Evaporation_Record(const TrajectoryBincount& bincount, int mpi_rank, unsigned long int trajectory_id, double completion_wall_time_sec, EvaporationRecord& rec)
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
	rec.completion_wall_time_sec = completion_wall_time_sec;
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
	// Normal evaporation output only accepts completed unbinding events. A
	// captured trajectory that has not completed is invalid for this estimator,
	// not a right-censored observation.
	rec.censored = false;
	rec.truncated = false;
	rec.termination_reason = bincount.termination_reason;
	rec.max_free_energy_drift_eV = bincount.max_free_energy_drift_eV;
	rec.max_free_energy_drift_rel = bincount.max_free_energy_drift_rel;
	rec.number_of_scatterings = bincount.number_of_scatterings;
	return true;
}

bool Make_Compact_Evaporation_Event(const EvaporationRecord& rec, CompactEvaporationEvent& event)
{
	if(!Is_Completed_Evaporation_Record(rec))
		return false;

	event.rank = rec.rank;
	event.trajectory_id = rec.trajectory_id;
	event.completion_wall_time_sec = rec.completion_wall_time_sec;
	event.lifetime_unbinding = rec.lifetime_unbinding;
	return true;
}

struct SnapshotEvaporationEntry
{
	uint64_t trajectory_id = 0;
	double completion_wall_time_sec = 0.0;
	double lifetime_unbinding_sec = -1.0;
};

struct RankedSnapshotEvaporationEntry
{
	int rank = -1;
	SnapshotEvaporationEntry entry;
};

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
	std::vector<SnapshotEvaporationEntry> new_evaporation_events;
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
	std::vector<RankedSnapshotEvaporationEntry> new_evaporation_events;
};

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

void Write_Snapshot_Evaporation_Entry_Binary(std::ofstream& file, const SnapshotEvaporationEntry& entry)
{
	Write_Binary_Value(file, entry.trajectory_id);
	Write_Binary_Value(file, entry.completion_wall_time_sec);
	Write_Binary_Value(file, entry.lifetime_unbinding_sec);
}

void Read_Snapshot_Evaporation_Entry_Binary(std::ifstream& file, SnapshotEvaporationEntry& entry)
{
	Read_Binary_Value(file, entry.trajectory_id);
	Read_Binary_Value(file, entry.completion_wall_time_sec);
	Read_Binary_Value(file, entry.lifetime_unbinding_sec);
}

SnapshotEvaporationEntry Make_Snapshot_Evaporation_Entry(const CompactEvaporationEvent& event)
{
	SnapshotEvaporationEntry entry;
	entry.trajectory_id = static_cast<uint64_t>(event.trajectory_id);
	entry.completion_wall_time_sec = event.completion_wall_time_sec;
	entry.lifetime_unbinding_sec = event.lifetime_unbinding;
	return entry;
}

CompactEvaporationEvent Make_Log_Event(int rank, const SnapshotEvaporationEntry& entry)
{
	CompactEvaporationEvent event;
	event.rank = rank;
	event.trajectory_id = entry.trajectory_id;
	event.completion_wall_time_sec = entry.completion_wall_time_sec;
	event.lifetime_unbinding = entry.lifetime_unbinding_sec;
	return event;
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

std::string Join_Path(const std::string& directory, const std::string& name)
{
	if(directory.empty() || directory.back() == '/')
		return directory + name;
	return directory + "/" + name;
}

std::string Evaporation_Log_Path_From_Output_Dir(const std::string& output_dir)
{
	return Join_Path(output_dir, "evaporation_times.txt");
}

long long Snapshot_Time_Label_Seconds(int snapshot_index, double interval_seconds);

void Write_Evaporation_Log_File_Header(std::ofstream& file, double mass_gev, double sigma_cm2)
{
	file << "# DaMaSCUS-SUN evaporation times\n";
	file << "# format_version = 3\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# sorted_by = lifetime_unbinding_sec rank trajectory_id\n";
	file << "# rank trajectory_id lifetime_unbinding_sec\n";
}

void Write_Evaporation_Log_Event(std::ostream& file, const CompactEvaporationEvent& event)
{
	file << event.rank << "\t" << event.trajectory_id << "\t" << std::scientific << std::setprecision(10)
	     << event.lifetime_unbinding << "\n";
}

bool Evaporation_Event_Order(const CompactEvaporationEvent& lhs, const CompactEvaporationEvent& rhs)
{
	if(lhs.lifetime_unbinding != rhs.lifetime_unbinding)
		return lhs.lifetime_unbinding < rhs.lifetime_unbinding;
	if(lhs.rank != rhs.rank)
		return lhs.rank < rhs.rank;
	return lhs.trajectory_id < rhs.trajectory_id;
}

void Write_Evaporation_Log_Events(std::ostream& file, const std::vector<CompactEvaporationEvent>& events)
{
	std::vector<CompactEvaporationEvent> sorted_events = events;
	std::sort(sorted_events.begin(), sorted_events.end(), Evaporation_Event_Order);
	for(const auto& event : sorted_events)
	{
		if(Is_Completed_Evaporation_Event(event))
			Write_Evaporation_Log_Event(file, event);
	}
}

bool Write_Final_Evaporation_Time_File(const std::string& path, double mass_gev, double sigma_cm2, const std::vector<CompactEvaporationEvent>& events)
{
	const std::string tmp_path = path + ".final.tmp." + std::to_string(getpid());
	std::ofstream file(tmp_path, std::ios::out | std::ios::trunc);
	if(!file.is_open())
		return false;
	Write_Evaporation_Log_File_Header(file, mass_gev, sigma_cm2);
	Write_Evaporation_Log_Events(file, events);
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

bool Ensure_Directory_Exists(const std::string& directory)
{
	if(directory.empty())
		return false;

	std::string normalized = directory;
	while(normalized.size() > 1 && normalized.back() == '/')
		normalized.pop_back();

	struct stat info;
	if(stat(normalized.c_str(), &info) == 0)
		return S_ISDIR(info.st_mode);

	const std::size_t separator = normalized.find_last_of('/');
	if(separator != std::string::npos)
	{
		std::string parent = normalized.substr(0, separator);
		if(parent.empty())
			parent = "/";
		if(parent != normalized && !Ensure_Directory_Exists(parent))
			return false;
	}

	if(mkdir(normalized.c_str(), 0755) != 0 && errno != EEXIST)
		return false;
	if(stat(normalized.c_str(), &info) != 0)
		return false;
	return S_ISDIR(info.st_mode);
}

long long Snapshot_Time_Label_Seconds(int snapshot_index, double interval_seconds)
{
	return static_cast<long long>(std::llround(snapshot_index * interval_seconds));
}

std::string Snapshot_Text_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s.txt";
}

std::string Snapshot_Evaporation_Time_File_Path(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds)) + "s_evaporation_times.txt";
}

// Retired block/manifest implementation. Snapshots are now standalone atomic
// bincount and evaporation-time files, so this path is deliberately inactive.
#if 0
std::string Snapshot_Evaporation_Block_Dir(const std::string& snapshot_root)
{
	return Join_Path(snapshot_root, "evaporation_blocks") + "/";
}

std::string Snapshot_Evaporation_Block_Path(const std::string& snapshot_root, int snapshot_index)
{
	std::ostringstream stream;
	stream << "block_" << std::setw(6) << std::setfill('0') << snapshot_index << ".txt";
	return Snapshot_Evaporation_Block_Dir(snapshot_root) + stream.str();
}

std::string Snapshot_Evaporation_Block_Dir_From_Output_Dir(const std::string& output_dir)
{
	return Snapshot_Evaporation_Block_Dir(Join_Path(output_dir, "snapshot"));
}

std::string Snapshot_Run_Manifest_Path(const std::string& snapshot_root)
{
	return Join_Path(snapshot_root, "run_manifest.txt");
}

std::string Evaporation_Partial_Log_Path_From_Output_Dir(const std::string& output_dir)
{
	return Join_Path(output_dir, "evaporation_times.partial.txt");
}

bool Is_Evaporation_Block_File_Name(const std::string& name)
{
	return name.size() > 10
	    && name.compare(0, 6, "block_") == 0
	    && name.compare(name.size() - 4, 4, ".txt") == 0;
}

struct EvaporationBlockFile
{
	int snapshot_index = -1;
	std::string path;
};

bool List_Evaporation_Block_Files(const std::string& snapshot_root, std::vector<EvaporationBlockFile>& block_files)
{
	block_files.clear();
	const std::string block_dir = Snapshot_Evaporation_Block_Dir(snapshot_root);
	DIR* dir = opendir(block_dir.c_str());
	if(dir == NULL)
		return errno == ENOENT;

	bool success = true;
	struct dirent* entry = NULL;
	while((entry = readdir(dir)) != NULL)
	{
		const std::string name = entry->d_name;
		if(!Is_Evaporation_Block_File_Name(name))
			continue;

		const std::string path = Join_Path(block_dir, name);
		EvaporationBlockMetadata metadata;
		if(!Read_Evaporation_Block_Metadata(path, metadata))
		{
			success = false;
			continue;
		}

		EvaporationBlockFile block_file;
		block_file.snapshot_index = metadata.snapshot_index;
		block_file.path = path;
		block_files.push_back(block_file);
	}
	closedir(dir);

	std::sort(block_files.begin(), block_files.end(), [](const EvaporationBlockFile& lhs, const EvaporationBlockFile& rhs)
	{
		return lhs.snapshot_index < rhs.snapshot_index;
	});

	return success;
}

bool Snapshot_Evaporation_Block_Dir_Has_Blocks(const std::string& snapshot_root)
{
	const std::string block_dir = Snapshot_Evaporation_Block_Dir(snapshot_root);
	DIR* dir = opendir(block_dir.c_str());
	if(dir == NULL)
		return false;

	bool has_blocks = false;
	struct dirent* entry = NULL;
	while((entry = readdir(dir)) != NULL)
	{
		if(Is_Evaporation_Block_File_Name(entry->d_name))
		{
			has_blocks = true;
			break;
		}
	}
	closedir(dir);
	return has_blocks;
}

bool Write_Snapshot_Run_Manifest(const std::string& snapshot_root, uint64_t run_id, double snapshot_interval_sec, double mass_gev, double sigma_cm2)
{
	const std::string path = Snapshot_Run_Manifest_Path(snapshot_root);
	const std::string tmp_path = path + ".tmp." + std::to_string(getpid());
	std::ofstream file(tmp_path, std::ios::out | std::ios::trunc);
	if(!file.is_open())
		return false;

	file << "# format_version = " << SNAPSHOT_RUN_MANIFEST_FORMAT_VERSION << "\n";
	file << "# run_id = " << run_id << "\n";
	file << "# snapshot_interval_sec = " << std::scientific << std::setprecision(17) << snapshot_interval_sec << "\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(17) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(17) << sigma_cm2 << "\n";
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

bool Read_Snapshot_Run_Manifest(const std::string& snapshot_root, SnapshotRunManifest& manifest)
{
	std::ifstream file(Snapshot_Run_Manifest_Path(snapshot_root));
	if(!file.is_open())
		return false;

	int format_version = 0;
	bool has_format_version = false;
	bool has_run_id = false;
	bool has_snapshot_interval = false;
	bool has_mass = false;
	bool has_sigma = false;

	std::string line;
	while(std::getline(file, line))
	{
		std::string key;
		std::string value;
		if(!Parse_Metadata_Line(line, key, value))
			continue;

		if(key == "format_version")
			has_format_version = Parse_Metadata_Int(value, format_version);
		else if(key == "run_id")
			has_run_id = Parse_Metadata_UInt64(value, manifest.run_id);
		else if(key == "snapshot_interval_sec")
			has_snapshot_interval = Parse_Metadata_Double(value, manifest.snapshot_interval_sec);
		else if(key == "DM_mass_GeV")
			has_mass = Parse_Metadata_Double(value, manifest.mass_gev);
		else if(key == "DM_sigma_cm2")
			has_sigma = Parse_Metadata_Double(value, manifest.sigma_cm2);
	}

	return has_format_version
	    && format_version == SNAPSHOT_RUN_MANIFEST_FORMAT_VERSION
	    && has_run_id
	    && has_snapshot_interval
	    && has_mass
	    && has_sigma;
}

void Remove_Stale_Snapshot_Evaporation_Block(const std::string& snapshot_root, int snapshot_index)
{
	std::string path = Snapshot_Evaporation_Block_Path(snapshot_root, snapshot_index);
	std::remove(path.c_str());
}

bool Write_Evaporation_Block_File(const std::string& snapshot_root, int snapshot_index, double snapshot_interval_sec, uint64_t run_id, double mass_gev, double sigma_cm2, const std::vector<CompactEvaporationEvent>& events, EvaporationLogState& state)
{
	const std::string path = Snapshot_Evaporation_Block_Path(snapshot_root, snapshot_index);
	if(!File_Is_Empty_Or_Missing(path))
	{
		if(Evaporation_Block_Metadata_Matches(path, run_id, snapshot_index, snapshot_interval_sec, mass_gev, sigma_cm2))
		{
			Recover_Evaporation_Log_State(path, state);
			return true;
		}
		if(std::remove(path.c_str()) != 0 && errno != ENOENT)
			return false;
	}

	std::vector<CompactEvaporationEvent> sorted_events = events;
	std::sort(sorted_events.begin(), sorted_events.end(), Evaporation_Event_Order);

	std::vector<CompactEvaporationEvent> pending_events;
	std::vector<EvaporationRecordKey> newly_written_keys;
	std::unordered_set<EvaporationRecordKey, EvaporationRecordKeyHash> pending_record_keys;
	for(const auto& event : sorted_events)
	{
		if(!Is_Completed_Evaporation_Event(event))
			continue;
		const EvaporationRecordKey key = Make_Evaporation_Record_Key(event.rank, event.trajectory_id);
		if(state.written_record_keys.count(key) != 0 || !pending_record_keys.insert(key).second)
			continue;
		pending_events.push_back(event);
		newly_written_keys.push_back(key);
	}

	Remove_Stale_Snapshot_Evaporation_Block(snapshot_root, snapshot_index);
	if(pending_events.empty())
		return true;

	const std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(snapshot_index);
	{
		std::ofstream file(tmp_path, std::ios::out | std::ios::trunc);
		if(!file.is_open())
			return false;
		Write_Evaporation_Block_Header(file, run_id, snapshot_index, snapshot_interval_sec, mass_gev, sigma_cm2);
		for(const auto& event : pending_events)
			Write_Evaporation_Log_Event(file, event);
		file.close();
		if(!file.good())
		{
			std::remove(tmp_path.c_str());
			return false;
		}
	}

	if(std::rename(tmp_path.c_str(), path.c_str()) != 0)
	{
		std::remove(tmp_path.c_str());
		return false;
	}

	state.written_record_keys.insert(newly_written_keys.begin(), newly_written_keys.end());
	return true;
}
#endif

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

bool Write_Snapshot_Evaporation_Time_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds, double mass_gev, double sigma_cm2, const std::vector<RankedSnapshotEvaporationEntry>& entries)
{
	std::vector<CompactEvaporationEvent> events;
	events.reserve(entries.size());
	for(const auto& entry : entries)
		events.push_back(Make_Log_Event(entry.rank, entry.entry));

	const std::string path = Snapshot_Evaporation_Time_File_Path(snapshot_root, snapshot_index, interval_seconds);
	return Write_Text_File_Atomically(path, snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = merged\n";
		file << "# snapshot_target_wall_time_s = " << Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# completed_evaporation_events_in_interval_only = 1\n";
		file << "# NOT_FOR_FINAL_SURVIVAL_ANALYSIS = 1\n";
		Write_Evaporation_Log_File_Header(file, mass_gev, sigma_cm2);
		Write_Evaporation_Log_Events(file, events);
	});
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
	const uint64_t new_evaporation_event_count = static_cast<uint64_t>(state.new_evaporation_events.size());
	Write_Binary_Value(file, new_evaporation_event_count);
	for(const auto& entry : state.new_evaporation_events)
		Write_Snapshot_Evaporation_Entry_Binary(file, entry);
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
	uint64_t new_evaporation_event_count = 0;
	Read_Binary_Value(file, new_evaporation_event_count);
	state.new_evaporation_events.clear();
	state.new_evaporation_events.resize(static_cast<size_t>(new_evaporation_event_count));
	for(size_t i = 0; i < state.new_evaporation_events.size(); i++)
		Read_Snapshot_Evaporation_Entry_Binary(file, state.new_evaporation_events[i]);

	if(!file)
		return false;

	return state.run_id == expected_run_id;
}

void Accumulate_Snapshot_Report_State(Snapshot_Report_State& report, const Rank_Snapshot_State& state)
{
	report.total_trajectories += state.local_total;
	report.captured_particles += state.local_captured;
	report.snapshot_bincount_captured_samples += state.local_captured;
	report.snapshot_bincount_not_captured_samples += (state.local_total - state.local_captured);
	const double lower_wall_time = (report.snapshot_time_label > report.snapshot_interval_seconds) ? (report.snapshot_time_label - report.snapshot_interval_seconds) : 0.0;
	const double upper_wall_time = static_cast<double>(report.snapshot_time_label);
	for(const auto& entry : state.new_evaporation_events)
	{
		if(entry.completion_wall_time_sec > lower_wall_time && entry.completion_wall_time_sec <= upper_wall_time)
		{
			RankedSnapshotEvaporationEntry ranked_entry;
			ranked_entry.rank = state.rank;
			ranked_entry.entry = entry;
			report.new_evaporation_events.push_back(ranked_entry);
		}
	}

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
		report.total_trajectories++;
		if(state.current_trajectory_captured)
		{
			report.captured_particles++;
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


}

bool Load_Snapshot_Report_State(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, Snapshot_Report_State& report)
{
	report = Snapshot_Report_State();
	report.snapshot_time_label = Snapshot_Time_Label_Seconds(snapshot_index, interval_seconds);
	report.snapshot_interval_seconds = interval_seconds;

	bool all_ranks_ready = true;

	for(int rank = 0; rank < mpi_processes; rank++)
	{
		Rank_Snapshot_State state;
		const std::string checkpoint_path = Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		if(Read_Rank_Snapshot_State(checkpoint_path, run_id, state))
		{
			Accumulate_Snapshot_Report_State(report, state);
			continue;
		}

		const std::string final_path = Rank_Snapshot_Final_Path(rank_snapshot_dir, rank);
		if(Read_Rank_Snapshot_State(final_path, run_id, state) && state.done && state.rank_elapsed_wall_sec <= report.snapshot_time_label)
		{
			Accumulate_Snapshot_Report_State(report, state);
			continue;
		}

		all_ranks_ready = false;
	}

	return all_ranks_ready;
}

bool Write_Snapshot_Report_File(const std::string& snapshot_root, int snapshot_index, double interval_seconds, double mass_gev, double sigma_cm2, const Snapshot_Report_State& report)
{
	if(!Write_Snapshot_Evaporation_Time_File(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report.new_evaporation_events))
		return false;

	if(!Write_Text_File_Atomically(Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds), snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = merged\n";
		file << "# Cumulative snapshot report\n";
		Write_Report_Header(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, false, report.snapshot_time_label, report.snapshot_interval_seconds);

		double snapshot_captured_samples = static_cast<double>(report.snapshot_bincount_captured_samples);
		double snapshot_not_captured_samples = static_cast<double>(report.snapshot_bincount_not_captured_samples);
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

bool Try_Write_Merged_Snapshot(const std::string& snapshot_root, const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes, uint64_t run_id, double mass_gev, double sigma_cm2, int caller_rank)
{
	if(caller_rank != 0)
		return false;

	const std::string snapshot_text_path = Snapshot_Text_File_Path(snapshot_root, snapshot_index, interval_seconds);
	const std::string snapshot_evaporation_path = Snapshot_Evaporation_Time_File_Path(snapshot_root, snapshot_index, interval_seconds);
	if(Snapshot_Text_File_Is_Merged(snapshot_text_path) && Snapshot_Text_File_Is_Merged(snapshot_evaporation_path))
	{
		Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
		return true;
	}

	Snapshot_Report_State report;
	if(!Load_Snapshot_Report_State(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes, run_id, report))
		return false;

	if(!Write_Snapshot_Report_File(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report))
		return false;

	Cleanup_Snapshot_Checkpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
	return true;
}

void Build_MPI_Gatherv_Layout(const std::vector<int>& item_counts, int fields_per_item, std::vector<int>& recv_counts, std::vector<int>& displacements, int& total_items)
{
	recv_counts.resize(item_counts.size());
	displacements.resize(item_counts.size());
	for(size_t i = 0; i < item_counts.size(); i++)
	{
		recv_counts[i] = item_counts[i] * fields_per_item;
		displacements[i] = (i == 0) ? 0 : displacements[i-1] + recv_counts[i-1];
	}
	total_items = std::accumulate(item_counts.begin(), item_counts.end(), 0);
}
}

// Retired with the block/manifest snapshot workflow.
#if 0
bool Recover_Evaporation_Time_File_From_Blocks(const std::string& snapshot_root, const std::string& output_path, uint64_t expected_run_id, double mass_gev, double sigma_cm2)
{
	std::vector<EvaporationBlockFile> block_files;
	if(!List_Evaporation_Block_Files(snapshot_root, block_files))
		return false;

	std::vector<CompactEvaporationEvent> recovered_events;
	EvaporationLogState recovered_state;
	bool found_matching_block = false;
	for(const auto& block_file : block_files)
	{
		EvaporationBlockMetadata metadata;
		if(!Read_Evaporation_Block_Metadata(block_file.path, metadata))
			return false;
		if(metadata.run_id != expected_run_id)
			continue;
		if(!Metadata_Double_Matches(metadata.mass_gev, mass_gev) || !Metadata_Double_Matches(metadata.sigma_cm2, sigma_cm2))
			continue;

		found_matching_block = true;
		if(!Read_Evaporation_Block_Events(block_file.path, recovered_events, recovered_state))
			return false;
	}

	if(!found_matching_block)
		return false;

	return Write_Final_Evaporation_Time_File(output_path, mass_gev, sigma_cm2, recovered_events);
}
#endif

Simulation_Data::Simulation_Data(unsigned int sample_size, unsigned int max_trajectories, double u_min, unsigned int iso_rings)
: requested_captured_particles(sample_size),
  normal_mode_mpi_sync_interval(NORMAL_MODE_MPI_SYNC_INTERVAL_FALLBACK),
  number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0),
  number_of_complete_evaporation_particles(0), number_of_censored_captured_particles(0),
  number_of_invalid_survival_captured_particles(0),
  number_of_initial_shift_failures(0), number_of_final_reflection_shift_failures(0), number_of_numerical_failures(0),
  average_number_of_scatterings(0.0), computing_time(0.0), early_stopped(false),
  mpi_rank(0), mpi_processes(1), isoreflection_rings(iso_rings), minimum_speed_threshold(u_min),
  number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)),
  data(iso_rings, std::vector<libphysica::DataPoint>())
{
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

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
}

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, unsigned long int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius      = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps       = max_free_steps;
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, obscura::DM_Distribution& halo_model, SnapshotConfig snapshot_cfg, unsigned int fixed_seed, bool capture_mode)
{
	if(capture_mode)
		snapshot_cfg.enabled = false;
	normal_mode_mpi_sync_interval = capture_mode ? 0UL : Normal_Mode_MPI_Sync_Interval(In_Units(DM.Sigma_Proton(), cm * cm));

	auto time_start = std::chrono::system_clock::now();
	unsigned long int local_captured = 0;
	unsigned long int local_total = 0;
	bool initial_shift_failure_warning_emitted = false;

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	simulator.max_trajectory_wall_time_sec = snapshot_cfg.max_trajectory_wall_time_sec;
	simulator.Enable_Capture_Mode(capture_mode);
	if(fixed_seed != 0)
	{
		const uint64_t rank_seed = static_cast<uint64_t>(fixed_seed) + 1000003ULL * static_cast<uint64_t>(mpi_rank);
		simulator.Fix_PRNG_Seed(static_cast<unsigned int>(rank_seed));
	}

	// Snapshot configuration
	const double snapshot_interval = (snapshot_cfg.interval_seconds > 0.0) ? snapshot_cfg.interval_seconds : 60.0;
	const double snapshot_mass_gev = In_Units(DM.mass, GeV);
	const double snapshot_sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);
	const std::string output_root = g_top_level_dir + "results_" + std::to_string(log10(snapshot_mass_gev)) + "_" + std::to_string(log10(snapshot_sigma_cm2)) + "/";
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
		snapshot_root = output_root + "snapshot/";
		rank_snapshot_dir = snapshot_root + "rank_snapshot/";
		int snapshot_init_ok = 1;
		if(mpi_rank == 0)
		{
			if(!Ensure_Directory_Exists(output_root))
				snapshot_init_ok = 0;

			if(snapshot_init_ok
			   && (!Ensure_Directory_Exists(snapshot_root)
			       || !Clear_Directory_Contents(snapshot_root)
			       || !Ensure_Directory_Exists(snapshot_root)
			       || !Ensure_Directory_Exists(rank_snapshot_dir)))
			{
				std::cerr << "Warning in Generate_Data(): failed to initialize snapshot directory; snapshots disabled for this run." << std::endl;
				snapshot_init_ok = 0;
			}
			if(snapshot_init_ok)
			{
				// Do not leave a previous run's final reports beside the new snapshots.
				std::remove(Join_Path(output_root, "bincount.txt").c_str());
				std::remove(Evaporation_Log_Path_From_Output_Dir(output_root).c_str());
				std::remove(Join_Path(output_root, "evaporation_diagnostics.txt").c_str());
			}

			if(snapshot_init_ok)
				snapshot_run_id = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
		}
		MPI_Bcast(&snapshot_init_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(!snapshot_init_ok)
		{
			snapshot_cfg.enabled = false;
		}
		else
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int rank_snapshot_dirs_ok =
			    (Ensure_Directory_Exists(snapshot_root)
			     && Ensure_Directory_Exists(rank_snapshot_dir)) ? 1 : 0;
			MPI_Allreduce(MPI_IN_PLACE, &rank_snapshot_dirs_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			if(!rank_snapshot_dirs_ok)
			{
				if(mpi_rank == 0)
					std::cerr << "Warning in Generate_Data(): failed to initialize rank snapshot directories; snapshots disabled for this run." << std::endl;
				snapshot_cfg.enabled = false;
			}
			else
				MPI_Bcast(&snapshot_run_id, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
		}
	}

	auto elapsed_since_start = [&]()
	{
		return 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - time_start).count();
	};

	size_t first_uncommitted_evaporation_entry = 0;

	auto build_rank_snapshot_state = [&](bool done, double snapshot_upper_wall_time, size_t evaporation_entry_begin, size_t evaporation_entry_end)
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
		(void)snapshot_upper_wall_time;
		state.new_evaporation_events.reserve(evaporation_entry_end - evaporation_entry_begin);
		for(size_t entry_index = evaporation_entry_begin; entry_index < evaporation_entry_end; entry_index++)
			state.new_evaporation_events.push_back(Make_Snapshot_Evaporation_Entry(compact_evaporation_events[entry_index]));
		return state;
	};

	int last_committed_snapshot_index = 0;
	auto try_write_ready_snapshots = [&](int max_snapshot_index)
	{
		if(mpi_rank != 0)
			return false;

		bool all_snapshots_merged = true;
		for(int snapshot_index = last_committed_snapshot_index + 1; snapshot_index <= max_snapshot_index; snapshot_index++)
		{
			if(!Try_Write_Merged_Snapshot(snapshot_root, rank_snapshot_dir, snapshot_index, snapshot_interval, mpi_processes, snapshot_run_id, snapshot_mass_gev, snapshot_sigma_cm2, mpi_rank))
			{
				all_snapshots_merged = false;
				break;
			}
			last_committed_snapshot_index = snapshot_index;
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
			const double snapshot_upper_wall_time = next_snapshot_index * snapshot_interval;
			size_t evaporation_entry_end = first_uncommitted_evaporation_entry;
			while(evaporation_entry_end < compact_evaporation_events.size()
			      && compact_evaporation_events[evaporation_entry_end].completion_wall_time_sec <= snapshot_upper_wall_time)
				evaporation_entry_end++;
			Rank_Snapshot_State state = build_rank_snapshot_state(false, snapshot_upper_wall_time, first_uncommitted_evaporation_entry, evaporation_entry_end);
			state.snapshot_index = next_snapshot_index;
			if(Write_Rank_Snapshot_State(Rank_Snapshot_Checkpoint_Path(rank_snapshot_dir, mpi_rank, next_snapshot_index, snapshot_interval), state))
				first_uncommitted_evaporation_entry = evaporation_entry_end;
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

	unsigned long int global_captured = 0;
	const unsigned long int mpi_sync_interval = capture_mode ? CAPTURE_MODE_MPI_SYNC_INTERVAL : normal_mode_mpi_sync_interval;
	auto select_trajectory_batch = [&](unsigned long int remaining_captures, unsigned long int& total_assigned)
	{
		total_assigned = 0;
		unsigned long int local_round_capacity = 0;
		if(local_total < max_trajectories_per_rank)
		{
			const unsigned long int local_capacity = max_trajectories_per_rank - local_total;
			local_round_capacity = std::min(local_capacity, mpi_sync_interval);
		}

		std::vector<unsigned long int> round_capacities(mpi_processes, 0);
		MPI_Allgather(&local_round_capacity, 1, MPI_UNSIGNED_LONG, round_capacities.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

		unsigned long int total_round_capacity = 0;
		for(const auto capacity : round_capacities)
			total_round_capacity += capacity;
		if(total_round_capacity == 0)
			return 0UL;

		const unsigned long int attempt_budget = capture_mode ? std::min(remaining_captures, total_round_capacity) : total_round_capacity;
		unsigned long int local_batch = 0;
		for(int rank = 0; rank < mpi_processes; rank++)
		{
			if(total_assigned >= attempt_budget)
				break;
			const unsigned long int assigned = std::min(round_capacities[rank], attempt_budget - total_assigned);
			if(rank == mpi_rank)
				local_batch = assigned;
			total_assigned += assigned;
		}
		return local_batch;
	};

	while(global_captured < requested_captured_particles)
	{
		const unsigned long int remaining_captures = requested_captured_particles - global_captured;
		unsigned long int assigned_trajectories = 0;
		const unsigned long int local_trajectories_this_round =
		    select_trajectory_batch(remaining_captures, assigned_trajectories);
		if(assigned_trajectories == 0)
		{
			early_stopped = true;
			break;
		}

		for(unsigned long int trajectory_in_round = 0; trajectory_in_round < local_trajectories_this_round; trajectory_in_round++)
		{
			Event IC = Initial_Conditions(halo_model, solar_model, simulator.PRNG);
			TrajectoryBincount failed_bincount;
			failed_bincount.termination_reason = TrajectoryTerminationReason::NumericalFailure;
			failed_bincount.survival_valid = false;
			const bool initial_shift_ok = Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
			if(!initial_shift_ok)
			{
				number_of_initial_shift_failures++;
				number_of_numerical_failures++;
			}
			Trajectory_Result trajectory = initial_shift_ok
			                                   ? simulator.Simulate(IC, DM, mpi_rank)
			                                   : Trajectory_Result(IC, IC, 0, failed_bincount);
			const double trajectory_completion_wall_time_sec = elapsed_since_start();

			local_total++;
			number_of_trajectories++;
			average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);
			const bool completed_outward_escape = Completed_Outward_Escape(trajectory.bincount.termination_reason);

			if(trajectory.bincount.is_captured)
			{
				number_of_captured_particles++;
				local_captured++;
				if(!capture_mode)
				{
					if(trajectory.bincount.event_observed)
						number_of_complete_evaporation_particles++;
					else
						number_of_invalid_survival_captured_particles++;
				}

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

					EvaporationRecord rec;
					if(Build_Evaporation_Record(trajectory.bincount, mpi_rank, number_of_trajectories, trajectory_completion_wall_time_sec, rec))
					{
						if(evaporation_diagnostics_enabled)
							evaporation_records.push_back(rec);
						CompactEvaporationEvent event;
						if(Make_Compact_Evaporation_Event(rec, event))
							compact_evaporation_events.push_back(event);
					}
				}

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

					if(trajectory.Particle_Free())
						number_of_free_particles++;
					else if(trajectory.Particle_Reflected())
					{
						number_of_reflected_particles++;
						if(Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU))
						{
							const double v_final = trajectory.final_event.Speed();
							if(trajectory.number_of_scatterings >= minimum_number_of_scatterings
							   && v_final > KDE_boundary_correction_factor * minimum_speed_threshold)
							{
								const unsigned int isoreflection_ring =
								    (isoreflection_rings == 1)
								    ? 0
								    : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
								data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
							}
						}
						else
						{
							number_of_final_reflection_shift_failures++;
							number_of_numerical_failures++;
						}
					}
				}

			}
		}

		MPI_Allreduce(&local_captured, &global_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		unsigned long int global_attempted_trajectories = 0;
		unsigned long int global_initial_shift_failures = 0;
		MPI_Allreduce(&number_of_trajectories, &global_attempted_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&number_of_initial_shift_failures, &global_initial_shift_failures, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

		if(global_attempted_trajectories > 0)
		{
			const double initial_shift_failure_fraction =
			    static_cast<double>(global_initial_shift_failures) / static_cast<double>(global_attempted_trajectories);
			if(initial_shift_failure_fraction > NUMERICAL_FAILURE_ABORT_FRACTION)
			{
				if(mpi_rank == 0)
					std::cerr << "Error in Generate_Data(): initial Kepler shift failure fraction "
					          << initial_shift_failure_fraction
					          << " exceeds abort threshold "
					          << NUMERICAL_FAILURE_ABORT_FRACTION
					          << ". Stopping this run to avoid biased capture statistics." << std::endl;
				early_stopped = true;
				break;
			}
			if(!initial_shift_failure_warning_emitted
			   && initial_shift_failure_fraction > NUMERICAL_FAILURE_WARNING_FRACTION)
			{
				if(mpi_rank == 0)
					std::cerr << "Warning in Generate_Data(): initial Kepler shift failure fraction "
					          << initial_shift_failure_fraction
					          << " exceeds warning threshold "
					          << NUMERICAL_FAILURE_WARNING_FRACTION
					          << "." << std::endl;
				initial_shift_failure_warning_emitted = true;
			}
		}

		// Progress bar (every 20%)
		if(mpi_rank == 0)
		{
			const double denominator = std::max(1u, requested_captured_particles);
			double progress = std::min(1.0, static_cast<double>(global_captured) / denominator);
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

	if(global_captured < requested_captured_particles)
		early_stopped = true;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();

	publish_checkpoint_snapshots();

	// Final snapshot state lets rank 0 merge only intervals that the run actually reached.
	if(snapshot_cfg.enabled)
	{
		Rank_Snapshot_State final_state = build_rank_snapshot_state(true, computing_time, first_uncommitted_evaporation_entry, compact_evaporation_events.size());
		final_state.snapshot_index = next_snapshot_index - 1;
		Write_Rank_Snapshot_State(Rank_Snapshot_Final_Path(rank_snapshot_dir, mpi_rank), final_state);

		MPI_Barrier(MPI_COMM_WORLD);
		double final_snapshot_elapsed = computing_time;
		MPI_Allreduce(MPI_IN_PLACE, &final_snapshot_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		int final_snapshot_index = std::max(0, static_cast<int>(std::floor(final_snapshot_elapsed / snapshot_interval)));
		if(final_snapshot_index == 0 || try_write_ready_snapshots(final_snapshot_index))
			Cleanup_Final_Snapshot_States(rank_snapshot_dir, mpi_processes);
	}

	if(mpi_rank == 0)
	{
		libphysica::Print_Progress_Bar(1.0, 0, 44, computing_time);
		std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Perform_MPI_Reductions(capture_mode);
}

void Simulation_Data::Perform_MPI_Reductions(bool capture_mode)
{
	MPI_Trace_Point(mpi_rank, "enter Perform_MPI_Reductions");
	average_number_of_scatterings *= number_of_trajectories;
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_trajectories");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_captured_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_initial_shift_failures");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_initial_shift_failures, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_final_reflection_shift_failures");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_final_reflection_shift_failures, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_numerical_failures");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_numerical_failures, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce average_number_of_scatterings");
	MPI_Allreduce(MPI_IN_PLACE, &average_number_of_scatterings, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings = (number_of_trajectories > 0) ? average_number_of_scatterings / number_of_trajectories : 0.0;

	// Reduce early_stopped flag (any rank early stopped => global flag)
	int local_es = early_stopped ? 1 : 0;
	int global_es = 0;
	MPI_Trace_Point(mpi_rank, "before allreduce early_stopped");
	MPI_Allreduce(&local_es, &global_es, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	early_stopped = (global_es > 0);

	if(capture_mode)
	{
		MPI_Trace_Point(mpi_rank, "before allreduce computing_time capture");
		MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Trace_Point(mpi_rank, "leave Perform_MPI_Reductions capture");
		return;
	}

	MPI_Trace_Point(mpi_rank, "before allreduce number_of_free_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_free_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_reflected_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_reflected_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_complete_evaporation_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_complete_evaporation_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_censored_captured_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_censored_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce number_of_invalid_survival_captured_particles");
	MPI_Allreduce(MPI_IN_PLACE, &number_of_invalid_survival_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

	// Reduce bincount histograms across all ranks
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_v2dt_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce captured_v2dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "before allreduce not_captured_v2dt_sq_hist");
	MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "after histogram allreduces");
	if(evaporation_diagnostics_enabled)
	{
		const int local_evap_count = static_cast<int>(evaporation_records.size());
		std::vector<int> evap_counts(mpi_processes, 0);
		MPI_Trace_Point(mpi_rank, "before gather diagnostic evap counts");
		MPI_Gather(&local_evap_count, 1, MPI_INT, mpi_rank == 0 ? evap_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		constexpr int EVAPORATION_MPI_INT_FIELDS = 8;
		constexpr int EVAPORATION_MPI_UINT_FIELDS = 2;
		constexpr int EVAPORATION_MPI_DOUBLE_FIELDS = 14;
		std::vector<int> local_evap_ints(local_evap_count * EVAPORATION_MPI_INT_FIELDS);
		std::vector<unsigned long long> local_evap_uints(local_evap_count * EVAPORATION_MPI_UINT_FIELDS);
		std::vector<double> local_evap_doubles(local_evap_count * EVAPORATION_MPI_DOUBLE_FIELDS);
		for(int i = 0; i < local_evap_count; i++)
		{
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i] = evaporation_records[i].rank;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 1] = evaporation_records[i].event_observed ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 2] = evaporation_records[i].boundary_escape_observed ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 3] = evaporation_records[i].survival_valid ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 4] = evaporation_records[i].numerically_invalid_escape ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 5] = evaporation_records[i].censored ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] = evaporation_records[i].truncated ? 1 : 0;
			local_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7] = TerminationReason_Index(evaporation_records[i].termination_reason);

			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i] = static_cast<unsigned long long>(evaporation_records[i].trajectory_id);
			local_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1] = static_cast<unsigned long long>(evaporation_records[i].number_of_scatterings);

			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i] = evaporation_records[i].completion_wall_time_sec;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 1] = evaporation_records[i].t_evap;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 2] = evaporation_records[i].t_capture;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 3] = evaporation_records[i].t_final_unbinding_scatter;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 4] = evaporation_records[i].t_boundary_escape;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 5] = evaporation_records[i].t_termination;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 6] = evaporation_records[i].observed_lifetime;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 7] = evaporation_records[i].lifetime_unbinding;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 8] = evaporation_records[i].lifetime_boundary;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 9] = evaporation_records[i].r_first_negative_km;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 10] = evaporation_records[i].E_first_negative_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 11] = evaporation_records[i].dE_first_negative_from_prev_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 12] = evaporation_records[i].max_free_energy_drift_eV;
			local_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 13] = evaporation_records[i].max_free_energy_drift_rel;
		}

		std::vector<int> recv_counts, displacements;
		std::vector<int> global_evap_ints;
		std::vector<unsigned long long> global_evap_uints;
		std::vector<double> global_evap_doubles;
		int total_evap = 0;
		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_INT_FIELDS, recv_counts, displacements, total_evap);
			global_evap_ints.resize(total_evap * EVAPORATION_MPI_INT_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap ints");
		MPI_Gatherv(local_evap_ints.data(), local_evap_count * EVAPORATION_MPI_INT_FIELDS, MPI_INT,
		            mpi_rank == 0 ? global_evap_ints.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_INT, 0, MPI_COMM_WORLD);

		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_UINT_FIELDS, recv_counts, displacements, total_evap);
			global_evap_uints.resize(total_evap * EVAPORATION_MPI_UINT_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap uints");
		MPI_Gatherv(local_evap_uints.data(), local_evap_count * EVAPORATION_MPI_UINT_FIELDS, MPI_UNSIGNED_LONG_LONG,
		            mpi_rank == 0 ? global_evap_uints.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

		if(mpi_rank == 0)
		{
			Build_MPI_Gatherv_Layout(evap_counts, EVAPORATION_MPI_DOUBLE_FIELDS, recv_counts, displacements, total_evap);
			global_evap_doubles.resize(total_evap * EVAPORATION_MPI_DOUBLE_FIELDS);
		}
		MPI_Trace_Point(mpi_rank, "before gatherv diagnostic evap doubles");
		MPI_Gatherv(local_evap_doubles.data(), local_evap_count * EVAPORATION_MPI_DOUBLE_FIELDS, MPI_DOUBLE,
		            mpi_rank == 0 ? global_evap_doubles.data() : nullptr,
		            mpi_rank == 0 ? recv_counts.data() : nullptr,
		            mpi_rank == 0 ? displacements.data() : nullptr,
		            MPI_DOUBLE, 0, MPI_COMM_WORLD);

		evaporation_records.clear();
		if(mpi_rank == 0)
		{
			evaporation_records.resize(total_evap);
			for(int i = 0; i < total_evap; i++)
			{
				evaporation_records[i].rank = global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i];
				evaporation_records[i].trajectory_id = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i]);
				evaporation_records[i].completion_wall_time_sec = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i];
				evaporation_records[i].t_evap = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 1];
				evaporation_records[i].t_capture = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 2];
				evaporation_records[i].t_final_unbinding_scatter = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 3];
				evaporation_records[i].t_boundary_escape = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 4];
				evaporation_records[i].t_termination = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 5];
				evaporation_records[i].observed_lifetime = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 6];
				evaporation_records[i].lifetime_unbinding = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 7];
				evaporation_records[i].lifetime_boundary = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 8];
				evaporation_records[i].r_first_negative_km = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 9];
				evaporation_records[i].E_first_negative_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 10];
				evaporation_records[i].dE_first_negative_from_prev_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 11];
				evaporation_records[i].event_observed = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 1] != 0);
				evaporation_records[i].boundary_escape_observed = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 2] != 0);
				evaporation_records[i].survival_valid = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 3] != 0);
				evaporation_records[i].numerically_invalid_escape = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 4] != 0);
				evaporation_records[i].censored = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 5] != 0);
				evaporation_records[i].truncated = (global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 6] != 0);
				evaporation_records[i].termination_reason = static_cast<TrajectoryTerminationReason>(TerminationReason_Index(static_cast<TrajectoryTerminationReason>(global_evap_ints[EVAPORATION_MPI_INT_FIELDS*i + 7])));
				evaporation_records[i].max_free_energy_drift_eV = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 12];
				evaporation_records[i].max_free_energy_drift_rel = global_evap_doubles[EVAPORATION_MPI_DOUBLE_FIELDS*i + 13];
				evaporation_records[i].number_of_scatterings = static_cast<unsigned long int>(global_evap_uints[EVAPORATION_MPI_UINT_FIELDS*i + 1]);
			}
		}
	}
	else
	{
		const int local_evap_count = static_cast<int>(compact_evaporation_events.size());
		std::vector<int> evap_counts(mpi_processes, 0);
		MPI_Trace_Point(mpi_rank, "before gather compact evap counts");
		MPI_Gather(&local_evap_count, 1, MPI_INT, mpi_rank == 0 ? evap_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int max_evap_count = 0;
		MPI_Trace_Point(mpi_rank, "before allreduce compact evap max count");
		MPI_Allreduce(&local_evap_count, &max_evap_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		std::vector<CompactEvaporationEvent> local_evap_padded(max_evap_count);
		for(int i = 0; i < local_evap_count; i++)
			local_evap_padded[i] = compact_evaporation_events[i];

		std::vector<CompactEvaporationEvent> global_evap_padded;
		if(mpi_rank == 0)
			global_evap_padded.resize(static_cast<size_t>(mpi_processes) * static_cast<size_t>(max_evap_count));

		const int padded_evap_bytes = max_evap_count * static_cast<int>(sizeof(CompactEvaporationEvent));
		MPI_Trace_Point(mpi_rank, "before gather compact evap padded bytes");
		MPI_Gather(max_evap_count == 0 ? nullptr : static_cast<void*>(local_evap_padded.data()),
		            padded_evap_bytes,
		            MPI_BYTE,
		            mpi_rank == 0 && max_evap_count > 0 ? static_cast<void*>(global_evap_padded.data()) : nullptr,
		            padded_evap_bytes,
		            MPI_BYTE, 0, MPI_COMM_WORLD);

		compact_evaporation_events.clear();
		if(mpi_rank == 0)
		{
			int total_evap = std::accumulate(evap_counts.begin(), evap_counts.end(), 0);
			compact_evaporation_events.reserve(total_evap);
			for(int rank = 0; rank < mpi_processes; rank++)
			{
				const size_t rank_offset = static_cast<size_t>(rank) * static_cast<size_t>(max_evap_count);
				for(int i = 0; i < evap_counts[rank]; i++)
					compact_evaporation_events.push_back(global_evap_padded[rank_offset + static_cast<size_t>(i)]);
			}
		}
	}

	MPI_Trace_Point(mpi_rank, "before allreduce computing_time final");
	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Trace_Point(mpi_rank, "leave Perform_MPI_Reductions");
}

void Simulation_Data::Write_Output_Files(const std::string& output_dir, obscura::DM_Particle& DM)
{
	if(mpi_rank != 0)
		return;

	if(!Ensure_Directory_Exists(output_dir))
	{
		std::cerr << "Warning in Write_Output_Files(): failed to create output directory " << output_dir << std::endl;
		return;
	}

	double mass_gev = In_Units(DM.mass, GeV);
	double sigma_cm2 = In_Units(DM.Sigma_Proton(), cm * cm);

	auto write_header = [&](std::ofstream& f) {
		Write_Report_Header(f, mass_gev, sigma_cm2, number_of_trajectories, number_of_captured_particles, early_stopped);
		f << "# valid_trajectories = " << Valid_Trajectories() << "\n";
		f << "# numerical_failures = " << number_of_numerical_failures << "\n";
		f << "# initial_shift_failures = " << number_of_initial_shift_failures << "\n";
		f << "# final_reflection_shift_failures = " << number_of_final_reflection_shift_failures << "\n";
		f << "# normal_mode_mpi_sync_interval = " << normal_mode_mpi_sync_interval << "\n";
		f << "# capture_rate_valid = " << std::fixed << std::setprecision(8) << Capture_Ratio_Valid() << "\n";
		f << "# numerical_failure_rate = " << std::fixed << std::setprecision(8) << Numerical_Failure_Ratio() << "\n";
	};

	std::remove((output_dir + "/captured_bincount.txt").c_str());
	std::remove((output_dir + "/not_captured_bincount.txt").c_str());
	std::remove((output_dir + "/evaporation_diagnostics.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_summary.txt").c_str());
	std::remove((output_dir + "/evaporation_" + "mode_" + "bincount.txt").c_str());
	std::remove((output_dir + "/computation_" + "time_summary.txt").c_str());

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

	// 2. Final evaporation-time list.  This is intentionally the only final
	// evaporation report; snapshot files are intermediate progress reports.
	bool evaporation_times_ok = Write_Final_Evaporation_Time_File(Evaporation_Log_Path_From_Output_Dir(output_dir), mass_gev, sigma_cm2, compact_evaporation_events);
	if(!evaporation_times_ok)
		std::cerr << "Warning in Write_Output_Files(): failed to write evaporation_times.txt" << std::endl;


}

double Simulation_Data::Free_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_free_particles / number_of_trajectories : 0.0;
}
double Simulation_Data::Capture_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_captured_particles / number_of_trajectories : 0.0;
}
double Simulation_Data::Reflection_Ratio(int isoreflection_ring) const
{
	if(isoreflection_ring < 0)
		return (number_of_trajectories > 0) ? 1.0 * number_of_reflected_particles / number_of_trajectories : 0.0;
	else
		return (number_of_trajectories > 0) ? 1.0 * data[isoreflection_ring].size() / number_of_trajectories : 0.0;
}

unsigned long int Simulation_Data::Valid_Trajectories() const
{
	return (number_of_trajectories > number_of_initial_shift_failures) ? (number_of_trajectories - number_of_initial_shift_failures) : 0UL;
}

double Simulation_Data::Free_Ratio_Valid() const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	return (valid_trajectories > 0) ? 1.0 * number_of_free_particles / valid_trajectories : 0.0;
}

double Simulation_Data::Capture_Ratio_Valid() const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	return (valid_trajectories > 0) ? 1.0 * number_of_captured_particles / valid_trajectories : 0.0;
}

double Simulation_Data::Reflection_Ratio_Valid(int isoreflection_ring) const
{
	const unsigned long int valid_trajectories = Valid_Trajectories();
	if(valid_trajectories == 0)
		return 0.0;
	if(isoreflection_ring < 0)
		return 1.0 * number_of_reflected_particles / valid_trajectories;
	else
		return 1.0 * data[isoreflection_ring].size() / valid_trajectories;
}

double Simulation_Data::Numerical_Failure_Ratio() const
{
	return (number_of_trajectories > 0) ? 1.0 * number_of_numerical_failures / number_of_trajectories : 0.0;
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
		double p = Capture_Ratio();
		double p_valid = Capture_Ratio_Valid();
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
		          << "Valid trajectories:\t\t" << Valid_Trajectories() << std::endl
		          << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
		          << "Capture rate raw:\t\t" << std::fixed << std::setprecision(8) << p << std::endl
		          << "Capture rate valid:\t\t" << std::fixed << std::setprecision(8) << p_valid << std::endl
		          << "Capture rate -error (95%):\t" << std::fixed << std::setprecision(8) << (p - ci_lower) << std::endl
		          << "Capture rate +error (95%):\t" << std::fixed << std::setprecision(8) << (ci_upper - p) << std::endl
		          << "Capture rate 95% CI raw:\t[" << std::fixed << std::setprecision(8) << ci_lower << ", " << ci_upper << "]" << std::endl
		          << "Numerical failure count:\t" << number_of_numerical_failures << std::endl
		          << "Initial shift failures:\t\t" << number_of_initial_shift_failures << std::endl
		          << "Final reflection shift failures:\t" << number_of_final_reflection_shift_failures << std::endl
		          << "Numerical failure rate:\t\t" << std::fixed << std::setprecision(8) << Numerical_Failure_Ratio() << std::endl;

		if(early_stopped)
			std::cout << "*** EARLY STOP: max_trajectories reached ***" << std::endl;
		if(Numerical_Failure_Ratio() > NUMERICAL_FAILURE_WARNING_FRACTION)
			std::cout << "*** WARNING: numerical failure rate exceeded "
			          << NUMERICAL_FAILURE_WARNING_FRACTION << " ***" << std::endl;

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
				  << "Valid trajectories:\t\t" << Valid_Trajectories() << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles raw [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles raw [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles raw [%]:\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl
				  << "Free particles valid [%]:\t" << libphysica::Round(100.0 * Free_Ratio_Valid()) << std::endl
				  << "Reflected particles valid [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio_Valid()) << std::endl
				  << "Captured particles valid [%]:\t" << libphysica::Round(100.0 * Capture_Ratio_Valid()) << std::endl
				  << "Captured count:\t\t\t" << number_of_captured_particles << std::endl
				  << "Normal-mode MPI sync interval:\t" << normal_mode_mpi_sync_interval << std::endl
				  << "Numerical failure count:\t" << number_of_numerical_failures << std::endl
				  << "Initial shift failures:\t\t" << number_of_initial_shift_failures << std::endl
				  << "Final reflection shift failures:\t" << number_of_final_reflection_shift_failures << std::endl
				  << "Numerical failure rate:\t\t" << std::fixed << std::setprecision(6) << Numerical_Failure_Ratio() << std::endl
				  << "Complete evaporation count:\t" << number_of_complete_evaporation_particles << std::endl
				  << "Censored captured count:\t" << number_of_censored_captured_particles << std::endl
				  << "Invalid survival count:\t\t" << number_of_invalid_survival_captured_particles << std::endl;

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
		if(Numerical_Failure_Ratio() > NUMERICAL_FAILURE_WARNING_FRACTION)
			std::cout << "*** WARNING: numerical failure rate exceeded "
			          << NUMERICAL_FAILURE_WARNING_FRACTION << " ***" << std::endl;

		// Median for observed unbinding events only; censored records belong in survival analysis.
		std::vector<double> observed_unbinding_lifetimes;
		if(evaporation_diagnostics_enabled)
		{
			for(const auto& rec : evaporation_records)
			{
				if(rec.survival_valid && rec.event_observed && Has_Positive_Evaporation_Time(rec.lifetime_unbinding))
					observed_unbinding_lifetimes.push_back(rec.lifetime_unbinding);
			}
		}
		else
		{
			for(const auto& event : compact_evaporation_events)
			{
				if(Has_Positive_Evaporation_Time(event.lifetime_unbinding))
					observed_unbinding_lifetimes.push_back(event.lifetime_unbinding);
			}
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


		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DaMaSCUS_SUN

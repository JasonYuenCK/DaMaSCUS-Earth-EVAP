#include "Snapshot_IO.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

namespace DaMaSCUS_SUN
{
namespace
{
template<typename T>
void WriteBinaryValue(std::ofstream& file, const T& value)
{
	file.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void WriteBinaryArray(std::ofstream& file, const std::array<T, N>& values)
{
	file.write(reinterpret_cast<const char*>(values.data()), N * sizeof(T));
}

template<typename T>
void ReadBinaryValue(std::ifstream& file, T& value)
{
	file.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T, size_t N>
void ReadBinaryArray(std::ifstream& file, std::array<T, N>& values)
{
	file.read(reinterpret_cast<char*>(values.data()), N * sizeof(T));
}

void WriteSnapshotEvaporationEntryBinary(std::ofstream& file, const SnapshotEvaporationProgressEntry& entry)
{
	WriteBinaryValue(file, entry.trajectory_id);
	WriteBinaryValue(file, entry.completion_wall_time_sec);
	WriteBinaryValue(file, entry.lifetime_unbinding_sec);
}

void ReadSnapshotEvaporationEntryBinary(std::ifstream& file, SnapshotEvaporationProgressEntry& entry)
{
	ReadBinaryValue(file, entry.trajectory_id);
	ReadBinaryValue(file, entry.completion_wall_time_sec);
	ReadBinaryValue(file, entry.lifetime_unbinding_sec);
}

CompactEvaporationEvent MakeLogEvent(int rank, const SnapshotEvaporationProgressEntry& entry)
{
	CompactEvaporationEvent event;
	event.rank = rank;
	event.trajectory_id = entry.trajectory_id;
	event.completion_wall_time_sec = entry.completion_wall_time_sec;
	event.lifetime_unbinding = entry.lifetime_unbinding_sec;
	return event;
}

bool HasBincountContribution(
	const std::array<double, NUM_BINS>& dt_hist,
	const std::array<double, NUM_BINS>& v2dt_hist)
{
	return std::any_of(dt_hist.begin(), dt_hist.end(), [](double value) { return value != 0.0; })
	    || std::any_of(v2dt_hist.begin(), v2dt_hist.end(), [](double value) { return value != 0.0; });
}

double SnapshotBinError(double sum, double sum_sq, double count)
{
	if(count <= 1.0)
		return 0.0;

	const double mean = sum / count;
	double variance = sum_sq / count - mean * mean;
	if(variance < 0.0)
		variance = 0.0;
	return std::sqrt(count * variance);
}

double CaptureRatioFromCounts(uint64_t total_trajectories, uint64_t captured_particles)
{
	return (total_trajectories > 0) ? static_cast<double>(captured_particles) / static_cast<double>(total_trajectories) : 0.0;
}

bool IsCompletedEvaporationEvent(const CompactEvaporationEvent& event)
{
	return std::isfinite(event.lifetime_unbinding) && event.lifetime_unbinding >= 0.0;
}

bool EvaporationEventOrder(const CompactEvaporationEvent& lhs, const CompactEvaporationEvent& rhs)
{
	if(lhs.lifetime_unbinding != rhs.lifetime_unbinding)
		return lhs.lifetime_unbinding < rhs.lifetime_unbinding;
	if(lhs.rank != rhs.rank)
		return lhs.rank < rhs.rank;
	return lhs.trajectory_id < rhs.trajectory_id;
}

void WriteEvaporationLogEvent(std::ostream& file, const CompactEvaporationEvent& event)
{
	file << event.rank << "\t" << event.trajectory_id << "\t" << std::scientific << std::setprecision(10)
	     << event.lifetime_unbinding << "\n";
}

void WriteEvaporationLogFileHeader(std::ofstream& file, double mass_gev, double sigma_cm2)
{
	file << "# DaMaSCUS-SUN snapshot evaporation times\n";
	file << "# format_version = 3\n";
	file << "# DIAGNOSTIC_ONLY = 1\n";
	file << "# NOT_FOR_FINAL_SURVIVAL_ANALYSIS = 1\n";
	file << "# completion_time_selected = 1\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# sorted_by = lifetime_unbinding_sec rank trajectory_id\n";
	file << "# rank trajectory_id lifetime_unbinding_sec\n";
}

void WriteEvaporationLogEvents(std::ostream& file, const std::vector<CompactEvaporationEvent>& events)
{
	std::vector<CompactEvaporationEvent> sorted_events = events;
	std::sort(sorted_events.begin(), sorted_events.end(), EvaporationEventOrder);
	for(const auto& event : sorted_events)
	{
		if(IsCompletedEvaporationEvent(event))
			WriteEvaporationLogEvent(file, event);
	}
}

bool WriteTextFileAtomically(const std::string& path, int unique_tag, const std::function<void(std::ofstream&)>& writer)
{
	const std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(unique_tag);
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

std::string JoinInts(const std::vector<int>& values)
{
	if(values.empty())
		return "";

	std::ostringstream stream;
	for(size_t i = 0; i < values.size(); i++)
	{
		if(i != 0)
			stream << ",";
		stream << values[i];
	}
	return stream.str();
}

bool SnapshotTextFileIsMerged(const std::string& path)
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

struct SnapshotReportState
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
	std::vector<SnapshotRankedEvaporationEntry> new_evaporation_events;
};

void AccumulateSnapshotReportState(SnapshotReportState& report, const SnapshotRankState& state)
{
	report.total_trajectories += state.local_total;
	report.captured_particles += state.local_captured;
	report.snapshot_bincount_captured_samples += state.bincount_captured_samples;
	report.snapshot_bincount_not_captured_samples += state.bincount_not_captured_samples;

	const double lower_wall_time = (report.snapshot_time_label > report.snapshot_interval_seconds) ? (report.snapshot_time_label - report.snapshot_interval_seconds) : 0.0;
	const double upper_wall_time = static_cast<double>(report.snapshot_time_label);
	for(const auto& entry : state.new_evaporation_events)
	{
		if(entry.completion_wall_time_sec > lower_wall_time && entry.completion_wall_time_sec <= upper_wall_time)
		{
			SnapshotRankedEvaporationEntry ranked_entry;
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

	if(state.trajectory_in_progress && HasBincountContribution(state.current_trajectory_dt_hist, state.current_trajectory_v2dt_hist))
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

SnapshotMergeResult LoadSnapshotReportState(
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	SnapshotReportState& report)
{
	SnapshotMergeResult result;
	report = SnapshotReportState();
	report.snapshot_time_label = SnapshotTimeLabelSeconds(snapshot_index, interval_seconds);
	report.snapshot_interval_seconds = interval_seconds;

	for(int rank = 0; rank < mpi_processes; rank++)
	{
		SnapshotRankState state;
		const std::string checkpoint_path = SnapshotRankCheckpointPath(rank_snapshot_dir, rank, snapshot_index, interval_seconds);
		if(ReadSnapshotRankState(checkpoint_path, run_id, state))
		{
			AccumulateSnapshotReportState(report, state);
			result.ready_ranks.push_back(rank);
			continue;
		}

		const std::string final_path = SnapshotRankFinalPath(rank_snapshot_dir, rank);
		if(ReadSnapshotRankState(final_path, run_id, state) && state.done && state.rank_elapsed_wall_sec <= report.snapshot_time_label)
		{
			AccumulateSnapshotReportState(report, state);
			result.ready_ranks.push_back(rank);
			continue;
		}

		result.missing_ranks.push_back(rank);
	}

	if(result.ready_ranks.empty())
		result.status = SnapshotMergeStatus::NoRanksReady;
	else if(result.missing_ranks.empty())
		result.status = SnapshotMergeStatus::Merged;
	else
		result.status = SnapshotMergeStatus::Partial;

	return result;
}

void WriteReportHeader(
	std::ofstream& file,
	double mass_gev,
	double sigma_cm2,
	uint64_t total_trajectories,
	uint64_t captured_particles,
	long long snapshot_time_label,
	double snapshot_interval_seconds)
{
	file << "# snapshot_target_wall_time_s = " << snapshot_time_label << "\n";
	file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << snapshot_interval_seconds << "\n";
	file << "# DM_mass_GeV = " << std::scientific << std::setprecision(6) << mass_gev << "\n";
	file << "# DM_sigma_cm2 = " << std::scientific << std::setprecision(6) << sigma_cm2 << "\n";
	file << "# total_trajectories = " << total_trajectories << "\n";
	file << "# captured_particles = " << captured_particles << "\n";

	const double p_cap = CaptureRatioFromCounts(total_trajectories, captured_particles);
	file << "# capture_rate = " << std::fixed << std::setprecision(8) << p_cap << "\n";

	const double N = static_cast<double>(total_trajectories);
	const double z = 1.96;
	const double sigma_p = (N > 0.0) ? std::sqrt(p_cap * (1.0 - p_cap) / N) : 0.0;
	file << "# capture_rate_err = " << std::fixed << std::setprecision(8) << sigma_p << "\n";
	if(N > 0.0)
	{
		const double denom = 1.0 + z * z / N;
		const double center = p_cap + z * z / (2.0 * N);
		const double spread = z * std::sqrt(p_cap * (1.0 - p_cap) / N + z * z / (4.0 * N * N));
		double ci_lower = (center - spread) / denom;
		double ci_upper = (center + spread) / denom;
		if(ci_lower < 0.0) ci_lower = 0.0;
		if(ci_upper > 1.0) ci_upper = 1.0;
		file << "# capture_rate_CI_95_lower = " << std::fixed << std::setprecision(8) << ci_lower << "\n";
		file << "# capture_rate_CI_95_upper = " << std::fixed << std::setprecision(8) << ci_upper << "\n";
	}
}

bool WriteSnapshotEvaporationTimeFile(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	double mass_gev,
	double sigma_cm2,
	const std::vector<SnapshotRankedEvaporationEntry>& entries,
	const SnapshotMergeResult& merge_result)
{
	std::vector<CompactEvaporationEvent> events;
	events.reserve(entries.size());
	for(const auto& entry : entries)
		events.push_back(MakeLogEvent(entry.rank, entry.entry));

	const std::string status = (merge_result.status == SnapshotMergeStatus::Merged) ? "merged" : "partial";
	const std::string path = SnapshotEvaporationTimeFilePath(snapshot_root, snapshot_index, interval_seconds);
	return WriteTextFileAtomically(path, snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = " << status << "\n";
		file << "# snapshot_target_wall_time_s = " << SnapshotTimeLabelSeconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# ready_ranks = " << JoinInts(merge_result.ready_ranks) << "\n";
		file << "# missing_ranks = " << JoinInts(merge_result.missing_ranks) << "\n";
		file << "# completed_evaporation_events_in_interval_only = 1\n";
		WriteEvaporationLogFileHeader(file, mass_gev, sigma_cm2);
		WriteEvaporationLogEvents(file, events);
	});
}

bool WriteSnapshotReportFile(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	double mass_gev,
	double sigma_cm2,
	const SnapshotReportState& report,
	const SnapshotMergeResult& merge_result)
{
	if(merge_result.status != SnapshotMergeStatus::Merged && SnapshotTextFileIsMerged(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds)))
		return true;

	if(!WriteSnapshotEvaporationTimeFile(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report.new_evaporation_events, merge_result))
		return false;

	const std::string status = (merge_result.status == SnapshotMergeStatus::Merged) ? "merged" : "partial";
	return WriteTextFileAtomically(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds), snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = " << status << "\n";
		file << "# Cumulative snapshot report\n";
		file << "# ready_ranks = " << JoinInts(merge_result.ready_ranks) << "\n";
		file << "# missing_ranks = " << JoinInts(merge_result.missing_ranks) << "\n";
		if(merge_result.status == SnapshotMergeStatus::Partial)
			file << "# WARNING = partial snapshot; not all ranks reached this checkpoint\n";
		WriteReportHeader(file, mass_gev, sigma_cm2, report.total_trajectories, report.captured_particles, report.snapshot_time_label, report.snapshot_interval_seconds);

		const double snapshot_captured_samples = static_cast<double>(report.snapshot_bincount_captured_samples);
		const double snapshot_not_captured_samples = static_cast<double>(report.snapshot_bincount_not_captured_samples);
		file << "#\n";
		file << "# [Bincount histogram]\n";
		file << "# bin_index  cap_dt[s]  cap_v2dt[km2/s]  cap_err_dt[s]  cap_err_v2dt[km2/s]  not_cap_dt[s]  not_cap_v2dt[km2/s]  not_cap_err_dt[s]  not_cap_err_v2dt[km2/s]\n";
		for(int bin = 0; bin < NUM_BINS; bin++)
		{
			const double cap_err_dt = SnapshotBinError(report.captured_dt_hist[bin], report.captured_dt_sq_hist[bin], snapshot_captured_samples);
			const double cap_err_v2dt = SnapshotBinError(report.captured_v2dt_hist[bin], report.captured_v2dt_sq_hist[bin], snapshot_captured_samples);
			const double not_cap_err_dt = SnapshotBinError(report.not_captured_dt_hist[bin], report.not_captured_dt_sq_hist[bin], snapshot_not_captured_samples);
			const double not_cap_err_v2dt = SnapshotBinError(report.not_captured_v2dt_hist[bin], report.not_captured_v2dt_sq_hist[bin], snapshot_not_captured_samples);

			file << bin << "\t" << std::scientific << std::setprecision(10)
			     << report.captured_dt_hist[bin] << "\t" << report.captured_v2dt_hist[bin]
			     << "\t" << cap_err_dt << "\t" << cap_err_v2dt
			     << "\t" << report.not_captured_dt_hist[bin] << "\t" << report.not_captured_v2dt_hist[bin]
			     << "\t" << not_cap_err_dt << "\t" << not_cap_err_v2dt << "\n";
		}
	});
}
}

long long SnapshotTimeLabelSeconds(int snapshot_index, double interval_seconds)
{
	return static_cast<long long>(std::llround(snapshot_index * interval_seconds));
}

std::string SnapshotTextFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s.txt";
}

std::string SnapshotEvaporationTimeFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds)
{
	return snapshot_root + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s_evaporation_times.txt";
}

std::string SnapshotRankCheckpointPath(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds)
{
	return rank_snapshot_dir + "snapshot_" + std::to_string(SnapshotTimeLabelSeconds(snapshot_index, interval_seconds)) + "s_rank" + std::to_string(rank) + ".bin";
}

std::string SnapshotRankFinalPath(const std::string& rank_snapshot_dir, int rank)
{
	return rank_snapshot_dir + "rank" + std::to_string(rank) + "_final.bin";
}

SnapshotEvaporationProgressEntry MakeSnapshotEvaporationProgressEntry(const CompactEvaporationEvent& event)
{
	SnapshotEvaporationProgressEntry entry;
	entry.trajectory_id = static_cast<uint64_t>(event.trajectory_id);
	entry.completion_wall_time_sec = event.completion_wall_time_sec;
	entry.lifetime_unbinding_sec = event.lifetime_unbinding;
	return entry;
}

bool WriteSnapshotRankState(const std::string& path, const SnapshotRankState& state)
{
	const std::string tmp_path = path + ".tmp." + std::to_string(getpid()) + "." + std::to_string(state.rank);
	std::ofstream file(tmp_path, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
		return false;

	WriteBinaryValue(file, state.run_id);
	WriteBinaryValue(file, state.snapshot_index);
	WriteBinaryValue(file, state.rank);
	WriteBinaryValue(file, state.done);
	WriteBinaryValue(file, state.trajectory_in_progress);
	WriteBinaryValue(file, state.local_captured);
	WriteBinaryValue(file, state.local_total);
	WriteBinaryValue(file, state.bincount_captured_samples);
	WriteBinaryValue(file, state.bincount_not_captured_samples);
	WriteBinaryValue(file, state.current_trajectory_id);
	WriteBinaryValue(file, state.rank_elapsed_wall_sec);
	WriteBinaryValue(file, state.current_trajectory_captured);
	WriteBinaryArray(file, state.current_trajectory_dt_hist);
	WriteBinaryArray(file, state.current_trajectory_v2dt_hist);
	WriteBinaryArray(file, state.captured_dt_hist);
	WriteBinaryArray(file, state.captured_v2dt_hist);
	WriteBinaryArray(file, state.captured_dt_sq_hist);
	WriteBinaryArray(file, state.captured_v2dt_sq_hist);
	WriteBinaryArray(file, state.not_captured_dt_hist);
	WriteBinaryArray(file, state.not_captured_v2dt_hist);
	WriteBinaryArray(file, state.not_captured_dt_sq_hist);
	WriteBinaryArray(file, state.not_captured_v2dt_sq_hist);
	const uint64_t event_count = static_cast<uint64_t>(state.new_evaporation_events.size());
	WriteBinaryValue(file, event_count);
	for(const auto& entry : state.new_evaporation_events)
		WriteSnapshotEvaporationEntryBinary(file, entry);
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

bool ReadSnapshotRankState(const std::string& path, uint64_t expected_run_id, SnapshotRankState& state)
{
	std::ifstream file(path, std::ios::binary);
	if(!file.is_open())
		return false;

	ReadBinaryValue(file, state.run_id);
	ReadBinaryValue(file, state.snapshot_index);
	ReadBinaryValue(file, state.rank);
	ReadBinaryValue(file, state.done);
	ReadBinaryValue(file, state.trajectory_in_progress);
	ReadBinaryValue(file, state.local_captured);
	ReadBinaryValue(file, state.local_total);
	ReadBinaryValue(file, state.bincount_captured_samples);
	ReadBinaryValue(file, state.bincount_not_captured_samples);
	ReadBinaryValue(file, state.current_trajectory_id);
	ReadBinaryValue(file, state.rank_elapsed_wall_sec);
	ReadBinaryValue(file, state.current_trajectory_captured);
	ReadBinaryArray(file, state.current_trajectory_dt_hist);
	ReadBinaryArray(file, state.current_trajectory_v2dt_hist);
	ReadBinaryArray(file, state.captured_dt_hist);
	ReadBinaryArray(file, state.captured_v2dt_hist);
	ReadBinaryArray(file, state.captured_dt_sq_hist);
	ReadBinaryArray(file, state.captured_v2dt_sq_hist);
	ReadBinaryArray(file, state.not_captured_dt_hist);
	ReadBinaryArray(file, state.not_captured_v2dt_hist);
	ReadBinaryArray(file, state.not_captured_dt_sq_hist);
	ReadBinaryArray(file, state.not_captured_v2dt_sq_hist);

	uint64_t event_count = 0;
	ReadBinaryValue(file, event_count);
	state.new_evaporation_events.clear();
	state.new_evaporation_events.resize(static_cast<size_t>(event_count));
	for(size_t i = 0; i < state.new_evaporation_events.size(); i++)
		ReadSnapshotEvaporationEntryBinary(file, state.new_evaporation_events[i]);

	if(!file)
		return false;
	return state.run_id == expected_run_id;
}

void CleanupSnapshotCheckpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(SnapshotRankCheckpointPath(rank_snapshot_dir, rank, snapshot_index, interval_seconds).c_str());
}

void CleanupFinalSnapshotStates(const std::string& rank_snapshot_dir, int mpi_processes)
{
	for(int rank = 0; rank < mpi_processes; rank++)
		std::remove(SnapshotRankFinalPath(rank_snapshot_dir, rank).c_str());
	rmdir(rank_snapshot_dir.c_str());
}

SnapshotMergeResult TryWriteSnapshot(
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	bool allow_partial)
{
	SnapshotMergeResult merged_result;
	merged_result.status = SnapshotMergeStatus::Merged;
	if(SnapshotTextFileIsMerged(SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds))
	   && SnapshotTextFileIsMerged(SnapshotEvaporationTimeFilePath(snapshot_root, snapshot_index, interval_seconds)))
	{
		CleanupSnapshotCheckpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
		return merged_result;
	}

	SnapshotReportState report;
	SnapshotMergeResult result = LoadSnapshotReportState(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes, run_id, report);
	if(result.status == SnapshotMergeStatus::NoRanksReady)
		return result;
	if(result.status == SnapshotMergeStatus::Partial && !allow_partial)
		return result;

	if(!WriteSnapshotReportFile(snapshot_root, snapshot_index, interval_seconds, mass_gev, sigma_cm2, report, result))
	{
		result.status = SnapshotMergeStatus::NoRanksReady;
		return result;
	}

	if(result.status == SnapshotMergeStatus::Merged)
		CleanupSnapshotCheckpoints(rank_snapshot_dir, snapshot_index, interval_seconds, mpi_processes);
	return result;
}

bool WriteMissedSnapshotMarker(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	double actual_write_wall_time_sec)
{
	const std::string path = SnapshotTextFilePath(snapshot_root, snapshot_index, interval_seconds);
	if(SnapshotTextFileIsMerged(path))
		return true;

	return WriteTextFileAtomically(path, snapshot_index, [&](std::ofstream& file)
	{
		file << "# snapshot_status = missed_writer_delay\n";
		file << "# snapshot_target_wall_time_s = " << SnapshotTimeLabelSeconds(snapshot_index, interval_seconds) << "\n";
		file << "# snapshot_interval_s = " << std::fixed << std::setprecision(3) << interval_seconds << "\n";
		file << "# actual_write_wall_time_s = " << std::fixed << std::setprecision(3) << actual_write_wall_time_sec << "\n";
		file << "# note = heartbeat woke too late; no state rewind is available\n";
	});
}

}	// namespace DaMaSCUS_SUN

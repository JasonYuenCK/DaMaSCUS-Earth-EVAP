#ifndef __Snapshot_IO_hpp_
#define __Snapshot_IO_hpp_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "Data_Generation.hpp"

namespace DaMaSCUS_SUN
{

struct SnapshotEvaporationProgressEntry
{
	uint64_t trajectory_id = 0;
	double completion_wall_time_sec = 0.0;
	double lifetime_unbinding_sec = -1.0;
};

struct SnapshotRankedEvaporationEntry
{
	int rank = -1;
	SnapshotEvaporationProgressEntry entry;
};

struct SnapshotRankState
{
	uint64_t run_id = 0;
	int32_t snapshot_index = 0;
	int32_t rank = 0;
	int32_t done = 0;
	int32_t trajectory_in_progress = 0;
	uint64_t local_captured = 0;
	uint64_t local_total = 0;
	uint64_t local_classified = 0;
	uint64_t local_numerical_failures = 0;
	uint64_t bincount_captured_samples = 0;
	uint64_t bincount_not_captured_samples = 0;
	uint64_t current_trajectory_id = 0;
	double rank_elapsed_wall_sec = 0.0;
	double current_trajectory_wall_sec = 0.0;
	double current_trajectory_simulated_elapsed_sec = 0.0;
	uint64_t current_trajectory_scatterings = 0;
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
	std::vector<SnapshotEvaporationProgressEntry> new_evaporation_events;
};

enum class SnapshotMergeStatus
{
	NoRanksReady,
	Partial,
	Merged
};

struct SnapshotMergeResult
{
	SnapshotMergeStatus status = SnapshotMergeStatus::NoRanksReady;
	bool cleanup_succeeded = true;
	std::vector<int> ready_ranks;
	std::vector<int> missing_ranks;
};

long long SnapshotTimeLabelSeconds(int snapshot_index, double interval_seconds);
std::string SnapshotTextFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds);
std::string SnapshotEvaporationTimeFilePath(const std::string& snapshot_root, int snapshot_index, double interval_seconds);
std::string SnapshotRankCheckpointPath(const std::string& rank_snapshot_dir, int rank, int snapshot_index, double interval_seconds);
std::string SnapshotRankFinalPath(const std::string& rank_snapshot_dir, int rank);

SnapshotEvaporationProgressEntry MakeSnapshotEvaporationProgressEntry(const CompactEvaporationEvent& event);

bool WriteSnapshotRankState(const std::string& path, const SnapshotRankState& state);
bool ReadSnapshotRankState(const std::string& path, uint64_t expected_run_id, SnapshotRankState& state);
bool CleanupSnapshotCheckpoints(const std::string& rank_snapshot_dir, int snapshot_index, double interval_seconds, int mpi_processes);
bool CleanupFinalSnapshotStates(const std::string& rank_snapshot_dir, int mpi_processes);

SnapshotMergeResult TryWriteSnapshot(
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	int snapshot_index,
	double interval_seconds,
	int mpi_processes,
	uint64_t run_id,
	double mass_gev,
	double sigma_cm2,
	bool allow_partial);

bool WriteMissedSnapshotMarker(
	const std::string& snapshot_root,
	int snapshot_index,
	double interval_seconds,
	uint64_t run_id,
	double actual_write_wall_time_sec);

}	// namespace DaMaSCUS_SUN

#endif

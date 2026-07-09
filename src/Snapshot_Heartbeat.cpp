#include "Snapshot_Heartbeat.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace DaMaSCUS_SUN
{

SnapshotHeartbeat::SnapshotHeartbeat(
	SnapshotSharedState& shared_state,
	int mpi_rank,
	int mpi_processes,
	uint64_t run_id,
	const std::string& snapshot_root,
	const std::string& rank_snapshot_dir,
	double snapshot_interval,
	double mass_gev,
	double sigma_cm2)
: shared_state_(shared_state),
  mpi_rank_(mpi_rank),
  mpi_processes_(mpi_processes),
  run_id_(run_id),
  snapshot_root_(snapshot_root),
  rank_snapshot_dir_(rank_snapshot_dir),
  snapshot_interval_(snapshot_interval),
  mass_gev_(mass_gev),
  sigma_cm2_(sigma_cm2),
  missed_tolerance_sec_(std::max(1.0, 0.1 * snapshot_interval))
{
}

SnapshotHeartbeat::~SnapshotHeartbeat()
{
	Stop();
}

void SnapshotHeartbeat::Start()
{
	std::lock_guard<std::mutex> lock(control_mutex_);
	if(started_)
		return;
	start_time_ = std::chrono::steady_clock::now();
	stop_requested_ = false;
	started_ = true;
	worker_ = std::thread(&SnapshotHeartbeat::ThreadMain, this);
}

void SnapshotHeartbeat::MarkDoneAndWriteFinal(double computing_time_sec)
{
	Stop();
	SnapshotRankState final_state = shared_state_.CopyFinal(computing_time_sec, first_uncommitted_evaporation_entry_);
	final_state.snapshot_index = last_rank_snapshot_written_;
	if(!WriteSnapshotRankState(SnapshotRankFinalPath(rank_snapshot_dir_, mpi_rank_), final_state))
	{
		std::cerr << "Warning in SnapshotHeartbeat::MarkDoneAndWriteFinal(): rank "
		          << mpi_rank_ << " failed to write final snapshot state." << std::endl;
	}
}

void SnapshotHeartbeat::FinalizeAfterAllRanksDone(double final_elapsed_sec)
{
	const int final_snapshot_index = std::max(0, static_cast<int>(std::floor(final_elapsed_sec / snapshot_interval_)));
	if(final_snapshot_index > 0)
		TryMergeUpTo(final_snapshot_index, false);

	for(int snapshot_index = 1; snapshot_index <= final_snapshot_index; snapshot_index++)
	{
		TryWriteSnapshot(
			snapshot_root_,
			rank_snapshot_dir_,
			snapshot_index,
			snapshot_interval_,
			mpi_processes_,
			run_id_,
			mass_gev_,
			sigma_cm2_,
			false);
	}
	CleanupFinalSnapshotStates(rank_snapshot_dir_, mpi_processes_);
}

void SnapshotHeartbeat::Stop()
{
	{
		std::lock_guard<std::mutex> lock(control_mutex_);
		if(!started_)
			return;
		stop_requested_ = true;
	}
	cv_.notify_all();
	if(worker_.joinable())
		worker_.join();
	std::lock_guard<std::mutex> lock(control_mutex_);
	started_ = false;
}

void SnapshotHeartbeat::ThreadMain()
{
	int next_index = 1;
	while(true)
	{
		const auto deadline = start_time_ + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
			std::chrono::duration<double>(next_index * snapshot_interval_));

		{
			std::unique_lock<std::mutex> lock(control_mutex_);
			if(cv_.wait_until(lock, deadline, [&] { return stop_requested_; }))
				break;
		}

		const double target_elapsed = next_index * snapshot_interval_;
		const double actual_elapsed = ElapsedSinceStart();
		if(actual_elapsed > target_elapsed + missed_tolerance_sec_)
		{
			if(mpi_rank_ == 0)
				WriteMissedSnapshotMarker(snapshot_root_, next_index, snapshot_interval_, actual_elapsed);
		}
		else
		{
			WriteRankCheckpoint(next_index, target_elapsed, actual_elapsed);
		}

		if(mpi_rank_ == 0)
			TryMergeUpTo(next_index, true);
		next_index++;
	}
}

void SnapshotHeartbeat::WriteRankCheckpoint(int snapshot_index, double target_wall_sec, double actual_elapsed_sec)
{
	size_t evaporation_entry_end = first_uncommitted_evaporation_entry_;
	SnapshotRankState state = shared_state_.CopyForSnapshot(
		snapshot_index,
		actual_elapsed_sec,
		target_wall_sec,
		first_uncommitted_evaporation_entry_,
		evaporation_entry_end);
	state.snapshot_index = snapshot_index;

	const std::string checkpoint_path = SnapshotRankCheckpointPath(rank_snapshot_dir_, mpi_rank_, snapshot_index, snapshot_interval_);
	if(WriteSnapshotRankState(checkpoint_path, state))
	{
		first_uncommitted_evaporation_entry_ = evaporation_entry_end;
		last_rank_snapshot_written_ = snapshot_index;
	}
	else
	{
		std::cerr << "Warning in SnapshotHeartbeat::WriteRankCheckpoint(): rank "
		          << mpi_rank_ << " failed to write " << checkpoint_path << std::endl;
	}
}

void SnapshotHeartbeat::TryMergeUpTo(int snapshot_index, bool allow_partial)
{
	for(int index = 1; index <= snapshot_index; index++)
	{
		SnapshotMergeResult result = TryWriteSnapshot(
			snapshot_root_,
			rank_snapshot_dir_,
			index,
			snapshot_interval_,
			mpi_processes_,
			run_id_,
			mass_gev_,
			sigma_cm2_,
			allow_partial);
		if(result.status == SnapshotMergeStatus::Merged)
			last_merge_attempted_ = std::max(last_merge_attempted_, index);
		else if(result.status == SnapshotMergeStatus::Partial && !allow_partial)
			break;
	}
}

double SnapshotHeartbeat::ElapsedSinceStart() const
{
	return 1.0e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(
		std::chrono::steady_clock::now() - start_time_).count();
}

}	// namespace DaMaSCUS_SUN

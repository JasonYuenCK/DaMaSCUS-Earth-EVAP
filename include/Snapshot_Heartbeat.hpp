#ifndef __Snapshot_Heartbeat_hpp_
#define __Snapshot_Heartbeat_hpp_

#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <string>
#include <thread>

#include "Snapshot_Shared_State.hpp"

namespace DaMaSCUS_SUN
{

class SnapshotHeartbeat
{
  public:
	SnapshotHeartbeat(
		SnapshotSharedState& shared_state,
		int mpi_rank,
		int mpi_processes,
		uint64_t run_id,
		const std::string& snapshot_root,
		const std::string& rank_snapshot_dir,
		double snapshot_interval,
		double mass_gev,
		double sigma_cm2);
	~SnapshotHeartbeat();

	void Start();
	void MarkDoneAndWriteFinal(double computing_time_sec);
	void FinalizeAfterAllRanksDone(double final_elapsed_sec);
	void Stop();

  private:
	void ThreadMain();
	void WriteRankCheckpoint(int snapshot_index, double target_wall_sec, double actual_elapsed_sec);
	void TryMergeUpTo(int snapshot_index, bool allow_partial);
	double ElapsedSinceStart() const;

	SnapshotSharedState& shared_state_;
	int mpi_rank_;
	int mpi_processes_;
	uint64_t run_id_;
	std::string snapshot_root_;
	std::string rank_snapshot_dir_;
	double snapshot_interval_;
	double mass_gev_;
	double sigma_cm2_;
	double missed_tolerance_sec_;

	std::thread worker_;
	std::mutex control_mutex_;
	std::condition_variable cv_;
	bool started_ = false;
	bool stop_requested_ = false;
	std::chrono::steady_clock::time_point start_time_;

	size_t first_uncommitted_evaporation_entry_ = 0;
	int last_rank_snapshot_written_ = 0;
	int last_merge_attempted_ = 0;
};

}	// namespace DaMaSCUS_SUN

#endif

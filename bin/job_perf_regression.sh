#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=/lustre/project/kennyng/backup_DM/bin/job-%x-%j.out
#SBATCH --error=/lustre/project/kennyng/backup_DM/bin/job-%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lingyuxia@link.cuhk.edu.hk

set -euo pipefail

PROJECT_ROOT=/lustre/project/kennyng/backup_DM
SCRIPT_DIR="$PROJECT_ROOT/bin"
VARIANT="${1:?Usage: job_perf_regression.sh <branch-variant> <config-file> <repetitions>}"
CONFIG_FILE="${2:?Usage: job_perf_regression.sh <branch-variant> <config-file> <repetitions>}"
REPETITIONS="${3:?Usage: job_perf_regression.sh <branch-variant> <config-file> <repetitions>}"

if ! [[ "$REPETITIONS" =~ ^[1-9][0-9]*$ ]]; then
	echo "Repetitions must be a positive integer: $REPETITIONS" >&2
	exit 2
fi

if [[ -n "${BENCHMARK_WORKTREE_OVERRIDE:-}" ]]; then
	WORKTREE="$BENCHMARK_WORKTREE_OVERRIDE"
elif [[ "$VARIANT" == "main" ]]; then
	WORKTREE="$PROJECT_ROOT"
else
	WORKTREE="$PROJECT_ROOT/damascus_worktrees/$VARIANT"
fi

if [[ "$CONFIG_FILE" != /* ]]; then
	CONFIG_FILE="$SCRIPT_DIR/$CONFIG_FILE"
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
	echo "Configuration file not found: $CONFIG_FILE" >&2
	exit 2
fi
if [[ ! -d "$WORKTREE" ]]; then
	echo "Worktree not found: $WORKTREE" >&2
	exit 2
fi
if [[ -n "$(git -C "$WORKTREE" status --porcelain)" ]]; then
	echo "Worktree is dirty; refusing to run an untraceable benchmark: $WORKTREE" >&2
	exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
	echo "python3 is required by benchmark_metrics.py." >&2
	exit 2
fi

module load intel/2022.2 boost/1.77.0

CMAKE_BIN="$PROJECT_ROOT/damascus_worktrees/tools/cmake-3.29.6/bin/cmake"
LIBCONFIG_PREFIX="$PROJECT_ROOT/damascus_worktrees/tools/libconfig-1.7.3"
if [[ ! -x "$CMAKE_BIN" ]]; then
	echo "CMake not found: $CMAKE_BIN" >&2
	exit 2
fi
if [[ ! -f "$LIBCONFIG_PREFIX/lib64/libconfig++.so" ]]; then
	echo "libconfig++ not found under: $LIBCONFIG_PREFIX" >&2
	exit 2
fi
export LD_LIBRARY_PATH="$LIBCONFIG_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

GNU_TIME_BIN="${GNU_TIME_BIN:-}"
if [[ -z "$GNU_TIME_BIN" && -x /usr/bin/time ]] && /usr/bin/time --version 2>&1 | grep -q "GNU time"; then
	GNU_TIME_BIN=/usr/bin/time
elif [[ -z "$GNU_TIME_BIN" ]] && command -v gtime >/dev/null 2>&1 && gtime --version 2>&1 | grep -q "GNU time"; then
	GNU_TIME_BIN="$(command -v gtime)"
fi
TIME_BACKEND="${GNU_TIME_BIN:-date-nanoseconds}"

COMMIT="$(git -C "$WORKTREE" rev-parse HEAD)"
HOST_NAME="$(hostname)"
CASE_ID="$(basename "$CONFIG_FILE" .cfg)"
SAFE_VARIANT="${VARIANT//\//_}"
BENCHMARK_ROOT="${BENCHMARK_ROOT:-$PROJECT_ROOT/bin/benchmark-results}"
JOB_ROOT="$BENCHMARK_ROOT/job-${SLURM_JOB_ID}-${SAFE_VARIANT}-${CASE_ID}-r${SLURM_NTASKS}"
BUILD_DIR="$JOB_ROOT/build"
METRICS_FILE="$JOB_ROOT/metrics.tsv"
mkdir -p "$JOB_ROOT"

{
	echo "variant=$VARIANT"
	echo "commit=$COMMIT"
	echo "worktree=$WORKTREE"
	echo "source_config=$CONFIG_FILE"
	echo "ranks=$SLURM_NTASKS"
	echo "repetitions=$REPETITIONS"
	echo "host=$HOST_NAME"
	echo "slurm_job_id=$SLURM_JOB_ID"
	echo "start_time=$(date --iso-8601=seconds)"
	echo "time_backend=$TIME_BACKEND"
} > "$JOB_ROOT/manifest.txt"
cp "$CONFIG_FILE" "$JOB_ROOT/source_config.cfg"

echo "BENCHMARK_ROOT=$JOB_ROOT"
echo "VARIANT=$VARIANT"
echo "COMMIT=$COMMIT"
echo "CONFIG=$CONFIG_FILE"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "REPETITIONS=$REPETITIONS"

build_start=$SECONDS
"$CMAKE_BIN" -S "$WORKTREE" -B "$BUILD_DIR" \
	-DCMAKE_BUILD_TYPE=Release \
	-DBUILD_TESTING=OFF \
	-DCODE_COVERAGE=OFF \
	-DLIBCONFIG_INCLUDE_DIRs="$LIBCONFIG_PREFIX/include" \
	-DLIBCONFIGPP_LIBRARY="$LIBCONFIG_PREFIX/lib64/libconfig++.so"
"$CMAKE_BIN" --build "$BUILD_DIR" --parallel "$SLURM_NTASKS"
BUILD_SECONDS=$((SECONDS - build_start))
echo "BUILD_SECONDS=$BUILD_SECONDS"
echo "build_seconds=$BUILD_SECONDS" >> "$JOB_ROOT/manifest.txt"

BINARY="$BUILD_DIR/src/DaMaSCUS-SUN"
if [[ ! -x "$BINARY" ]]; then
	echo "DaMaSCUS-SUN executable not found: $BINARY" >&2
	exit 2
fi

overall_status=0
for ((repeat = 1; repeat <= REPETITIONS; repeat++)); do
	RUN_DIR="$JOB_ROOT/run-$repeat"
	OUTPUT_DIR="$RUN_DIR/output/"
	RUN_CONFIG="$RUN_DIR/config.cfg"
	RUN_LOG="$RUN_DIR/run.log"
	TIME_LOG="$RUN_DIR/time.txt"
	mkdir -p "$OUTPUT_DIR"

	# Preserve the source configuration byte-for-byte except for the two fields
	# that must be isolated to prevent concurrent/repeated runs from colliding.
	awk -v output_dir="$OUTPUT_DIR" -v run_id="${CASE_ID}_${SAFE_VARIANT}_r${SLURM_NTASKS}_repeat${repeat}" '
		BEGIN { found_output = 0; found_id = 0 }
		/^[[:space:]]*output_dir[[:space:]]*=/ {
			printf "\toutput_dir\t=\t\"%s\";\n", output_dir
			found_output = 1
			next
		}
		/^[[:space:]]*ID[[:space:]]*=/ {
			printf "\tID\t\t=\t\"%s\";\n", run_id
			found_id = 1
			next
		}
		{ print }
		END { if(!found_output || !found_id) exit 3 }
	' "$CONFIG_FILE" > "$RUN_CONFIG"

	echo "RUN_REPEAT=$repeat"
	set +e
	if [[ -n "$GNU_TIME_BIN" ]]; then
		"$GNU_TIME_BIN" -v -o "$TIME_LOG" \
			mpirun -n "$SLURM_NTASKS" "$BINARY" "$RUN_CONFIG" > "$RUN_LOG" 2>&1
		run_status=$?
	else
		run_start_ns="$(date +%s%N)"
		mpirun -n "$SLURM_NTASKS" "$BINARY" "$RUN_CONFIG" > "$RUN_LOG" 2>&1
		run_status=$?
		run_end_ns="$(date +%s%N)"
		awk -v start="$run_start_ns" -v end="$run_end_ns" \
			'BEGIN { printf "Elapsed (wall clock) time (h:mm:ss or m:ss): %.6f\n", (end - start) / 1000000000.0 }' \
			> "$TIME_LOG"
	fi
	set -e

	python3 "$SCRIPT_DIR/benchmark_metrics.py" extract \
		--log "$RUN_LOG" \
		--time "$TIME_LOG" \
		--output "$METRICS_FILE" \
		--variant "$VARIANT" \
		--case "$CASE_ID" \
		--ranks "$SLURM_NTASKS" \
		--repeat "$repeat" \
		--commit "$COMMIT" \
		--host "$HOST_NAME" \
		--run-status "$run_status"

	echo "RUN_STATUS=$run_status"
	tail -n 1 "$METRICS_FILE"
	if [[ "$run_status" -ne 0 ]]; then
		overall_status=1
	fi
done

echo "end_time=$(date --iso-8601=seconds)" >> "$JOB_ROOT/manifest.txt"
echo "METRICS_FILE=$METRICS_FILE"
exit "$overall_status"

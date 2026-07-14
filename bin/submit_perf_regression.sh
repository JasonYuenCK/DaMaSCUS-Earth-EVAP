#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	cat <<'EOF'
Usage:
  ./submit_perf_regression.sh [options] <config.cfg> [config.cfg ...]

Options:
  --variants LIST      Comma-separated worktree names (default: perf-initial-validation)
  --ranks LIST         Comma-separated MPI ranks (default: 1,4,16,32)
  --repetitions N      Repeats per branch/config/rank group (default: 5)
  --nodelist NODE      Pin every job to one node for branch comparisons
  --exclusive          Request exclusive use of the allocated node
  --parallel           Submit independent jobs without dependencies
  --help               Show this message

By default every job is chained with an afterany dependency. This prevents
different benchmark cases from competing for the same node and still lets the
suite continue if one case fails.
EOF
}

VARIANTS="${BENCHMARK_VARIANTS:-perf-initial-validation}"
RANKS="${BENCHMARK_RANKS:-1,4,16,32}"
REPETITIONS="${BENCHMARK_REPETITIONS:-5}"
NODELIST="${BENCHMARK_NODELIST:-}"
EXCLUSIVE="${BENCHMARK_EXCLUSIVE:-0}"
SEQUENTIAL=1

while [[ $# -gt 0 ]]; do
	case "$1" in
		--variants)
			VARIANTS="${2:?--variants requires a value}"
			shift 2
			;;
		--ranks)
			RANKS="${2:?--ranks requires a value}"
			shift 2
			;;
		--repetitions)
			REPETITIONS="${2:?--repetitions requires a value}"
			shift 2
			;;
		--nodelist)
			NODELIST="${2:?--nodelist requires a value}"
			shift 2
			;;
		--exclusive)
			EXCLUSIVE=1
			shift
			;;
		--parallel)
			SEQUENTIAL=0
			shift
			;;
		--help|-h)
			usage
			exit 0
			;;
		--*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 2
			;;
		*)
			break
			;;
	esac
done

if [[ $# -lt 1 ]]; then
	usage >&2
	exit 2
fi
if ! [[ "$REPETITIONS" =~ ^[1-9][0-9]*$ ]]; then
	echo "Repetitions must be a positive integer: $REPETITIONS" >&2
	exit 2
fi

IFS=',' read -r -a variant_list <<< "$VARIANTS"
IFS=',' read -r -a rank_list <<< "$RANKS"
config_list=("$@")

for config_index in "${!config_list[@]}"; do
	config="${config_list[$config_index]}"
	if [[ "$config" != /* ]]; then
		config="$SCRIPT_DIR/$config"
	fi
	if [[ ! -f "$config" ]]; then
		echo "Configuration file not found: $config" >&2
		exit 2
	fi
	config_list[$config_index]="$config"
done

previous_job_id=""
submitted_job_ids=()
for variant in "${variant_list[@]}"; do
	if [[ -z "$variant" ]]; then
		echo "Variant list contains an empty entry: $VARIANTS" >&2
		exit 2
	fi
	for config in "${config_list[@]}"; do
		case_id="$(basename "$config" .cfg)"
		for ranks in "${rank_list[@]}"; do
			if ! [[ "$ranks" =~ ^[1-9][0-9]*$ ]]; then
				echo "MPI rank count must be a positive integer: $ranks" >&2
				exit 2
			fi
			safe_variant="${variant//\//_}"
			job_name="reg-${safe_variant:0:12}-${case_id:0:12}-r${ranks}"
			submit_args=(
				--parsable
				--ntasks="$ranks"
				--job-name="$job_name"
			)
			if [[ "$SEQUENTIAL" -eq 1 && -n "$previous_job_id" ]]; then
				submit_args+=(--dependency="afterany:$previous_job_id")
			fi
			if [[ -n "$NODELIST" ]]; then
				submit_args+=(--nodelist="$NODELIST")
			fi
			if [[ "$EXCLUSIVE" -eq 1 ]]; then
				submit_args+=(--exclusive)
			fi

			job_id="$(sbatch "${submit_args[@]}" "$SCRIPT_DIR/job_perf_regression.sh" "$variant" "$config" "$REPETITIONS")"
			job_id="${job_id%%;*}"
			submitted_job_ids+=("$job_id")
			if [[ "$SEQUENTIAL" -eq 1 ]]; then
				previous_job_id="$job_id"
			fi
			echo "$variant | $case_id | $ranks ranks | $REPETITIONS repeats: submitted $job_id"
		done
	done
done

printf 'Submitted job IDs:'
printf ' %s' "${submitted_job_ids[@]}"
printf '\n'

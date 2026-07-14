# Performance and physics regression suite

This suite builds each selected branch once per Slurm job, then performs repeated timed runs on one node. Build time is recorded separately and is never included in the per-run metrics.

## Submit the low-cross-section matrix

Use the two production configurations already stored in `bin/` on CHPC:

```bash
./submit_perf_regression.sh \
  --variants perf-initial-validation \
  --ranks 1,4,16,32 \
  --repetitions 5 \
  config_perf_m0p1_s40_n10000_r32.cfg \
  config_perf_m0p5_s40_n10000_r32.cfg
```

Jobs are sequential by default so benchmark cases do not contend with each other. Use `--parallel` only when node-to-node throughput is more important than a controlled comparison.

For a branch regression comparison, provide matching configuration files and both worktrees:

```bash
./submit_perf_regression.sh \
  --variants main,perf-initial-validation \
  --ranks 4,32 \
  --repetitions 5 \
  --nodelist chpc-cn066 \
  --exclusive \
  config_benchmark_main.cfg
```

When comparing branches, the physics values in each supplied configuration must be identical. Only `output_dir` and `ID` are rewritten at runtime, into a unique directory per repeat, so concurrent runs cannot overwrite each other. `--exclusive` may be expensive on a large node; use it for strict performance gates, not routine physics checks.

## Artifacts and comparison

Each job writes under `bin/benchmark-results/job-<id>-.../`:

- `manifest.txt`: branch, full commit, source configuration, host, rank count and build time;
- `source_config.cfg`: unmodified input configuration;
- `run-N/config.cfg`: isolated runtime configuration;
- `run-N/run.log` and `run-N/time.txt`: complete program and wall-time/resource output;
- `metrics.tsv`: one normalized row per repeat.

Compare completed baseline and candidate metrics with:

```bash
python3 benchmark_metrics.py compare \
  --baseline main \
  --candidate perf-initial-validation \
  benchmark-results/job-*/metrics.tsv
```

The default gate fails when the two variants ran on different hosts, when the capture-rate difference exceeds 3 standard deviations, or when the candidate numerical-failure rate exceeds `1e-4`. Runtime repetitions use the same seed and are therefore treated as timing repetitions, not independent Monte Carlo samples. An optional performance gate can be added with `--min-speedup 1.05`; use `--allow-host-mismatch` only for physics-only comparisons.

The metrics also include MPI synchronization rounds, last-round trajectory count and capture overshoot. These fields are the evidence used to decide whether the current fixed low-cross-section batch size should be replaced by an adaptive policy.

When GNU time is unavailable, the job automatically falls back to nanosecond wall-clock timing. In that mode `max_rss_kb` is intentionally left empty instead of failing the benchmark.

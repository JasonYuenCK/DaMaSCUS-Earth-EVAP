# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun — Evaporation Branch

## Overview

This project extends [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) to study dark matter (DM) capture and evaporation in the Sun. It simulates DM particle trajectories inside the Sun using adaptive RK45 integration, tracks gravitational capture via energy monitoring, and computes time-weighted radial distribution histograms (bincount) on-the-fly.

## Key Features

- **Online bincount accumulation**: Replaces trajectory file output with in-memory histogram accumulation (r-histogram and v²-histogram, 2000 bins from 0 to 2R☉)
- **Capture detection**: Marks first capture after a scattering event leaves the particle with negative total energy; later free-propagation checks update the last bound time
- **Evaporation survival records**: Records every captured trajectory with observed/censored lifetime fields and free-flight energy-drift diagnostics
- **MPI parallelization**: Each rank independently targets `ceil(sample_size / N_ranks)` captured particles
- **Safety valve**: `max_trajectories` parameter prevents infinite runtime for low capture-rate scenarios

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The executable is placed in `bin/DaMaSCUS-SUN`.

### Dependencies

- C++17 compiler
- MPI (OpenMPI or MPICH)
- CMake ≥ 3.14
- libconfig++

External dependencies (`libphysica`, `obscura`) are fetched automatically by CMake.

## Usage

```bash
mpirun -np 8 bin/DaMaSCUS-SUN bin/config.cfg
```

### Configuration File

Key parameters in the `.cfg` file:

| Parameter | Description |
|-----------|-------------|
| `sample_size` | Target number of captured particles |
| `max_trajectories` | Maximum total trajectories before forced stop (default: sample_size × 1000) |
| `output_dir` | Directory for output files |
| `DM_mass` | Dark matter mass in GeV |
| `DM_cross_section_nucleon` | DM-nucleon cross section in cm² |
| `run_mode` | "Parameter point" or "Parameter scan" |
| `evaporation_mode_bincount_enabled` | Optional split of captured bincounts by evaporation-time peaks |
| `snapshot_evaporation_log_enabled` | Optional snapshot append log for evaporation/survival deltas; defaults to `true` when `snapshot_enabled = true` |
| `evaporation_diagnostics_enabled` | Optional full evaporation survival/diagnostic output; default `false` |

Optional evaporation-mode settings:

```cfg
evaporation_mode_bincount_enabled = true;
evaporation_mode_boundaries_log10_s = (4.5, 11.1);
evaporation_mode_labels = ("P1_fast", "P2_theory", "P3_tail");
evaporation_mode_include_truncated = false;
```

`evaporation_mode_include_truncated` is retained for config compatibility, but mode assignment uses only valid observed unbinding events.

For reproducible MPI runs, a nonzero fixed seed is treated as a base seed; each rank uses `base_seed + 1000003 * mpi_rank` to avoid duplicate trajectories. Computational cutoffs (`wall_time_limit`, `max_free_steps`, `max_scatterings`) are marked invalid for survival analysis rather than normal right-censoring. Before production evaporation runs, compare `interpolation_points = 0`, `1000`, and `2000` for median lifetime, tail fraction, complete-event fraction, and mean scatterings.

### Output Files

For each parameter point, the main generated files are:

- `bincount.txt` — Combined captured/not-captured time-weighted radial histogram with error estimates; incomplete non-captured trajectories are excluded from the not-captured histogram
- `evaporation_times.txt` — Without snapshot logging, the default final evaporation-time sequence is a compact table of complete valid evaporation events: `rank trajectory_id t_evap_s`. With `snapshot_evaporation_log_enabled = true` or `evaporation_diagnostics_enabled = true`, it is an append log with ordered snapshot/final blocks and right-censoring flags.
- `evaporation_summary.txt` — Only written when `evaporation_diagnostics_enabled = true`; one full survival-analysis record for every captured trajectory, including right-censoring and numerical diagnostics.
- `evaporation_mode_summary.txt` — Counts and log10(lifetime_unbinding/s) boundaries for each configured evaporation mode, using only valid observed unbinding events, when enabled
- `evaporation_mode_bincount.txt` — Per-mode captured radial bincounts (`dt`, `v2dt`, and errors), when enabled
- `computation_time_summary.txt` — Wall-clock time, RK45 step-count statistics, captured/complete/censored/invalid-survival counts, and trajectory termination-reason counts

When `snapshot_enabled = true`, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt` — Cumulative snapshot report and bincount histogram. The rank diagnostic table combines the checkpoint source with the current `running`/`done` status.

Snapshot evaporation deltas are appended only to the single `evaporation_times.txt` file. Full 23+ column evaporation diagnostics are still controlled separately by `evaporation_diagnostics_enabled` and written once to `evaporation_summary.txt`.

## References

- T. Emken, C. Kouvaris, [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani & Palomares-Ruiz (2017), Evaporation of dark matter in the Sun

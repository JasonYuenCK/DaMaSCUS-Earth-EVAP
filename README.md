# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun — Evaporation Branch

## Overview

This project extends [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) to study dark matter (DM) capture and evaporation in the Sun. It simulates DM particle trajectories inside the Sun using adaptive RK45 integration, tracks gravitational capture via energy monitoring, and computes time-weighted radial distribution histograms (bincount) on-the-fly.

## Key Features

- **Online bincount accumulation**: Replaces trajectory file output with in-memory histogram accumulation (r-histogram and v²-histogram, 2000 bins from 0 to 2R☉)
- **Capture detection**: Marks first capture after a scattering event leaves the particle with negative total energy; later free-propagation checks update the last bound time
- **Evaporation survival records**: Records every captured trajectory with observed/censored lifetime fields and free-flight energy-drift diagnostics
- **MPI parallelization**: Each rank independently targets `ceil(sample_size / N_ranks)` captured particles
- **Optional safety valve**: set `max_trajectories` only when an explicit hard cap is required; unset means no trajectory-count limit

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
| `max_trajectories` | Optional maximum total trajectories before forced stop (default: no limit) |
| `output_dir` | Directory for output files |
| `DM_mass` | Dark matter mass in GeV |
| `DM_cross_section_nucleon` | DM-nucleon cross section in cm² |
| `run_mode` | "Parameter point" or "Parameter scan" |
| `snapshot_evaporation_log_enabled` | Optional snapshot append log for evaporation/survival deltas; defaults to `true` when `snapshot_enabled = true` |
| `evaporation_diagnostics_enabled` | Optional full evaporation survival/diagnostic output; default `false` |

For reproducible MPI runs, a nonzero fixed seed is treated as a base seed; each rank uses `base_seed + 1000003 * mpi_rank` to avoid duplicate trajectories. Computational cutoffs (`wall_time_limit`, `max_free_steps`, `max_scatterings`) are marked invalid for survival analysis rather than normal right-censoring. Before production evaporation runs, compare `interpolation_points = 0`, `1000`, and `2000` for median lifetime, tail fraction, complete-event fraction, and mean scatterings.

### Output Files

For each parameter point, the main generated files are:

- `bincount.txt` — Combined captured/not-captured time-weighted radial histogram with error estimates; incomplete non-captured trajectories are excluded from the not-captured histogram
- `evaporation_times.txt` — Compact table of complete valid evaporation events only: `rank trajectory_id lifetime_unbinding_sec`.
- `evaporation_diagnostics.txt` — Only written when `evaporation_diagnostics_enabled = true`; one full survival-analysis record for every captured trajectory, including right-censoring and numerical diagnostics.

When `snapshot_enabled = true`, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt` — Cumulative snapshot report and bincount histogram.
- `snapshot_{time}s_evaporation_times.txt` — Complete, valid evaporation events newly finished in that snapshot interval.

At each merged snapshot, `evaporation_times.txt` is atomically refreshed with all complete, valid evaporation events observed so far. Snapshot reports also record the minimum and maximum rank checkpoint wall times; these describe an asynchronous checkpoint range rather than an exact globally synchronized instant. Full evaporation diagnostics are controlled separately by `evaporation_diagnostics_enabled` and written once to `evaporation_diagnostics.txt`.

## References

- T. Emken, C. Kouvaris, [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani & Palomares-Ruiz (2017), Evaporation of dark matter in the Sun

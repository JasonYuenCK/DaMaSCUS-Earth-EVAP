# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun — Evaporation Branch

## Overview

This project extends [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) to study dark matter (DM) capture and evaporation in the Sun. It simulates DM particle trajectories inside the Sun using adaptive RK45 integration, tracks gravitational capture via energy monitoring, and computes time-weighted radial distribution histograms (bincount) on-the-fly.

## Key Features

- **Online bincount accumulation**: Replaces trajectory file output with in-memory histogram accumulation (r-histogram and v²-histogram, 2000 bins from 0 to 2R☉)
- **Capture detection**: Marks first capture after a scattering event leaves the particle with negative total energy; later free-propagation checks update the last bound time
- **Evaporation time**: Records positive gravitationally bound durations for captured particles, with truncation flags for censored data
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

Optional evaporation-mode settings:

```cfg
evaporation_mode_bincount_enabled = true;
evaporation_mode_boundaries_log10_s = (4.5, 11.1);
evaporation_mode_labels = ("P1_fast", "P2_theory", "P3_tail");
evaporation_mode_include_truncated = false;
```

### Output Files

For each parameter point, the main generated files are:

- `bincount.txt` — Combined captured/not-captured time-weighted radial histogram with error estimates; incomplete non-captured trajectories are excluded from the not-captured histogram
- `evaporation_summary.txt` — Positive evaporation durations, keyed by `rank` and rank-local `trajectory_id`, with the post-scatter first-negative radius, energy, previous-step energy difference, truncation flag, and termination reason; zero-duration captures are treated as non-evaporation and are not included in the statistics
- `evaporation_mode_summary.txt` — Counts and log10(t_evap/s) boundaries for each configured evaporation mode, when enabled
- `evaporation_mode_bincount.txt` — Per-mode captured radial bincounts (`dt`, `v2dt`, and errors), when enabled
- `computation_time_summary.txt` — Wall-clock time, RK45 step-count statistics, and trajectory termination-reason counts

When `snapshot_enabled = true`, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt` — Cumulative snapshot report and bincount histogram. The rank diagnostic table combines the checkpoint source with the current `running`/`done` status.
- `snapshot_{time}s_evaporation.txt` — Rank-aware positive evaporation-duration list for captured trajectories completed by that snapshot time, with the same first-negative and termination diagnostics as `evaporation_summary.txt`.

## References

- T. Emken, C. Kouvaris, [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani & Palomares-Ruiz (2017), Evaporation of dark matter in the Sun

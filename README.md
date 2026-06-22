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
- `evaporation_summary.txt` — One survival-analysis record for every captured trajectory, keyed by `rank` and rank-local `trajectory_id`; includes capture time, final unbinding-scatter time, boundary escape time, termination time, observed/censored flags, unbinding and boundary lifetimes, termination reason, and maximum free-flight energy drift diagnostics
- `evaporation_mode_summary.txt` — Counts and log10(observed_lifetime/s) boundaries for each configured evaporation mode, when enabled
- `evaporation_mode_bincount.txt` — Per-mode captured radial bincounts (`dt`, `v2dt`, and errors), when enabled
- `computation_time_summary.txt` — Wall-clock time, RK45 step-count statistics, captured/complete/censored counts, and trajectory termination-reason counts

When `snapshot_enabled = true`, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt` — Cumulative snapshot report and bincount histogram. The rank diagnostic table combines the checkpoint source with the current `running`/`done` status.
- `snapshot_{time}s_evaporation.txt` — Rank-aware survival-analysis records for captured trajectories completed by that snapshot time, with the same lifetime, censoring, termination, and free-flight energy-drift diagnostics as `evaporation_summary.txt`.

## References

- T. Emken, C. Kouvaris, [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani & Palomares-Ruiz (2017), Evaporation of dark matter in the Sun

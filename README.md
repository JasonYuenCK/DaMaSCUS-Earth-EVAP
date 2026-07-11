# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun, with capture- and evaporation-focused
extensions.

## Overview

DaMaSCUS-SUN-EVAP builds on
[DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) and keeps its central
Monte Carlo picture: dark matter particles are propagated through the solar
potential, scatter on solar targets, and are classified from their trajectory
history. This branch is organized around production workflows for low-mass dark
matter capture and evaporation studies rather than broad direct-detection scans.

The current code path is centered on two practical modes:

- **Capture mode**: a fast capture-rate workflow. It terminates a trajectory once
  a post-scatter bound state is identified and avoids the full evaporation and
  histogram output path.
- **Parameter-point simulation**: the main evaporation workflow for one mass and
  cross section. It accumulates time-weighted radial histograms, records complete
  valid evaporation events, and can emit wall-clock snapshot progress files for
  long MPI runs.

The older parameter-scan machinery is still present, but the most actively
maintained outputs in this branch are the capture summary and the single
parameter-point evaporation products.

## Main Changes

Compared with the upstream DaMaSCUS-SUN workflow, this branch emphasizes:

- in-memory radial bincount accumulation instead of writing full trajectory
  files;
- capture detection from the first post-scatter negative-energy state;
- compact final evaporation-time output for complete valid unbinding events;
- optional richer survival diagnostics behind explicit diagnostic paths;
- MPI-aware snapshot output for long parameter-point jobs;
- safeguards for pathological trajectories and low-capture-rate runs;
- server-friendly local configuration conventions, with generated binaries,
  job scripts, and run configs kept outside version control under `bin/`.

## Build And Deployment

### Dependencies

- CMake 3.12 or newer
- C++11-capable compiler
- MPI implementation such as OpenMPI or MPICH
- Boost
- libconfig++

`obscura` and `libphysica` are fetched by CMake under `external/`.

### Local Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF
cmake --build build --config Release --parallel
cmake --install build --config Release
```

The installed executable is expected at:

```bash
bin/DaMaSCUS-SUN
```

The repository ignores `bin/`, so local run scripts, cluster submission scripts,
and private configuration files can live there without being committed.

The installed binary is currently checkout-bound: its build records the source
directory used to locate `data/model_agss09.dat`. Keep the checkout in place and
rebuild after moving it; copying the executable alone is not supported.

### Cluster Deployment

A typical cluster checkout follows the same pattern:

```bash
git clone git@github.com:Funyday-k/DaMaSCUS-SUN-EVAP.git
cd DaMaSCUS-SUN-EVAP
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF
cmake --build build --config Release --parallel
cmake --install build --config Release
```

Keep machine-specific files such as `bin/config_Lingyu.cfg`,
`bin/job_Lingyu.sh`, and scheduler output files in `bin/`. They are local
deployment state, not source files.

Run a parameter point with MPI:

```bash
mpirun -np 8 bin/DaMaSCUS-SUN bin/config.cfg
```

For a scheduler, wrap the same executable/config pair in the local batch script
used by that machine.

## Configuration

Configuration files use libconfig syntax. The most important controls are:

| Setting | Meaning |
| --- | --- |
| `run_mode` | `"Parameter point"` for the main evaporation workflow, `"Capture"` for capture-rate runs, or `"Parameter scan"` for the older scan path. |
| `capture_mode` | Boolean override for capture-only behavior. `run_mode = "Capture"` also enables capture mode. |
| `sample_size` | Target number of captured particles. Normal-mode MPI runs batch trajectories between reductions, so the final captured count may slightly exceed this target. |
| `fixed_seed` | Optional non-negative PRNG seed. `0` or an omitted setting uses nondeterministic seeding; a nonzero value is expanded independently by MPI rank. |
| `max_trajectories` | Optional hard cap on generated trajectories. `0` or unset means no trajectory-count cap. |
| `interpolation_points` | Scattering-rate interpolation grid size. `0` disables interpolation; production runs should compare representative values before fixing this. |
| `output_dir` | Root directory for generated result folders. |
| `DM_mass` | Dark matter mass in GeV. |
| `DM_cross_section_nucleon` | DM-nucleon cross section in cm^2. |
| `DM_cross_section_electron` | DM-electron cross section in cm^2 where relevant. |
| `maximum_number_of_scatterings` | Per-trajectory computational cutoff. Cutoff-terminated captures are not treated as clean physical evaporation events. |

Normal-mode MPI synchronization uses an automatic batch size selected from the
DM-nucleon cross section. The selected value is written to output headers as
`normal_mode_mpi_sync_interval`.

| `DM_cross_section_nucleon` range [cm^2] | MPI sync interval per rank |
| --- | ---: |
| `>= 1e-35` | 64 |
| `[1e-36, 1e-35)` | 128 |
| `[1e-37, 1e-36)` | 1024 |
| `[1e-38, 1e-37)` | 8192 |
| `[1e-39, 1e-38)` | 65536 |
| `< 1e-39` including `1e-40` and below | 1048576 |
| `snapshot_enabled` | Enables intermediate wall-clock progress reports for parameter-point runs. Disabled automatically in capture mode. |
| `snapshot_interval` | Positive integer wall-clock spacing, in seconds, for snapshot reports. Defaults to 60 seconds when snapshots are enabled. |
| `max_trajectory_wall_time_sec` | Optional per-trajectory wall-time guard. Snapshot recorder overhead is excluded from this budget. |

For reproducible MPI runs, a nonzero fixed seed is expanded by rank as
`base_seed + 1000003 * mpi_rank`. Computational cutoffs are tracked separately
from physical right-censoring so that final evaporation-time files contain only
complete valid unbinding events.

## Outputs

For non-capture parameter-point runs, the final files are written after MPI
reduction:

- `bincount.txt`: captured and not-captured time-weighted radial histograms with
  error estimates.
- `evaporation_times.txt`: compact complete-event table with
  `rank trajectory_id lifetime_unbinding_sec`, sorted by
  `lifetime_unbinding_sec` with `rank trajectory_id` tie-breakers.

The `bincount.txt` and snapshot report headers expose both `capture_rate_raw`
(captured / all attempted) and `capture_rate_valid` (captured / physically
classified), with separate standard errors and Wilson intervals. The legacy
`capture_rate*` fields remain aliases for the raw definition.

When snapshots are enabled, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt`: cumulative progress report at the snapshot wall time.
- `snapshot_{time}s_evaporation_times.txt`: complete valid evaporation events
  first published by that checkpoint, sorted by `lifetime_unbinding_sec`.
  An event committed concurrently with a snapshot boundary is assigned once to
  the next checkpoint rather than being dropped.

Snapshot files are progress diagnostics. They do not replace the final
post-reduction `bincount.txt` and `evaporation_times.txt` products, and they are
not restart checkpoints. A report can temporarily have `snapshot_status =
partial` while ranks publish their state. If a rank misses the deadline, the
report remains incomplete rather than reconstructing state that was not
captured. Final rank states are retained when the last merge is incomplete.
Snapshot reports likewise expose attempted, classified, and unresolved counts
plus raw and valid capture-rate intervals.

The executable requests `MPI_THREAD_FUNNELED`; only the main thread calls MPI,
while one heartbeat thread per rank performs local state copies and file I/O.
Concurrent jobs must use distinct output directories.

Snapshot checkpoint I/O is supported on homogeneous POSIX (Linux/macOS) MPI
nodes sharing the output filesystem. The binary checkpoint representation is
local to a run and is not a portable interchange or restart format.

Capture-mode runs skip the full output path and print the capture summary
instead.

## Citation

If this branch is used in analysis, cite this repository using `CITATION.cff`
and cite the original DaMaSCUS-SUN work where appropriate.

Useful references:

- T. Emken and C. Kouvaris,
  [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani and Palomares-Ruiz, evaporation of dark matter in the Sun

## License

This project is distributed under the MIT License. See `LICENSE` for the
upstream and modification copyright notices.

# miniWeather — MPI / OpenMP / OpenACC Port and Performance Study

Distributed, GPU-offloaded port of the **miniWeather** mini-app, with parallel NetCDF output, an
MPI-wide timing framework, CI, and a multi-GPU scaling study on **Leonardo Booster** (CINECA).

Coursework for *P1.8 — Best Practices in Scientific Software Development*, Master in High Performance
Computing (ICTP / SISSA, Trieste), Nov–Dec 2025.

![Simulation](code/results/movie0100.jpeg)

> **This is a presentation copy.** The project was developed by three people in
> [`prabhkodes/fightClub`](https://github.com/prabhkodes/fightClub) — that repository holds the full
> commit history, all 40 reviewed pull requests, and the raw Nsight Systems traces (~83 MB, omitted here
> to keep this clone small). Everything I claim below links to the pull request that carries it.

## Provenance

**The solver is not ours.** The 2-D compressible Euler dynamical core — finite-volume discretisation,
dimensional splitting, RK3, hyperviscosity — is [**miniWeather**](https://github.com/mrnorman/miniWeather)
by Matthew R. Norman (ORNL 2018, NVIDIA 2021), used under its BSD licence
([`code/serial/LICENSE`](code/serial/LICENSE)). It is written as a parallel-programming training app.

The project adds the parallelisation, the engineering around it, and the measurement: modular refactor ·
MPI domain decomposition with halo exchange · OpenMP · OpenACC offload with multi-GPU rank→device
mapping · parallel NetCDF I/O · CMake, Docker, Singularity · GitHub Actions CI with NetCDF regression
tests · Doxygen · MPI-wide timing, NVTX ranges, Nsight traces, `perf` counters.

Three developers, 155 commits, 40 reviewed pull requests on a `feature → dev → main` model:
[@RaionG18](https://github.com/RaionG18) (Emilio Gordillo) ·
[@formidablefrank](https://github.com/formidablefrank) (J. Franco Ray) ·
[@prabhkodes](https://github.com/prabhkodes) (Prabhsharan Singh).

## My contribution — @prabhkodes

| What | Where | PRs |
|---|---|---|
| **Parallel timing framework** — written from scratch. Instruments named regions across all ranks, reduces to per-routine max / exclusive / average / call count, reports the slowest rank per routine. Every performance number in this repo was measured with it. | [`parallel_timer.f90`](code/serial/parallel_timer.f90) (198 ln) | [#3](https://github.com/prabhkodes/fightClub/pull/3) [#14](https://github.com/prabhkodes/fightClub/pull/14) [#17](https://github.com/prabhkodes/fightClub/pull/17) |
| **Parallel NetCDF output** — authored the module. Per-rank hyperslab writes into one shared file, `start`/`count` from the domain decomposition, unlimited time dimension. | [`module_output.F90`](code/serial/module_output.F90) (200 ln) | [#30](https://github.com/prabhkodes/fightClub/pull/30) [#36](https://github.com/prabhkodes/fightClub/pull/36) |
| **OpenMP threading of the tendency stencils** — the *x*/*z* tendency loops, innermost hot loops of the RK stages. | `module_types.F90` | [#20](https://github.com/prabhkodes/fightClub/pull/20) |
| **MPI ↔ OpenACC integration** — reconciled the distributed and GPU branches into one source tree building CPU-only, MPI+OpenMP and MPI+OpenACC. | 6 files | [#34](https://github.com/prabhkodes/fightClub/pull/34) |
| **Benchmark harness and the scaling campaign** — I/O-free benchmark mode, SLURM sweeps, every CPU and multi-GPU run below, `perf` collection. | `submit_sweep.sh`, `nvtx_batch.sh` | [#39](https://github.com/prabhkodes/fightClub/pull/39) [#21](https://github.com/prabhkodes/fightClub/pull/21) [#25](https://github.com/prabhkodes/fightClub/pull/25) [#43](https://github.com/prabhkodes/fightClub/pull/43) [#46](https://github.com/prabhkodes/fightClub/pull/46) |
| **Integration and release** — merged the team's PRs, cut both releases to `main`. | — | [#52](https://github.com/prabhkodes/fightClub/pull/52) [#54](https://github.com/prabhkodes/fightClub/pull/54) |

Not mine: the **OpenACC port** is Emilio's ([#32](https://github.com/prabhkodes/fightClub/pull/32)
[#37](https://github.com/prabhkodes/fightClub/pull/37) [#42](https://github.com/prabhkodes/fightClub/pull/42)
[#49](https://github.com/prabhkodes/fightClub/pull/49)); **CMake, Doxygen and the regression test** are
Franco's ([#26](https://github.com/prabhkodes/fightClub/pull/26)
[#40](https://github.com/prabhkodes/fightClub/pull/40) [#44](https://github.com/prabhkodes/fightClub/pull/44)).

## Results — Leonardo Booster

Grid 2000 × 1000, `dt` = 0.0333 s. A100 (`cc80`), NVHPC 24.5, HPC-X MPI 2.19, netcdf-fortran 4.6.1,
built `-fast -acc -gpu=cc80`.

| GPUs | Nodes × ranks | Runtime | vs 1 GPU | Efficiency |
|---:|---|---:|---:|---:|
| 1 | 1 × 1 | 49.9 s | 1.00× | 100% |
| 4 | 1 × 4 | 21.1 s | **2.37×** | 59% |
| 8 | 2 × 4 | 16.4 s | 3.04× | 38% |
| 12 | 3 × 4 | 16.7 s | 2.99× | 25% |
| 16 | 4 × 4 | 17.1 s | 2.92× | 18% |

Saturates past 8 GPUs — the per-device subdomain gets small enough that halo exchange and launch
overhead dominate.

Against the best CPU configuration (8 nodes, 128 MPI × 2 OMP, 256 cores, **~67.4 s**): 1 × A100 is
1.35×, 4 × A100 is **3.2×**.

> **Normalisation:** the CPU baseline ran to `final time = 1000`, the GPU runs to `500`. Its raw 134.7 s
> is halved above for a like-for-like comparison. Logs: [`code/results/gpu/`](code/results/gpu/) —
> `best_cpu.init`, `multigpu_scaling_*.init`.

![CPU vs GPU](code/results/gpu/cpu_vs_gpu_comparison.png)
![Scaling](code/results/gpu/scaling_analysis.png)

`perf stat` on the CPU build: 2.40 IPC, 3.01% L1-d miss, 3.06% LLC miss, **74.6% bad speculation** — the
dominant stall. Full counters in [`perf_results/base_perf.txt`](perf_results/base_perf.txt). Raw Nsight
traces with NVTX ranges live in
[fightClub's `code/results/gpu/`](https://github.com/prabhkodes/fightClub/tree/dev/code/results/gpu).

## Build and run

Needs `cmake`, `gfortran` (or `nvfortran`), MPI, `netcdf-fortran`; `doxygen`/`graphviz` for docs.
`./runenv.sh` gives you all of it in Docker instead.

```bash
cd code/serial && mkdir build && cd build
cmake ..                      # MPI + OpenMP
cmake .. -DUSE_OPENACC=ON     # MPI + OpenACC, on a GPU system
make -j4
make test                     # NetCDF regression comparison
mpirun -n 4 ./model 100 1000 10   # nx, timesteps, output frequency
```

Writes `output.nc` — view with `ncview` or VisIt. SLURM scripts for Leonardo (native, containerised,
InfiniBand-forced, Nsight-profiling) are in [`code/serial/`](code/serial/).

---


## Source code structure
- `model.F90` contains the main routine of the program
- `module_output.F90` contains the routines relevant for output file generation in parallel
- `module_parameters.F90` contains parameters for domain decomposition and solvers, as well as physical constants
- `module_physics.F90` takes care of initial and boundary conditions, numerical solution and calculation of mass and energy budgets
- `module_types.F90` calculates the atmospheric state type, flux and tendency, as well as handle halo exchange
- `parallel_timer.F90` is used to time the execution of several routines among all ranks and generates a summary on total time, max time, average time and number of calls.


## Running on Leonardo
You can see the slurm scripts to see how the program was built and run on the cluster:
- `cpu_model.sh`, `gpu_model.sh`
- `serial/batch.sh`, `serial/gpu_batch.sh`, `serial/gpu.sh`

Sample command to upload your files to the cluster. This assumes that you have SSH key acquired from `step`:
```bash
rsync -arvzP code leo:/leonardo_scratch/large/userexternal/jrayo000
```

### Project Path
```bash
/leonardo/pub/userexternal/jgordill/fightClub
```

### Modules
#### MPI+OpenMP
```bash
module purge

# Compiler
module load cmake/3.27.9
module load gcc/12.2.0

# MPI
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2

# NetCDF Fortran
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22

```

### MPI+OpenACC
```bash
module purge

# Compiler
module load cmake/3.27.9
module load nvhpc/24.5

# MPI
module load hpcx-mpi/2.19

# NetCDF Fortran
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load binutils/2.42
```

### Python (for running output comparison test)
Load packages
```bash
module purge

module load python/3.11
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module load netcdf-c/4.9.2--openmpi--4.1.6--gcc--12.2.0-spack0.22
module load parallel-netcdf/1.12.3--openmpi--4.1.6--gcc--12.2.0-spack0.22
```

Create virtual environment and install packages. You may also have the option not to explicitly activate/deactivate the virtual env by directly using the executable files inside.
```bash
python3 -m venv pyenv
pyenv/bin/pip install numpy netCDF4
```

Run the program
```bash
pyenv/bin/python nccmp3.py output-serial.nc output-serial-optimized.nc output.nc
```

---

## Licence

miniWeather is BSD-licensed by ORNL and NVIDIA ([`code/serial/LICENSE`](code/serial/LICENSE)).
Modifications are released under the same terms.

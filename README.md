# miniWeather — MPI / OpenMP / OpenACC port

Taking an existing serial Fortran weather model and making it run on 256 CPU cores and on 16 GPUs.
190 s down to 2.1 s on CPU, and another 8.2× on GPU.

**Stack:** Fortran · MPI · OpenMP · OpenACC · NetCDF · CMake · Docker · GitHub Actions · Nsight
Systems / NVTX · `perf` · SLURM · Leonardo Booster (A100)

**Team of three.** Coursework for *P1.8 — Best Practices in Scientific Software Development*, Master in
High Performance Computing (ICTP / SISSA, Trieste), Nov–Dec 2025.

![Simulation](results/plots/movie0100.jpeg)

## Where the source code comes from

**We did not write the solver.** It is
[**miniWeather**](https://github.com/mrnorman/miniWeather) by Matthew R. Norman, Oak Ridge National
Laboratory (2018) and NVIDIA (2021). BSD licensed — the licence is in [`LICENSE`](LICENSE) and stays
there.

Norman wrote the physics: 2-D compressible Euler equations for a stratified atmosphere, finite-volume
discretisation with fourth-order flux interpolation, dimensional splitting, three-stage Runge-Kutta,
hyperviscosity. miniWeather is written to be a teaching code for parallel programming, which is exactly
what we used it for.

Everything below is what we added on top of it.

## What we did

| | |
|---|---|
| **Refactor** | Split the single-file solver into typed Fortran modules |
| **OpenMP** | Threaded the tendency stencils — the hot loops of the RK stages |
| **MPI** | 1-D domain decomposition in *x*, non-blocking halo exchange, `Allreduce` for conservation |
| **OpenACC** | GPU offload, explicit device data residency, one GPU per rank, GPU-aware MPI |
| **Parallel I/O** | NetCDF, per-rank hyperslab writes into one shared file |
| **Build** | CMake, one tree, three configurations (CPU-only, MPI+OpenMP, MPI+OpenACC) |
| **CI** | GitHub Actions — builds the Docker image, runs the test suite on every push and PR |
| **Tests** | Mass/energy conservation thresholds + NetCDF output comparison, wired into CTest |
| **Profiling** | `perf` counters, an MPI-wide timing framework, NVTX ranges for Nsight Systems |
| **Docs** | Doxygen with call graphs |

Full write-up in [`docs/presentation.pdf`](docs/presentation.pdf).

## Results

Leonardo Booster (CINECA) — Intel Ice Lake Xeon 8358, 32 cores/node, 4× A100 64 GB per node.

### What is actually being measured

Everything is **double precision** (`wp = real64`, 8 bytes). Four prognostic variables — density,
*x*-momentum, *z*-momentum, ρθ. Halo of 2 cells, fourth-order stencil. Domain is fixed at
20 km × 10 km, so `nz = nx/2` and `dt = 1.5·dx/450` from the CFL condition.

| | NX = 400 (CPU runs) | nx = 2000, nz = 1000 (GPU runs) |
|---|---:|---:|
| Cells | 80,000 | 2,000,000 |
| dx, dt | 50 m, 0.1667 s | 10 m, 0.0333 s |
| Simulated time / timesteps | 1000 s / 6,000 | 500 s / 15,000 |
| 2 × state, `(nx+4)(nz+4)·4` | 5.0 MiB | 122.8 MiB |
| flux `(nx+1)(nz+1)·4` | 2.5 MiB | 61.1 MiB |
| tendency `nx·nz·4` | 2.4 MiB | 61.0 MiB |
| **Working set** | **10.0 MiB** | **245.2 MiB** |
| Cell-updates | 4.8 × 10⁸ | 3.0 × 10¹⁰ |
| FLOPs (≈) | 274 GFLOP | 17.1 TFLOP |

At NX = 400 the whole working set is 10 MiB, which fits in L3 on these nodes. That's why the serial
profile shows only a 3.01% L1-d miss rate but 74.6% bad speculation — the bottleneck is branching, not
memory.

FLOP counts are counted by hand from the stencil body: about 87 flops per cell per directional sweep,
six sweeps per timestep (3 RK stages × 2 directions), plus the state update. Divisions count as one and
the `**cdocv` power is left out, so take it as ±20% and a lower bound.

### CPU, NX = 400

| Configuration | Time | Speedup | Mcell-upd/s | GFLOP/s |
|---|---:|---:|---:|---:|
| 1 rank, 1 thread (baseline) | 190.1 s | 1× | 2.5 | 1.4 |
| 1 rank, 32 threads | 18.2 s | 10.4× | 26.4 | 15.0 |
| 4 ranks × 8 threads, 1 node | 9.5 s | 20.0× | 50.5 | 28.8 |
| **16 ranks × 2 threads, 8 nodes** | **2.14 s** | **88.8×** | **224.3** | **127.9** |

MPI+OpenMP beat pure OpenMP at the same core count on one node, 9.5 s against 18.2 s. Four ranks each
pinned to their own NUMA domain do better than one rank spread across the whole socket.

### Multi-GPU, nx = 2000, nz = 1000

| GPUs | Time | Mcell-upd/s | GFLOP/s | Communication | Working set per GPU |
|---:|---:|---:|---:|---:|---:|
| 1 | 58.2 s | 516 | 294 | 0% | 245 MiB |
| 4 | 23.3 s | 1288 | 734 | 4% | 61 MiB |
| 8 | 17.6 s | 1705 | 972 | 10% | 31 MiB |
| 12 | 17.5 s | 1714 | 977 | 21% | 20 MiB |
| 16 | 18.2 s | 1648 | 939 | 24% | 15 MiB |

**Scaling stops at 8 GPUs.** The timer shows why: communication goes from 0% to 24% of the runtime as
each GPU gets a smaller piece of the grid. At 16 GPUs a card holds 15 MiB and barely does any work per
timestep, but the halo exchange still costs the same. The fix is to overlap communication with
computation. We didn't get to it.

294 GFLOP/s on one GPU is about 3% of the A100's 9.7 TFLOP/s fp64 peak. That's normal for a stencil
like this: not much arithmetic per byte loaded, plus divisions and a power call.

### CPU vs GPU, per component

128 MPI × 2 OMP (256 cores) against 8 GPUs on 2 nodes, nx = 2000, nz = 1000, 1000 s simulated.

| Component | CPU | GPU | Speedup |
|---|---:|---:|---:|
| Step (main loop) | 128.0 s | 14.7 s | **8.68×** |
| Communication | 32.0 s | 1.7 s | 18.52× |
| Init / thermal / hydrostatic | — | — | 0.03–0.07× |
| **Total** | | | **8.20×** |

The setup routines are 25-30× *slower* on GPU. They run once, they're not offloaded, and they pay for
device init. Over a real run length it doesn't matter, but it's there. Communication looks 18.5× better
mostly because 8 ranks have far less to exchange than 128.

![CPU vs GPU](results/plots/cpu_vs_gpu_comparison.png)
![Scaling](results/plots/scaling_analysis.png)

## My part — @prabhkodes

| What | Where | PRs |
|---|---|---|
| **Parallel timing framework**, written from scratch. Scoped timer using a Fortran `final` binding. Reports per-routine max / exclusive / average / call count across all ranks, and names the slowest rank per routine. Every number on this page was measured with it. | [`src/parallel_timer.f90`](src/parallel_timer.f90) | [#3](https://github.com/prabhkodes/fightClub/pull/3) [#14](https://github.com/prabhkodes/fightClub/pull/14) [#17](https://github.com/prabhkodes/fightClub/pull/17) |
| **Parallel NetCDF output.** Wrote the module — hyperslab writes with per-rank offsets from the decomposition, unlimited time dimension. | [`src/module_output.F90`](src/module_output.F90) | [#30](https://github.com/prabhkodes/fightClub/pull/30) [#36](https://github.com/prabhkodes/fightClub/pull/36) |
| **OpenMP threading** of the *x*/*z* tendency stencils. | `src/module_types.F90` | [#20](https://github.com/prabhkodes/fightClub/pull/20) |
| **Merged the MPI and OpenACC branches** into one source tree building all three configurations. | 6 files | [#34](https://github.com/prabhkodes/fightClub/pull/34) |
| **Benchmark harness and the whole scaling campaign** — I/O-free benchmark mode, SLURM sweeps, every run in the tables above, `perf` collection. | `scripts/slurm/` | [#39](https://github.com/prabhkodes/fightClub/pull/39) [#43](https://github.com/prabhkodes/fightClub/pull/43) [#46](https://github.com/prabhkodes/fightClub/pull/46) |
| **Integration and releases** — merged the team's PRs, cut both releases to `main`. | — | [#52](https://github.com/prabhkodes/fightClub/pull/52) [#54](https://github.com/prabhkodes/fightClub/pull/54) |

**Not mine.** The OpenACC port is Emilio's ([#32](https://github.com/prabhkodes/fightClub/pull/32),
[#37](https://github.com/prabhkodes/fightClub/pull/37),
[#42](https://github.com/prabhkodes/fightClub/pull/42),
[#49](https://github.com/prabhkodes/fightClub/pull/49)). CMake, Doxygen and the NetCDF regression test
are Franco's ([#26](https://github.com/prabhkodes/fightClub/pull/26),
[#40](https://github.com/prabhkodes/fightClub/pull/40),
[#44](https://github.com/prabhkodes/fightClub/pull/44)).

[@RaionG18](https://github.com/RaionG18) Emilio Gordillo ·
[@formidablefrank](https://github.com/formidablefrank) J. Franco Ray ·
[@prabhkodes](https://github.com/prabhkodes) Prabhsharan Singh

This is a presentation copy. The project was built in
[`prabhkodes/fightClub`](https://github.com/prabhkodes/fightClub), which has the full history, all 40
pull requests, and the raw Nsight traces (~83 MB, left out here).

## Build and run

Needs `cmake`, `gfortran` or `nvfortran`, MPI, `netcdf-fortran`. `doxygen` and `graphviz` for docs.
`scripts/runenv.sh` gives you all of it in Docker.

```bash
cmake -S . -B build                      # MPI + OpenMP
cmake -S . -B build -DUSE_OPENACC=ON     # MPI + OpenACC
cmake --build build -j
ctest --test-dir build --output-on-failure
mpirun -n 4 ./build/model 100 1000 10    # nx, timesteps, output frequency
```

Writes `output.nc` — open with `ncview` or VisIt. `cmake --build build --target doc` builds the Doxygen
docs into `build/doc/html/`.

Verified on macOS/arm64 (gfortran 14, Open MPI 5.0, netcdf-fortran 4.6) and in the Ubuntu 22.04
container CI uses.

## How the work was run

**Git.** Feature branch → PR onto `dev`. A PR merged only with a peer approval *and* a green CI run.
Releases went out as `dev → main` PRs. 40 PRs across three people.

**CI.** [`.github/workflows/ci.yml`](.github/workflows/ci.yml) builds the Docker image (Buildx, with the
Actions cache as backend) and runs the tests inside it. About 50 seconds, so review never waited on it.
Building the image instead of compiling directly is what keeps local, CI and cluster the same.

**Tests.** Two checks. First, conservation: the run passes only if fractional mass drift is under 1e-13
and energy drift under 1e-3. A broken halo exchange shows up straight away as mass appearing or
disappearing. Second, comparing the output NetCDF against reference files with `nccmp3.py`. Both run on
every configuration, so if a GPU result drifted from the CPU one we'd know.

*Note:* the reference `.nc` files were never committed, so the comparison currently no-ops and only the
conservation check runs. Worth fixing.

**Profiling.** `perf stat` on the serial build first: 256 s, 2.02 trillion instructions, 2.40 IPC, 3.01%
L1-d miss, 3.06% LLC miss — and **74.6% bad speculation**. Not memory-starved, branch-bound.

Then the timer, which found the hot spot immediately: `Computation: step` was 190.08 s of a 190.10 s
run. 99.99%. Everything after that went into one routine. The timer covers computation and
communication and deliberately skips file I/O, which would otherwise dominate and hide what we were
trying to fix.

On GPU we added NVTX ranges around the transfers and compute regions, so Nsight traces show
`Halo Exchange X` and `Runge-Kutta Integration` as named spans per rank instead of unnamed kernels.

**Then the loop:** profile → find the hot spot → change it → re-check scaling → repeat.

## Some implementation notes

**Domain decomposition.** `setup_domain_decomposition` splits the global grid in *x*, spreads the
remainder across the first ranks so load stays even when `nx` doesn't divide cleanly, and wires up
periodic neighbours.

**Halo exchange.** Non-blocking `MPI_Irecv`/`MPI_Isend` with a single `MPI_Waitall`. Mass and energy
reduce with `MPI_Allreduce`.

**GPU data residency.** `enter data create/copyin/attach` at setup, `update self` only when output is
actually needed, `exit data delete` at teardown. Compute regions use `present(...)` instead of implicit
copies, so nothing quietly copies back and forth every timestep.

**GPU-aware MPI.** Halo exchange uses `!$acc host_data use_device(...)`, handing device pointers straight
to MPI so buffers move device-to-device without staging through the host. Halo buffers are allocated
once and stay resident.

## What we'd do differently

Overlap communication with computation. The multi-GPU table shows exactly what that's worth. And run
the tests on the cluster from CI, not just in a container on GitHub's runners, so toolchain failures
that only happen on the cluster get caught by the pipeline instead of by a person.

Smaller lessons: module load order matters, and gives you confusing link errors when it's wrong. One
GPU per rank. Check your results every time, not just when something looks off. Profile even when the
code seems fine — nobody would have guessed the 74.6% bad speculation without measuring it.

## Layout

```
CMakeLists.txt          single build system; -DUSE_OPENACC=ON switches toolchain
Dockerfile
LICENSE                 miniWeather BSD licence (ORNL / NVIDIA)
src/                    Fortran sources
tests/                  NetCDF comparator, Python requirements
scripts/
  runenv.sh             build the image and shell into it
  slurm/                Leonardo batch scripts and the scaling sweep
  profiling/            Nsight / NVTX job scripts
docs/                   Doxyfile, CSS, physics notes, presentation
results/
  plots/                figures used above
  cpu/                  CPU run logs and perf counters
  gpu/                  multi-GPU run logs, Nsight session view
  analysis/             the plotting scripts that made the figures
```

| Source file | Contents |
|---|---|
| `model.F90` | Main driver — init, RK time stepping, diagnostics |
| `module_physics.f90` | Initial and boundary conditions, solution, mass/energy budgets |
| `module_types.F90` | State, flux and tendency types; halo exchange |
| `module_parameters.f90` | Decomposition and solver parameters, physical constants |
| `module_output.F90` | Parallel NetCDF output |
| `parallel_timer.f90` | Per-routine timing across all ranks |
| `module_nvtx.F90` | NVTX ranges; no-ops when built without NVTX |

## Running on Leonardo

Batch scripts are in [`scripts/slurm/`](scripts/slurm/) — `cpu_omp.sh`, `gpu.sh`, `gpu_batch.sh`,
`submit_sweep.sh` — and [`scripts/profiling/nvtx_batch.sh`](scripts/profiling/nvtx_batch.sh) for Nsight
runs.

**MPI + OpenMP**

```bash
module purge
module load cmake/3.27.9 gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22
```

**MPI + OpenACC**

```bash
module purge
module load cmake/3.27.9 nvhpc/24.5 hpcx-mpi/2.19
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load binutils/2.42
```

Rank pinning used `OMP_PROC_BIND=close` and `OMP_PLACES=cores`.

## Licence

miniWeather is BSD licensed by ORNL and NVIDIA — see [`LICENSE`](LICENSE). Our changes are under the
same terms.

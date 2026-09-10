# miniWeather — MPI / OpenMP / OpenACC port

Taking an existing serial Fortran weather model and making it run on 256 CPU cores and on 16 GPUs.
190 s down to 2.1 s on CPUs; another 8.2× on GPUs.

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

Norman wrote the physics: the 2-D compressible Euler equations for a stratified atmosphere,
finite-volume discretisation with fourth-order flux interpolation, dimensional splitting, three-stage
Runge–Kutta, hyperviscosity. miniWeather exists specifically as a teaching code for parallel
programming, which is what it was used for here.

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

Leonardo Booster (CINECA) — Intel Ice Lake, 4× A100 per node.

### CPU, NX = 400

| Configuration | Time | Speedup |
|---|---:|---:|
| 1 rank, 1 thread (baseline) | 190.1 s | 1× |
| 1 rank, 32 threads | 18.2 s | 10.4× |
| 4 ranks × 8 threads, 1 node | 9.5 s | 20.0× |
| **16 ranks × 2 threads, 8 nodes** | **2.14 s** | **88.8×** |

Hybrid MPI+OpenMP beat pure OpenMP at the same core count on one node — 9.5 s against 18.2 s. Four
ranks pinned to their own NUMA domains beat one rank spanning the socket.

### Multi-GPU, nx = 2000, nz = 1000

| GPUs | Time | Communication |
|---:|---:|---:|
| 1 | 58.2 s | 0% |
| 4 | 23.3 s | 4% |
| 8 | 17.6 s | 10% |
| 12 | 17.5 s | 21% |
| 16 | 18.2 s | 24% |

**Scaling stops at 8 GPUs.** The timer says why: communication climbs from 0% to 24% of runtime as the
per-device subdomain shrinks. Past 8 GPUs the halo exchange costs more than the compute it saves.
Overlapping communication with computation is the obvious next step and we didn't get to it.

### CPU vs GPU, per component

128 MPI × 2 OMP (256 cores) against 8 GPUs on 2 nodes.

| Component | CPU | GPU | Speedup |
|---|---:|---:|---:|
| Step (main loop) | 128.0 s | 14.7 s | **8.68×** |
| Communication | 32.0 s | 1.7 s | 18.52× |
| Init / thermal / hydrostatic | — | — | 0.03–0.07× |
| **Total** | | | **8.20×** |

The setup routines are 25–30× *slower* on GPU. They run once, they aren't offloaded, and they pay
device initialisation. It doesn't matter over a real run length, but it's there. Communication improves
18.5× mostly because 8 ranks exchange far less than 128 do.

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
Building the image rather than compiling directly is what keeps local, CI and cluster identical.

**Tests.** Two checks. First, physical conservation — the run passes only if fractional mass drift is
below 1e-13 and energy drift below 1e-3. A broken halo exchange shows up immediately as mass appearing
or disappearing. Second, an output comparison against reference NetCDF files via `nccmp3.py`. Both run
against every configuration, so a GPU result that diverged from the CPU result got caught.

*Note:* the reference `.nc` files were never committed, so the comparison currently no-ops and only the
conservation check runs. Worth fixing.

**Profiling.** `perf stat` on the serial build first: 256 s, 2.02 trillion instructions, 2.40 IPC, 3.01%
L1-d miss, 3.06% LLC miss — and **74.6% bad speculation**. Not memory-starved, branch-bound.

Then the timer, which localised the cost immediately: `Computation: step` was 190.08 s of a 190.10 s
run. 99.99%. Everything after that targeted one routine. Timing covers computation and communication
and deliberately excludes file I/O, which would otherwise dominate and hide what we were optimising.

On GPU, NVTX ranges around the transfers and compute regions, so Nsight traces show `Halo Exchange X`
and `Runge-Kutta Integration` as named spans per rank instead of anonymous kernels.

**Then the loop:** profile → find the hot spot → change it → re-check scaling → repeat.

## Some implementation notes

**Domain decomposition.** `setup_domain_decomposition` splits the global grid in *x*, spreads the
remainder across the first ranks so load stays even when `nx` doesn't divide cleanly, and wires up
periodic neighbours.

**Halo exchange.** Non-blocking `MPI_Irecv`/`MPI_Isend` with a single `MPI_Waitall`. Mass and energy
reduce with `MPI_Allreduce`.

**GPU data residency.** `enter data create/copyin/attach` at setup, `update self` only when output is
actually needed, `exit data delete` at teardown. Compute regions use `present(...)` rather than implicit
copies, so nothing silently round-trips every timestep.

**GPU-aware MPI.** Halo exchange uses `!$acc host_data use_device(...)`, handing device pointers straight
to MPI so buffers move device-to-device without staging through the host. Halo buffers are allocated
once and stay resident.

## What we'd do differently

Overlap communication with computation — the multi-GPU table shows exactly what that's worth. And run
the tests on the cluster from CI, not just in a container on GitHub's runners, so cluster-only toolchain
failures get caught by the pipeline instead of by a person.

Smaller lessons: module load order matters and produces baffling link errors when wrong. One GPU per
rank. Check results every time, not just when something looks off. Profile even when the code looks
fine — the 74.6% bad-speculation number was invisible until measured.

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

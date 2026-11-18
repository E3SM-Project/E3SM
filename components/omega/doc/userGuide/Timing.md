(omega-user-timing)=

# Timing

Omega uses the Pacer library for timing the code and incorporates timers around
various parts of the code.

By default, the timing output is written to two files: `omega.summary` and `omega.timing0`.
The `omega.summary` file presents accumulated timing statistics across all MPI ranks.
The `omega.timing.0` show timing result only from the first rank.

There are four parameters that are set by the user in the input configuration
file that control the timing behavior. These are:
```yaml
Timing:
   Level: 2
   AutoFence: true
   TimingBarriers: false
   PrintAllRanks: false
```
The `Level` parameter is a non-negative integer that determines the granularity of timers.
Increasing it will turn on more timers.
Having more timers provides more detailed information, but it also comes with increased overhead,
and may be counter-productive if a high-level look at model performance is sufficient.

The `Autofence` Boolean option determines if Kokkos fences are automatically added before every timer call.
This option **needs** to be true for accurate timing using Omega timers on GPU-based systems.
However, there are circumstances when turning off automatic fences is useful.
The main use case is using external profiling tools.
Another one is measuring the overhead of automatic synchronization for very high timing levels.

The `TimingLevel` Boolean option determines if MPI barriers are added before or after certain timers.
Adding barriers may be necessary to properly measure communication time, but it can add top much overhead in
production runs.

The `PrintAllRanks` Boolean option determines if all ranks should print their timing information. If this
option is set to `true` the output will include additional files `omega.timing.i` with the
timing data from rank `i`.

## Integration with external profiling tools

External profilers often include APIs to mark regions of code for detailed profiling.
On some platforms, Omega timers automatically add these annotations.
Currently, this is only implemented on systems with NVIDIA GPUs using NVTX.

This allows, for example, to use the Nsight Compute kernel profiler to obtain
detailed kernel information for all kernels enclosed in the `Tend:computeVelocityTendencies`
Omega timer.
```bash
mpirun -np 1 ncu --nvtx --nvtx-include "Tend:computeVelocityTendencies/" omega.exe
```

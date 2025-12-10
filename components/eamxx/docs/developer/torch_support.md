# Support for Torch-based emulators

In EAMxx, it is possible to wrap a ML emulator inside an atmosphere process (AP) relatively easily.
Currently, there are two ways to do so:

- via python: EAMxx can call python functions via pybind11 interfaces, so one can write a small
  py module that wraps a pytorch model
- via LAPIS: [LAPIS](https://github.com/sandialabs/lapis) is a lightweight library
  that can be used to translate torch-mlir into pure c++/kokkos code, so one can
  write a small py script that dumps a torch model as MLIR, and then use LAPIS to
  generate a kokkos equivalent version of that model.

In the folder `src/physics/cld_fraction/cld_frac_net` one can find example files for how to
wrap an emulator inside an AP. These files implement an AP called "CldFracNet",
which wraps a torch-based ML emulator. We emphasize that the emulator is for illustration purposes only,
and the quality of the emulator is beyond the scope of the example.

## Python

In order for EAMxx to call python, the code must be compiled with the cmake option `EAMXX_ENABLE_PYTHON=ON`.
For a CIME case, one can do `xmlchange --append SCREAM_CMAKE_OPTIONS='EAMXX_ENABLE PYTHON ON'` from the case dir.

When python support is enabled, the AP base class will take care of creating python-compatible
arrays (from now on, "PyFields") for each of the process input/output `Field`'s,
provided that the process input parameters contain a non-null string for the option `py_module_name`.

At runtime, the derived process implementation is in charge of gathering all PyField's, and call the
desired function from the python module. Here's an example from the `cld_frac_net` process:

```c++
    auto py_qi  = get_py_field_dev("qi");
    auto py_liq = get_py_field_dev("cldfrac_liq");
    auto py_ice = get_py_field_dev("cldfrac_ice");
    auto py_tot = get_py_field_dev("cldfrac_tot");

    double ice_threshold = m_params.get<double>("ice_cloud_threshold");

    py_module_call("forward",ice_threshold,py_qi,py_liq,py_ice,py_tot);
```

The function `get_py_field_dev` is implemented in AP base class, and retrieves the PyField
corresponding to the given input or output field. The function `py_module_call` is also implemented
in the AP base class, and wraps the call to the python function, with some exception
handling to ensure errors are nicely printed in case something goes wrong.

Some comments:

- Since Python does not have the concept of "const", it is technically possible to pass a const
  field as an output when calling python, so care is needed to ensure arguments are passed in the correct
  order.
- By default, the AP class looks for the module provided via `py_module_name` in the current working
  directory. The user can override this by specifying the optional `py_module_path` parameter.
- The AP base class also provides `get_py_field_host`. When the model is running on GPU, this allows
  to obtain a PyField that shares the same pointer of the Host view in the Field. There is NO
  automatic synchronization of device and host views. Like for regular fields, it is the process
  responsibility to copy to host at the beginning of `run_impl`, and make sure that the field is
  sync-ed back to device upon return.
- The file `cld_frac_net.py` provides the definition of the torch model, as well as the initialization
  and run routines that are called at runtime. The names of these routines are arbitrary, but they must
  of course match what the process interface calls in the C++ interface.

## C++/Kokkos

In order to generate C++/Kokkos, we can make use of the LAPIS library. LAPIS is a lightweight C++
library built on top of LLVM, which can interpret MLIR code and translate it to C++ with Kokkos.
At the time of this writing (Dec 2025), LAPIS relies on a particular version of LLVM (as well as
`torch_mlir`, which is optional). Instructions on how to install these LAPIS dependencies
can be found on the LAPIS documentation on GitHub.

For the CldFracNet example, we used LAPIS built with torch-mlir support. The file `gen_cpp_cld_frac_net.py`
shows how to load a torch based model, dump it as MLIR, and then use lapis tools to convert this
to C++ and Kokkos. LAPIS can convert a model that is meant to be called as a "top-level" kernel, as well
as a function inside a kernel. In our example, we opted for the second strategy, since it is the most
versatile one (e.g., one could replace a small portion of a device kernel). In the process implementation,
we provide the syntactic sugar necessary to then wrap the converted model inside a `Kokkos::parallel_for`.

The general structure of the code that calls the emulator is

```c++
    auto policy = ...;
    // Allow up to this much level 0 scratch (aka shared mem) per team
    // when deciding which temporary views to spill into level 1 scratch (global pool).
    constexpr int max_shared = 4096;
    constexpr int scratch0_required = forward_L0_scratch_required(max_shared);
    constexpr int scratch1_required = forward_L1_scratch_required(max_shared);
    GlobalViews_forward globalViews;
    policy.set_scratch_size(0, Kokkos::PerTeam(scratch0_required));
    policy.set_scratch_size(1, Kokkos::PerTeam(scratch1_required));

    // Execute the appropriate specialization based on shared amount
    auto lambda = KOKKOS_LAMBDA (const MemberType& team) {
      // If needed, prepare inputs for the forward function. E.g., subview input/output
      // views at a particular index corresponding to the league_rank

      char* scratch0 = (char*)(team.team_scratch(0).get_shmem(scratch0_required));
      char* scratch1 = (char*)(team.team_scratch(1).get_shmem(scratch1_required));

      forward<ExeSpace, max_shared>(team, globalViews, <output views>, <input views>, scratch0, scratch1);
    };
    Kokkos::parallel_for(policy,lambda);
```

where the function names `forward_L0_scratch_required`, `forward_L1_scratch_required`, and `forward` are
hard-coded by LAPIS, while `MemberType` and `ExeSpace` are typedefs similar to what is done in other
areas of EAMxx.

Some comments:

- The code needed to wrap the C++/Kokkos version of CldFracNet is larger. This is because we cannot provide
  any "common" implementation in the AP base class (like it was done for the Python version), and all
  the necessary plumbing must happen in the derived process implementation.
- The script `gen_cpp_cld_frac_net.py` illustrates how to convert the torch model to C++/Kokkos. The script
  first creates and initializes the torch model (defined in `cld_frac_net.py`), then dumps MLIR code using
  the torch-mlir package, and finally calls a couple of LAPIS executables to generate the hpp/cpp files.
  These files must be included, built, and linked by the host project.

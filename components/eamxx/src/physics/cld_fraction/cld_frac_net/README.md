# CldFracNet: a dummy emulator for EAMxx's CldFraction atmosphere process

This emulator serves the goal of illustrating how to wrap a torch model in EAMxx.
There are two ways to wrap a torch model in an atm process:

 - via pytorch: EAMxx can call an arbitrary python module, so one can write a small
   py module that wraps a pytorch model
 - via LAPIS: [LAPIS](https://github.com/sandialabs/lapis) is a lightweight library
   that can be used to translate torch-mlir into pure c++/kokkos code. So one can
   write a small py script that dumps a torch model as MLIR, and then use LAPIS to
   generate a kokkos equivalent version of that model

The files in this folder illustrate this process:

 - `cld_frac_net.py`: this file contains the model definition and a free function `forward`,
   that can be called from C++ (via python bindings)
 - `cld_frac_net_weights.pth`: contains the weights for the model
 - `gen_cpp_cld_frac_net.py`: this script shows how one can load the model, dump it as MLIR,
   and finally call lapis to generate the hpp/cpp file to build in the project
 - `cld_frac_net.*pp`: these are the generated files

We point out that, at the time of this writing (Dec 2025), LAPIS relies on a particular version
of LLVM as well as `torch_mlir`. Instructions on how to install these LAPIS dependencies
can be found on the LAPIS repo documentation.

Current EAMxx testing uses all the `cld_frac_net*` files, while the `gen_cpp_fld_frac_net.py` file
is shipped just as an example for how to generate the hpp/cpp files with LAPIS

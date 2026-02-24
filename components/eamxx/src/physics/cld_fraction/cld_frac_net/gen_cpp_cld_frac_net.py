#!/usr/bin/env python3

from cld_frac_net import create_cld_frac_net
import os
import sys
import torch
import torch_mlir.torchscript

CURR_FILE_PATH = os.path.dirname(os.path.abspath(__file__))
EAMXX_SCRIPTS_PATH = os.path.join(os.path.abspath(CURR_FILE_PATH), '../../../../scripts')

sys.path.append(EAMXX_SCRIPTS_PATH)
from utils import run_cmd_no_fail

def main ():

    model = create_cld_frac_net()
    col = torch.rand(model.nlevs)
    mlir_module = torch_mlir.torchscript.compile(model, (col,col), output_type='linalg-on-tensors')
    with open("cld_frac_net.mlir",'w') as fd:
        fd.write(str(mlir_module))

    # Lower to kokkos dialect
    run_cmd_no_fail("lapis-opt --team-compiler-kokkos cld_frac_net.mlir -o cld_frac_net_lowered.mlir", from_dir=os.getcwd())

    # Generate c++ code from kokkos mlir
    run_cmd_no_fail("lapis-translate cld_frac_net_lowered.mlir --team-level -o cld_frac_net.cpp -hpp cld_frac_net.hpp", from_dir=os.getcwd())

if __name__ == "__main__":
    main()

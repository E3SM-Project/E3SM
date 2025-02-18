
from machines_specs import assert_machine_supported, get_machine, get_mach_env_setup_command
from utils import expect

CHOICES = (
    "cxx_compiler",
    "c_compiler",
    "f90_compiler",
    "batch",
    "env",
    "baseline_root",
    "cuda",
    "comp_j",
    "test_j",
)

###############################################################################
def query_eamxx(machine, param):
###############################################################################
    assert_machine_supported(machine)
    expect(param in CHOICES, f"Unknown param {param}")

    mach = get_machine(machine)
    if param == "cxx_compiler":
        return mach.cxx_compiler
    elif param == "c_compiler":
        return mach.c_compiler
    elif param == "f90_compiler":
        return mach.ftn_compiler
    elif param == "batch":
        return mach.batch
    elif param == "env":
        return get_mach_env_setup_command(machine)
    elif param == "baseline_root":
        return mach.baselines_dir
    elif param == "cuda":
        return str(mach.gpu_arch == "cuda")
    elif param == "comp_j":
        return num_bld_res
    elif param == "test_j":
        return gnum_run_res
    else:
        expect(False, f"Unhandled param {param}")

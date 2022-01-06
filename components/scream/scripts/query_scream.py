
from machines_specs import assert_machine_supported, \
    get_mach_cxx_compiler, get_mach_c_compiler, get_mach_f90_compiler, \
    get_mach_batch_command, get_mach_env_setup_command, \
    get_mach_baseline_root_dir, is_cuda_machine, \
    get_mach_compilation_resources, get_mach_testing_resources
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
def query_scream(machine, param):
###############################################################################
    assert_machine_supported(machine)
    expect(param in CHOICES, f"Unknown param {param}")

    if param == "cxx_compiler":
        return get_mach_cxx_compiler(machine)
    elif param == "c_compiler":
        return get_mach_c_compiler(machine)
    elif param == "f90_compiler":
        return get_mach_f90_compiler(machine)
    elif param == "batch":
        return get_mach_batch_command(machine)
    elif param == "env":
        return get_mach_env_setup_command(machine)
    elif param == "baseline_root":
        return get_mach_baseline_root_dir(machine)
    elif param == "cuda":
        return str(is_cuda_machine(machine))
    elif param == "comp_j":
        return get_mach_compilation_resources()
    elif param == "test_j":
        return get_mach_testing_resources(machine)
    else:
        expect(False, f"Unhandled param {param}")

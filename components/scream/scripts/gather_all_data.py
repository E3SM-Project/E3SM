from utils import run_cmd_no_fail

import os
import concurrent.futures as threading3

# MACHINE -> (env_setup, compiler, batch submit prefix)
MACHINE_METADATA = {
    "melvin"   : (["module purge && module load sems-env && module load sems-gcc/7.3.0 sems-openmpi/1.10.1 sems-gcc/7.3.0 sems-git/2.10.1 sems-cmake/3.10.3 sems-python/3.5.2", "export CTEST_PARALLEL_LEVEL=48"],
                  "$(which mpicxx)",
                  ""),
    "bowman"   : (["module purge && module load openmpi/1.10.6/intel/17.2.174 git/2.8.2 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-bowman/bin:$PATH", "export CTEST_PARALLEL_LEVEL=48"],
                  "$(which mpicxx)",
                  "srun"),
    "blake"    : (["module purge && module load openmpi/2.1.5/intel/19.1.144 git/2.9.4 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-blake/bin:$PATH", "export CTEST_PARALLEL_LEVEL=48"],
                  "$(which mpicxx)",
                  "srun"),
    "waterman" : (["module purge && module load devpack/latest/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 git/2.10.1", "module switch cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-waterman/bin:$PATH"],
                  "$(which mpicxx)",
                  "bsub -I -q rhel7W"),
    "white"    : (["module purge && module load devpack/20181011/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 git/2.10.1 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-white/bin:$PATH"],
                  "$(which mpicxx)",
                  "bsub -I -q rhel7G"),
    "lassen" : (["module purge && module load gcc/7.3.1 cuda/10.1.243 cmake/3.14.5 spectrum-mpi netcdf/4.7.0 python/3.7.2 && export LLNL_USE_OMPI_VARS='y'"],
                  "$(which mpicxx)",
                  "bsub -Ip"),
    "quartz" : (["module purge &&  module load StdEnv cmake/3.14.5 mkl netcdf/4.4.1.1 mvapich2/2.3 intel/19.0.4 python/3.7.2"],
                  "$(which mpicxx)",
                  "srun"),
    "syrah" : (["module purge && module load StdEnv cmake/3.14.5 mkl netcdf/4.4.1.1 mvapich2/2.3 intel/19.0.4 python/3.7.2"],
                  "$(which mpicxx)",
                  "srun"),
}

###############################################################################
class GatherAllData(object):
###############################################################################

    ###########################################################################
    def __init__(self, run, commit, machines, scream_docs, local, kokkos):
    ###########################################################################
        self._run         = run
        self._commit      = commit
        self._machines    = machines
        self._scream_docs = scream_docs
        self._local       = local
        self._kokkos      = kokkos

    ###########################################################################
    def formulate_command(self, machine):
    ###########################################################################
        env_setup, compiler, batch = MACHINE_METADATA[machine]

        # Compute locations of key repos
        if self._local:
            scream_docs_repo = os.path.abspath("./scream-docs/micro-apps")
            scream_repo      = os.path.abspath("./scream/components/scream")
        else:
            scream_docs_repo = "~/scream-docs-perf-{}/micro-apps".format(machine)
            scream_repo      = "~/scream-perf-{}/components/scream".format(machine)

        repo       = scream_docs_repo if self._scream_docs else scream_repo

        # Need to know kokkos location in order to set up OMPI_CXX
        if self._kokkos:
            kokkos_loc = self._kokkos.replace("$machine", machine).replace("$compiler", compiler)
        else:
            kokkos_loc = os.path.join(os.path.dirname(os.path.dirname(scream_repo)), "externals", "kokkos")

        local_cmd = self._run
        # Do magic replacements here
        local_cmd = local_cmd.\
            replace("$compiler", compiler).\
            replace("$kokkos", kokkos_loc).\
            replace("$machine", machine).\
            replace("$scream_docs", scream_docs_repo).\
            replace("$scream", scream_repo)

        # Scream-docs tests may depend on scream's Kokkos and scripts, so update scream repo too if we're doing
        # a scream-docs test.
        setup = ""
        if (not self._local and self._scream_docs):
            setup = "cd {} && git fetch && git reset --hard origin/master && git submodule update --init && "\
                .format(scream_repo)

        extra_env = ""
        if machine in ["waterman", "white", "lassen"] and "OMPI_CXX" not in " ".join(env_setup):
            extra_env = "OMPI_CXX={}/bin/nvcc_wrapper ".format(kokkos_loc)
        else:
            extra_env = "OMP_PROC_BIND=FALSE "

        repo_setup = "true" if (self._local) else "git fetch && git checkout {} && git submodule update --init".format(self._commit)

        cmd = "{}cd {} && {} && {} && {}{} {}".format(setup, repo, " && ".join(env_setup), repo_setup, extra_env, batch, local_cmd)

        return cmd

    ###########################################################################
    def run_on_machine(self, machine):
    ###########################################################################
        cmd = self.formulate_command(machine)
        print("Starting analysis on {} with cmd: {}".format(machine, cmd))

        if self._local:
            run_cmd_no_fail(cmd, arg_stdout=None, arg_stderr=None, verbose=True, exc_type=RuntimeError)
        else:
            try:
                output = run_cmd_no_fail("ssh -o StrictHostKeyChecking=no {} '{}'".format(machine, cmd), exc_type=RuntimeError, combine_output=True)
            except RuntimeError as e:
                output = str(e)
                raise
            finally:
                result_path = os.path.join("gather-all-results", self._commit, machine)
                if os.path.exists(result_path):
                    old_path = result_path + ".old"
                    while (os.path.exists(old_path)):
                        old_path += ".old"
                    print("Warning moving old results to {}".format(old_path))
                    os.rename(result_path, old_path)

                with open(result_path, "w") as fd:
                    fd.write(output)

        print("Completed analysis on {}".format(machine))

    ###########################################################################
    def gather_all_data(self):
    ###########################################################################
        if not self._local:
            os.makedirs(os.path.join("gather-all-results", self._commit), exist_ok=True)

        success = True

        with threading3.ThreadPoolExecutor(max_workers=len(self._machines)) as executor:
            future_to_machine = {executor.submit(self.run_on_machine, machine): machine for machine in self._machines}
            for future in threading3.as_completed(future_to_machine):
                machine = future_to_machine[future]
                try:
                    future.result()
                except RuntimeError:
                    print('{} failed'.format(machine))
                    success = False

        return success

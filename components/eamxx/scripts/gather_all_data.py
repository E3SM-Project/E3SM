from utils import run_cmd_no_fail

import os, pathlib
import concurrent.futures as threading3
from machines_specs import get_mach_env_setup_command, get_mach_batch_command, \
                           get_mach_cxx_compiler, get_mach_f90_compiler, get_mach_c_compiler

###############################################################################
class GatherAllData(object):
###############################################################################

    ###########################################################################
    def __init__(self, run, commit, machines, scream_docs, local, kokkos, root_dir, dry_run):
    ###########################################################################
        self._run         = run
        self._commit      = commit
        self._machines    = machines
        self._scream_docs = scream_docs
        self._local       = local
        self._kokkos      = kokkos
        self._root_dir    = root_dir
        self._dry_run     = dry_run

    ###########################################################################
    def formulate_command(self, machine):
    ###########################################################################
        # gather-all runs on login, not compute nodes, so we cannot rely on
        # probed values to be accurate, so we should not set CTEST_PARALLEL_LEVEL
        # based on these potentially wrong values.
        env_setup    = get_mach_env_setup_command(machine, ctest_j=-1)

        cxx_compiler = get_mach_cxx_compiler(machine)
        f90_compiler = get_mach_f90_compiler(machine)
        c_compiler   = get_mach_c_compiler(machine)
        batch        = get_mach_batch_command(machine)

        env_setup_str = " && ".join(env_setup)

        root_dir = pathlib.Path(str(self._root_dir).replace("$machine", machine))

        # Compute locations of key repos
        if self._local:
            if self._scream_docs:
                scream_docs_repo = root_dir
                scream_repo      = pathlib.Path(__file__).resolve().parent.parent
            else:
                scream_docs_repo = None
                scream_repo      = root_dir
        else:
            if self._scream_docs:
                scream_docs_repo = root_dir
                scream_repo      = pathlib.Path("~/scream-perf-{}/components/eamxx".format(machine))
            else:
                scream_docs_repo = None
                scream_repo      = root_dir

        repo = scream_docs_repo if self._scream_docs else scream_repo

        # Compute kokkos location. scream itself won't allow this, and test-all-scream will error out,
        # but one can use gather-all-data for plenty of other apps, so we need to support this here.
        local_cmd = self._run
        if self._kokkos:
            # Note: kokkos_loc can use other magic strings, like $cxx_compiler, or $machine
            kokkos_loc = self._kokkos
        else:
            kokkos_loc = pathlib.Path(scream_repo.parent.parent, "externals", "ekat", "extern", "kokkos")

        # Do magic replacements here
        local_cmd = local_cmd.\
            replace("$kokkos", str(kokkos_loc)).\
            replace("$cxx_compiler", cxx_compiler).\
            replace("$f90_compiler", f90_compiler).\
            replace("$c_compiler", c_compiler).\
            replace("$machine", machine).\
            replace("$scream_docs", str(scream_docs_repo)).\
            replace("$scream", str(scream_repo))

        # Scream-docs tests may depend on scream's Kokkos and scripts, so update scream repo too if we're doing
        # a scream-docs test.
        setup = ""
        if (not self._local and self._scream_docs):
            setup = "cd {} && git fetch && git reset --hard origin/master && git submodule update --init --recursive && "\
                .format(scream_repo)

        repo_setup = "true" if (self._local) else "git fetch && git checkout {} && git submodule update --init --recursive".format(self._commit)

        cmd = "{setup}cd {repo} && {env_setup} && {repo_setup} && {batch} {local_cmd}".\
              format(setup=setup, repo=repo, env_setup=env_setup_str, repo_setup=repo_setup, batch=batch, local_cmd=local_cmd)

        return cmd

    ###########################################################################
    def run_on_machine(self, machine):
    ###########################################################################
        cmd = self.formulate_command(machine)
        print("Starting analysis on {} with cmd: {}".format(machine, cmd))

        if self._local:
            run_cmd_no_fail(cmd, arg_stdout=None, arg_stderr=None, verbose=True, dry_run=self._dry_run, exc_type=RuntimeError)
        else:
            try:
                ssh_cmd = "ssh -o StrictHostKeyChecking=no {} '{}'".format(machine, cmd)
                output = run_cmd_no_fail(ssh_cmd, dry_run=self._dry_run, exc_type=RuntimeError, combine_output=True)
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

                with open(result_path, "w", encoding="utf-8") as fd:
                    fd.write(output)

        print("Completed analysis on {}".format(machine))

    ###########################################################################
    def gather_all_data(self):
    ###########################################################################
        if not self._local:
            os.makedirs(os.path.join("gather-all-results", self._commit), exist_ok=True)

        success = True

        if len(self._machines) > 1:
            with threading3.ThreadPoolExecutor(max_workers=len(self._machines)) as executor:
                future_to_machine = {executor.submit(self.run_on_machine, machine): machine for machine in self._machines}
                for future in threading3.as_completed(future_to_machine):
                    machine = future_to_machine[future]
                    try:
                        future.result()
                    except RuntimeError:
                        print('{} failed'.format(machine))
                        success = False

        else:
            machine = self._machines[0]
            try:
                self.run_on_machine(machine)
            except RuntimeError:
                print('{} failed'.format(machine))
                success = False

        return success

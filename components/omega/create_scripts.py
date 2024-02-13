""" Collect compiler and machine info from CIME

This script should be invoked by Omega cmake build system
"""

import argparse
import getpass
import os
import re
import stat
import subprocess
import sys
import typing

pat_envvar = re.compile(r'^([_\d\w]+)=(.*)$', flags=re.MULTILINE)


def parse_cmdline():

    here, progname = os.path.split(__file__)

    parser = argparse.ArgumentParser(
        prog=progname,
        description='generate machine specific infomation',
        epilog='Contact: <T.B.D.>')

    parser.add_argument('-p', '--cimepath', default=os.path.realpath(
                        os.path.join(here, "..", "..", "cime")))  # CIME root
    parser.add_argument('-o', '--outpath', default="_Omega.cmake")  # outfile
    parser.add_argument('-m', '--machine')  # machine
    parser.add_argument('-c', '--compiler')  # compiler
    parser.add_argument('-d', '--debug', action='store_const',
                        const="TRUE", default="FALSE")  # debug mode
    parser.add_argument('-v', '--verbose')  # verbose output

    return parser.parse_args()


args = parse_cmdline()

# import Machines class from CIME
sys.path.insert(0, args.cimepath)
from CIME.XML.machines import Machines  # noqa: E402

sys.path.pop(0)


# run a shell command
def run_cmd_no_fail(cmd, env):
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, env=env)
    retval = str(out.stdout, 'UTF-8').strip()
    return retval


# main class that extends Machines class
class OmegaMachines(Machines):

    def __init__(self, cimepath, outpath, machine, compiler, debug):

        if isinstance(machine, str):
            machine = machine.strip()

        super(OmegaMachines, self).__init__(machine=machine)

        self.cimepath = cimepath
        self.outpath = outpath
        self.debug = debug
        self.mpilib = ""
        self.mpiexec = ""

        if compiler is None:
            self.compiler = self.get_default_compiler().strip()
        else:
            self.compiler = compiler.strip()

        self.machpath = os.path.join(self.cimepath, "..", "cime_config",
                                     "machines")
        self.macrospath = os.path.join(self.machpath, "cmake_macros",
                                       "Macros.cmake")

        self.machname = self.get_machine_name()
        self.machos = self.get_value("OS")
        self.mpilibs = self.get_value("MPILIBS").split(",")

    # modifed based on generic_xml.py in CIME
    def get_processed_value(self, raw_value, outvar):

        reference_re = re.compile(r"\${?(\w+)}?")
        env_ref_re = re.compile(r"\$ENV\{(\w+)\}")
        shell_ref_re = re.compile(r"\$SHELL\{([^}]+)\}")
        math_re = re.compile(r"\s[+-/*]\s")
        item_data = raw_value

        if item_data is None:
            return None

        if not isinstance(item_data, str):
            return item_data

        for m in env_ref_re.finditer(item_data):
            env_var = m.groups()[0]

            if env_var in outvar:
                item_data = item_data.replace(m.group(), outvar[env_var])

            elif env_var in os.environ:
                item_data = item_data.replace(m.group(), os.environ[env_var])

            else:
                raise Exception(f"Undefined env var '{env_var}'")

        for s in shell_ref_re.finditer(item_data):
            shell_cmd = s.groups()[0]
            item_data = item_data.replace(s.group(),
                                          run_cmd_no_fail(shell_cmd, outvar))

        for m in reference_re.finditer(item_data):
            var = m.groups()[0]
            ref = self.get_value(var)

            if ref is not None:
                procval = self.get_processed_value(str(ref), outvar)
                item_data = item_data.replace(m.group(), procval)
            elif var == "CIMEROOT":
                item_data = item_data.replace(m.group(), self.cimepath)
            elif var == "SRCROOT":
                cpath = os.path.join(self.cimepath, "..")
                item_data = item_data.replace(m.group(), cpath)
            elif var == "USER":
                item_data = item_data.replace(m.group(), getpass.getuser())

        if math_re.search(item_data):
            try:
                tmp = eval(item_data)
            except Exception:
                tmp = item_data
            item_data = str(tmp)

        return item_data

    # get module info
    def get_modules(self, outvar, exclude_envs):

        modcmds: typing.List[str] = []
        self.__OMEGA_MODULE_COMMANDS__ = modcmds

        module_system_node = self.get_child("module_system")
        module_system_type = self.get(module_system_node, "type")

        if module_system_type != "module":
            print((f"ERROR: '{module_system_type}' "
                   "module system is not supported."))
            exit(-1)

        # get env. variables *before* applying module configurations
        # specified in config_machines.xml

        out1 = subprocess.check_output("env", shell=True)
        env1 = str(out1, 'UTF-8')

        module_nodes = self.get_children(
            "modules", root=module_system_node
        )

        modcmd = "module"

        shcmds = []

        # read module commands for the specified compiler and mpi library

        for module_node in module_nodes:
            compiler = self.get(module_node, "compiler")
            mpilib = self.get(module_node, "mpilib")
            debug = self.get(module_node, "DEBUG")

            if not (compiler is None or
                    re.match("^" + compiler + "$", self.compiler)):
                continue

            if not (mpilib is None or
                    re.match("^" + mpilib + "$", self.mpilib)):
                continue

            if not (debug is None or
                    re.match("^" + debug + "$", self.debug)):
                continue

            command_nodes = self.get_children("command",
                                              root=module_node)
            for command_node in command_nodes:
                name = self.get(command_node, "name")
                module = self.text(command_node)
                if module is None:
                    shcmd = f"{modcmd} {name}"
                else:
                    shcmd = f"{modcmd} {name} {module}"

                modcmds.append(shcmd)
                shcmds.append(shcmd)

        shcmds.append("env")

        # get env. variables *after* applying module configurations
        # specified in config_machines.xml

        out2 = subprocess.check_output(";".join(shcmds), shell=True)
        env2 = str(out2, 'UTF-8')

        parsed1 = {}
        parsed2 = {}

        for (name, value) in pat_envvar.findall(env1):
            parsed1[name] = value

        for (name, value) in pat_envvar.findall(env2):
            parsed2[name] = value

        for (name, value) in parsed2.items():
            if name in parsed1.keys():
                if parsed1[name] != value:
                    outvar[name] = value
            else:
                outvar[name] = value

        for (name, value) in parsed1.items():
            if name not in parsed2.keys():
                exclude_envs.append(name)

    # get environmental variables info
    def get_envs(self, outvar):

        exports: typing.Dict[str, str] = {}
        self.__OMEGA_SCRIPT_EXPORTS__ = exports

        envvar_nodes = self.get_children("environment_variables",
                                         root=self.machine_node)

        # possible attribs : compiler, DEBUG, SMP_PRESENT, mpilib
        for envvar_node in envvar_nodes:
            compiler = self.get(envvar_node, "compiler")
            mpilib = self.get(envvar_node, "mpilib")
            debug = self.get(envvar_node, "DEBUG")

            if not (compiler is None or
                    re.match("^" + compiler + "$", self.compiler)):
                continue

            if not (mpilib is None or
                    re.match("^" + mpilib + "$", self.mpilib)):
                continue

            if not (debug is None or
                    re.match("^" + debug + "$", self.debug)):
                continue

            env_nodes = self.get_children("env", root=envvar_node)
            for env_node in env_nodes:
                name = self.get(env_node, "name")
                value = self.get_processed_value(self.text(env_node).strip(),
                                                 outvar)
                if value is None:
                    outvar[name] = ""
                    exports[name] = outvar[name]
                elif value.startswith("$ENV{") and value.endswith("}"):
                    vname = value[5:-1]
                    if vname in outvar:
                        outvar[name] = outvar[vname]
                    elif vname in os.environ:
                        outvar[name] = os.environ[vname]
                    else:
                        outvar[name] = ""
                    exports[name] = outvar[name]
                elif value.startswith("$SHELL{") and value.endswith("}"):
                    print("Warning: SHELL evaluation is not supported.")

                else:
                    outvar[name] = value
                    exports[name] = outvar[name]

    # get mpirun info
    def get_mpirun(self, outvar):

        mpirun_node = self.get_child("mpirun")
        mpirun_mpilib = self.get(mpirun_node, "mpilib")

        exec_nodes = self.get_children(
            "executable", root=mpirun_node
        )

        if len(exec_nodes) != 1:
            print("ERROR: 'more than one executable nodes in mpirun node")
            exit(-1)

        self.mpilib = mpirun_mpilib
        if self.mpilib == "default":
            self.mpilib = self.mpilibs[0]

        self.mpiexec = self.text(exec_nodes[0])

    # collect machine info
    def gen_machinfo(self):

        outvar: typing.Dict[str, str] = {}
        exclude_envs: typing.List[str] = []

        self.get_mpirun(outvar)
        self.get_modules(outvar, exclude_envs)
        self.get_envs(outvar)
        self.write_output(outvar, exclude_envs)
        self.generate_scripts(outvar)

    # create a temporary cmake script to be included
    # in the main Omega cmake build system
    def write_output(self, outvar, exclude_envs):

        with open(self.outpath, "w") as f:
            f.write("message(STATUS \"Reading E3SM machine info\")\n")

            for key, value in outvar.items():
                if not key.startswith("__OMEGA_"):
                    f.write("set(ENV{%s} \"%s\")\n" % (key, value))

            for name in exclude_envs:
                f.write("unset(ENV{%s})\n" % name)

            f.write(f"set(MACH {self.machname})\n")
            f.write(f"set(OS {self.machos})\n")
            f.write(f"set(COMPILER {self.compiler})\n")
            f.write(f"set(MPI_EXEC {self.mpiexec})\n")
            f.write(f"set(CASEROOT {self.machpath})\n")
            f.write(f"include({self.macrospath})\n")

    # create scripts
    def generate_scripts(self, outvar):

        omega_env = os.path.join(os.path.dirname(self.outpath),
                                 "omega_env.sh")
        omega_build = os.path.join(os.path.dirname(self.outpath),
                                   "omega_build.sh")
        omega_run = os.path.join(os.path.dirname(self.outpath),
                                 "omega_run.sh")
        omega_ctest = os.path.join(os.path.dirname(self.outpath),
                                   "omega_ctest.sh")

        with open(omega_env, "w") as f:
            f.write("#!/usr/bin/env bash\n\n")

            f.write("# module commands\n")
            for cmd in self.__OMEGA_MODULE_COMMANDS__:
                f.write(cmd + "\n")

            f.write("\n# env. variables\n")
            for key, value in self.__OMEGA_SCRIPT_EXPORTS__.items():
                f.write(f"export {key}=\"{value}\"\n")

        with open(omega_build, "w") as f:
            f.write("#!/usr/bin/env bash\n\n")

            f.write("source ./omega_env.sh\n")
            nthreads_build = self.get_value("GMAKE_J")
            f.write(f"make -j{nthreads_build}\n")

        with open(omega_run, "w") as f:
            f.write("#!/usr/bin/env bash\n\n")

            f.write("source ./omega_env.sh\n")
            f.write("./src/omega.exe\n")

        with open(omega_ctest, "w") as f:
            f.write("#!/usr/bin/env bash\n\n")

            f.write("source ./omega_env.sh\n")
            f.write("ctest $* # --rerun-failed --output-on-failure\n")

        st = os.stat(omega_env)
        os.chmod(omega_env, st.st_mode | stat.S_IEXEC)

        st = os.stat(omega_build)
        os.chmod(omega_build, st.st_mode | stat.S_IEXEC)

        st = os.stat(omega_run)
        os.chmod(omega_run, st.st_mode | stat.S_IEXEC)

        st = os.stat(omega_ctest)
        os.chmod(omega_ctest, st.st_mode | stat.S_IEXEC)


def main():

    mach = OmegaMachines(args.cimepath, args.outpath, args.machine,
                         args.compiler, args.debug)
    mach.gen_machinfo()


if __name__ == "__main__":
    main()

import argparse
import os
import re
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
                        os.path.join(here, "..", "..", "cime")))    # CIME root
    parser.add_argument('-o', '--outpath', default="_Omega.cmake")  # outfile
    parser.add_argument('-c', '--compiler')     # compiler
    parser.add_argument('-v', '--verbose')      # verbose output

    # cimepath, outpath, compiler, verbose
    return parser.parse_args()


args = parse_cmdline()

sys.path.insert(0, args.cimepath)
from CIME.XML.machines import Machines  # noqa: E402

sys.path.pop(0)


class OmegaMachines(Machines):

    def __init__(self, cimepath, outpath, compiler):

        super(OmegaMachines, self).__init__()

        self.cimepath = cimepath
        self.outpath = outpath
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

    def get_modules(self, outvar):

        module_system_node = self.get_child("module_system")
        module_system_type = self.get(module_system_node, "type")

        # check  = (True if self.get(module_system_node,
        #           "allow_error") == "true" else False)

        if module_system_type != "module":
            print((f"ERROR: '{module_system_type}' "
                   "module system is not supported."))
            exit(-1)

        out1 = subprocess.check_output("env", shell=True)
        env1 = str(out1, 'UTF-8')

        module_nodes = self.get_children(
            "modules", root=module_system_node
        )

        modcmd = "module"

        shcmds = []

        for module_node in module_nodes:
            compiler = self.get(module_node, "compiler")
            if compiler is None or re.match(compiler, self.compiler):
                command_nodes = self.get_children("command",
                                                  root=module_node)
                for command_node in command_nodes:
                    name = self.get(command_node, "name")
                    module = self.text(command_node)
                    if module is None:
                        shcmds.append(f"{modcmd} {name}")
                    else:
                        shcmds.append(f"{modcmd} {name} {module}")

        shcmds.append("env")

        out2 = subprocess.check_output(";".join(shcmds), shell=True)
        env2 = str(out2, 'UTF-8')

        parsed1 = {}

        for (name, value) in pat_envvar.findall(env1):
            parsed1[name] = value

        for (name, value) in pat_envvar.findall(env2):
            if name in parsed1:
                if parsed1[name] != value:
                    outvar[name] = value
            else:
                outvar[name] = value

    def get_envs(self, outvar):

        envvar_nodes = self.get_children(
            root=self.get_child("environment_variables"))

        for envvar_node in envvar_nodes:
            compiler = self.get(envvar_node, "compiler")
            if compiler is None or re.match(compiler, self.compiler):
                env_nodes = self.get_children("env", root=envvar_node)
                for env_node in env_nodes:
                    name = self.get(env_node, "name")
                    value = self.text(env_node).strip()
                    if value is None:
                        outvar[name] = ""
                    elif value.startswith("$ENV{") and value.endswith("}"):
                        vname = value[5:-1]
                        if vname in outvar:
                            outvar[name] = outvar[vname]
                        elif vname in os.environ:
                            outvar[name] = os.environ[vname]
                        else:
                            outvar[name] = ""
                    else:
                        outvar[name] = value

    def get_mpirun(self, outvar):

        mpirun_node = self.get_child("mpirun")
        mpirun_mpilib = self.get(mpirun_node, "mpilib")

        if mpirun_mpilib != "default":
            print(f"ERROR: '{mpirun_mpilib}' mpilib is not supported.")
            exit(-1)

        exec_nodes = self.get_children(
            "executable", root=mpirun_node
        )

        if len(exec_nodes) != 1:
            print("ERROR: 'more than one executable nodes in mpirun node")
            exit(-1)

        self.mpiexec = self.text(exec_nodes[0])

    def gen_machinfo(self):

        outvar: typing.Dict[str, int] = {}

        self.get_modules(outvar)
        self.get_envs(outvar)
        self.get_mpirun(outvar)
        self.write_output(outvar)

    def write_output(self, outvar):

        with open(self.outpath, "w") as f:
            f.write("message(STATUS \"Reading E3SM machine info\")\n")

            for key, value in outvar.items():
                f.write("set(ENV{%s} \"%s\")\n" % (key, value))

            f.write(f"set(MACH {self.machname})\n")
            f.write(f"set(OS {self.machos})\n")
            f.write(f"set(COMPILER {self.compiler})\n")
            f.write(f"set(MPI_EXEC {self.mpiexec})\n")
            f.write(f"set(CASEROOT {self.machpath})\n")
            f.write(f"include({self.macrospath})\n")

            f.write(("message(STATUS \"End of reading E3SM "
                     "machine info\")\n"))


def main():

    mach = OmegaMachines(args.cimepath, args.outpath, args.compiler)
    mach.gen_machinfo()


if __name__ == "__main__":
    main()

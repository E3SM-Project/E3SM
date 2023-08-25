import os, sys, shutil, argparse, re, subprocess

# TODO: print basic info on screen and verbose info too

pat_envvar = re.compile(r'^([_\d\w]+)=(.*)$', flags=re.MULTILINE)

def parse_cmdline():

    here, progname = os.path.split(__file__)

    parser = argparse.ArgumentParser(
                    prog=progname,
                    description='generate machine specific infomation',
                    epilog='Contact: <T.B.D.>')

    parser.add_argument('-p', '--cimepath', default=os.path.realpath(os.path.join(
                            here, "..", "..", "cime")))             # CIME root path
    parser.add_argument('-o', '--outpath', default="_Omega.cmake")  # outfile path
    parser.add_argument('-c', '--compiler') # compiler
    parser.add_argument('-v', '--verbose') # verbose output

    args = parser.parse_args()
   
    return (args.cimepath, args.outpath, args.compiler, args.verbose)

def main():

    cimepath, outpath, compiler, verbose = parse_cmdline()

    sys.path.insert(0, cimepath)

    from CIME.XML.machines import Machines

    class OmegaMachines(Machines):

        def __init__(self, cimepath, outpath, compiler):

            super(OmegaMachines, self).__init__()

            self.cimepath = cimepath
            self.outpath = outpath
 
            if compiler is None:
                self.compiler = self.get_default_compiler().strip()
            else:
                self.compiler = compiler.strip()

            self.machpath = os.path.join(self.cimepath, "..", "cime_config", "machines")
            self.macrospath = os.path.join(self.machpath, "cmake_macros", "Macros.cmake")

            self.machname = self.get_machine_name()
            self.machos = self.get_value("OS")

        def get_modules(self, outvar):

            module_system_node = self.get_child("module_system") 
            module_system_type = self.get(module_system_node, "type")
            check  = (True if self.get(module_system_node,
                            "allow_error")=="true" else False)

            if module_system_type != "module":
                print("ERROR: '%s' module system is not supported." % module_system_type)
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
                    command_nodes = self.get_children("command", root=module_node)
                    for command_node in command_nodes:
                        name = self.get(command_node, "name")
                        module = self.text(command_node)
                        if module is None:
                            shcmds.append("%s %s" % (modcmd, name))
                        else:
                            shcmds.append("%s %s %s" % (modcmd, name, module))

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

            envvar_nodes = self.get_children(root=self.get_child("environment_variables"))

            for envvar_node in envvar_nodes:
                compiler = self.get(envvar_node, "compiler")
                if compiler is None or re.match(compiler, self.compiler):
                    env_nodes = self.get_children("env", root=envvar_node)
                    for env_node in env_nodes:
                        name = self.get(envnode, "name")
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

        def gen_machinfo(self):

            outvar = {}

            self.get_modules(outvar)
            self.get_envs(outvar)
            self.write_output(outvar)

        def write_output(self, outvar):

            with open(self.outpath, "w") as f:
                f.write("message(STATUS \"Reading E3SM machine info\")\n")

                for key, value in outvar.items():
                    f.write("set(ENV{%s} \"%s\")\n" % (key, value))
                    
                f.write("set(MACH %s)\n" % self.machname)
                f.write("set(OS %s)\n" % self.machos)
                f.write("set(COMPILER %s)\n" % self.compiler)
                f.write("set(CASEROOT %s)\n" % os.path.realpath(self.machpath))
                f.write("include(%s)\n" % os.path.realpath(self.macrospath))

                f.write("message(STATUS \"End of reading E3SM machine info\")\n")

    mach = OmegaMachines(cimepath, outpath, compiler)

    mach.gen_machinfo()


if __name__ == "__main__":
    main()

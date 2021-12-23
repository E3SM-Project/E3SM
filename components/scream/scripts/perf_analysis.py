from utils import run_cmd_no_fail, run_cmd, expect, median

import os, tempfile, re

###############################################################################
class ScalingExp(object):
###############################################################################

    def __init__(self, argmap, threads, arg_str, machine):
        try:
            self.varname, self.scale_factor, self.upper_limit = arg_str.split(":")
            self.scale_factor = float(self.scale_factor)
            self.upper_limit  = int(self.upper_limit)
        except Exception:
            expect(False, "Scaling experiment needs to be in format VARNAME:SCALE_FACTOR:MAX")

        self.args = list(argmap.keys())
        for name, val in argmap.items():
            setattr(self, name, val)

        self.threads = threads
        self.machine = machine

        expect(self.varname in dir(self), "Unknown varname '{}'".format(self.varname))

    def should_continue(self):
        return self.get_scaling_var() <= self.upper_limit

    def update_values(self):
        """
        >>> se = ScalingExp({'ni':10, 'nk:1'}, 1, 'ni:2:1000')
        >>> se.update_values()
        >>> se.ni
        20
        >>> se.nk
        1
        """
        setattr(self, self.varname, int(self.get_scaling_var() * self.scale_factor))

    def get_scaling_var(self):
        return getattr(self, self.varname)

    def values(self, incl_threads=True):
        results = [getattr(self, name) for name in self.args]
        if incl_threads:
            results.append(self.threads)
        return tuple(results)

    def plot(self, results):
        prov_msg = "Provenance: "
        item_list = [(name, getattr(self, name)) for name in self.args]
        item_list.append(("threads", self.threads))
        for name, val in item_list:
            if name != self.varname:
                prov_msg += " {}={}".format(name, val)

        st, out, _ = run_cmd("git rev-parse --short HEAD")
        git_commit = out if st == 0 else "Unknown"
        print("{} machine={} commit={}".format(prov_msg, self.machine, git_commit))

        for test_name, test_results in results.items():
            print(test_name, self.varname)
            for cols, med_time, scaling_var in test_results:
                cols_sec = float(cols) / med_time
                print("{}, {:.2f}".format(scaling_var, cols_sec))

###############################################################################
class PerfAnalysis(object):
###############################################################################

    ###########################################################################
    def __init__(self, argmap, force_threads, num_runs, tests, cmake_options, use_existing, scaling_exp, plot_friendly, machine, scream_docs, verbose, cd):
    ###########################################################################
        self._argmap        = argmap
        self._force_threads = force_threads
        self._num_runs      = num_runs
        self._tests         = tests
        self._cmake_options = cmake_options
        self._use_existing  = use_existing
        self._scaling_exp   = scaling_exp
        self._plot_friendly = plot_friendly
        self._machine       = machine
        self._scream_docs   = scream_docs
        self._verbose       = verbose
        self._cd            = cd

    ###############################################################################
    def build(self):
    ###############################################################################
        cmake_cmd = "cmake {} ..".format(self._cmake_options)

        if self._verbose:
            print("In dir {}, building with cmake command: {}\nOutput will be stored in build.perf.log".\
                  format(os.getcwd(), cmake_cmd))

        with open("build.perf.log", "w", encoding="utf-8") as fd:

            make_cmd  = "make -j16 VERBOSE=1"
            fd.write(cmake_cmd + "\n")
            fd.write(run_cmd_no_fail(cmake_cmd, combine_output=True) + "\n\n")
            fd.write(make_cmd + "\n")
            fd.write(run_cmd_no_fail(make_cmd, combine_output=True) + "\n")

    ###############################################################################
    def get_time(self, output):
    ###############################################################################
        r"""
        >>> output = 'Foo\nTime = 0.047 seconds.\nbar'
        >>> get_time(output)
        0.047
        >>> output = 'Foo\nTime = 1.732e+01 seconds.\nbar'
        >>> get_time(output)
        17.32
        """
        regex = re.compile(r'Time\s*=\s*([^\s]+)\s*seconds')
        the_time = None
        for line in output.splitlines():
            m = regex.match(line)
            if m:
                expect(the_time is None, "Multiple matches!")
                the_time = float(m.groups()[0])

        return the_time

    ###############################################################################
    def get_threads(self, output):
    ###############################################################################
        r"""
        >>> output = 'Foo\nARCH: dp 1 avx  FPE 0 nthread 48\nTime = 0.047 seconds.\nbar'
        >>> get_threads(output)
        48
        """
        for line in output.splitlines():
            if "nthread" in line:
                items = line.split()
                threads = int(items[items.index("nthread") + 1])
                return threads

        if "OMP_NUM_THREADS" in os.environ:
            return int(os.environ["OMP_NUM_THREADS"])
        else:
            return 1

    ###############################################################################
    def formulate_cmd(self, test_exe):
    ###############################################################################
        prefix = "" if "NUMA_PREFIX" not in os.environ else "{} ".format(os.environ["NUMA_PREFIX"])
        replaced = []
        for name, val in zip(self._argmap.keys(), self._scaling_exp.values(incl_threads=False)):
            if name.upper() in test_exe:
                test_exe = test_exe.replace(name.upper(), str(val))
                replaced.append(name)

        for name, val in zip(self._argmap.keys(), self._scaling_exp.values(incl_threads=False)):
            if name not in replaced:
                test_exe += " {}".format(val)

        cmd = "{}./{}".format(prefix, test_exe)
        return cmd

    ###############################################################################
    def run_test(self, test_cmd):
    ###############################################################################
        if self._cd:
            test_path, test_exe = os.path.split(test_cmd)
            test_path = None if not test_path else test_path
        else:
            test_exe = test_cmd
            test_path = None

        self.machine_specific_init(self._scaling_exp.threads)
        self.test_specific_init(test_exe, self._scaling_exp.threads)

        cmd = self.formulate_cmd(test_exe)
        results = []
        with open("{}.perf.log".format(os.path.split(test_exe)[1].split(" ")[0]), "w", encoding="utf-8") as fd:
            fd.write(cmd + "\n\n")
            fd.write("ENV: \n{}\n\n".format(run_cmd_no_fail("env")))
            for _ in range(self._num_runs):
                output = run_cmd_no_fail(cmd, from_dir=test_path, verbose=(not self._plot_friendly or self._verbose))
                fd.write(output + "\n\n")
                results.append(self.get_time(output))

            threads = self.get_threads(output)

        return median(results), threads

    ###############################################################################
    def user_explain(self, test, cols, med_time, reference, threads):
    ###############################################################################
        msg = "{} ran in {} seconds with {} threads, {:.2f} cols/sec".format(test, med_time, threads, float(cols)/med_time)
        if reference:
            speedup = (1.0 - (med_time / reference)) * 100
            msg += ", speedup={:.2f}%".format(speedup)

        print(msg)

    ###############################################################################
    def perf_analysis(self):
    ###############################################################################
        if self._use_existing:
            expect(os.path.exists("CMakeCache.txt"),
                   "{} doesn't look like a build directory".format(os.getcwd()))

        else:
            if self._scream_docs:
                expect(os.path.basename(os.getcwd()) == "micro-apps", "Please run from micro-apps directory")

            tmpdir = tempfile.mkdtemp(prefix="build", dir=os.getcwd())
            os.chdir(tmpdir)

            if not self._plot_friendly:
                print("BUILDING")

            self.build()

        results = {}
        while (self._scaling_exp.should_continue()):
            if not self._plot_friendly:
                print()
                print("RUNNING {}".format(" ".join(["{}={}".format(name, val) for name, val in zip(self._argmap.keys(), self._scaling_exp.values(incl_threads=False))])))

            reference = None
            for test, test_cmd in self._tests.items():
                med_time, threads = self.run_test(test_cmd)
                self._scaling_exp.threads = threads

                if self._plot_friendly:
                    results.setdefault(test, []).append((self._scaling_exp.values()[0], med_time, self._scaling_exp.get_scaling_var()))
                else:
                    self.user_explain(test, self._scaling_exp.values()[0], med_time, reference, threads)

                reference = med_time if reference is None else reference

            self._scaling_exp.update_values()

        if self._plot_friendly:
            self._scaling_exp.plot(results)

        return True

    ###############################################################################
    def test_specific_init(self, exename, force_threads):
    ###############################################################################
        # This appears to be slower with 48
        #if self._machine == "blake" and exename == "p3_ref":
        #    os.environ["OMP_NUM_THREADS"] = "48"

        if force_threads:
            os.environ["OMP_NUM_THREADS"] = str(force_threads)

    ###############################################################################
    def machine_specific_init(self, force_threads=None):
    ###############################################################################
        force_threads = self._force_threads if force_threads is None else force_threads

        if self._machine == "bowman":
            os.environ["OMP_NUM_THREADS"] = "272"
            os.environ["NUMA_PREFIX"] = "numactl -i 1"
        elif self._machine == "blake":
            os.environ["OMP_NUM_THREADS"] = "96"
        elif self._machine in ["white", "weaver", "melvin"]:
            pass
        elif self._machine == "mappy":
            os.environ["OMP_NUM_THREADS"] = "46"
        else:
            print("WARNING: Unrecognized machine {}".format(self._machine))

        if force_threads:
            os.environ["OMP_NUM_THREADS"] = str(force_threads)

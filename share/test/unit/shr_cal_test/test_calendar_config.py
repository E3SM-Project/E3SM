"""Run with python3; validate CIME calendar configurations without model cases."""

import ast
import importlib.machinery
import importlib.util
import datetime
import gzip
import inspect
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "cime"))
sys.path.insert(0, str(ROOT / "cime" / "CIME" / "Tools"))
os.environ.setdefault("CIME_MODEL", "e3sm")

from CIME.nmlgen import NamelistGenerator
from CIME.get_tests import get_test_data
from CIME.SystemTests.test_mods import find_test_mods
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import find_system_test, get_full_test_name, parse_test_name


class CalendarConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.drivers = {}
        for driver in ("driver-mct", "driver-moab"):
            path = ROOT / driver / "cime_config" / "buildnml"
            loader = importlib.machinery.SourceFileLoader(driver, str(path))
            spec = importlib.util.spec_from_loader(loader.name, loader)
            module = importlib.util.module_from_spec(spec)
            loader.exec_module(module)
            cls.drivers[driver] = module

    def generator(self, driver, calendar, profile=None, user_namelist=None):
        case = Mock()
        case.get_value.side_effect = {"DIN_LOC_ROOT": str(ROOT)}.get
        definition = ROOT / driver / "cime_config" / "namelist_definition_drv.xml"
        generator = NamelistGenerator(case, [str(definition)])
        infiles = []
        if profile is not None:
            infiles.append(str(ROOT / "cime_config" / "testmods_dirs" / "allactive" /
                               ("calendar_" + profile) / "user_nl_cpl"))
        if user_namelist is not None:
            infiles.append(str(user_namelist))
        generator.init_defaults(infiles, {"calendar": calendar}, skip_entry_loop=True)
        for name in ("orb_mode", "orb_iyear", "orb_iyear_align", "orb_use_equation_of_time",
                 "fudge_noleap_for_gregorian", "fudge_gregorian_to_noleap"):
            generator.add_default(name)
        return generator

    def test_epoch_defaults_and_eot_default(self):
        policies = {ast.dump(ast.parse(inspect.getsource(module._set_orbital_defaults)))
                    for module in self.drivers.values()}
        self.assertEqual(len(policies), 1, "MCT/MOAB orbital policy drift; update both buildnml helpers")
        for driver, module in self.drivers.items():
            for calendar in ("NO_LEAP", "GREGORIAN"):
                for token, expected in (("1850", "1850"), ("2010", "2010"), ("20TR", "1850")):
                    with self.subTest(driver=driver, calendar=calendar, token=token):
                        generator = self.generator(driver, calendar)
                        module._set_orbital_defaults(generator, token)
                        self.assertEqual(generator.get_value("orb_iyear"), expected)
                        self.assertEqual(generator.get_value("orb_iyear_align"), expected)
                        self.assertEqual(generator.get_value("orb_mode"), "fixed_year")
                        self.assertEqual(generator.get_value("orb_use_equation_of_time"), ".true.")

    def test_profiles(self):
        for driver, module in self.drivers.items():
            for profile in ("nf0", "nf1", "nv0", "nv1", "gf0", "gf1", "gv0", "gv1"):
                with self.subTest(driver=driver, profile=profile):
                    calendar = "GREGORIAN" if profile.startswith("g") else "NO_LEAP"
                    generator = self.generator(driver, calendar, profile)
                    module._set_orbital_defaults(generator, "1850")
                    variable = profile[1] == "v"
                    self.assertEqual(generator.get_value("orb_mode"),
                                     "variable_orbit_year" if variable else "fixed_year")
                    self.assertEqual(generator.get_value("orb_iyear"), "1999" if variable else "1850")
                    self.assertEqual(generator.get_value("orb_iyear_align"), "1999" if variable else "1850")
                    self.assertEqual(generator.get_value("orb_use_equation_of_time"),
                                     ".true." if profile.endswith("1") else ".false.")

    def test_explicit_epoch_and_legacy_alias(self):
        for driver, module in self.drivers.items():
            for calendar in ("NO_LEAP", "GREGORIAN"):
                with self.subTest(driver=driver, calendar=calendar):
                    generator = self.generator(driver, calendar)
                    generator.set_value("orb_mode", "variable_year")
                    generator.set_value("orb_iyear", "1950")
                    generator.set_value("orb_iyear_align", "2000")
                    module._set_orbital_defaults(generator, "1850")
                    self.assertEqual(generator.get_value("orb_mode"), "variable_year")
                    self.assertEqual(generator.get_value("orb_iyear"), "1950")
                    self.assertEqual(generator.get_value("orb_iyear_align"), "2000")

    def test_suite_modifications(self):
        for suffix, expected in (("smoke", 16), ("short", 128), ("cycle", 64)):
            tests = get_test_data("e3sm_calendar_" + suffix)[-1]
            self.assertEqual(len(tests), expected)
            self.assertEqual(len(set(tests)), expected)
            for test in tests:
                with self.subTest(suite=suffix, test=test):
                    test_name, modifiers = test.rsplit(".", 1)
                    full_name = get_full_test_name(
                        test_name, machine="pm-cpu", compiler="gnu", testmods_string=modifiers)
                    parsed = parse_test_name(full_name)
                    paths = find_test_mods("mct", parsed[6])
                    self.assertEqual(len(paths), 2)
                    for path in paths:
                        self.assertTrue((Path(path) / "shell_commands").is_file())

    def test_restart_pointers(self):
        for calendar in ("GREGORIAN", "NO_LEAP"):
            cycle_restart = ".2000-12-31-00000" if calendar == "GREGORIAN" else ".2001-01-01-00000"
            for start, stop, interval, expected in (
                    ("1999-01-01", 1461 if calendar == "GREGORIAN" else 1460, 730, cycle_restart),
                    ("2000-02-27", 6, 1, ".2000-02-28-00000"),
                    ("1999-12-29", 6, 1, ".1999-12-30-00000")):
                with self.subTest(calendar=calendar, start=start):
                    values = dict(STOP_N=stop, STOP_OPTION="ndays", COMP_ATM="datm",
                                  ATM_NCPL=24, COMP_INTERFACE="mct", NCPL_BASE_PERIOD="day",
                                  CALENDAR=calendar, START_TOD=0, RUN_STARTDATE=start)
                    case = Mock()
                    case.get_value.side_effect = values.get
                    case.get_values.return_value = ["CPL", "ATM"]
                    context = SimpleNamespace(
                        _case=case, _rest_n=interval,
                        _leap_year_correction=SystemTestsCommon._leap_year_correction)
                    SystemTestsCommon._set_restart_interval(context)
                    self.assertEqual(context._rest_time, expected)

    def compile_orbit_probe(self, program=None, extra_sources=(), link_args=(), expected_returncode=0):
        stubs = """module shr_kind_mod
implicit none
integer,parameter::SHR_KIND_R8=kind(1d0),SHR_KIND_IN=kind(1),SHR_KIND_I8=selected_int_kind(18)
end module
module shr_const_mod
implicit none
real(8),parameter::SHR_CONST_PI=acos(-1d0)
end module
module shr_log_mod
implicit none
integer::shr_log_Level=0,shr_log_Unit=6
end module
module shr_sys_mod
implicit none
contains
subroutine shr_sys_abort(message)
character(*),optional::message
if(present(message)) print *,message
stop 99
end subroutine
end module
"""
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            flags = [shutil.which("gfortran"), "-cpp", "-ffree-line-length-none", "-fcheck=all",
                     "-ffpe-trap=invalid,zero,overflow", "-Werror=implicit-interface", "-J", directory, "-I", directory]
            subprocess.run(flags + ["-x", "f95", "-c", "-", "-o", str(folder / "stubs.o")],
                           input=stubs, text=True, cwd=folder, capture_output=True, check=True)
            sources = [str(ROOT / "share/util/shr_orb_mod.F90")] + [str(ROOT / path) for path in extra_sources]
            if program is None:
                sources.append(str(ROOT / "share/test/old_unit_testers/test_shr_orb.F90"))
            else:
                sources += ["-x", "f95", "-", "-x", "none"]
            compiled = subprocess.run(flags + sources + [str(folder / "stubs.o"), *link_args, "-o", str(folder / "probe")],
                                      input=program, text=True, cwd=folder, capture_output=True)
            self.assertEqual(compiled.returncode, 0, compiled.stderr)
            executed = subprocess.run([str(folder / "probe")], cwd=folder, capture_output=True, text=True)
            self.assertEqual(executed.returncode, expected_returncode, executed.stdout + executed.stderr)
            return executed


    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for stub-backed reader checks")
    def test_restart_character_readers(self):
        header = """module reader_probe
use shr_sys_mod, only: shr_sys_abort
implicit none
integer,parameter::CL=256,logunit=6
character(*),parameter::prefix='seq_io_'
type file_desc_t
integer::value
end type
type var_desc_t
integer::value
end type
contains
"""
        stubs = """
subroutine seq_io_read_openfile(filename,pioid,addprefix)
character(*)::filename
type(file_desc_t)::pioid
logical::addprefix
addprefix=.false.
end subroutine
integer function pio_inq_varid(pioid,name,varid)
type(file_desc_t)::pioid
character(*)::name
type(var_desc_t)::varid
pio_inq_varid=0
end function
integer function pio_get_var(pioid,varid,text)
type(file_desc_t)::pioid
type(var_desc_t)::varid
character(*)::text
text=repeat('x',len(text))
pio_get_var=0
end function
end module
program probe
use reader_probe
implicit none
type(file_desc_t)::file
character(16)::short
character(128*(1024+66))::long
call seq_io_read_char('',file,short,'')
call seq_io_read_char('',file,long,'')
if(short/=repeat('x',len(short)).or.long/=repeat('x',len(long))) stop 1
end program
"""
        for driver in ("driver-mct", "driver-moab"):
            with self.subTest(driver=driver):
                source = (ROOT / driver / "shr/seq_io_read_mod.F90").read_text()
                routine = re.search(r"(?is)  subroutine seq_io_read_char\(.*?end subroutine seq_io_read_char", source)
                self.assertIsNotNone(routine)
                opening = source.lower().split("subroutine seq_io_read_openfile(", 1)[1].split(
                    "end subroutine seq_io_read_openfile", 1)[0]
                self.assertNotIn("len(rdata)", opening)
                self.compile_orbit_probe(header + routine.group() + stubs)


    @unittest.skipUnless(shutil.which("gfortran") and shutil.which("clang++"),
                         "gfortran and clang++ are required for linked callback checks")
    def test_input_calendar_cpp_callback(self):
        header = (ROOT / "components/eamxx/src/share/core/eamxx_config.hpp").read_text()
        source = (ROOT / "components/eamxx/src/share/core/eamxx_config.cpp").read_text()
        resolver_type = re.search(r"using input_calendar_resolver = .*?;", header).group()
        seconds_type = re.search(r"using input_calendar_seconds = .*?;", header).group()
        seconds = re.search(r"double gregorian_input_seconds\(.*?\n\}", source, re.S).group()
        check_range = re.search(r"void check_input_calendar_range\(.*?\n\}", source, re.S).group()
        resolver = re.search(r"std::string resolve_input_calendar\(.*?\n\}", source, re.S).group()
        program = """#include <string>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdexcept>
#define EKAT_REQUIRE_MSG(condition,message) if (condition) {} else throw std::runtime_error(message)
#define EKAT_ERROR_MSG(message) throw std::runtime_error(message)
extern "C" int input_calendar_resolve_c(const char*,int,const char*,int,int,char*,char*,char*,char*);
bool is_scream_standalone() {return true;}
""" + resolver_type + seconds_type + """
input_calendar_resolver calendar_resolver=input_calendar_resolve_c;
input_calendar_seconds calendar_seconds=nullptr;
""" + seconds + check_range + resolver + """
int main() {
    std::string provenance,flavor;
    for (int callback=0; callback<2; ++callback) {
        calendar_resolver=callback ? input_calendar_resolve_c : nullptr;
        if (resolve_input_calendar("a.nc","proleptic_gregorian",true,provenance,flavor)!="GREGORIAN" ||
                flavor!="PROLEPTIC_GREGORIAN" || provenance!="metadata") return 1;
        if (resolve_input_calendar("a.nc","standard",true,provenance,flavor)!="GREGORIAN" ||
                flavor!="STANDARD") return 2;
        check_input_calendar_range(flavor,20000101,0,0,"modern");
        check_input_calendar_range("PROLEPTIC_GREGORIAN",10101,0,730181.*86400,"historical");
        int rejected=0;
        try {check_input_calendar_range(flavor,10101,0,730181.*86400,"ancient epoch");}
        catch (const std::runtime_error&) {++rejected;}
        try {check_input_calendar_range(flavor,15821015,0,-1,"reform crossing");}
        catch (const std::runtime_error&) {++rejected;}
        if (rejected!=2) return 3;
    }
}
"""
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            subprocess.run([shutil.which("clang++"), "-std=c++17", "-Wall", "-Wextra", "-Werror", "-c",
                            "-x", "c++", "-", "-o", str(folder / "callback.o")], input=program,
                           text=True, capture_output=True, cwd=folder, check=True)
            runtime = ["-lc++"] if sys.platform == "darwin" else ["-lstdc++"]
            subprocess.run([shutil.which("gfortran"), "-ffree-line-length-none", "-fcheck=all", "-J", directory,
                            str(ROOT / "share/util/shr_input_calendar_mod.F90"), str(folder / "callback.o"),
                            *runtime, "-o", str(folder / "probe")],
                           text=True, capture_output=True, cwd=folder, check=True)
            subprocess.run([str(folder / "probe")], text=True, capture_output=True, cwd=folder, check=True)


    def test_noleap_input_defaults(self):
        for driver in self.drivers:
            for calendar in ("NO_LEAP", "GREGORIAN"):
                with self.subTest(driver=driver, calendar=calendar):
                    generator = self.generator(driver, calendar)
                    self.assertEqual(generator.get_value("fudge_noleap_for_gregorian"),
                                     ".false.")
                    self.assertEqual(generator.get_value("fudge_gregorian_to_noleap"),
                                     ".false.")


    def test_gregorian_input_user_overrides(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "user_nl_cpl"
            for driver in self.drivers:
                for enabled in (".true.", ".false."):
                    with self.subTest(driver=driver, enabled=enabled):
                        path.write_text("fudge_gregorian_to_noleap=" + enabled + "\n")
                        generator = self.generator(driver, "NO_LEAP", user_namelist=path)
                        self.assertEqual(generator.get_value("fudge_gregorian_to_noleap"), enabled)
                        self.assertEqual(generator.get_value("fudge_noleap_for_gregorian"), ".false.")


    def test_noleap_input_user_overrides(self):
        with tempfile.TemporaryDirectory() as directory:
            user_namelist = Path(directory) / "user_nl_cpl"
            for driver in self.drivers:
                for enabled in (".true.", ".false."):
                    with self.subTest(driver=driver, enabled=enabled):
                        user_namelist.write_text(
                            "fudge_noleap_for_gregorian = " + enabled + "\n")
                        generator = self.generator(
                            driver, "GREGORIAN", user_namelist=user_namelist)
                        self.assertEqual(generator.get_value("fudge_noleap_for_gregorian"),
                                         enabled)


    def test_explicit_source_calendar_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "user_nl_cpl"
            path.write_text("fudge_noleap_for_gregorian = .true.\n"
                            "input_calendar_override_selectors = 'file:/forcing/a.nc', 'stream:forcing.txt'\n"
                            "input_calendar_override_calendars = 'noleap', 'gregorian'\n")
            for driver in self.drivers:
                with self.subTest(driver=driver):
                    generator = self.generator(driver, "GREGORIAN", user_namelist=path)
                    self.assertEqual(generator.get_value("input_calendar_override_selectors"),
                                     ["file:/forcing/a.nc", "stream:forcing.txt"])
                    self.assertEqual(generator.get_value("input_calendar_override_calendars"),
                                     ["noleap", "gregorian"])


    def test_noleap_input_suites(self):
        for suffix, expected in (("smoke", 2), ("short", 56), ("negative", 11), ("cycle", 2)):
            tests = get_test_data("e3sm_calendar_noleap_" + suffix)[-1]
            self.assertEqual(len(tests), expected)
            self.assertEqual(len(set(tests)), expected)
            for test in tests:
                with self.subTest(suite=suffix, test=test):
                    test_name, modifiers = test.rsplit(".", 1)
                    full_name = get_full_test_name(
                        test_name, machine="pm-cpu", compiler="gnu", testmods_string=modifiers)
                    parsed = parse_test_name(full_name)
                    paths = find_test_mods("mct", parsed[6])
                    self.assertEqual(len(paths), 3)
                    for path in paths:
                        self.assertTrue((Path(path) / "shell_commands").is_file())
                    self.assertTrue(test.startswith("NLN_" if suffix == "negative" else "NLF_"))


    def test_gregorian_input_suites(self):
        for suffix, expected in (("smoke", 2), ("short", 56), ("negative", 18), ("cycle", 2)):
            tests = get_test_data("e3sm_calendar_gregorian_to_noleap_" + suffix)[-1]
            self.assertEqual(len(tests), expected)
            self.assertEqual(len(set(tests)), expected)
            for test in tests:
                with self.subTest(suite=suffix, test=test):
                    test_name, modifiers = test.rsplit(".", 1)
                    full_name = get_full_test_name(test_name, machine="pm-cpu", compiler="gnu", testmods_string=modifiers)
                    paths = find_test_mods("mct", parse_test_name(full_name)[6])
                    self.assertEqual(len(paths), 3)
                    self.assertIn("calendar_nf1", modifiers)
                    for path in paths:
                        self.assertTrue((Path(path) / "shell_commands").is_file())
                    self.assertTrue(test.startswith("NLN_" if suffix == "negative" else "NLF_"))
                    if suffix == "cycle":
                        self.assertTrue(test.startswith("NLF_Ld1460."))


    def test_noleap_system_test_loader(self):
        case = Mock()
        case.get_compset_components.return_value = ["datm", "docn", "dice"]
        directories = {"any": str(ROOT / "cime/CIME/SystemTests"),
                       "allactive": str(ROOT / "cime_config/SystemTests")}
        case.get_value.side_effect = lambda name, attribute: directories.get(attribute["component"])
        for name in ("NLF", "NLN"):
            with self.subTest(test=name):
                test_class = find_system_test(name, case)
                self.assertEqual(test_class.__name__, name)
                self.assertTrue(issubclass(test_class, SystemTestsCommon))



if __name__ == "__main__":
    unittest.main()
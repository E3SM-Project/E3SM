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

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for stub-backed numerical checks")
    def test_fixed_year_bridge_parity(self):
        program = """program probe
use iso_c_binding
use shr_orb_mod, only: shr_orb_decl_fixed_year
use shr_orb_mod_c2f, only: shr_orb_decl_fixed_year_c2f
implicit none
real(c_double)::delta,factor,eot,expected_delta,expected_factor,expected_eot,day
integer::current,previous,sample,enabled
do current=365,366
do previous=365,366
do sample=1,1464
day=1._c_double+real(sample-1,c_double)/4._c_double
do enabled=0,1
call shr_orb_decl_fixed_year_c2f(day,real(current,c_double),real(previous,c_double), &
    0.0167_c_double,1.8_c_double,0._c_double,0.409_c_double,delta,factor,eot,logical(enabled==1,c_bool))
call shr_orb_decl_fixed_year(day,real(current,c_double),real(previous,c_double), &
     0.0167_c_double,1.8_c_double,0._c_double,0.409_c_double,expected_delta,expected_factor,expected_eot)
if(enabled==0) expected_eot=0._c_double
if(delta/=expected_delta.or.factor/=expected_factor.or.eot/=expected_eot) stop 1
end do
end do
end do
end do
end program
"""
        self.compile_orbit_probe(program, ["components/eamxx/src/physics/rrtmgp/shr_orb_mod_c2f.F90"])


    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for stub-backed numerical checks")
    def test_ocean_albedo_scheduled_eot(self):
        for driver in ("driver-mct", "driver-moab"):
            with self.subTest(driver=driver):
                source = (ROOT / driver / "main/seq_flux_mct.F90").read_text()
                block = re.search(r"          calday = nextsw_cday\n.*?          ! Compute albedos", source, re.S)
                self.assertIsNotNone(block)
                cosine = re.search(r"^\s*cosz = shr_orb_cosz\([^\n]*eot_seconds[^\n]*\)", source, re.M)
                self.assertIsNotNone(cosine)
                program = """program probe
use shr_kind_mod, only: r8=>SHR_KIND_R8
use shr_const_mod, only: pi=>SHR_CONST_PI
use shr_sys_mod, only: shr_sys_abort
use shr_orb_mod
implicit none
character(*),parameter::seq_timemgr_gregorian='GREGORIAN'
character(16)::calendar
integer,parameter::infodata=0,nextsw_ymd=20001103
real(r8)::calday,nextsw_cday,nextsw_tod,delta,eccf,eot_seconds,cosz,reference,nextsw_dt_avg,sample_tod
real(r8)::rlat=0._r8,rlon=0._r8,eccen=0._r8,mvelpp=0._r8,lambm0=0._r8,obliqr=0._r8
logical::nextsw_valid,use_eot,context_called
integer::model,enabled,sign,edge,window,sample
do model=1,2
calendar='NO_LEAP'
if(model==2) calendar='GREGORIAN'
do enabled=0,1
use_eot=enabled==1
do sign=-1,1,2
do edge=1,2
do window=0,1
nextsw_dt_avg=real(window,r8)*1800._r8
nextsw_tod=21300._r8
if(edge==2) nextsw_tod=65100._r8
nextsw_cday=200._r8+nextsw_tod/86400._r8
nextsw_valid=model==2.or.use_eot
context_called=.false.
""" + block.group().replace("          ! Compute albedos", "") + cosine.group() + """
if(context_called.neqv.nextsw_valid) stop 1
reference=-cos(2._r8*pi*(nextsw_tod+real(enabled*sign,r8)*960._r8)/86400._r8)
if(window==1) then
reference=0._r8
do sample=1,12000
sample_tod=nextsw_tod+(real(sample,r8)-0.5_r8)*nextsw_dt_avg/12000._r8
reference=reference+max(0._r8,-cos(2._r8*pi*(sample_tod+real(enabled*sign,r8)*960._r8)/86400._r8))/12000._r8
end do
end if
if(abs(cosz-reference)>2.e-8_r8) stop 2
if(.not.nextsw_valid.and.calday/=nextsw_cday) stop 3
end do
end do
end do
end do
end do
contains
function seq_timemgr_GetCalendar() result(value)
character(16)::value
value=calendar
end function
logical function seq_infodata_UseEquationOfTime(context)
integer,intent(in)::context
seq_infodata_UseEquationOfTime=use_eot
end function
subroutine seq_infodata_SolarDecl(context,date,tod,cal,declination,factor,eot)
integer,intent(in)::context,date
real(r8),intent(in)::tod
character(*),intent(in)::cal
real(r8),intent(out)::declination,factor,eot
context_called=.true.
if(date/=nextsw_ymd.or.tod/=nextsw_tod.or.cal/=calendar) stop 4
declination=0._r8
factor=1._r8
eot=real(enabled*sign,r8)*960._r8
end subroutine
end program
"""
                self.compile_orbit_probe(program)


    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for compiled ocean checks")
    def test_ocean_interval_flux_closure(self):
        for driver in ("driver-mct", "driver-moab"):
            with self.subTest(driver=driver):
                source = (ROOT / driver / "main/seq_flux_mct.F90").read_text()
                if driver == "driver-mct":
                    flux = re.search(r"       do n=1,nloc_o\n          avsdr =.*?       end do", source, re.S)
                    declarations = """integer,parameter::nloc_o=1
integer,parameter::index_xao_So_anidr=1,index_xao_So_anidf=2,index_xao_So_avsdr=3,index_xao_So_avsdf=4
integer,parameter::index_xao_Faox_swdn=5,index_xao_Faox_swup=6
integer,parameter::index_a2x_Faxa_swndr=1,index_a2x_Faxa_swndf=2,index_a2x_Faxa_swvdr=3,index_a2x_Faxa_swvdf=4
type vector
real(r8)::rAttr(6,1)
end type
type(vector)::xao_o,a2x_o
real(r8)::anidr,anidf,avsdr,avsdf,swupc,swdnc
"""
                    setup = "xao_o%rAttr(1:4,1)=saved\na2x_o%rAttr(1:4,1)=incoming\n"
                    result = "down=xao_o%rAttr(5,1)\nup=xao_o%rAttr(6,1)\n"
                    albedo_declaration = "real(r8)::anidr,anidf,avsdr,avsdf\n"
                    albedo_result = "values=[anidr,anidf,avsdr,avsdf]\n"
                else:
                    flux = re.search(r"          do n=1,nloc_of\n             avsdr1 =.*?          end do", source, re.S)
                    declarations = """integer,parameter::nloc_of=1
real(r8)::anidr(1),anidf(1),avsdr(1),avsdf(1),swdnc(1),swupc(1)
real(r8)::swndr(1),swndf(1),swvdr(1),swvdf(1),anidr1,anidf1,avsdr1,avsdf1,swdnc1,swupc1
"""
                    setup = "anidr=saved(1)\nanidf=saved(2)\navsdr=saved(3)\navsdf=saved(4)\n"
                    setup += "swndr=incoming(1)\nswndf=incoming(2)\nswvdr=incoming(3)\nswvdf=incoming(4)\n"
                    result = "down=swdnc(1)\nup=swupc(1)\n"
                    albedo_declaration = "real(r8)::anidr(1),anidf(1),avsdr(1),avsdf(1)\ninteger,parameter::n=1\n"
                    albedo_result = "values=[anidr(1),anidf(1),avsdr(1),avsdf(1)]\n"
                albedo = re.search(r"             if \(cosz  >  0.0_r8\).*?             end if", source, re.S)
                self.assertIsNotNone(flux)
                self.assertIsNotNone(albedo)
                self.assertLess(flux.start(), albedo.start())
                program = """program probe
use shr_kind_mod,only:r8=>SHR_KIND_R8
use shr_orb_mod,only:shr_orb_cosz
implicit none
real(r8)::day,cosz,next_cosz,saved(4),following(4),incoming(4),down,up,absorbed
integer::edge,band
do edge=1,3
day=80._r8+(5._r8+50._r8/60._r8)/24._r8
if(edge==2) day=80._r8+(17._r8+50._r8/60._r8)/24._r8
if(edge==3) day=80._r8
cosz=shr_orb_cosz(day,0._r8,0._r8,0._r8,1800._r8)
call ocean_albedo(cosz,saved)
if(edge==1.and.shr_orb_cosz(day,0._r8,0._r8,0._r8)>=0._r8) stop 1
if(edge<=2.and.any(saved>=1._r8)) stop 2
if(edge==3.and.any(saved/=1._r8)) stop 3
next_cosz=shr_orb_cosz(day+1800._r8/86400._r8,0._r8,0._r8,0._r8,1800._r8)
call ocean_albedo(next_cosz,following)
if(edge==1.and.all(saved==following)) stop 11
do band=1,4
incoming=0._r8
incoming(band)=100._r8*max(0._r8,cosz)
call ocean_exchange(incoming,saved,down,up)
if(abs(down-sum(incoming))>1.e-12_r8) stop 4
if(abs(up+sum(incoming*saved))>1.e-12_r8) stop 5
absorbed=sum(incoming*(1._r8-saved))
if(abs(down+up-absorbed)>1.e-12_r8) stop 6
if(edge<=2.and.absorbed<=0._r8) stop 7
if(edge==3.and.(down/=0._r8.or.up/=0._r8)) stop 8
end do
incoming=[10._r8,20._r8,30._r8,40._r8]*max(0._r8,next_cosz)
call ocean_exchange(incoming,following,down,up)
if(abs(down-sum(incoming))>1.e-12_r8) stop 12
if(abs(up+sum(incoming*following))>1.e-12_r8) stop 13
if(abs(down+up-sum(incoming*(1._r8-following)))>1.e-12_r8) stop 14
end do
saved=1._r8
incoming=[10._r8,20._r8,30._r8,40._r8]
call ocean_exchange(incoming,saved,down,up)
if(down/=100._r8.or.up/=-100._r8) stop 9
saved=[0.1_r8,0.2_r8,0.3_r8,0.4_r8]
call ocean_exchange(incoming,saved,down,up)
if(abs(down-100._r8)>1.e-12_r8.or.abs(up+30._r8)>1.e-12_r8) stop 10
contains
subroutine ocean_albedo(cosz,values)
real(r8),intent(in)::cosz
real(r8),intent(out)::values(4)
real(r8),parameter::seq_flux_mct_albdif=0.06_r8
""" + albedo_declaration + albedo.group() + "\n" + albedo_result + """end subroutine
subroutine ocean_exchange(incoming,saved,down,up)
real(r8),intent(in)::incoming(4),saved(4)
real(r8),intent(out)::down,up
integer::n
""" + declarations + setup + flux.group() + "\n" + result + "end subroutine\nend program\n"
                self.compile_orbit_probe(program)
                controller = (ROOT / driver / "main/cime_comp_mod.F90").read_text()
                setup_routine = controller.split("subroutine cime_run_atmocn_setup(hashint)", 1)[1].split(
                    "end subroutine cime_run_atmocn_setup", 1)[0]
                self.assertLess(setup_routine.index("call prep_ocn_mrg"),
                                setup_routine.index("call cime_run_ocn_albedos"))


    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for compiled restart checks")
    def test_ocean_interval_metadata(self):
        header = """module interval_probe
use, intrinsic :: ieee_arithmetic
implicit none
integer,parameter::SHR_KIND_R8=kind(1d0)
integer::aborts=0
type seq_infodata_type
logical::ocn_present=.true.,flux_albav=.false.,read_restart=.false.
integer::ocean_solar_version=0
real(8)::nextsw_dt_avg=1800d0,restart_solar_dt_avg=1800d0
end type
contains
subroutine shr_sys_abort(message)
character(*)::message
aborts=aborts+1
end subroutine
"""
        program = """
end module
program probe
use interval_probe
implicit none
type(seq_infodata_type)::data
integer::restart,present_ocean,constant_albedo,version,width
logical::expected
do restart=0,1
do present_ocean=0,1
do constant_albedo=0,1
do version=-1,2
do width=1,4
data%read_restart=restart==1
data%ocn_present=present_ocean==1
data%flux_albav=constant_albedo==1
data%ocean_solar_version=version
select case(width)
case(1)
data%nextsw_dt_avg=-1d0
case(2)
data%nextsw_dt_avg=0d0
case(3)
data%nextsw_dt_avg=1800d0
case(4)
data%nextsw_dt_avg=ieee_value(0d0,ieee_quiet_nan)
end select
expected=width==2.or.width==3
if(restart==1) expected=width==3.and.version==1
if(present_ocean==0.or.constant_albedo==1) expected=.true.
aborts=0
call validate_ocean_solar(data)
if(expected.neqv.(aborts==0)) stop 1
end do
end do
end do
end do
end do
data=seq_infodata_type()
data%read_restart=.true.
data%ocean_solar_version=1
aborts=0
call publish_interval(data,80d0,900d0)
if(aborts==0) stop 2
aborts=0
call publish_interval(data,80d0,1800d0)
if(aborts/=0) stop 3
data%restart_solar_dt_avg=-1d0
aborts=0
call publish_interval(data,80d0,1800d0)
if(aborts==0) stop 4
data%restart_solar_dt_avg=ieee_value(0d0,ieee_quiet_nan)
aborts=0
call publish_interval(data,80d0,1800d0)
if(aborts==0) stop 5
end program
"""
        for driver, writers in (("driver-mct", 1), ("driver-moab", 2)):
            with self.subTest(driver=driver):
                metadata = (ROOT / driver / "shr/seq_infodata_mod.F90").read_text()
                routine = re.search(r"  subroutine validate_ocean_solar\(.*?end subroutine validate_ocean_solar",
                                    metadata, re.S)
                self.assertIsNotNone(routine)
                setter = "    if (present(nextsw_dt_avg)) then" + metadata.split(
                    "    if (present(nextsw_dt_avg)) then", 1)[1].split("    if ( present(precip_fact)", 1)[0]
                publisher = """
subroutine publish_interval(infodata,nextsw_cday,nextsw_dt_avg)
type(seq_infodata_type),intent(inout)::infodata
real(8),optional,intent(in)::nextsw_cday,nextsw_dt_avg
""" + setter + "end subroutine\n"
                self.compile_orbit_probe(header + routine.group() + publisher + program)
                self.assertEqual(metadata.count("call shr_mpi_bcast(infodata%nextsw_dt_avg,"), 4)
                restart_source = (ROOT / driver / "main/seq_rest_mod.F90").read_text()
                for field in ("seq_infodata_nextsw_dt_avg", "seq_infodata_ocean_solar_version"):
                    self.assertEqual(restart_source.count("'" + field + "'"), writers)


    @unittest.skipUnless(shutil.which("gfortran") and shutil.which("clang++"),
                         "gfortran and clang++ are required for linked scheduler checks")
    def test_eamxx_shortwave_interval_abi(self):
        source = (ROOT / "components/eamxx/src/mct_coupling/eamxx_cxx_f90_interface.cpp").read_text()
        scheduler = re.search(r"void scream_nextsw_time\(.*?\n\}", source, re.S)
        utilities = (ROOT / "components/eamxx/src/physics/rrtmgp/rrtmgp_utils.hpp").read_text()
        predicate = re.search(r"inline bool radiation_do\(.*?\n\}", utilities, re.S)
        geometry_source = (ROOT / "components/eamxx/src/share/physics/physics_share.cpp").read_text()
        geometry = re.search(r"  const bool coupled_orbit =.*?(?=  if \(!coupled_orbit\))",
                             geometry_source, re.S)
        for extracted in (scheduler, predicate, geometry):
            self.assertIsNotNone(extracted)
        cpp = """#include <cmath>
#include <stdexcept>
#define EKAT_REQUIRE_MSG(condition,message) if(condition) {} else throw std::runtime_error(message)
int frequency=1,step=0;
struct Stamp {
    int get_year() const{return 2000;}
    int get_month() const{return 3;}
    int get_day() const{return 20;}
    int get_num_steps() const{return step;}
    double frac_of_year_in_days() const{return 79.25;}
};
struct Driver {
    const int* get_atm_processes() const{return &frequency;}
    Stamp get_atm_time_stamp() const{return {};}
};
Driver driver;
Driver& get_ad_nonconst(){return driver;}
int radiation_frequency(const int& value){return value;}
template<class Function> void fpe_guard_wrapper(Function function){function();}
namespace scream {namespace rrtmgp {
""" + predicate.group() + """
}}
bool leap=false,geometry_available=true;
int geometry_calls=0,declination_calls=0;
bool use_leap_year(){return leap;}
bool shr_orb_coupled_geometry_c2f(int,double,double* delta,double* eccf,double* eot) {
    ++geometry_calls; *delta=0.25; *eccf=1; *eot=0; return geometry_available;
}
bool shr_orb_coupled_decl_c2f(int,double,double* delta,double* eccf) {
    ++declination_calls; *delta=0.5; *eccf=1; return true;
}
bool select_geometry(bool has_coupled_geometry) {
    int solar_ymd=20000320;
    double solar_tod=21600,delta=0,eccf=0,eot_seconds=0;
""" + geometry.group() + """
    return coupled_orbit;
}
extern "C" {
void set_probe(int value,int index){frequency=value;step=index;}
""" + scheduler.group() + """
int invalid_dt() {
    int date; double tod,day,width;
    try {scream_nextsw_time(0,&date,&tod,&day,&width);}
    catch(const std::runtime_error&) {return 1;}
    return 0;
}
int geometry_probe() {
    if(!select_geometry(true) || geometry_calls!=1 || declination_calls!=0) return 1;
    if(select_geometry(false) || geometry_calls!=1 || declination_calls!=0) return 2;
    leap=true;
    if(!select_geometry(false) || declination_calls!=1) return 3;
    geometry_available=false;
    try {select_geometry(true);}
    catch(const std::runtime_error&) {return 0;}
    return 4;
}
}
"""
        program = """program probe
use iso_c_binding
use eamxx_f2c_mod,only:scream_nextsw_time
implicit none
interface
subroutine set_probe(value,index) bind(c)
import c_int
integer(c_int),value::value,index
end subroutine
integer(c_int) function invalid_dt() bind(c)
import c_int
end function
integer(c_int) function geometry_probe() bind(c)
import c_int
end function
end interface
integer(c_int)::frequency,step,date
real(c_double)::tod,day,width
logical::active
do frequency=0,3
do step=0,5
call set_probe(frequency,step)
call scream_nextsw_time(900_c_int,date,tod,day,width)
active=.false.
if(available.and.frequency>0) active=step==0.or.mod(step,frequency)==0
if(available) then
if(width/=real(frequency,c_double)*900._c_double) stop 1
else
if(width/=0._c_double) stop 2
end if
if(active.neqv.(day>=0._c_double)) stop 3
if(active.and.date/=20000320) stop 4
if(active.and.tod/=21600._c_double) stop 5
end do
end do
if(invalid_dt()/=1) stop 6
if(geometry_probe()/=0) stop 7
end program
"""
        for enabled in (False, True):
            with self.subTest(rrtmgp=enabled), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                flags = ["-DEAMXX_HAS_RRTMGP"] if enabled else []
                subprocess.run([shutil.which("clang++"), "-std=c++17", "-Wall", "-Wextra", "-Werror",
                                *flags, "-x", "c++", "-c", "-", "-o", str(folder / "scheduler.o")],
                               input=cpp, text=True, capture_output=True, cwd=folder, check=True)
                runtime = ["-lc++"] if sys.platform == "darwin" else ["-lstdc++"]
                declaration = "logical,parameter::available=" + (".true." if enabled else ".false.") + "\n"
                fortran = program.replace("logical::active\n", "logical::active\n" + declaration)
                subprocess.run([shutil.which("gfortran"), "-cpp", "-ffree-line-length-none", "-fcheck=all",
                                "-ffpe-trap=invalid,zero,overflow", "-J", directory, "-I", directory,
                                str(ROOT / "components/eamxx/src/mct_coupling/eamxx_f2c_mod.F90"),
                                "-x", "f95", "-", "-x", "none", str(folder / "scheduler.o"),
                                *runtime, "-o", str(folder / "probe")],
                               input=fortran, text=True, capture_output=True, cwd=folder, check=True)
                subprocess.run([str(folder / "probe")], text=True, capture_output=True, cwd=folder, check=True)
        coupling = (ROOT / "components/eamxx/src/mct_coupling/atm_comp_mct.F90").read_text()
        self.assertEqual(coupling.count("call publish_solar_time(infodata, dt_scream)"), 2)

        with self.subTest("shared solar geometry"):
            header = (ROOT / "components/eamxx/src/share/physics/physics_share.hpp").read_text()
            api = header.split("namespace physics {", 1)[1].split('extern "C"', 1)[0].rsplit("}", 1)[0]
            implementation = geometry_source.split("namespace physics {", 1)[1].split(
                "//\n// A C++ interface", 1)[0]
            self.assertIn("SolarGeometry solar_geometry", api)
            self.assertIn("double SolarGeometry::cosine_zenith", implementation)
            shared_cpp = """#include <cmath>
#include <ctime>
#include <limits>
#include <stdexcept>
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"
#define EKAT_REQUIRE_MSG(condition,message) if(condition) {} else throw std::runtime_error(message)
namespace ekat {
struct ParameterList {
    bool present=false,value=true;
    bool isParameter(const char*) const{return present;}
    template<class Value> Value get(const char*) const{return value;}
};
}
namespace scream {
bool leap_enabled=false;
bool use_leap_year(){return leap_enabled;}
namespace util {
bool is_leap_year(int year){return year%4==0 && (year%100!=0 || year%400==0);}
struct TimeStamp {
    int year,month,day,hour,minute,second;
    TimeStamp(int year_in,int month_in,int day_in,int hour_in,int minute_in,int second_in)
        :year(year_in),month(month_in),day(day_in),hour(hour_in),minute(minute_in),second(second_in) {}
    int get_year() const{return year;}
    int get_month() const{return month;}
    int get_day() const{return day;}
    double midnight_seconds() const {
        std::tm date={}; date.tm_year=year-1900; date.tm_mon=month-1; date.tm_mday=day;
        return static_cast<double>(timegm(&date));
    }
    double days_from(const TimeStamp& other) const {
        return (midnight_seconds()-other.midnight_seconds())/86400;
    }
    double frac_of_year_in_days() const {
        double day_index=days_from(TimeStamp(year,1,1,0,0,0));
        if(!use_leap_year() && is_leap_year(year) && month>2) --day_index;
        return day_index+(hour*3600+minute*60+second)/86400.;
    }
    int days_in_curr_year() const{return use_leap_year() && is_leap_year(year) ? 366 : 365;}
};
}
namespace physics {
""" + api + implementation + """
}
}
bool callback_available=true,callback_eot=false;
int callback_date=0;
double callback_tod=0;
bool solar_callback(int date,double tod,double* delta,double* eccf,double* eot) {
    callback_date=date; callback_tod=tod;
    *delta=0; *eccf=1.023; *eot=callback_eot ? 120 : 0;
    return callback_available;
}
extern "C" int shared_solar_probe() {
    using namespace scream;
    const int dates[]={19991231,20000101,20000229,20000301,20001231,20010101,19000301,21000301};
    const double widths[]={0,600,3600};
    for(int calendar=0;calendar<2;++calendar) {
        leap_enabled=calendar;
        for(int setting=-1;setting<=1;++setting) {
            const bool use_eot=setting!=0;
            ekat::ParameterList params{setting>=0,use_eot};
            for(int date:dates) {
                if(!leap_enabled && date==20000229) continue;
                util::TimeStamp timestamp(date/10000,date/100%100,date%100,5,50,0);
                for(int mode=0;mode<3;++mode) {
                    const int configured_year=mode==0 ? 1850 : -1;
                    const double configured_eccen=mode==2 ? 0.017 : -1;
                    const double configured_obliq=mode==2 ? 23.4 : -1;
                    const double configured_mvelp=mode==2 ? 102.7 : -1;
                    const auto solar=physics::solar_geometry(timestamp,params,configured_year,
                            configured_eccen,configured_obliq,configured_mvelp);
                    int orbital_year=mode==2 ? shr_orb_undef_int_c2f : mode==0 ? 1850 : timestamp.get_year();
                    double eccen=configured_eccen,obliq=configured_obliq,mvelp=configured_mvelp;
                    double obliqr,lambm0,mvelpp,delta,eccf,eot=0;
                    const double day=timestamp.frac_of_year_in_days()+1;
                    if(calendar && mode==1) {
                        shr_orb_params_continuous_c2f(orbital_year,(day-1)/timestamp.days_in_curr_year(),
                                &eccen,&obliq,&mvelp,&obliqr,&lambm0,&mvelpp);
                    } else {
                        shr_orb_params_c2f(&orbital_year,&eccen,&obliq,&mvelp,&obliqr,&lambm0,&mvelpp);
                    }
                    if(calendar && mode==0) {
                        shr_orb_decl_fixed_year_c2f(day,timestamp.days_in_curr_year(),
                                util::is_leap_year(timestamp.get_year()-1) ? 366 : 365,
                                eccen,mvelpp,lambm0,obliqr,&delta,&eccf,&eot,use_eot);
                    } else if(calendar) {
                        const auto elapsed=timestamp.days_from(util::TimeStamp(2000,1,1,0,0,0))+day-std::floor(day);
                        if(use_eot) shr_orb_decl_gregorian_eot_c2f(elapsed,eccen,mvelpp,lambm0,obliqr,&delta,&eccf,&eot);
                        else shr_orb_decl_gregorian_c2f(elapsed,eccen,mvelpp,lambm0,obliqr,&delta,&eccf);
                    } else {
                        if(use_eot) shr_orb_decl_eot_c2f(day,eccen,mvelpp,lambm0,obliqr,&delta,&eccf,&eot);
                        else shr_orb_decl_c2f(day,eccen,mvelpp,lambm0,obliqr,&delta,&eccf);
                    }
                    if(std::abs(solar.eccentricity_factor()-eccf)>1.e-14) return 1;
                    for(double width:widths) {
                        const double expected=use_eot ? shr_orb_cosz_eot_c2f(day,0.2,-0.1,delta,width,eot)
                                                                                 : shr_orb_cosz_c2f(day,0.2,-0.1,delta,width);
                        const auto actual=solar.cosine_zenith(0.2,-0.1,width);
                        if(!std::isfinite(actual) || std::abs(actual-expected)>1.e-13) return 2;
                    }
                }
            }
        }
        for(int enabled=0;enabled<2;++enabled) {
            callback_eot=enabled;
            scream_set_solar_geometry_callback(solar_callback,enabled);
            const util::TimeStamp timestamp(2000,3,20,5,50,0);
            const auto solar=physics::solar_geometry(timestamp,ekat::ParameterList{},1850,-1,-1,-1);
            if(solar.eccentricity_factor()!=1.023 || callback_date!=20000320 ||
                    std::abs(callback_tod-21000)>1.e-6) return 3;
            for(double width:widths) {
                const auto day=timestamp.frac_of_year_in_days()+1;
                const auto expected=enabled ? shr_orb_cosz_eot_c2f(day,0,0,0,width,120)
                                                                        : shr_orb_cosz_c2f(day,0,0,0,width);
                if(std::abs(solar.cosine_zenith(0,0,width)-expected)>1.e-14) return 4;
            }
            int rejected=0;
            try {physics::solar_geometry(timestamp,ekat::ParameterList{true,!enabled},1850,-1,-1,-1);}
            catch(const std::runtime_error&) {++rejected;}
            callback_available=false;
            try {physics::solar_geometry(timestamp,ekat::ParameterList{},1850,-1,-1,-1);}
            catch(const std::runtime_error&) {++rejected;}
            callback_available=true;
            const double invalid_widths[]={-1,std::numeric_limits<double>::quiet_NaN(),
                                                                        std::numeric_limits<double>::infinity()};
            for(double width:invalid_widths) {
                try {solar.cosine_zenith(0,0,width);}
                catch(const std::runtime_error&) {++rejected;}
            }
            if(rejected!=5) return 5;
            if(!enabled && !(solar.cosine_zenith(0,0,600)<=0 && solar.cosine_zenith(0,0,3600)>0)) return 6;
            scream_set_solar_geometry_callback(nullptr,false);
        }
    }
    return 0;
}
"""
            registrations = re.findall(r"void scream_set_solar_(?:geometry_)?callback\(.*?\n\}", source, re.S)
            self.assertEqual(len(registrations), 2)
            shared_cpp = shared_cpp.replace("#define EKAT_REQUIRE_MSG", "\n".join(registrations) +
                                            "\n#define EKAT_REQUIRE_MSG", 1)
            program = """program probe
use iso_c_binding
implicit none
interface
integer(c_int) function shared_solar_probe() bind(c)
import c_int
end function
end interface
integer::status
status=shared_solar_probe()
if(status/=0) stop status
end program
"""
            with tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                subprocess.run([shutil.which("clang++"), "-std=c++17", "-Wall", "-Wextra", "-Werror",
                        "-I", str(ROOT / "components/eamxx/src"), "-x", "c++", "-c", "-",
                        "-o", str(folder / "solar.o")], input=shared_cpp, text=True,
                       capture_output=True, cwd=folder, check=True)
                runtime = ["-lc++"] if sys.platform == "darwin" else ["-lstdc++"]
                self.compile_orbit_probe(program, ["components/eamxx/src/physics/rrtmgp/shr_orb_mod_c2f.F90",
                                  folder / "solar.o"], link_args=runtime)
            for process in ("rrtmgp/eamxx_rrtmgp_process_interface.cpp", "mam/eamxx_mam_microphysics_process_interface.cpp"):
                process_source = (ROOT / "components/eamxx/src/physics" / process).read_text()
                self.assertIn("physics::solar_geometry(", process_source)
                self.assertIn("solar.cosine_zenith(", process_source)
                self.assertNotIn("shr_orb_", process_source)
            shared_build = (ROOT / "components/eamxx/src/share/physics/CMakeLists.txt").read_text()
            radiation_build = (ROOT / "components/eamxx/src/physics/rrtmgp/CMakeLists.txt").read_text()
            chemistry_build = (ROOT / "components/eamxx/src/physics/mam/CMakeLists.txt").read_text()
            self.assertIn("../../physics/rrtmgp/shr_orb_mod_c2f.F90", shared_build)
            self.assertNotIn("shr_orb_mod_c2f.F90", radiation_build)
            self.assertNotIn("shr_orb_mod_c2f.F90", chemistry_build)


    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for interval getter checks")
    def test_atmosphere_interval_producers(self):
        for variant in ("rrtmg", "rrtmgp", "crm/rrtmgp"):
            with self.subTest(variant=variant):
                source = (ROOT / "components/eam/src/physics" / variant / "radiation.F90").read_text()
                getter = re.search(r"(?:real\(r8\) function radiation_sw_dt_avg\(\)).*?"
                                   r"end function radiation_sw_dt_avg", source, re.S)
                self.assertIsNotNone(getter)
                program = """module interval_probe
use shr_kind_mod,only:r8=>SHR_KIND_R8
implicit none
real(r8)::dt_avg=0._r8
contains
""" + getter.group() + """
end module
program probe
use interval_probe
implicit none
if(radiation_sw_dt_avg()/=0._r8) stop 1
dt_avg=1800._r8
if(radiation_sw_dt_avg()/=1800._r8) stop 2
end program
"""
                self.compile_orbit_probe(program)
        eam = (ROOT / "components/eam/src/cpl/atm_comp_mct.F90").read_text()
        datm = (ROOT / "components/data_comps/datm/src/atm_comp_mct.F90").read_text()
        self.assertEqual(eam.count("nextsw_dt_avg=radiation_sw_dt_avg()"), 3)
        self.assertEqual(datm.count("nextsw_dt_avg=0._r8"), 3)
        first_phase = eam.split("if (nsrest == 0) then", 1)[1].split("else ! so here first_time", 1)[0]
        self.assertLess(first_phase.index("       end if"), first_phase.index("nextsw_dt_avg="))




if __name__ == "__main__":
    unittest.main()
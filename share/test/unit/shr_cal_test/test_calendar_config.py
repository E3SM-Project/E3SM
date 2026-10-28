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

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for stub-backed numerical checks")
    def test_compiled_shared_orbit_numerics(self):
        self.compile_orbit_probe()
        with self.subTest(geometry="forcing geometric horizon"):
            source = (ROOT / "share/streams/shr_tInterp_mod.F90").read_text()
            routine = re.search(r"   subroutine shr_tInterp_getCosz\(.*?  end subroutine shr_tInterp_getCosz", source, re.S)
            averaging = re.search(r"  subroutine shr_tInterp_getAvgCosz\(.*?  end subroutine shr_tInterp_getAvgCosz", source, re.S)
            scaling_source = (ROOT / "share/streams/shr_strdata_mod.F90").read_text()
            scaling = re.search(
                r"             do i = 1,lsizeF\n                if \(cosz\(i\).*?             enddo", scaling_source, re.S)
            self.assertIsNotNone(routine)
            self.assertIsNotNone(averaging)
            self.assertIsNotNone(scaling)
            elm_source = (ROOT / "components/elm/src/utils/elm_varorb.F90").read_text()
            elm_decl = re.search(r"  subroutine elm_orbital_decl\(.*?  end subroutine elm_orbital_decl", elm_source, re.S)
            elm_cosine = re.search(r"  function elm_solar_cosine\(.*?  end function elm_solar_cosine", elm_source, re.S)
            self.assertIsNotNone(elm_decl)
            self.assertIsNotNone(elm_cosine)
            bypass_source = (ROOT / "components/elm/src/cpl/lnd_import_export.F90").read_text()
            bypass_block = bypass_source.split("        !Shortwave radiation (cosine zenith angle interpolation)", 1)[1].split(
                "        if (atm2lnd_vars%metsource == 5)", 1)[0]
            bypass_cases = []
            def arguments(values):
                return ",".join(str(value) for value in values)

            for date, tod, step, ratio, calendar in (
                    ((2000, 1, 1), 0, 600, 6, "GREGORIAN"),
                    ((2000, 3, 1), 0, 600, 6, "GREGORIAN"),
                    ((2000, 3, 1), 0, 600, 6, "NO_LEAP"),
                    ((1900, 3, 1), 0, 600, 6, "GREGORIAN"),
                    ((2100, 3, 1), 0, 600, 6, "NO_LEAP"),
                    ((2000, 1, 1), 0, 601, 6, "GREGORIAN"),
                    ((2000, 3, 20), 21602, 1, 1, "NO_LEAP")):
                midpoint = -step / 2
                period_start = math.floor((tod + midpoint) / (step * ratio)) * step * ratio
                offsets = [midpoint]
                if ratio > 1:
                    offsets += [period_start + (sample + 0.5) * step - tod for sample in range(ratio)]
                dates = [datetime.datetime(*date) + datetime.timedelta(seconds=tod + math.floor(offset)) for offset in offsets]
                if calendar == "NO_LEAP":
                    dates = [date.replace(day=28) if date.month == 2 and date.day == 29 else date for date in dates]
                bypass_cases.append(
                    f"model_step={step}\nmodel_calendar='{calendar}'\nexpected_calls={len(offsets)}\n"
                    "clock_calls=0\nrecord_count=0\nrecord_geometry=.true.\n"
                    f"expected_offsets(1:expected_calls)=[{arguments(math.floor(offset) for offset in offsets)}]\n"
                    f"expected_dates(1:expected_calls)=[{arguments(date.year*10000+date.month*100+date.day for date in dates)}]\n"
                    f"expected_seconds(1:expected_calls)=[{arguments(date.hour*3600+date.minute*60+date.second for date in dates)}]\n"
                    f"call bypass_weight({tod},{ratio},flux)\n"
                    "if(clock_calls/=expected_calls.or.record_count/=expected_calls) stop 10\n"
                    "if(any(record_dates(1:record_count)/=expected_dates(1:record_count))) stop 11\n"
                    f"if(any(abs(record_tods(1:record_count)-real(expected_seconds(1:record_count),r8) &\n"
                    f"     -{abs(midpoint-math.floor(midpoint))}_r8)>1.e-12_r8)) stop 12\n"
                    "if(.not.ieee_is_finite(flux).or.flux<0._r8.or.flux>10._r8) stop 13\n"
                    + ("if(flux/=1._r8) stop 14\n" if ratio == 1 else ""))
            elm_probe = """module shr_cal_mod
use forcing_probe,only:shr_cal_date2julian
end module
module elm_probe
use forcing_probe,only:r8,shr_orb_context_type,seq_infodata_SolarDecl=>shr_cal_solarDecl
use shr_sys_mod,only:shr_sys_abort
implicit none
type(shr_orb_context_type)::solar_context
logical::solar_context_initialized=.true.
contains
""" + elm_decl.group() + "\n" + elm_cosine.group() + """
end module
module bypass_probe
use elm_probe,only:elm_solar_cosine
use forcing_probe,only:r8
implicit none
integer::model_step,expected_calls,clock_calls,expected_offsets(32),expected_dates(32),expected_seconds(32)
character(16)::model_calendar
contains
integer function get_step_size()
get_step_size=model_step
end function
function get_calendar() result(calendar)
character(16)::calendar
calendar=model_calendar
end function
subroutine get_curr_date(year,month,day,seconds,offset)
integer,intent(out)::year,month,day,seconds
integer,intent(in)::offset
clock_calls=clock_calls+1
if(clock_calls>expected_calls) stop 20
if(offset/=expected_offsets(clock_calls)) stop 21
year=expected_dates(clock_calls)/10000
month=modulo(expected_dates(clock_calls)/100,100)
day=modulo(expected_dates(clock_calls),100)
seconds=expected_seconds(clock_calls)
end subroutine
subroutine bypass_weight(tod,ratio,weight)
integer,intent(in)::tod,ratio
real(r8),intent(out)::weight
real(r8),parameter::SHR_CONST_PI=acos(-1._r8)
integer,parameter::g=1
integer::tm,swrad_period_len,swrad_period_start,solar_year,solar_month,solar_day,solar_seconds
real(r8)::thiscosz,avgcosz,solar_offset,solar_tod,solar_latitude,solar_longitude,wt2(14)
type forcing_type
real(r8)::npf(14)
end type
type domain_type
real(r8)::latc(1),lonc(1)
end type
type(forcing_type)::atm2lnd_vars
type(domain_type)::ldomain
atm2lnd_vars%npf=real(ratio,r8)
ldomain%latc=0._r8
ldomain%lonc=0._r8
""" + bypass_block + """
weight=wt2(4)
end subroutine
end module
"""
            program = """module forcing_probe
use shr_kind_mod
use shr_orb_mod,only:shr_orb_cosz,shr_orb_decl,shr_orb_params
use shr_sys_mod,only:shr_sys_abort
implicit none
integer,parameter::r8=SHR_KIND_R8
character(*),parameter::shr_cal_gregorian='GREGORIAN'
type shr_orb_context_type
logical::use_equation_of_time=.false.
real(r8)::eccen,mvelpp,lambm0,obliqr,eot
end type
logical::record_geometry=.false.
integer::record_count=0,record_dates(32)
real(r8)::record_tods(32)
type ESMF_Time
integer(SHR_KIND_I8)::seconds
end type
type ESMF_TimeInterval
integer(SHR_KIND_I8)::seconds
end type
interface operator(>)
module procedure time_greater
end interface
interface operator(<)
module procedure time_less
end interface
interface operator(-)
module procedure time_difference
end interface
contains
logical function time_greater(left,right)
type(ESMF_Time),intent(in)::left,right
time_greater=left%seconds>right%seconds
end function
logical function time_less(left,right)
type(ESMF_Time),intent(in)::left,right
time_less=left%seconds<right%seconds
end function
function time_difference(left,right) result(difference)
type(ESMF_Time),intent(in)::left,right
type(ESMF_TimeInterval)::difference
difference%seconds=left%seconds-right%seconds
end function
subroutine ESMF_TimeIntervalGet(interval,s_i8)
type(ESMF_TimeInterval),intent(in)::interval
integer(SHR_KIND_I8),intent(out)::s_i8
s_i8=interval%seconds
end subroutine
subroutine shr_cal_timeSet(time,ymd,tod,calendar)
type(ESMF_Time),intent(out)::time
integer,intent(in)::ymd,tod
character(*),intent(in)::calendar
time%seconds=int(tod,SHR_KIND_I8)
end subroutine
subroutine shr_cal_advDateInt(seconds,unit,ymd0,tod0,ymd,tod,calendar)
integer,intent(in)::seconds,ymd0,tod0
character(*),intent(in)::unit,calendar
integer,intent(out)::ymd,tod
ymd=ymd0
tod=tod0+seconds
end subroutine
function shr_cal_calendarName(calendar) result(name)
character(*),intent(in)::calendar
character(len=len(calendar))::name
name=calendar
end function
subroutine shr_cal_date2julian(ymd,tod,calday,calendar)
integer,intent(in)::ymd,tod
character(*),intent(in)::calendar
real(r8),intent(out)::calday
calday=80._r8+real(tod,r8)/86400._r8
end subroutine
subroutine shr_cal_solarDecl(context,ymd,tod,calendar,delta,factor,eot)
type(shr_orb_context_type),intent(in)::context
integer,intent(in)::ymd
real(r8),intent(in)::tod
character(*),intent(in)::calendar
real(r8),intent(out)::delta,factor,eot
if(record_geometry) then
record_count=record_count+1
record_dates(record_count)=ymd
record_tods(record_count)=tod
end if
call shr_orb_decl(80._r8+tod/86400._r8,context%eccen,context%mvelpp,context%lambm0,context%obliqr,delta,factor)
eot=context%eot
end subroutine
""" + routine.group() + "\n" + averaging.group() + """
subroutine scale_flux(cosz,tavCosz,flux)
real(r8),intent(in)::cosz(1),tavCosz(1)
real(r8),intent(out)::flux
integer,parameter::lsizeF=1,n=1
integer::i
real(r8),parameter::solZenMin=.001_r8
type vector
real(r8)::rAttr(1,1)
end type
type data_type
type(vector)::avs(1),avFLB(1)
end type
type(data_type)::SDAT
SDAT%avFLB(1)%rAttr=100._r8
""" + scaling.group() + """
flux=SDAT%avs(1)%rAttr(1,1)
end subroutine
end module
""" + elm_probe + """
program probe
use forcing_probe
use elm_probe,only:elm_solar_cosine,solar_context
use bypass_probe,only:bypass_weight,model_step,model_calendar,expected_calls,clock_calls, &
    expected_offsets,expected_dates,expected_seconds
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
implicit none
type(shr_orb_context_type)::context
integer::orbital_year,calendar_index,eot_index,tod,sample,period
real(r8)::obliq,mvelp,delta,factor,expected(7),actual(7),average,flux,latitude(1),longitude(1),cosz(1)
real(r8)::reference_average,averaged(1)
integer::sample_count,sample_seconds(4),sample_weights(4)
character(16)::calendar
orbital_year=2000
context%eccen=-9999._r8
obliq=-9999._r8
mvelp=-9999._r8
call shr_orb_params(orbital_year,context%eccen,obliq,mvelp,context%obliqr,context%lambm0,context%mvelpp,.false.)
latitude=0._r8
longitude=0._r8
do calendar_index=1,2
calendar='NO_LEAP'
if(calendar_index==2) calendar='GREGORIAN'
do eot_index=-1,1
context%use_equation_of_time=eot_index/=0
context%eot=real(eot_index,r8)*960._r8
solar_context=context
sample_count=3
sample_seconds=[21601,21606,21611,0]
sample_weights=[5,5,1,0]
if(calendar=='NO_LEAP'.and..not.context%use_equation_of_time) then
sample_count=4
sample_seconds=[21601,21604,21607,21610]
sample_weights=[3,3,3,2]
end if
reference_average=0._r8
do sample=1,sample_count
tod=sample_seconds(sample)
call shr_orb_decl(80._r8+real(tod,r8)/86400._r8,context%eccen,context%mvelpp,context%lambm0,context%obliqr,delta,factor)
reference_average=reference_average+real(sample_weights(sample),r8)*max(0._r8,cos(delta) &
    *cos((real(tod,r8)/86400._r8-.5_r8)*2._r8*acos(-1._r8)+2._r8*acos(-1._r8)*context%eot/86400._r8))/11._r8
end do
call shr_tInterp_getAvgCosz(averaged,longitude,latitude,20000320,21601,20000320,21612, &
    context%eccen,context%mvelpp,context%lambm0,context%obliqr,5,calendar,context)
if(abs(averaged(1)-reference_average)>1.e-12_r8) stop 7
call shr_tInterp_getAvgCosz(averaged,longitude,latitude,20000320,0,20000320,11, &
    context%eccen,context%mvelpp,context%lambm0,context%obliqr,5,calendar,context)
if(averaged(1)/=0._r8) stop 8
do period=0,23
do sample=1,7
tod=period*3600+(sample-1)*600+1
call shr_orb_decl(80._r8+real(tod,r8)/86400._r8,context%eccen,context%mvelpp,context%lambm0,context%obliqr,delta,factor)
expected(sample)=max(0._r8,cos(delta)*cos((real(tod,r8)/86400._r8-.5_r8)*2._r8*acos(-1._r8) &
     +2._r8*acos(-1._r8)*context%eot/86400._r8))
call shr_tInterp_getCosz(cosz,longitude,latitude,20000320,tod,context%eccen,context%mvelpp, &
     context%lambm0,context%obliqr,calendar,context)
actual(sample)=cosz(1)
if(abs(actual(sample)-expected(sample))>1.e-12_r8) stop 1
if(abs(max(0._r8,elm_solar_cosine(20000320,real(tod,r8),calendar,latitude(1),longitude(1))) &
    -expected(sample))>1.e-12_r8) stop 6
end do
average=sum(expected)/real(size(expected),r8)
do sample=1,7
call scale_flux([actual(sample)],[average],flux)
if(.not.ieee_is_finite(flux)) stop 2
if(average>0._r8) then
if(abs(flux-100._r8*expected(sample)/average)>1.e-8_r8*max(1._r8,abs(flux))) stop 3
else
if(flux/=0._r8) stop 4
end if
end do
end do
end do
end do
call scale_flux([0._r8],[0._r8],flux)
if(flux/=0._r8) stop 5
context%use_equation_of_time=.false.
context%eot=0._r8
solar_context=context
""" + "\n".join(bypass_cases) + """
end program
"""
            self.compile_orbit_probe(program)
            for calendar in ("NO_LEAP", "GREGORIAN"):
                for step, end_seconds in ((0, 21612), (5, 21601)):
                    with self.subTest(calendar=calendar, invalid_average=(step, end_seconds)):
                        invalid = f"""calendar='{calendar}'
call shr_tInterp_getAvgCosz(averaged,longitude,latitude,20000320,21601,20000320,{end_seconds}, &
     context%eccen,context%mvelpp,context%lambm0,context%obliqr,{step},calendar,context)
stop 77
do calendar_index=1,2"""
                        self.compile_orbit_probe(program.replace("do calendar_index=1,2", invalid, 1), expected_returncode=99)
        with self.subTest(geometry="scheduled Icepack clock"):
            source = (ROOT / "components/mpas-seaice/src/icepack/columnphysics/icepack_orbital.F90").read_text()
            routine = re.search(r"      subroutine compute_coszen \(.*?      end subroutine compute_coszen", source, re.S)
            self.assertIsNotNone(routine)
            program = """#define CESMCOUPLED
    module icepack_parameters
    implicit none
    contains
    logical function icepack_chkoptargflag(first_call)
    logical,intent(in)::first_call
    icepack_chkoptargflag=first_call
    end function
    end module
    module orbital_probe
    use shr_kind_mod,only:r8=>SHR_KIND_R8
    use shr_orb_mod,only:shr_orb_params,shr_orb_decl,shr_orb_cosz
    implicit none
    integer,parameter::dbl_kind=r8,int_kind=kind(1),log_kind=kind(.true.),char_len=64
    real(r8),parameter::c0=0._r8,c2=2._r8,p5=0.5_r8,pi=acos(-1._r8),secday=86400._r8
    real(r8)::eccen,mvelpp,lambm0,obliqr,decln,eccf
    contains
    subroutine icepack_warnings_add(message)
    character(*),intent(in)::message
    print *,message
    end subroutine
    subroutine icepack_warnings_setabort(value,filename,line)
    logical,intent(in)::value
    character(*),intent(in)::filename
    integer,intent(in)::line
    if(value) stop 90
    end subroutine
    """ + routine.group() + """
    end module
    program probe
    use orbital_probe
    implicit none
    integer::orbital_year,calendar_index,season,tod,latitude_degrees,longitude_degrees,eot_index
    real(r8)::obliq,mvelp,calday,latitude,longitude,delta,factor,expected,actual,eot
    character(char_len)::calendar
    logical::use_eot
    orbital_year=2000
    eccen=-9999._r8
    obliq=-9999._r8
    mvelp=-9999._r8
    call shr_orb_params(orbital_year,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.false.)
    do calendar_index=1,2
    calendar='NO_LEAP'
    if(calendar_index==2) calendar='GREGORIAN'
    do season=1,365,91
    do tod=0,86400,600
    calday=real(season,r8)+real(tod,r8)/secday
    call shr_orb_decl(calday,eccen,mvelpp,lambm0,obliqr,delta,factor)
    do eot_index=-1,1
    eot=real(eot_index,r8)*960._r8
    use_eot=eot_index/=0
    do latitude_degrees=-90,90,30
    latitude=real(latitude_degrees,r8)*pi/180._r8
    do longitude_degrees=-180,180,60
    longitude=real(longitude_degrees,r8)*pi/180._r8
    expected=shr_orb_cosz(calday,latitude,longitude,delta,eot_seconds=eot)
    actual=-99._r8
    call compute_coszen(latitude,longitude,real(season,r8),modulo(tod+1800,86400),actual, &
        365,calday,calendar,20000320,real(modulo(tod,86400),r8),delta,.true.,eot,use_eot)
    if(abs(actual-expected)>1.e-12_r8) stop 1
    end do
    end do
    end do
    end do
    end do
    end do
    actual=0.123_r8
    call compute_coszen(0._r8,0._r8,80._r8,21000,actual,365,-1._r8,calendar, &
        0,-1._r8,0._r8,.false.,0._r8,.true.)
    if(actual/=0.123_r8) stop 2
    end program
    """
            self.compile_orbit_probe(program)

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

    def test_runtime_context_and_restart_phase_wiring(self):
        land = (ROOT / "components/elm/src/cpl/lnd_comp_mct.F90").read_text()
        land_run = land.split("subroutine lnd_run_mct(", 1)[1].split("end subroutine lnd_run_mct", 1)[0]
        self.assertEqual(land_run.count("call elm_set_orbital_context(infodata)"), 1)
        self.assertLess(land_run.index("call elm_set_orbital_context(infodata)"), land_run.index("call lnd_import("))
        bypass = (ROOT / "components/elm/src/cpl/lnd_import_export.F90").read_text()
        self.assertNotIn("szenith", bypass)
        self.assertEqual(bypass.count("offset=floor(solar_offset)"), 2)
        source = (ROOT / "components/eam/src/cpl/atm_comp_mct.F90").read_text()
        run = source.lower().split("subroutine atm_run_mct(", 1)[1].split("end subroutine atm_run_mct", 1)[0]
        self.assertIn("call orbit_set_context(infodata)", run)
        for driver, writers in (("driver-mct", 1), ("driver-moab", 2)):
            restart = (ROOT / driver / "main/seq_rest_mod.F90").read_text()
            self.assertEqual(restart.count("'seq_infodata_solar_phase_version'"), writers)
            metadata = (ROOT / driver / "shr/seq_infodata_mod.F90").read_text()
            self.assertIn("if (infodata%solar_phase_version /= 1)", metadata)
            self.assertIn("call shr_mpi_bcast(infodata%solar_phase_version, mpicom)", metadata)

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

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for shared solar diagnostics")
    def test_solar_diagnostic_module(self):
        program = """
    program probe
    use shr_kind_mod,only:r8=>shr_kind_r8
    use shr_solar_diag_mod
    implicit none
    real(r8) :: albedo(4,2)
    albedo=.2_r8
    call shr_solar_diag_configure(.true.)
    call shr_solar_diag_schedule(20000229,3600._r8,.true.)
    call shr_solar_diag_import(1,albedo)
    call shr_solar_diag_import(3,albedo)
    call shr_solar_diag_consume(1,20000229,3600._r8,albedo)
    call shr_solar_diag_consume(3,20000229,3600._r8,albedo)
    call shr_solar_diag_column(1,20000229,3600._r8,[300._r8],[200._r8],[100._r8])
    call shr_solar_diag_radiation(1,20000229,3600._r8,albedo,[300._r8,300._r8], &
     [200._r8,200._r8],reshape([.0098_r8,.0098_r8,.0098_r8,.0098_r8],[2,2]), &
    reshape([50000._r8,50000._r8,50000._r8,50000._r8],[2,2]),9.8_r8,.true., &
    top=reshape([310._r8,310._r8,10._r8,10._r8],[2,2]))
    call shr_solar_diag_buffer([20._r8],[30._r8],[40._r8],[10._r8])
    call shr_solar_diag_clock(20000229,3600._r8,1800._r8)
    if (shr_solar_diag_dt /= 1800._r8) stop 2
    call shr_solar_diag_configure(.false.)
    call shr_solar_diag_configure(.true.)
    call shr_solar_diag_schedule(20000229,3600._r8,.true.)
    call shr_solar_diag_import(17,albedo,19)
    call shr_solar_diag_import(19,albedo,19)
    call shr_solar_diag_consume(17,20000229,3600._r8,albedo)
    call shr_solar_diag_consume(19,20000229,3600._r8,albedo)
    end program
    """
        sources = ("share/util/shr_solar_diag_mod.F90",)
        result = self.compile_orbit_probe(program=program, extra_sources=sources)
        self.assertIn("SOLAR_CONSUME_V1", result.stdout)
        self.assertIn("SOLAR_COLUMN_V1", result.stdout)
        for original, replacement, message in (
            ("[100._r8]", "[90._r8]", "TOA/heating closure failed"),
            ("[10._r8]", "[20._r8]", "buffer energy closure failed"),
            ("310._r8,310._r8,10._r8,10._r8", "310._r8,310._r8,9._r8,9._r8", "upper layer closure failed"),
            ("consume(1,20000229,3600._r8", "consume(1,20000229,5400._r8", "stale albedo consumption time"),
            ("call shr_solar_diag_consume(1", "albedo=.3_r8\ncall shr_solar_diag_consume(1", "imported and consumed albedos differ")):
            rejected = self.compile_orbit_probe(program=program.replace(original, replacement),
                            extra_sources=sources, expected_returncode=99)
            self.assertIn(message, rejected.stdout)
            if message in ("TOA/heating closure failed", "upper layer closure failed"):
                self.assertIn("chunk=1", rejected.stdout)
                self.assertIn("column=1", rejected.stdout)
                self.assertIn("tolerance=", rejected.stdout)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for the ocean buffer probe")
    def test_solar_ocean_accumulator(self):
        source = (ROOT / "driver-mct/main/prep_ocn_mod.F90").read_text()
        routines = []
        for name in ("prep_ocn_accum", "prep_ocn_accum_avg", "prep_ocn_solar_record"):
            match = re.search(r"  subroutine " + name + r"\(.*?  end subroutine " + name + r"\b", source, re.DOTALL)
            self.assertIsNotNone(match)
            routines.append(match.group(0))
        program = """
module ocean_probe
use shr_kind_mod,only:r8=>shr_kind_r8
use shr_sys_mod,only:shr_sys_abort
use shr_const_mod,only:shr_const_pi
use shr_solar_diag_mod
implicit none
type mct_avect
real(r8),allocatable::rAttr(:,:)
end type
type mct_ggrid
type(mct_avect)::data
end type
type(mct_avect),target::held
type(mct_ggrid),target::grid
type(mct_avect),allocatable::x2oacc_ox(:)
integer::ocn(1)=[1],num_inst_ocn=1,mpicom_CPLID=0,CPLID=0,logunit=6,x2oacc_ox_cnt=0
contains
function component_get_x2c_cx(instance) result(value)
integer,intent(in)::instance
type(mct_avect),pointer::value
value=>held
end function
function component_get_dom_cx(instance) result(value)
integer,intent(in)::instance
type(mct_ggrid),pointer::value
value=>grid
end function
integer function mct_aVect_indexRA(vector,name)
type(mct_avect),intent(in)::vector
character(*),intent(in)::name
mct_aVect_indexRA=1
if(name=='ifrac') mct_aVect_indexRA=2
end function
subroutine mct_avect_copy(source,destination)
type(mct_avect),intent(in)::source
type(mct_avect),intent(inout)::destination
destination%rAttr=source%rAttr
end subroutine
subroutine mct_avect_accum(source,destination)
type(mct_avect),intent(in)::source
type(mct_avect),intent(inout)::destination
destination%rAttr=destination%rAttr+source%rAttr
end subroutine
subroutine mct_avect_avg(vector,count)
type(mct_avect),intent(inout)::vector
integer,intent(in)::count
vector%rAttr=vector%rAttr/count
end subroutine
subroutine t_drvstartf(timer,barrier)
character(*),intent(in)::timer
integer,optional,intent(in)::barrier
end subroutine
subroutine t_drvstopf(timer)
character(*),intent(in)::timer
end subroutine
subroutine shr_mpi_sum(local,global,comm,tag)
real(r8),intent(in)::local(:)
real(r8),intent(out)::global(:)
integer,intent(in)::comm
character(*),intent(in)::tag
global=local
end subroutine
subroutine seq_comm_getData(comm,iam)
integer,intent(in)::comm
integer,intent(out)::iam
iam=0
end subroutine
""" + "\n".join(routines) + """
end module
program probe
use ocean_probe
implicit none
type(mct_avect)::fractions(1)
integer::step
allocate(x2oacc_ox(1))
allocate(held%rAttr(1,2),x2oacc_ox(1)%rAttr(1,2),fractions(1)%rAttr(2,2),grid%data%rAttr(1,2))
x2oacc_ox(1)%rAttr=0
grid%data%rAttr=2._r8*shr_const_pi
fractions(1)%rAttr(1,:)=[.2_r8,.7_r8]
fractions(1)%rAttr(2,:)=[.8_r8,.3_r8]
call shr_solar_diag_configure(.true.)
do step=1,4
call shr_solar_diag_clock(20000229,real(step*1800,r8),1800._r8)
held%rAttr(1,:)=[10._r8,20._r8]*step
call prep_ocn_accum('probe',fractions)
enddo
if(x2oacc_ox_cnt/=4) stop 10
call prep_ocn_accum_avg('probe',fractions)
if(x2oacc_ox_cnt/=0) stop 11
if(any(abs(held%rAttr(1,:)-[25._r8,50._r8])>1.e-12_r8)) stop 12
end program
"""
        result = self.compile_orbit_probe(program=program, extra_sources=("share/util/shr_solar_diag_mod.F90",))
        rows = [line.split() for line in result.stdout.splitlines() if "SOLAR_BUFFER_V1" in line]
        self.assertEqual(len(rows), 5)
        self.assertEqual([row[4] for row in rows], ["ADD"] * 4 + ["SEND"])
        self.assertEqual(int(rows[-1][5]), 4)
        self.assertAlmostEqual(float(rows[-1][8]), 270000.)
        self.assertEqual(float(rows[-1][9]), 0.)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for atmosphere source probes")
    def test_solar_atmosphere_source(self):
        heating_source = (ROOT / "components/eam/src/physics/rrtmgp/radiation_utils.F90").read_text()
        heating = re.search(r"   subroutine calculate_heating_rate\(.*?   end subroutine calculate_heating_rate",
                            heating_source, re.DOTALL).group(0)
        import_source = (ROOT / "components/eam/src/cpl/atm_import_export.F90").read_text()
        capture = re.search(r"    if \(shr_solar_diag_enabled\) then.*?    endif", import_source, re.DOTALL).group(0)
        for implementation in ("rrtmg", "rrtmgp"):
            source = (ROOT / "components/eam/src/physics" / implementation / "radiation.F90").read_text()
            routine = re.search(r"\bsubroutine radiation_tend\b.*?end subroutine radiation_tend", source, re.DOTALL).group(0)
            self.assertIn("use shr_solar_diag_mod", routine)
            self.assertNotIn("solar_top(:,:) =>", routine)
            self.assertIn("pbuf_add_field('SOLAR_TOP','global'", source)
            block = re.search(r"(?m)^\s*if \(shr_solar_diag_enabled\) then.*?^\s*endif", routine, re.DOTALL).group(0)
            if implementation == "rrtmg":
                top_source = (ROOT / "components/eam/src/physics/rrtmg/radsw.F90").read_text()
                top_block = re.search(r"   if \(present\(solar_top\)\) then.*?   endif", top_source, re.DOTALL).group(0)
                top_setup = """
block
integer,parameter::Nday=2,rrtmg_levs=3
integer::IdxDay(2)=[1,2]
real(r8),parameter::cpair=1004._r8,dps=1._r8/86400._r8
real(r8)::fsntoa(2),swhr(2,3),pintmb(2,4)
fsntoa=310
swhr=0
swhr(:,3)=10._r8*gravit/(50._r8*cpair*dps)
pintmb=0
pintmb(:,3)=.5_r8
""" + top_block + "\nend block\n"
            else:
                self.assertIn("only: clip_values, calculate_heating_rate", routine)
                top_block = re.search(r"         if \(shr_solar_diag_requested\) then.*?         endif", routine, re.DOTALL).group(0)
                top_setup = """
block
integer,parameter::ktop=2
type flux_type
real(r8)::flux_up(2,4),flux_dn(2,4)
end type
type(flux_type)::fluxes_allsky
real(r8)::pint(2,4),diag_hr(2,ktop-1)
pint=reshape([0._r8,0._r8,50._r8,50._r8,50000._r8,50000._r8,100000._r8,100000._r8],[2,4])
fluxes_allsky%flux_up=0
fluxes_allsky%flux_dn=reshape([310._r8,310._r8,300._r8,300._r8,250._r8,250._r8,200._r8,200._r8],[2,4])
""" + top_block + "\nend block\n"
            top_routine = ("subroutine extra_layer(solar_top)\n"
                           "real(r8),optional,intent(out)::solar_top(:,:)\n" + top_setup + "\nend subroutine\n")
            top_setup = "call extra_layer(solar_top)\n"
            program = """
module physconst
use shr_kind_mod,only:r8=>shr_kind_r8
implicit none
real(r8),parameter::gravit=9.8_r8
end module
module heating_probe
use shr_kind_mod,only:r8=>shr_kind_r8
implicit none
contains
subroutine assert(condition,message)
use shr_sys_mod,only:shr_sys_abort
logical,intent(in)::condition
character(*),intent(in)::message
if(.not.condition) call shr_sys_abort(message)
end subroutine
""" + heating + """
end module
program probe
use shr_kind_mod,only:r8=>shr_kind_r8
use shr_const_mod,only:shr_const_pi
use shr_solar_diag_mod
use physconst,only:gravit
use heating_probe,only:calculate_heating_rate
implicit none
integer,parameter::pcols=2,pver=2,ncol=2,begchunk=1,endchunk=1
type input_type
real(r8)::asdir(2),asdif(2),aldir(2),aldif(2)
end type
type state_type
integer::lchnk=1
real(r8)::pdel(2,2)
end type
type(input_type)::cam_in
type(state_type)::state
integer::diag_year,diag_month,diag_day,diag_seconds
real(r8)::diag_area(2),fsnt(2),fsns(2),qrs(2,2),pint(2,3),up(2,3),down(2,3),solar_top(2,2)
cam_in%asdir=.1_r8
cam_in%asdif=.2_r8
cam_in%aldir=.3_r8
cam_in%aldif=.4_r8
call shr_solar_diag_configure(.true.)
call shr_solar_diag_schedule(20000229,3600._r8,.true.)
call capture_import([cam_in])
pint=reshape([0._r8,0._r8,50000._r8,50000._r8,100000._r8,100000._r8],[2,3])
up=0
down=reshape([300._r8,300._r8,250._r8,250._r8,200._r8,200._r8],[2,3])
call calculate_heating_rate(up,down,pint,qrs)
state%pdel=pint(:,2:3)-pint(:,1:2)
fsnt=300
fsns=200
""" + top_setup + block + """
contains
""" + top_routine + """
subroutine capture_import(cam_in)
type(input_type),intent(in)::cam_in(1)
integer::c,ncols
""" + capture + """
end subroutine
integer function get_ncols_p(chunk)
integer,intent(in)::chunk
get_ncols_p=2
end function
subroutine get_area_all_p(chunk,count,area)
integer,intent(in)::chunk,count
real(r8),intent(out)::area(count)
area=2._r8*shr_const_pi
end subroutine
subroutine get_curr_date(year,month,day,seconds)
integer,intent(out)::year,month,day,seconds
year=2000
month=2
day=29
seconds=3600
end subroutine
logical function radiation_do(operation)
character(*),intent(in)::operation
radiation_do=.true.
end function
end program
"""
            sources = ("share/util/shr_solar_diag_mod.F90",)
            result = self.compile_orbit_probe(program=program, extra_sources=sources)
            row = next(line.split() for line in result.stdout.splitlines() if "SOLAR_COLUMN_V1" in line)
            self.assertAlmostEqual(float(row[6]), 310.)
            self.assertAlmostEqual(float(row[7]), 200.)
            self.assertAlmostEqual(float(row[8]), 110.)
            self.assertAlmostEqual(float(row[10]), 1.)
            if implementation == "rrtmgp":
                no_extra = program.replace("integer,parameter::ktop=2", "integer,parameter::ktop=1")
                no_extra = no_extra.replace("fsnt=300", "fsnt=310").replace(
                    "down=reshape([300._r8,300._r8,", "down=reshape([310._r8,310._r8,")
                self.compile_orbit_probe(program=no_extra, extra_sources=sources)
            bad = program.replace("call capture_import([cam_in])", "call capture_import([cam_in])\ncam_in%asdir(2)=.9_r8")
            rejected = self.compile_orbit_probe(program=bad, extra_sources=sources, expected_returncode=99)
            self.assertIn("imported and consumed albedos differ", rejected.stdout)

    def test_solar_exchange_checker(self):
        spec = importlib.util.spec_from_file_location(
            "solar_conservation_check", ROOT / "cime_config/SystemTests/scv.py")
        checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(checker)
        for calendar in ("GREGORIAN", "NO_LEAP"):
            with self.subTest(calendar=calendar):
                start = checker.calendar_seconds(20000228, 84600, calendar)
                next_date = 20000229 if calendar == "GREGORIAN" else 20000301
                records = [(next_date, 0, [-100., 30., 50., 10., 10.]),
                           (next_date, 1800, [-120., 40., 60., 10., 10.])]
                result = checker.check_exchange(records, calendar, start, start + 3600, 1800)
                self.assertEqual(result["residual_j_m2"], 0)
                self.assertEqual(result["component_energy_j_m2"][0], -396000)
                for bad in (records[:1], records[::-1], records + records[-1:],
                            [(next_date, 0, [-100., 30., 60., 10., 10.]), records[1]],
                            [(next_date, 0, [100., -30., -50., -10., -10.]), records[1]],
                            [(next_date, 0, [0., 0., 0., 0., 0.]),
                             (next_date, 1800, [0., 0., 0., 0., 0.])],
                            [(next_date, 0, [-100., 30., math.nan, 10., 10.]), records[1]]):
                    with self.assertRaises(AssertionError):
                        checker.check_exchange(bad, calendar, start, start + 3600, 1800)
        self.assertEqual(checker.calendar_seconds(20000301, 0, "GREGORIAN") -
                         checker.calendar_seconds(20000228, 0, "GREGORIAN"), 172800)
        self.assertEqual(checker.calendar_seconds(21000301, 0, "GREGORIAN") -
                         checker.calendar_seconds(21000228, 0, "GREGORIAN"), 86400)
        with self.assertRaises(ValueError):
            checker.calendar_seconds(20000229, 0, "NO_LEAP")

        source = (ROOT / "driver-mct/main/seq_diag_mct.F90").read_text()
        block = re.search(r"          if \(ip == p_inst \.and\. plev >= 4\) then.*?          endif",
                  source, re.DOTALL).group(0)
        self.assertIn("dataGpr(f_hswnet,c_ocn_os,ip)", block)
        self.assertIn("dataGpr(f_hswnet,c_inh_ir,ip)", block)
        self.assertIn("dataGpr(f_hswnet,c_ish_ir,ip)", block)
        output = self.compile_orbit_probe(program="""
    program probe
    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_log_mod, only: logunit => shr_log_Unit
    implicit none
    integer, parameter :: p_inst=1, f_hswnet=1, c_atm_ar=1, c_lnd_lr=2
    integer, parameter :: c_ocn_os=3, c_inh_ir=4, c_ish_ir=5
    integer :: ip=1, plev=4, cdate=20000229, sec=1800
    type solar_context
    logical :: nextsw_valid=.true.
    integer :: nextsw_ymd=20000229, ocean_solar_version=0
    real(r8) :: nextsw_tod=5400, nextsw_cday=60._r8+5400._r8/86400, nextsw_dt_avg=1800
    end type
    type(solar_context) :: infodata
    real(r8) :: dataGpr(1,5,1)
    dataGpr(1,:,1)=[-100._r8,30._r8,50._r8,10._r8,10._r8]
    """ + block + "\nend program probe\n")
        self.assertEqual(checker.exchange_records(output.stdout),
                 [(20000229, 1800., [-100., 30., 50., 10., 10.])])
        self.assertEqual(checker.check_schedule(checker.schedule_records(output.stdout), "GREGORIAN",
                               checker.exchange_records(output.stdout), 1800, 1800), 1)
        for bad in (output.stdout.replace("T 20000229", "T 20000228"),
                output.stdout.replace("SOLAR_SCHEDULE_V1", "MISSING_SCHEDULE")):
            self.assertRaises(AssertionError, checker.check_schedule,
                              checker.schedule_records(bad), "GREGORIAN",
                              checker.exchange_records(output.stdout), 1800, 1800)
        with self.assertRaises(AssertionError):
            checker.exchange_records("SOLAR_EXCHANGE_V1 20000229 1800 1")

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for the Icepack source probe")
    def test_icepack_solar_partition(self):
        source = (ROOT / "components/mpas-seaice/src/icepack/columnphysics/icepack_shortwave.F90").read_text()
        routine = re.search(r"      subroutine icepack_prep_radiation\(.*?      end subroutine icepack_prep_radiation",
                            source, re.DOTALL)
        self.assertIsNotNone(routine)
        self.compile_orbit_probe(program="""
module icepack_probe
use shr_kind_mod, only: dbl_kind=>shr_kind_r8, int_kind=>shr_kind_in
implicit none
integer, parameter :: ncat=2,nilyr=2,nslyr=2
real(dbl_kind), parameter :: c0=0,c1=1,puny=1.e-11_dbl_kind
contains
""" + routine.group(0) + """
end module
program probe
use icepack_probe
implicit none
real(dbl_kind) :: area(2),surface(2),interior(2),transmitted(2),reference,netsw
real(dbl_kind) :: visdir(2),visdif(2),nirdir(2),nirdif(2)
real(dbl_kind) :: snow(2,2),ice(2,2),penetration(3,2),open_absorbed,ocean,retained,reflected
integer :: category
area=[.2_dbl_kind,.3_dbl_kind]
reference=112
surface=[.2_dbl_kind,.4_dbl_kind]*reference
interior=[.5_dbl_kind,.3_dbl_kind]*reference
transmitted=.3_dbl_kind*reference
visdir=.1_dbl_kind*transmitted
visdif=.2_dbl_kind*transmitted
nirdir=.3_dbl_kind*transmitted
nirdif=.4_dbl_kind*transmitted
do category=1,2
snow(:,category)=.1_dbl_kind*interior(category)
ice(:,category)=.4_dbl_kind*interior(category)
penetration(:,category)=transmitted(category)
enddo
call icepack_prep_radiation(.5_dbl_kind,area,60._dbl_kind,200._dbl_kind,40._dbl_kind,300._dbl_kind, &
 .7_dbl_kind,.6_dbl_kind,.4_dbl_kind,.3_dbl_kind,reference,surface,interior,transmitted, &
 visdir,visdif,nirdir,nirdif,penetration,snow,ice)
netsw=60*.3_dbl_kind+200*.4_dbl_kind+40*.6_dbl_kind+300*.7_dbl_kind
call require(maxval(abs(surface+interior+transmitted-netsw))<1.e-10_dbl_kind)
call require(abs(sum(area*(surface+interior+transmitted))-.5_dbl_kind*netsw)<1.e-10_dbl_kind)
call require(maxval(abs(visdir+visdif+nirdir+nirdif-transmitted))<1.e-10_dbl_kind)
do category=1,2
call require(abs(sum(snow(:,category))+sum(ice(:,category))-interior(category))<1.e-10_dbl_kind)
enddo
open_absorbed=60*.9_dbl_kind+200*.8_dbl_kind+40*.7_dbl_kind+300*.6_dbl_kind
ocean=.5_dbl_kind*open_absorbed+sum(area*transmitted)
retained=sum(area*(surface+interior))
reflected=.5_dbl_kind*(600-netsw)+.5_dbl_kind*(600-open_absorbed)
call require(abs(600-reflected-retained-ocean)<1.e-10_dbl_kind)
call require(abs(600-reflected-retained-ocean-sum(area*transmitted))>1._dbl_kind)
reference=netsw
call icepack_prep_radiation(.5_dbl_kind,area,0._dbl_kind,0._dbl_kind,0._dbl_kind,0._dbl_kind, &
 .7_dbl_kind,.6_dbl_kind,.4_dbl_kind,.3_dbl_kind,reference,surface,interior,transmitted, &
 visdir,visdif,nirdir,nirdif,penetration,snow,ice)
call require(maxval(abs(surface+interior+transmitted))==0._dbl_kind)
reference=0
call icepack_prep_radiation(.5_dbl_kind,area,60._dbl_kind,200._dbl_kind,40._dbl_kind,300._dbl_kind, &
 1._dbl_kind,1._dbl_kind,1._dbl_kind,1._dbl_kind,reference,surface,interior,transmitted, &
 visdir,visdif,nirdir,nirdif,penetration,snow,ice)
call require(reference==1._dbl_kind)
call require(maxval(abs(surface+interior+transmitted))==0._dbl_kind)
contains
subroutine require(condition)
logical,intent(in)::condition
if(.not.condition) stop 14
end subroutine
end program
""")

    def test_solar_full_oracles(self):
        spec = importlib.util.spec_from_file_location("solar_full_check", ROOT / "cime_config/SystemTests/scv.py")
        checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(checker)
        for calendar, next_date in (("GREGORIAN", 20000229), ("NO_LEAP", 20000301)):
            times = [(20000228, 82800), (20000228, 84600), (next_date, 0)]
            events = [(20000228, 84600), (next_date, 0), (next_date, 1800)]
            exchanges = [(date, seconds, [-100., 30., 50., 10., 10.]) for date, seconds in times]
            columns = [(date, seconds, 1, 1, 0., 300., 100., 200., 600., 1.) for date, seconds in events]
            consumed = [(date, seconds, 1, 1) for date, seconds in events]
            schedules = [(date, seconds, True, next_date, future, 60 + future/86400, 1800., 0.)
                         for (date, seconds), future in zip(times, (0, 1800, 3600))]
            result = checker.check_columns(columns, consumed, schedules, exchanges, calendar, 1800, 1800)
            self.assertEqual(result["toa_energy_j_m2"], 1620000)
            self.assertEqual(result["atmosphere_heating_j_m2"], 1080000)
            for bad_columns, bad_consumed in (
                    (columns[:-1], consumed), (columns, consumed[:-1]), (columns, consumed + consumed[:1]),
                    ([row[:7] + (199.,) + row[8:] for row in columns], consumed),
                    ([row[:-1] + (.5,) for row in columns], consumed),
                    (columns, [(next_date, 3600, 1, 1)] + consumed[1:])):
                self.assertRaises(AssertionError, checker.check_columns, bad_columns, bad_consumed,
                                  schedules, exchanges, calendar, 1800, 1800)
            buffers = [
                (*times[0], 1, "ADD", 2, 90000., 90000., 0., 180000.),
                (*times[0], 1, "SEND", 2, 180000., 0., 180000., 0.),
                (*times[1], 1, "ADD", 1, 0., 90000., 0., 90000.),
                (*times[2], 1, "ADD", 2, 90000., 90000., 0., 180000.)]
            ledger = checker.check_buffers(buffers, exchanges, calendar, 1800, 3600)
            self.assertEqual(ledger["ocean_pending_start_j_m2"], 90000)
            self.assertEqual(ledger["ocean_pending_end_j_m2"], 180000)
            self.assertEqual(ledger["ocean_delivered_j_m2"], 180000)
            self.assertEqual(ledger["ocean_incoming_j_m2"] - ledger["ocean_delivered_j_m2"],
                             ledger["ocean_pending_change_j_m2"])
            for bad in (buffers[:-1], buffers[:2] + buffers[3:], buffers + buffers[-1:],
                        [row[:5] + tuple(value*2 for value in row[5:]) for row in buffers]):
                self.assertRaises(AssertionError, checker.check_buffers, bad, exchanges, calendar, 1800, 3600)

    def test_solar_suite_definitions(self):
        for suffix, count in (("smoke", 4), ("boundaries", 28), ("cadence", 4), ("unequal", 12)):
            tests = get_test_data("e3sm_solar_conservation_" + suffix)[-1]
            self.assertEqual(len(tests), count)
            self.assertEqual(len(set(tests)), count)
            for test in tests:
                name, modifiers = test.rsplit(".", 1)
                full_name = get_full_test_name(name, machine="pm-cpu", compiler="gnu", testmods_string=modifiers)
                paths = find_test_mods("mct", parse_test_name(full_name)[6])
                self.assertEqual(len(paths), 3 + int("solar_conservation_rad2" in modifiers) +
                                 int("solar_conservation_ocn" in modifiers))
                self.assertTrue(test.startswith("SCV_Ld"))
                self.assertIn("WCYCL1850NS", name)
                for path in paths:
                    self.assertTrue(any((Path(path) / filename).is_file() for filename in
                                        ("shell_commands", "user_nl_cpl", "user_nl_eam")))
        case = Mock()
        case.get_compset_components.return_value = ["eam", "elm", "mpaso", "mpassi"]
        directories = {"any": str(ROOT / "cime/CIME/SystemTests"),
                       "allactive": str(ROOT / "cime_config/SystemTests")}
        case.get_value.side_effect = lambda name, attribute: directories.get(attribute["component"])
        test_class = find_system_test("SCV", case)
        self.assertEqual(test_class.__name__, "SCV")
        self.assertTrue(issubclass(test_class, SystemTestsCommon))
        generator = self.generator("driver-mct", "GREGORIAN", user_namelist=(
            ROOT / "cime_config/testmods_dirs/allactive/solar_conservation/user_nl_cpl"))
        self.assertEqual(generator.get_value("budget_inst"), "4")

    def test_solar_runner(self):
        spec = importlib.util.spec_from_file_location("solar_runner_check", ROOT / "cime_config/SystemTests/scv.py")
        checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(checker)
        with tempfile.TemporaryDirectory() as directory:
            values = dict(COMP_INTERFACE="mct", COMPSET="WCYCL1850NS", NINST=1,
                          NCPL_BASE_PERIOD="day", ATM_NCPL=48, LND_NCPL=48, ICE_NCPL=48,
                          OCN_NCPL=48, BUDGETS=True, STOP_OPTION="ndays", STOP_N=2,
                          CONTINUE_RUN=False, CALENDAR="GREGORIAN", RUN_STARTDATE="2000-02-28",
                          START_TOD=0, RUNDIR=directory, CASE="SCV.solar_conservation")
            runner = object.__new__(checker.SCV)
            runner._case = Mock()
            runner._case.get_value.side_effect = values.get
            runner._case.get_values.return_value = ["CPL", "ATM", "LND", "ICE", "OCN", "ROF"]
            runner._rest_time = ".2000-02-29-00000"
            log_numbers = iter(range(20))
            lose_pending = False

            def emit(*args, **kwargs):
                suffix = kwargs.get("suffix", "base")
                start = datetime.datetime(2000, 2, 29 if suffix == "rest" else 28)
                lines = []
                ratio = 48 // values["OCN_NCPL"]
                pending = 0 if lose_pending else ratio - 1
                for step in range(1, values["STOP_N"] * 48 + 1):
                    stamp = start + datetime.timedelta(seconds=step * 1800)
                    future = stamp + datetime.timedelta(seconds=3600)
                    date = stamp.strftime("%Y%m%d")
                    seconds = stamp.hour * 3600 + stamp.minute * 60
                    next_seconds = future.hour * 3600 + future.minute * 60
                    cday = future.timetuple().tm_yday + next_seconds / 86400
                    lines.append("SOLAR_EXCHANGE_V1 {} {} -100 30 50 10 10".format(date, seconds))
                    lines.append("SOLAR_SCHEDULE_V1 {} {} T {} {} {} 1800 {}".format(
                        date, seconds, future.strftime("%Y%m%d"), next_seconds, cday, int(suffix == "rest")))
                    event = stamp + datetime.timedelta(seconds=1800)
                    event_seconds = event.hour * 3600 + event.minute * 60
                    lines.append("SOLAR_COLUMN_V1 {} {} 1 1 0 300 100 200 600 1".format(
                        event.strftime("%Y%m%d"), event_seconds))
                    lines.append("SOLAR_CONSUME_V1 {} {} 1 1".format(event.strftime("%Y%m%d"), event_seconds))
                    count = (step - 1 + pending) % ratio + 1
                    lines.append("SOLAR_BUFFER_V1 {} {} 1 ADD {} {} 90000 0 {}".format(
                        date, seconds, count, (count - 1) * 90000, count * 90000))
                    if count == ratio:
                        lines.append("SOLAR_BUFFER_V1 {} {} 1 SEND {} {} 0 {} 0".format(
                            date, seconds, count, count * 90000, count * 90000))
                (Path(directory) / ("cpl.log.{}.{}".format(suffix, next(log_numbers)))).write_text("\n".join(lines) + "\n")

            with patch.object(checker.ERS, "run_indv", side_effect=emit):
                for frequency in (48, 24, 12):
                    values.update(CONTINUE_RUN=False, STOP_N=2, OCN_NCPL=frequency)
                    runner.run_indv()
                    values.update(CONTINUE_RUN=True, STOP_N=1)
                    runner.run_indv(suffix="rest")
                    for suffix in ("base", "rest"):
                        report = json.loads((Path(directory) / ("solar-conservation." + suffix + ".json")).read_text())
                        self.assertEqual(report["ocean_pending_start_j_m2"], (48 // frequency - 1) * 90000)
                        self.assertEqual(report["ocean_pending_end_j_m2"], report["ocean_pending_start_j_m2"])
                lose_pending = True
                self.assertRaisesRegex(AssertionError, "Restart solar energy", runner.run_indv, suffix="rest")
            self.assertTrue((Path(directory) / "solar-conservation.base.json").is_file())
            self.assertTrue((Path(directory) / "solar-conservation.rest.json").is_file())
            with patch.object(checker.ERS, "run_indv"):
                self.assertRaises(AssertionError, runner.run_indv, suffix="rest")
            values["OCN_NCPL"] = 25
            self.assertRaises(ValueError, runner.run_indv)
            values.update(OCN_NCPL=12, ROF_NCPL=96)
            self.assertRaisesRegex(ValueError, "faster than the atmosphere", runner.run_indv)
            log_path = Path(directory) / "cpl.log.stream"
            log_path.write_text("before\n")
            before = checker.coupler_logs(directory)
            self.assertIsInstance(before[log_path][1], bytes)
            log_path.write_text("before\nafter\n")
            self.assertEqual(list(checker.fresh_log_lines(directory, before)), ["after\n"])
            compressed_path = log_path.with_suffix(log_path.suffix + ".gz")
            with gzip.open(compressed_path, "wt") as stream:
                stream.write("before\nafter\n")
            self.assertEqual(list(checker.fresh_log_lines(directory, before)), ["after\n"])
            log_path.unlink()
            self.assertEqual(list(checker.fresh_log_lines(directory, before)), ["after\n"])
            log_path.write_text("before\nconflict\n")
            self.assertRaisesRegex(AssertionError, "Conflicting compressed", checker.coupler_logs, directory)
            self.assertRaisesRegex(AssertionError, "Conflicting compressed", list,
                                   checker.fresh_log_lines(directory, before))
            log_path.unlink()
            with gzip.open(compressed_path, "wt") as stream:
                stream.write("changed\n")
            self.assertRaisesRegex(AssertionError, "replaced during SCV", list,
                                   checker.fresh_log_lines(directory, before))
            with gzip.open(compressed_path, "wt") as stream:
                stream.write("short")
            self.assertRaisesRegex(AssertionError, "replaced during SCV", list,
                                   checker.fresh_log_lines(directory, before))
            self.assertEqual(checker.read_solar_records(iter(["unrelated " * 10000])), [[], [], [], [], []])

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
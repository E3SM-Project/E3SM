#!/usr/bin/env python3
"""Compile focused irrigation regressions without a CIME case.

Run with FC=gfortran python3 test_irrigation_coupling.py. The harness compiles
the production allocation module and extracts the production clock expressions,
alignment checks, and MOSART direct-transfer block. It does not replace a coupled
smoke/budget/restart test on a supported E3SM machine.
"""

import itertools
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[4]


def source(path):
    return (ROOT / path).read_text()


def phase_body(driver):
    text = source(f"{driver}/shr/seq_timemgr_mod.F90")
    start = text.index("    irrigation_land_offset =")
    end = text.index("\n\n", start)
    return text[start:end]


class IrrigationCouplingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shlex.split(os.environ.get("FC", "gfortran"))
        if not shutil.which(compiler[0]):
            raise unittest.SkipTest("Set FC to an available GNU Fortran compiler")
        cls.temp = tempfile.TemporaryDirectory(prefix="irrigation-regression-")
        cls.addClassCleanup(cls.temp.cleanup)
        cls.work = Path(cls.temp.name)
        (cls.work / "shr_kind_mod.F90").write_text(
            "module shr_kind_mod\n"
            "! !DESCRIPTION: Real kind stub for the standalone irrigation test.\n"
            "integer, parameter :: shr_kind_r8 = kind(1.d0)\n"
            "end module\n"
        )
        canopy = source("components/elm/src/biogeophys/CanopyHydrologyMod.F90")
        guards = canopy[canopy.index("    if (river_dt <= 0") : canopy.index("    begp = bounds%begp")]
        routing = source("components/mosart/src/riverroute/RtmMod.F90")
        transfer = routing[
            routing.index("          !---- all dto water") :
            routing.index("       enddo\n    enddo\n    call mct_avect_zero(avdst_direct)")
        ]
        template = r"""
module harness
! !DESCRIPTION: Execute production clock and direct-transfer code in isolation.
use shr_kind_mod, only: r8 => shr_kind_r8
use IrrigationAllocationMod
implicit none
integer, parameter :: SHR_KIND_R8=r8, seq_timemgr_nclock_lnd=1, seq_timemgr_nclock_rof=2
integer, parameter :: nt_nliq=1
integer :: step_dt
logical :: rejected
type routing_type
  real(r8) :: qdto(1,2), qgwl(1,2), qsub(1,2), qsur(1,2), wr(1,2)
  real(r8) :: area(1), qirr_return(1)
  integer :: mask(1)
  logical :: euler_calc(2)
end type
type vector_type
  real(r8) :: rAttr(2,1)
end type
contains
integer function get_step_size()
! !DESCRIPTION: Supply the configured land step to the production guards.
get_step_size=step_dt
end function
subroutine endrun(message)
! !DESCRIPTION: Record a rejected interval configuration without exiting.
character(*), intent(in) :: message ! production error message [-]
rejected=.true.
end subroutine
subroutine check_alignment(land_time,land_dt,river_dt,river_phase)
! !DESCRIPTION: Execute the production interval guards.
real(r8), intent(in) :: land_time ! land clock time since reference date [days]
integer, intent(in) :: land_dt ! land coupling interval [s]
integer, intent(in) :: river_dt ! river coupling interval [s]
integer, intent(in) :: river_phase ! river-run phase on the land clock [s]
@GUARDS@
end subroutine
subroutine phase_mct(irrigation_clock_time,lnd_cpl_dt,rof_cpl_dt,offset,seq_timemgr_irr_phase)
! !DESCRIPTION: Execute the MCT production phase calculation.
real(r8), intent(in) :: irrigation_clock_time ! initialization time since reference date [days]
integer, intent(in) :: lnd_cpl_dt, rof_cpl_dt ! land and river intervals [s]
integer, intent(in) :: offset(2) ! driver alarm offsets including driver step [s]
integer, intent(out) :: seq_timemgr_irr_phase ! river-run phase on the land clock [s]
integer :: irrigation_land_offset
@MCT@
end subroutine
subroutine phase_moab(irrigation_clock_time,lnd_cpl_dt,rof_cpl_dt,offset,seq_timemgr_irr_phase)
! !DESCRIPTION: Execute the MOAB production phase calculation.
real(r8), intent(in) :: irrigation_clock_time ! initialization time since reference date [days]
integer, intent(in) :: lnd_cpl_dt, rof_cpl_dt ! land and river intervals [s]
integer, intent(in) :: offset(2) ! driver alarm offsets including driver step [s]
integer, intent(out) :: seq_timemgr_irr_phase ! river-run phase on the land clock [s]
integer :: irrigation_land_offset
@MOAB@
end subroutine
subroutine route(gwl,returned,mask,inundflag,direct,channel)
! !DESCRIPTION: Execute the production MOSART direct-transfer block.
real(r8), intent(in) :: gwl ! ordinary glacier/wetland/lake runoff [m3/s]
real(r8), intent(in) :: returned ! unused irrigation delivery [m3/s]
integer, intent(in) :: mask ! river routing cell mask [-]
logical, intent(in) :: inundflag ! inundation scheme flag [-]
real(r8), intent(out) :: direct ! liquid flow sent directly to outlet [m3/s]
real(r8), intent(out) :: channel ! liquid flow retained for channel routing [m3/s]
type(routing_type) :: rtmCTL, TRunoff, TUnit
type(vector_type) :: avsrc_direct
integer :: nr,nt,cnt
real(r8) :: qgwl_volume,river_volume_minimum
real(r8), parameter :: delt_mosart=900._r8,river_depth_minimum=1.e-6_r8
rtmCTL%qdto=0._r8
rtmCTL%qirr_return=returned
rtmCTL%mask=mask
rtmCTL%area=1.e6_r8
TRunoff%qgwl=0._r8
TRunoff%qgwl(1,1)=gwl
TRunoff%qsur=0._r8
TRunoff%qsub=0._r8
TRunoff%wr=1.e9_r8
TUnit%mask=mask
TUnit%euler_calc=.true.
avsrc_direct%rAttr=0._r8
nr=1
cnt=1
do nt=1,2
@TRANSFER@
enddo
direct=avsrc_direct%rAttr(1,1)
channel=TRunoff%qgwl(1,1)
if (TRunoff%qgwl(1,2)/=0._r8) error stop 'liquid return entered the ice tracer'
if (any(TRunoff%qsur/=0._r8).or.any(TRunoff%qsub/=0._r8)) &
  error stop 'return entered hillslope or tributary forcing'
end subroutine
end module
program test
use harness
implicit none
real(r8) :: start_seconds,land_time,gwl,returned,direct,channel
integer :: mode,ios,ld,rd,dd,lo,ro,phase,other,i,mask,wet
integer :: closes(18),offset(2)
do
  read(*,*,iostat=ios) mode
  if (ios/=0) exit
  if (mode==1) then
    read(*,*) start_seconds,dd,ld,rd,lo,ro
    offset=[dd+lo,dd+ro]
    call phase_mct(start_seconds/86400._r8,ld,rd,offset,phase)
    call phase_moab(start_seconds/86400._r8,ld,rd,offset,other)
    if (phase/=other) error stop 'MCT and MOAB disagree'
    step_dt=ld
    rejected=.false.
    do i=1,18
      land_time=(start_seconds+real(i,r8)*ld)/86400._r8
      call check_alignment(land_time,ld,rd,phase)
      closes(i)=merge(1,0,IrrigationWindowCloses(land_time,rd,phase))
    enddo
    write(*,*) merge(1,0,rejected),closes
  else
    read(*,*) gwl,returned,mask,wet
    call route(gwl,returned,mask,wet==1,direct,channel)
    write(*,*) direct,channel
  endif
enddo
end program
"""
        for key, value in {
            "GUARDS": guards,
            "MCT": phase_body("driver-mct"),
            "MOAB": phase_body("driver-moab"),
            "TRANSFER": transfer,
        }.items():
            template = template.replace(f"@{key}@", value)
        (cls.work / "test.F90").write_text(template)
        subprocess.run(
            compiler + ["-ffree-line-length-none", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                        "shr_kind_mod.F90", str(ROOT / "components/elm/src/biogeophys/IrrigationAllocationMod.F90"),
                        "test.F90", "-o", "test"],
            cwd=cls.work, check=True, capture_output=True, text=True,
        )

    def run_cases(self, cases):
        result = subprocess.run(
            [str(self.work / "test")], input="\n".join(cases) + "\n",
            capture_output=True, text=True, check=True,
        )
        return result.stdout.splitlines()

    def test_clock_matches_driver_alarm_schedule(self):
        cases, expected = [], []
        # Includes the reported 900/1800/10800 case, unequal user offsets,
        # rejected non-overlapping alarms, arbitrary restart seconds and long runs.
        for start, dd, ld, rd in itertools.product(
            (0, 1234, 370 * 365 * 86400 + 1234), (900, 1800), (1800, 3600), (1800, 10800)
        ):
            for lo, ro in itertools.product((-ld, -dd, 0, dd, ld), (-rd, 0, ld, rd)):
                def first_alarm(dt, offset):
                    alarm = offset + dd - dt
                    while alarm <= 0:
                        alarm += dt
                    return alarm
                first_land = first_alarm(ld, lo)
                first_river = first_alarm(rd, ro)
                rejected = rd % ld != 0 or (first_river - first_land) % ld != 0
                closes = [int((first_land + i * ld - first_river) % rd == 0) for i in range(18)]
                cases += ["1", f"{start} {dd} {ld} {rd} {lo} {ro}"]
                expected.append([int(rejected)] + closes)
        actual = [[int(x) for x in line.split()] for line in self.run_cases(cases)]
        self.assertEqual(actual, expected)

    def test_return_reaches_channel_without_changing_gwl_transfer(self):
        cases, expected = [], []
        for gwl, returned, mask, wet in itertools.product((-3, 0, 4), (0, 5), (0, 1, 3), (0, 1)):
            cases += ["2", f"{gwl} {returned} {mask} {wet}"]
            routed = returned if mask > 0 else 0
            expected.append((gwl + returned - routed, routed))
        actual = [tuple(float(x) for x in line.split()) for line in self.run_cases(cases)]
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()

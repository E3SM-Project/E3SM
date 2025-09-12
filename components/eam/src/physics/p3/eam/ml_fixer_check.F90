module ml_fixer_check

contains

subroutine ml_fixer_calc(dt,qc,nc,qr,nr,qctend,nctend,qrtend,nrtend,fixer,qc_fixer, nc_fixer, qr_fixer, nr_fixer)

use shr_kind_mod,   only: r8=>shr_kind_r8
use physconst,      only: pi, rhoh2o

real(r8), intent(in) :: dt
real(r8), intent(in) :: qc
real(r8), intent(in) :: nc
real(r8), intent(in) :: qr
real(r8), intent(in) :: nr
real(r8), intent(inout) :: qctend
real(r8), intent(inout) :: nctend
real(r8), intent(inout) :: qrtend
real(r8), intent(inout) :: nrtend

real(r8), intent(out) :: qc_fixer
real(r8), intent(out) :: nc_fixer
real(r8), intent(out) :: qr_fixer
real(r8), intent(out) :: nr_fixer

real(r8), intent(out) :: fixer

real(r8) :: qc_tmp, nc_tmp, qr_tmp, nr_tmp
integer :: i

fixer = 0._r8

qc_fixer = 0._r8
qr_fixer = 0._r8
nc_fixer = 0._r8
nr_fixer = 0._r8

qc_tmp = qc+qctend*dt
nc_tmp = nc+nctend*dt
qr_tmp = qr+qrtend*dt
nr_tmp = nr+nrtend*dt

if( qc_tmp.lt.0._r8 ) then
   fixer = 1._r8
   qctend = -qc/dt
   qrtend = qc/dt
   nctend = -nc/dt   
end if
if( qr_tmp.lt.0._r8 ) then
   fixer = 1._r8
   qrtend = -qr/dt
   qctend = qr/dt
   nrtend = -nr/dt   
end if
if( nc_tmp.lt.0._r8 ) then
   fixer = 1._r8
   if( qc_tmp.gt.0._r8 ) then
      nc_tmp = qc_tmp/(4._r8/3._r8*pi*(5.e-5_r8)**3._r8*rhoh2o)
      nctend = (nc_tmp-nc)/dt
   else
      nctend = -nc/dt
   end if
end if
if( nr_tmp.lt.0._r8 ) then
   fixer = 1._r8
   if(qr_tmp.gt.0._r8) then
      nr_tmp = qr_tmp/(4._r8/3._r8*pi*(5.e-5_r8)**3._r8*rhoh2o)
      nrtend = (nr_tmp-nr)/dt
   else
      nrtend = -nr/dt
   end if
end if

qc_fixer = qc+qctend*dt-qc_tmp
qr_fixer = qr+qrtend*dt-qr_tmp
nc_fixer = nc+nctend*dt-nc_tmp
nr_fixer = nr+nrtend*dt-nr_tmp

end subroutine ml_fixer_calc

end module ml_fixer_check

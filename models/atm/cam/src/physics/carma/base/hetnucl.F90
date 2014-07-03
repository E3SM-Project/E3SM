! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!!  heterogeneous deposition nucleation only. The parameters are adjusted
!! for mesospheric conditions, based upon the recommendations of Keesee.
!!
!!  Based on expressions from ...
!!    Keesee [JGR,1989]
!!    Pruppacher and Klett [2000]
!!    Rapp and Thomas [JASTP, 2006]
!!    Trainer et al. [2008]
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! To avoid nucleation into an evaporating bin, this subroutine must
!! be called after growp, which evaluates evaporation loss rates <evaplg>.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Oct-2000, Jan-2010
subroutine hetnucl(carma, cstate, iz, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: igas    ! gas index
  integer                              :: igroup  ! group index
  integer                              :: ibin    ! bin index
  integer                              :: iepart  ! element for condensing group index
  integer                              :: inuc    ! nucleating element index
  integer                              :: ienucto ! index of target nucleation element
  integer                              :: ignucto ! index of target nucleation group
  real(kind=f)                         :: rmw
  real(kind=f)                         :: R_H2O
  real(kind=f)                         :: rnh2o
  real(kind=f)                         :: rlogs
  real(kind=f)                         :: ag
  real(kind=f)                         :: contang
  real(kind=f)                         :: xh
  real(kind=f)                         :: phih
  real(kind=f)                         :: rath
  real(kind=f)                         :: fv3h
  real(kind=f)                         :: fv4h
  real(kind=f)                         :: fh
  real(kind=f)                         :: delfg
  real(kind=f)                         :: expon

  ! Heterogeneous nucleation factors
  real(kind=f), parameter              :: gdes    = 2.9e-13_f
  real(kind=f), parameter              :: gsd     = 2.9e-14_f
  real(kind=f), parameter              :: zeld    = 0.1_f
  real(kind=f), parameter              :: vibfreq = 1.e13_f
  real(kind=f), parameter              :: diflen  = 0.1e-7_f
  real(kind=f)                         :: rmiv
      
  rmiv    = 0.95_f

  ! rmiv - Eq. 2, Trainer et al. [2008]
!  rmiv = 0.94_f - (6005._f * exp(-0.065_f * max(150._f, t(iz))))
!  rmiv = max(0._f, 0.94_f - (6005._f * exp(-0.065_f * t(iz))))

  !  Loop over particle groups.
  do igroup = 1, NGROUP

    igas = inucgas(igroup)                ! condensing gas
 
    if (igas .ne. 0) then

      iepart = ienconc(igroup)              ! particle number density element

      rmw = gwtmol(igas) / AVG
      R_H2O = RGAS / gwtmol(igas)
      rnh2o = gc(iz,igas) * R_H2O / BK

      ! Calculate nucleation loss rates.  Do not allow nucleation into
      ! an evaporating bin.
      !
      ! <ienucto> is index of target nucleation element;
      ! <ignucto> is index of target nucleation group.
      do inuc = 1, nnuc2elem(iepart)

        ienucto = inuc2elem(inuc,iepart)
        
        if (ienucto .ne. 0) then
          ignucto = igelem(ienucto)
        else
          ignucto = 0
        endif
  
        ! Only compute nucleation rate for heterogenous nucleation
        if (inucproc(iepart,ienucto) .eq. I_HETNUC) then
  
          ! Loop over particle bins.  Loop from largest to smallest for 
          ! evaluation of index of smallest bin nucleated during time step <inucstep>.
          do ibin = NBIN, 1, -1
    
            ! Bypass calculation if few particles are present
            if (pconmax(iz,igroup) .gt. FEW_PC) then
    
              ! Only proceed if ice supersaturated
              !
              ! NOTE: We are only trying to model PMC partcles, so turn of nucleation
              ! where the CAM microphysics takes over (~1 mb = 1000 dyne).
              if ((p(iz) .lt. 1.e3_f) .and. (supsati(iz,igas) .gt. 0._f)) then
                rlogs = log(supsati(iz,igas) + 1._f)
      
                ! Critical ice germ radius formed in the sulfate solution
                !
                !   Eq. 2, Rapp & Thomas [2006]
                ag = 2._f * gwtmol(igas) * surfctia(iz) / rgas / t(iz) / RHO_I / rlogs
      
                ! Heterogeneous nucleation geometric factor
                !
                !  Eq. 9-22, Pruppacher & Klett [2000]
                contang = acos(rmiv)
                xh = r(ibin,igroup) / ag
                phih = sqrt(1._f - 2._f * rmiv * xh + xh**2 )
                rath = (xh-rmiv) / phih
                fv3h = xh**3 * (2._f - 3._f * rath + rath**3 )
                fv4h = 3._f * rmiv * xh**2 * (rath - 1._f)
                
                if (abs(rath) .gt. 1._f - 1.e-8_f)  fv3h = 0._f
                if (abs(rath) .gt. 1._f - 1.e-10_f) fv4h = 0._f
      
                fh = 0.5_f * (1._f + ((1._f - rmiv * xh) / phih)**3 + fv3h + fv4h)
      
                ! Gibbs free energy of ice germ formation in the ice/sulfate solution
                !
                !  Eq. 3, Rapp & Thomas [2006]
                delfg = 4._f * PI * ag**2 * surfctia(iz) - 4._f * PI * RHO_I * ag**3 *BK * t(iz) * rlogs / 3._f / rmw
      
                ! Ice nucleation rate in a 0.2 micron aerosol (/sec)
                expon = (2._f * gdes - gsd - fh*delfg) / BK / t(iz)

                ! NOTE: Excessive nucleation makes it difficult for the substepping to find a
                ! stable solution, so put a cap on really large nucleation values that can be produced.
                rnuclg(ibin,igroup,ignucto) = min(1e10_f, zeld * BK * t(iz) * diflen * ag * sin(contang) * &
                  4._f * PI * r(ibin,igroup)**2 * rnh2o**2 / (fh * rmw * vibfreq) * exp(expon))
              endif
            endif   ! pconmax(ixyz,igroup) .gt. FEW_PC 
          enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_DROPACT
      enddo        ! inuc = 1,nnuc2elem(iepart)
    endif        ! (igas = inucgas(igroup) .ne. 0)
  enddo         ! igroup = 1,NGROUP

  !  Return to caller with particle loss rates due to nucleation evaluated.
  return
end

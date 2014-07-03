!----------------------------------------------------------------------------------
! low level utility module for cloud aerosols
!
! Created by Francis Vitt
!----------------------------------------------------------------------------------
module cldaero_mod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver

  implicit none
  private

  public :: cldaero_uptakerate
  public :: cldaero_conc_t
  public :: cldaero_allocate
  public :: cldaero_deallocate

  type cldaero_conc_t
     real(r8), pointer :: so4c(:,:)
     real(r8), pointer :: nh4c(:,:)
     real(r8), pointer :: no3c(:,:)
     real(r8), pointer :: xlwc(:,:)
     real(r8) :: so4_fact
  end type cldaero_conc_t

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function cldaero_allocate( ) result( cldconc )
    type(cldaero_conc_t), pointer:: cldconc

    allocate( cldconc )
    allocate( cldconc%so4c(pcols,pver) )
    allocate( cldconc%nh4c(pcols,pver) )
    allocate( cldconc%no3c(pcols,pver) )
    allocate( cldconc%xlwc(pcols,pver) )

    cldconc%so4c(:,:) = 0._r8
    cldconc%nh4c(:,:) = 0._r8
    cldconc%no3c(:,:) = 0._r8
    cldconc%xlwc(:,:) = 0._r8
    cldconc%so4_fact  = 2._r8

  end function cldaero_allocate

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  subroutine cldaero_deallocate( cldconc )
    type(cldaero_conc_t), pointer :: cldconc

    if ( associated(cldconc%so4c) ) then
       deallocate(cldconc%so4c)
       nullify(cldconc%so4c)
    endif

    if ( associated(cldconc%nh4c) ) then
       deallocate(cldconc%nh4c)
       nullify(cldconc%nh4c)
    endif
    
    if ( associated(cldconc%no3c) ) then
       deallocate(cldconc%no3c)
       nullify(cldconc%no3c)
    endif

    if ( associated(cldconc%xlwc) ) then
       deallocate(cldconc%xlwc)
       nullify(cldconc%xlwc)
    endif

    deallocate( cldconc )
    nullify( cldconc )

  end subroutine cldaero_deallocate

!----------------------------------------------------------------------------------
! utility function for cloud-borne aerosols
!----------------------------------------------------------------------------------

  function cldaero_uptakerate( xl, cldnum, cfact, cldfrc, tfld,  press ) result( uptkrate )
    use mo_constants, only : pi

    real(r8), intent(in) :: xl, cldnum, cfact, cldfrc, tfld,  press

    real(r8) :: uptkrate

    real(r8) :: &
         rad_cd, radxnum_cd, num_cd, &
         gasdiffus, gasspeed, knudsen, &
         fuchs_sutugin, volx34pi_cd

!-----------------------------------------------------------------------
! compute uptake of h2so4 and msa to cloud water
!
! first-order uptake rate is
! 4*pi*(drop radius)*(drop number conc)
! *(gas diffusivity)*(fuchs sutugin correction)

! num_cd = (drop number conc in 1/cm^3)
        num_cd = 1.0e-3_r8*cldnum*cfact/cldfrc
        num_cd = max( num_cd, 0.0_r8 )

! rad_cd = (drop radius in cm), computed from liquid water and drop number,
! then bounded by 0.5 and 50.0 micrometers
! radxnum_cd = (drop radius)*(drop number conc)
! volx34pi_cd = (3/4*pi) * (liquid water volume in cm^3/cm^3)

        volx34pi_cd = xl*0.75_r8/pi

! following holds because volx34pi_cd = num_cd*(rad_cd**3)
        radxnum_cd = (volx34pi_cd*num_cd*num_cd)**0.3333333_r8

! apply bounds to rad_cd to avoid the occasional unphysical value
        if (radxnum_cd .le. volx34pi_cd*4.0e4_r8) then
            radxnum_cd = volx34pi_cd*4.0e4_r8
            rad_cd = 50.0e-4_r8
        else if (radxnum_cd .ge. volx34pi_cd*4.0e8_r8) then
            radxnum_cd = volx34pi_cd*4.0e8_r8
            rad_cd = 0.5e-4_r8
        else
            rad_cd = radxnum_cd/num_cd
        end if

! gasdiffus = h2so4 gas diffusivity from mosaic code (cm^2/s)
! (pmid must be Pa)
        gasdiffus = 0.557_r8 * (tfld**1.75_r8) / press

! gasspeed = h2so4 gas mean molecular speed from mosaic code (cm/s)
        gasspeed = 1.455e4_r8 * sqrt(tfld/98.0_r8)

! knudsen number
        knudsen = 3.0_r8*gasdiffus/(gasspeed*rad_cd)

! following assumes accomodation coefficient = 0.65
! (Adams & Seinfeld, 2002, JGR, and references therein)
! fuchs_sutugin = (0.75*accom*(1. + knudsen)) /
! (knudsen*(1.0 + knudsen + 0.283*accom) + 0.75*accom)
        fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) / &
                        (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)

! instantaneous uptake rate
        uptkrate = 12.56637_r8*radxnum_cd*gasdiffus*fuchs_sutugin

  end function cldaero_uptakerate

end module cldaero_mod

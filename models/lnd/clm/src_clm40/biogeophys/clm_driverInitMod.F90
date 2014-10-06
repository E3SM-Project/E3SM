
module clm_driverInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_driverInitMod
!
! !DESCRIPTION:
! Initialization of clm driver variables needed from previous timestep
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_driverInit
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_driverInit
!
! !INTERFACE:
  subroutine clm_driverInit(lbc, ubc, lbp, ubp, &
             num_nolakec, filter_nolakec, num_lakec, filter_lakec)
!
! !DESCRIPTION:
! Initialization of clm driver variables needed from previous timestep
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevsno
    use subgridAveMod, only : p2c
    use clm_varcon, only    : h2osno_max, rair, cpair, grav, istice_mec, lapse_glcmec
    use clm_atmlnd, only    : clm_a2l
    use domainMod, only     : ldomain
    use clmtype
    use QsatMod,    only    : Qsat

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_lakec                   ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1)     ! column filter for non-lake points
!
! !CALLED FROM:
! subroutine driver1
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in variables
!
    real(r8), pointer :: pwtgcell(:)           ! weight of pft wrt corresponding gridcell
    integer , pointer :: snl(:)                ! number of snow layers
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    integer , pointer :: frac_veg_nosno_alb(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer , pointer :: frac_veg_nosno(:)     ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (pft-level)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
!
! local pointers to original implicit out variables
!
    logical , pointer :: do_capsnow(:)         ! true => do snow capping
    real(r8), pointer :: h2osno_old(:)         ! snow water (mm H2O) at previous time step
    real(r8), pointer :: frac_iceold(:,:)      ! fraction of ice relative to the tot water
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer :: g, l, c, p, f, j, fc            ! indices

    real(r8), pointer :: qflx_glcice(:)     ! flux of new glacier ice (mm H2O/s) [+ = ice grows]
    real(r8), pointer :: eflx_bot(:)        ! heat flux from beneath soil/ice column (W/m**2)
    real(r8), pointer :: glc_topo(:)        ! sfc elevation for glacier_mec column (m)
    real(r8), pointer :: forc_t(:)          ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_th(:)         ! atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q(:)          ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_pbot(:)       ! atmospheric pressure (Pa)
    real(r8), pointer :: forc_rho(:)        ! atmospheric density (kg/m**3)
    integer , pointer :: cgridcell(:)       ! column's gridcell
    integer , pointer :: clandunit(:)       ! column's landunit
    integer , pointer :: plandunit(:)       ! pft's landunit
    integer , pointer :: ityplun(:)         ! landunit type

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_c,Hbot
    real(r8) :: zbot_g, tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g
    real(r8) :: zbot_c, tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c
    real(r8) :: egcm_c, rhos_c
    real(r8) :: dum1,   dum2

!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit-level)

    ityplun            => lun%itype

    ! Assign local pointers to derived type members (column-level)

    snl                => cps%snl
    h2osno             => cws%h2osno
    h2osno_old         => cws%h2osno_old
    do_capsnow         => cps%do_capsnow
    frac_iceold        => cps%frac_iceold
    h2osoi_ice         => cws%h2osoi_ice
    h2osoi_liq         => cws%h2osoi_liq
    frac_veg_nosno_alb => pps%frac_veg_nosno_alb
    frac_veg_nosno     => pps%frac_veg_nosno
    qflx_glcice        => cwf%qflx_glcice
    eflx_bot           => cef%eflx_bot
    glc_topo           => cps%glc_topo
    forc_t             => ces%forc_t
    forc_th            => ces%forc_th
    forc_q             => cws%forc_q
    forc_pbot          => cps%forc_pbot
    forc_rho           => cps%forc_rho
    clandunit          => col%landunit
    cgridcell          => col%gridcell

    ! Assign local pointers to derived type members (pft-level)

    pwtgcell           => pft%wtgcell
    plandunit          => pft%landunit

    do c = lbc, ubc

      l = clandunit(c)
      g = cgridcell(c)

      ! Initialize column forcing

      forc_t(c)    = clm_a2l%forc_t(g)
      forc_th(c)   = clm_a2l%forc_th(g)
      forc_q(c)    = clm_a2l%forc_q(g)
      forc_pbot(c) = clm_a2l%forc_pbot(g)
      forc_rho(c)  = clm_a2l%forc_rho(g)

      ! Save snow mass at previous time step
      h2osno_old(c) = h2osno(c)

      ! Decide whether to cap snow
      if (h2osno(c) > h2osno_max) then
         do_capsnow(c) = .true.
      else
         do_capsnow(c) = .false.
      end if
      eflx_bot(c)    = 0._r8
      
      ! Initialize qflx_glcice, but only over ice_mec landunits (elsewhere, it is spval)
      if (ityplun(l) == istice_mec) then
         qflx_glcice(c) = 0._r8
      end if

    end do

    ! Initialize fraction of vegetation not covered by snow (pft-level)

    do p = lbp,ubp
       l = plandunit(p)
       ! Note: Some glacier_mec points may have zero weight
       if (pwtgcell(p)>0._r8 .or. ityplun(l) == istice_mec) then
          frac_veg_nosno(p) = frac_veg_nosno_alb(p)
       else
          frac_veg_nosno(p) = 0._r8
       end if
    end do

    ! Initialize set of previous time-step variables
    ! Ice fraction of snow at previous time step

    do j = -nlevsno+1,0
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if (j >= snl(c) + 1) then
            frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
         end if
      end do
    end do

   ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
   ! For glacier_mec columns the downscaling is based on surface elevation.
   ! For other columns the downscaling is a simple copy.

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       l = clandunit(c)
       g = cgridcell(c)

       if (ityplun(l) == istice_mec) then   ! downscale to elevation classes

          ! This is a simple downscaling procedure taken from subroutine clm_mapa2l.
          ! Note that forc_hgt, forc_u, and forc_v are not downscaled.

          hsurf_g      = ldomain%topo(g)          ! gridcell sfc elevation
          hsurf_c      = glc_topo(c)              ! column sfc elevation

          tbot_g       = clm_a2l%forc_t(g)        ! atm sfc temp
          thbot_g      = clm_a2l%forc_th(g)       ! atm sfc pot temp
          qbot_g       = clm_a2l%forc_q(g)        ! atm sfc spec humid
          pbot_g       = clm_a2l%forc_pbot(g)     ! atm sfc pressure
          zbot_g       = clm_a2l%forc_hgt(g)      ! atm ref height

          zbot_c = zbot_g
          tbot_c = tbot_g-lapse_glcmec*(hsurf_c-hsurf_g)   ! sfc temp for column

          Hbot   = rair*0.5_r8*(tbot_g+tbot_c)/grav        ! scale ht at avg temp
          pbot_c = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)     ! column sfc press
          thbot_c= tbot_c*exp((zbot_c/Hbot)*(rair/cpair))  ! pot temp calc

          call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
          call Qsat(tbot_c,pbot_c,es_c,dum1,qs_c,dum2)

          qbot_c = qbot_g*(qs_c/qs_g)
          egcm_c = qbot_c*pbot_c/(0.622+0.378*qbot_c)
          rhos_c = (pbot_c-0.378*egcm_c) / (rair*tbot_c)

          forc_t(c)    = tbot_c
          forc_th(c)   = thbot_c
          forc_q(c)    = qbot_c
          forc_pbot(c) = pbot_c
          forc_rho(c)  = rhos_c

       endif

    enddo    ! num_nolakec

  end subroutine clm_driverInit

end module clm_driverInitMod

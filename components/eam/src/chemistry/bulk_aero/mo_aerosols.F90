module mo_aerosols
  !-----------------------------------------------------------------
  !
  ! this module computes the production of ammonium nitrate
  ! using the formulation by Seinfeld and Pandis (p531, 1998)
  ! with the simplification of activity coefficients and
  ! aerosol molality using the parameterizations
  ! from Metzger et al. (JGR, ACH-16, 107(D16), 2002) 
  !
  ! written by Jean-Francois Lamarque (April 2004)
  ! adapted for CAM (May 2004)
  !
  !-----------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pver
  use cam_logfile,  only: iulog

  private
  public :: aerosols_inti,aerosols_formation
  public :: has_aerosols

  save

  integer, target  :: spc_ndx(5)
  integer, pointer :: nh3_ndx, nh4no3_ndx, nh4_ndx
  integer, pointer :: so4_ndx, hno3_ndx
  integer  :: xhno3_ndx, xnh4no3_ndx
  integer  :: nu_i(2)
  real(r8)     :: zeta_inv
  real(r8)     :: z_i(2)
  logical  :: has_aerosols = .true.

contains

  subroutine aerosols_inti()

    use mo_chem_utls, only : get_spc_ndx
    use cam_history,  only : addfld
    use spmd_utils,   only : masterproc

    implicit none

    !-----------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------
    integer :: m

    nh3_ndx    => spc_ndx(1)
    nh4no3_ndx => spc_ndx(2)
    so4_ndx    => spc_ndx(3)
    hno3_ndx   => spc_ndx(4)
    nh4_ndx    => spc_ndx(5)

    !-----------------------------------------------------------------
    ! 	... set species index
    !-----------------------------------------------------------------
    nh3_ndx    = get_spc_ndx( 'NH3'    )
    nh4no3_ndx = get_spc_ndx( 'NH4NO3' )
    so4_ndx    = get_spc_ndx( 'SO4'    )
    hno3_ndx   = get_spc_ndx( 'HNO3'   )
    nh4_ndx    = get_spc_ndx( 'NH4'    )
    xnh4no3_ndx = get_spc_ndx( 'XNH4NO3' )
    xhno3_ndx   = get_spc_ndx( 'XHNO3'   )

    has_aerosols = all( spc_ndx(:) > 0 )
    if( .not. has_aerosols ) then
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) 'mozart will NOT do nh4no3'
          write(iulog,*) 'following species are missing'
          do m = 1,size(spc_ndx)
             if( spc_ndx(m) < 1 ) then
                write(iulog,*) m
             end if
          end do
          write(iulog,*) '-----------------------------------------'
       endif
       return
    else
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) 'mozart will do nh4no3'
          write(iulog,*) '-----------------------------------------'
       end if
    end if

    !
    ! define parameters 
    !
    ! ammonium nitrate (NH4NO3)
    !
    zeta_inv = 1._r8/4._r8
    nu_i(1)  = 4
    z_i (1)  = 1._r8
    !
    ! ammonium sulfate
    !
    nu_i(2)  = 4
    z_i (2)  = 0.5_r8

    call addfld ('TSO4_VMR',(/ 'lev' /), 'A', 'mol/mol','total sulfate in mo_aerosols')
    call addfld ('THNO3_VMR',(/ 'lev' /), 'A','mol/mol','total nitric acid in mo_aerosols')

    return
  end subroutine aerosols_inti
  
  subroutine aerosols_formation( ncol, lchnk, tfld, rh, qin)

    use ppgrid, only        : pcols, pver
    use chem_mods,  only    : gas_pcnst, adv_mass
    use cam_history, only   : outfld
    
    implicit none
    !
    ! input arguments
    !
    !
    ! input arguments
    !
    integer,  intent(in)    :: ncol        ! number columns in chunk
    integer,  intent(in)    :: lchnk       ! chunk index
    real(r8), intent(in)    :: tfld(:,:)   ! temperature
    real(r8), intent(in)    :: rh(:,:)     ! relative humidity
    real(r8), intent(inout) :: qin(:,:,:)  ! xported species ( vmr )

    !
    ! local variables
    !
    integer :: i,j,k,n
    integer                           :: domain_number   ! concentration domain
    real(r8)                          :: sulfate_state   ! fraction of sulfate neutralized by ammonia
    real(r8), dimension(ncol,pver)    :: tso4, &         ! total sulfate
                                         thno3,&         ! total nitric acid
                                         txhno3          ! total nitric acid ( XHNO3 )
    real(r8)                          :: tnh3, &         ! total ammonia
                                         fnh3, &         ! free ammonia
                                         rhd,  &         ! relative humidity of deliquescence
                                         gamma, &        ! activity coefficient
                                         ssm_nh4no3, &   ! single solute molality for NH4NO3
                                         ta, &           ! total ammonia
                                         tn, &           ! total nitrate
                                         kp, &           ! equilibrium constant
                                         nh4no3          ! ammonium nitrate produced
    real(r8) :: log_t
    real(r8) :: ti
    real(r8) :: xnh4no3

    do k=1,pver
       do i=1,ncol

          !
          ! compute total concentrations
          !
          tnh3      = (qin(i,k,  nh3_ndx)+qin(i,k,nh4_ndx))
          tso4(i,k) = qin(i,k,so4_ndx)
          !
          ! define concentration domain
          !
          if ( tnh3 < tso4(i,k) ) then
             domain_number = 4
             sulfate_state = 1.0_r8
          elseif ( tnh3 < 2._r8*tso4(i,k) ) then
             domain_number = 3
             sulfate_state = 1.5_r8
          else
             domain_number = 2
             sulfate_state = 2.0_r8
          endif
          !
          ! define free ammonia (ammonia available for ammonium nitrate production)
          !
          fnh3 = tnh3 - sulfate_state * tso4(i,k)
          fnh3 = max(0._r8,fnh3)
          !
          ! convert initial concentrations to ppbv
          !
          tso4(i,k)  = tso4(i,k) * 1.e9_r8
          tnh3       = tnh3 * 1.e9_r8
          fnh3       = fnh3 * 1.e9_r8
          thno3(i,k) = (qin(i,k,hno3_ndx)+qin(i,k,nh4no3_ndx)) * 1.e9_r8
          if ( xhno3_ndx > 0 .and. xnh4no3_ndx > 0 ) then
            txhno3(i,k) = (qin(i,k,xhno3_ndx)+qin(i,k,xnh4no3_ndx)) * 1.e9_r8
          endif
          !
          ! compute relative humidity of deliquescence (%) for NH4NO3
          ! (Seinfeld and Pandis, p532)
          !
          ti    = 1._r8/tfld(i,k)
          rhd   = 0.01_r8 * exp( 1.6954_r8 + 723.7_r8*ti )
          log_t = log( tfld(i,k)/298._r8 )
          if ( rh(i,k) < rhd ) then
             !
             ! crystalline ammonium nitrate
             !
             ! compute equilibrium constant
             !
             kp = exp( 84.6_r8 - 24220._r8*ti - 6.1_r8*log_t )
             !
          else
             !
             ! aqueous phase ammonium nitrate
             !
             ! compute activity coefficients (from Menzger et al.)
             !
             n = domain_number
             gamma = (rh(i,k)**n/(1000._r8/n*(1._r8-rh(i,k))+n))**zeta_inv
             !
             ! compute single solute molality for NH4NO3
             !
             ssm_nh4no3 = (1000._r8 * 0.81_r8 * nu_i(1) * (1._r8/rh(i,k)-1._r8)/80._r8)**z_i(1)
             !
             ! compute equilibrium constant
             !
             kp = (gamma*ssm_nh4no3)**2 * exp( 53.19_r8 - 15850.62_r8*ti + 11.51_r8*log_t )

          endif
          !
          ! calculate production of NH4NO3 (in ppbv) using Seinfeld and Pandis (p534, 1998)
          !
          ta = fnh3
          tn = thno3(i,k)
          nh4no3 = 0.5_r8 * (ta + tn - sqrt(max(0._r8,(ta+tn)**2 - 4._r8*(ta*tn-kp))))
          nh4no3 = max(0._r8,nh4no3)
          if ( xhno3_ndx > 0 .and. xnh4no3_ndx > 0 ) then
             tn = txhno3(i,k)
             xnh4no3 = 0.5_r8 * (ta + tn - sqrt(max(0._r8,(ta+tn)**2 - 4._r8*(ta*tn-kp))))
             xnh4no3 = max(0._r8,xnh4no3)
          endif
          !
          ! reset concentrations according to equilibrium calculation
          !
          qin(i,k,nh4no3_ndx)  = nh4no3
          if ( xhno3_ndx > 0 ) then
             qin(i,k,xnh4no3_ndx)  = xnh4no3
          endif
          qin(i,k,nh3_ndx   )  = max(0._r8,(fnh3-nh4no3))
          qin(i,k,nh4_ndx   )  = max(0._r8,(tnh3-(fnh3-nh4no3)))
          qin(i,k,hno3_ndx  )  = max(0._r8,(thno3(i,k)-nh4no3))
          if ( xhno3_ndx > 0 ) then
             qin(i,k,xhno3_ndx  )  = max(0._r8,(txhno3(i,k)-xnh4no3))
          endif
          qin(i,k,so4_ndx   )  = tso4(i,k)
          !
          ! convert from ppbv to vmr
          !
          qin(i,k,nh4no3_ndx)  = qin(i,k,nh4no3_ndx) * 1.e-9_r8
          qin(i,k,nh3_ndx   )  = qin(i,k,nh3_ndx   ) * 1.e-9_r8
          qin(i,k,nh4_ndx   )  = qin(i,k,nh4_ndx   ) * 1.e-9_r8
          qin(i,k,hno3_ndx  )  = qin(i,k,hno3_ndx  ) * 1.e-9_r8
          qin(i,k,so4_ndx   )  = qin(i,k,so4_ndx   ) * 1.e-9_r8
          if ( xhno3_ndx > 0 ) then
             qin(i,k,xnh4no3_ndx)  = qin(i,k,xnh4no3_ndx) * 1.e-9_r8
          endif
          if ( xhno3_ndx > 0 ) then
             qin(i,k,xhno3_ndx  )  = qin(i,k,xhno3_ndx  ) * 1.e-9_r8
          endif

       end do
    end do
    !
    ! outputs
    !
    call outfld ('TSO4_VMR' ,tso4 (:ncol,:), ncol, lchnk )
    call outfld ('THNO3_VMR',thno3(:ncol,:), ncol, lchnk )

    return
  end subroutine aerosols_formation

end module mo_aerosols

!===========================================================================
! Now used (also) for chemUCI usrrxt's because it has access to the photolysis rates
!   Michael Prather 4/7/2020
!
! Originally for LLNL SuperFast chemistry:
! Combine several reactions into one pseudo reaction to correct the
! photolysis rate J(O1D) to incorporate the effect of the other reactions.
!
! Creator: Philip Cameron-Smith
!===========================================================================

module llnl_O1D_to_2OH_adj

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: O1D_to_2OH_adj, O1D_to_2OH_adj_init

  save

!-LLNL-superfast
  integer :: jo1d_ndx

!-chemUCI
  integer :: jo1dU_ndx
  integer :: uci1_ndx, uci2_ndx, uci3_ndx, uci4_ndx, uci5_ndx, uci6_ndx, uci7_ndx, uci8_ndx, uci9_ndx
  integer :: HO2NO2f_ndx, N2O5f_ndx, PANf_ndx, ucih1_ndx, ucih2_ndx, ucih3_ndx


contains
!===========================================================================

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj_init
    use mo_chem_utls, only : get_rxt_ndx
    use cam_logfile,  only : iulog
    use spmd_utils,       only : masterproc

    implicit none

!-LLNL-superfast
    jo1d_ndx  = get_rxt_ndx( 'j2oh' )
!
!-chemUCI
!
    jo1dU_ndx        = get_rxt_ndx( 'jo1dU' )   ! O3 + hv -> O(1D) + O2
    uci1_ndx         = get_rxt_ndx( 'uci1' )    ! O3(O1D) + H2O -> 2 OH
    uci2_ndx         = get_rxt_ndx( 'uci2' )    ! O3(O1D) + H2 -> OH _ HO2
    uci3_ndx         = get_rxt_ndx( 'uci3' )    ! O3(O1D) + CH4 -> OH + CH3OO
    uci4_ndx         = get_rxt_ndx( 'uci4' )    ! HO2 + HO2 + M
    uci5_ndx         = get_rxt_ndx( 'uci5' )    ! HO2 + HO2
    uci6_ndx         = get_rxt_ndx( 'uci6' )    ! OH + HNO3
    uci7_ndx         = get_rxt_ndx( 'uci7' )    ! HO2NO2r thermal decomp
    uci8_ndx         = get_rxt_ndx( 'uci8' )    ! N2O5r thermal
    uci9_ndx         = get_rxt_ndx( 'uci9' )    ! PANr thermal
    HO2NO2f_ndx      = get_rxt_ndx( 'HO2NO2f' ) ! HO2NO2 formation
    N2O5f_ndx        = get_rxt_ndx( 'N2O5f' )   ! N2O5 formation
    PANf_ndx         = get_rxt_ndx( 'PANf' )    ! PAN formation
    ucih1_ndx        = get_rxt_ndx( 'ucih1' )   ! N2O5 + ASAD (het chem)
    ucih2_ndx        = get_rxt_ndx( 'ucih2' )   ! NO3 + ASAD
    ucih3_ndx        = get_rxt_ndx( 'ucih3' )   ! HO2 + ASAD

    if (masterproc) then
       write (iulog,*) 'O1D_to_2OH_adj_init: Found j2oh index in O1D_to_2OH_adj_init of   ', jo1d_ndx
       write (iulog,*) 'O1D_to_2OH_adj_init: O1D_to_2OH_adj is active'

!chemUCI
       write(iulog,*) 'O1D_to_2OH_adj_init: chemUCI diagnostics '
       write(iulog,'(10i5)') uci1_ndx, uci2_ndx, uci3_ndx, uci4_ndx, uci5_ndx, uci6_ndx, uci7_ndx, uci8_ndx, uci9_ndx
       write(iulog,'(10i5)') jo1dU_ndx, ucih1_ndx, ucih2_ndx, ucih3_ndx
    endif

  end subroutine O1D_to_2OH_adj_init

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj( rxt, invariants, m, ncol, temp )

    use chem_mods,    only : nfs, phtcnt, rxntot, nfs !PJC added rxntot, nfs
    use ppgrid,       only : pcols, pver              !PJC added pcols
    use mo_setinv,    only : n2_ndx, o2_ndx, h2o_ndx  !PJC + MJP
    use mo_usrrxt,    only : comp_exp

    implicit none

    !--------------------------------------------------------------------
    ! ... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: invariants(ncol,pver,nfs)   ! density (#/cm3) of species nfs
    real(r8), intent(in) :: m(ncol,pver)                ! air denisty (#/cm3)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)    ! reaction rates, including photolytic
    real(r8), intent(in)    :: temp(pcols,pver)         ! midpoint temperature (K)

    !--------------------------------------------------------------------
    ! ... local variables
    !--------------------------------------------------------------------
    integer  ::  i, k
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  ko(ncol)
    real(r8) ::  kinf(ncol)
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
    real(r8) ::  fc(ncol)
    real(r8) ::  prod_O1D(ncol)
!    real(r8) ::  c_n2o5, c_ho2, c_no3
!    real(r8) ::  gamma_n2o5, gamma_ho2, gamma_no3

    real(r8) :: n2_rate(ncol,pver)
    real(r8) :: o2_rate(ncol,pver)
    real(r8) :: h2o_rate(ncol,pver)

! maybe need these is calling aero_model_surfarea()
    real(r8), pointer :: sfc(:), dm_aer(:)
    integer :: ntot_amode
    real(r8), pointer :: sfc_array(:,:,:), dm_array(:,:,:)
!    real(r8) :: sad_total(pcols,pver)    ! total surface area density (cm2/cm3)
!    ! worry about strato_sad later, as CARMA is not used for now.
!    ! otherwise, it needs to be passed in.
!    real(r8) :: strato_sad(pcols,pver)     ! stratospheric aerosol sad (1/cm)

    real(r8), parameter :: x1 = 2.15e-11_r8
    real(r8), parameter :: x2 = 3.30e-11_r8
    real(r8), parameter :: x3 = 1.63e-10_r8
    real(r8), parameter :: y1 = 110.0_r8
    real(r8), parameter :: y2 =  55.0_r8
    real(r8), parameter :: y3 =  60.0_r8

!
!-chemUCI reset of all user rates, including those with the O(1D) fix
!     remember that the rxt()'s heref have been calculated rxt =k[X&Y] * [X] * [Y]
!

    level_loop : do k = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,k)
       tp(:)             = 300._r8 * tinv(:)
       sqrt_t(:)         = sqrt( temp(:ncol,k) )

       if( jo1dU_ndx > 0 ) then
          prod_O1D(:) = rxt(:,k,jo1dU_ndx)   ! 1/sec, = production of O(1D) when *[O3]
             call comp_exp( exp_fac, 110._r8*tinv, ncol )
          fc(:) =         2.15e-11_r8 * exp_fac(:) * invariants(:,k,n2_ndx)
             call comp_exp( exp_fac, 55._r8*tinv, ncol )
          fc(:) = fc(:) + 3.30e-11_r8 * exp_fac(:) * invariants(:,k,o2_ndx)
             call comp_exp( exp_fac, 60._r8*tinv, ncol )
          fc(:) = fc(:) + 1.63e-10_r8 * exp_fac(:) * invariants(:,k,h2o_ndx)

          if( uci1_ndx > 0 ) then
!  note RATE = prod_O1D*[O3] *[H2O]*rxt(,,uci1)/fc, last 3 terms give fraction of O1D to that reaction
             rxt(:,k,uci1_ndx) = prod_O1D(:) * rxt(:,k,uci1_ndx) / fc(:)
          endif

          if( uci2_ndx > 0 ) then
             rxt(:,k,uci2_ndx) = prod_O1D(:) * rxt(:,k,uci2_ndx) / fc(:)
          endif

          if( uci3_ndx > 0 ) then
             rxt(:,k,uci3_ndx) = prod_O1D(:) * rxt(:,k,uci3_ndx) / fc(:)
          endif
       endif

       if( uci4_ndx > 0 ) then
          if( h2o_ndx > 0 ) then
             call comp_exp( exp_fac, 2200._r8*tinv, ncol )
             fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,k,h2o_ndx) * exp_fac(:)
             rxt(:,k,uci4_ndx) = rxt(:,k,uci4_ndx) * fc(:)
             rxt(:,k,uci5_ndx) = rxt(:,k,uci5_ndx) * fc(:)
          end if
       endif

       if( uci6_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.50e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          kinf(:) = 2.70e-17_r8 * exp_fac(:)
          rxt(:,k,uci6_ndx) = ko(:)*kinf(:)/(ko(:) + kinf(:))
       endif

! thermal decomposition rates for HO2NO2(uci7), N2O5(uci8), and PAN(uci9)
       if( uci7_ndx > 0 ) then
          if ( HO2NO2f_ndx > 0 ) then
             rxt(:,k,uci7_ndx) = rxt(:,k,HO2NO2f_ndx)/rxt(:,k,uci7_ndx)
          endif
       endif

       if( uci8_ndx > 0 ) then
          if ( N2O5f_ndx > 0 ) then
             rxt(:,k,uci8_ndx) = rxt(:,k,N2O5f_ndx)/rxt(:,k,uci8_ndx)
          endif
       endif

       if( uci9_ndx > 0 ) then
          if ( PANf_ndx > 0 ) then
             rxt(:,k,uci9_ndx) = rxt(:,k,PANf_ndx)/rxt(:,k,uci9_ndx)
          endif
       endif

!       if( ucih1_ndx > 0 ) then
!           rxt(:,k,ucih1_ndx) = 0._r8
!       endif
!
!       if( ucih2_ndx > 0 ) then
!           rxt(:,k,ucih2_ndx) = 0._r8
!       endif
!
!       if( ucih3_ndx > 0 ) then
!           rxt(:,k,ucih3_ndx) = 0._r8
!       endif


!  hetrxtrate sums over all the aerosol types: hydrolysis reactions on wetted aerosols
!    NEEDS surface area(sfc) & diameter(dm_aer) for as many aerosols as there are.
!    SUM  rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
!         rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))
!       if( ucih1_ndx > 0 .or. ucih2_ndx > 0 .or. ucih3_ndx > 0 ) then
!          uci_loop : do i = 1,ncol
!             sfc    => sfc_array(i,k,:)
!             dm_aer => dm_array(i,k,:)
!             c_n2o5 = 1.40e3_r8 * sqrt_t(i)         ! mean molecular speed of n2o5
!             c_no3  = 1.85e3_r8 * sqrt_t(i)         ! mean molecular speed of no3
!             c_ho2  = 2.53e3_r8 * sqrt_t(i)         ! mean molecular speed of ho2
!             gamma_n2o5 = rxt(i,k,ucih1_ndx)
!             gamma_no3  = rxt(i,k,ucih2_ndx)
!             gamma_ho2  = rxt(i,k,ucih3_ndx)
!             if( ucih1_ndx > 0 ) then
!                rxt(i,k,ucih1_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
!             end if
!             if( ucih2_ndx > 0 ) then
!                rxt(i,k,ucih2_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3 , gamma_no3  )
!             end if
!             if( ucih3_ndx > 0 ) then
!                rxt(i,k,ucih3_ndx) = hetrxtrate( sfc, dm_aer, dg, c_ho2 , gamma_ho2  )
!             end if
!          enddo uci_loop
!      endif

    end do level_loop
!
!-chemUCI----------------------------------------------------------------

!---revert to superfast fixes here ONLY if not doing chemUCI
   if ( jo1d_ndx > 0 ) then

    n2_rate(:,:)  = x1 * Exp ( y1 / temp(:ncol,:)) * invariants(:,:,n2_ndx)
    o2_rate(:,:)  = x2 * Exp ( y2 / temp(:ncol,:)) * invariants(:,:,o2_ndx)
    h2o_rate(:,:) = x3 * Exp ( y3 / temp(:ncol,:)) * invariants(:,:,h2o_ndx)

    rxt(:,:,jo1d_ndx) = rxt(:,:,jo1d_ndx) *   &
                          (h2o_rate(:,:) / (h2o_rate(:,:) + n2_rate(:,:) + o2_rate(:,:)))

   endif


  end subroutine O1D_to_2OH_adj

end module llnl_O1D_to_2OH_adj

!===========================================================================
! Now used (also) for chemUCI usrrxt's because it has access to the photolysis rates
!   Michael Prather 4/4/2020
!
! Originally for LLNL SuperFast chemistry:
! Combine several reactions into one pseudo reaction to correct the
! photolysis rate J(O1D) to incorporate the effect of the other reactions.
!
! Creator: Philip Cameron-Smith
!===========================================================================

module llnl_O1D_to_2OH_adj

  save
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: O1D_to_2OH_adj, O1D_to_2OH_adj_init

  integer :: jo1d_ndx

!-chemUCI
  integer :: uci1_ndx, uci2_ndx, uci3_ndx, uci4_ndx, uci5_ndx, uci6_ndx, uci7_ndx, uci8_ndx, uci9_ndx
  integer :: jo1dU_ndx, HNO4_ndx, N2O5_ndx, PAN_ndx, ucih1_ndx, ucih2_ndx, ucih3_ndx


contains
!===========================================================================

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj_init
    use mo_chem_utls, only : get_rxt_ndx
    use cam_logfile,  only : iulog
    use spmd_utils,       only : masterproc

    implicit none

    jo1d_ndx  = get_rxt_ndx( 'j2oh' )
!
!-chemUCI
!
    uci1_ndx         = get_rxt_ndx( 'uci1' )  ! O3(O1D) + H2O -> 2 OH
    uci2_ndx         = get_rxt_ndx( 'uci2' )  ! O3(O1D) + H2 -> OH _ HO2
    uci3_ndx         = get_rxt_ndx( 'uci3' )  ! O3(O1D) + CH4 -> OH + CH3OO
    uci4_ndx         = get_rxt_ndx( 'uci4' )  ! HO2 + HO2 + M
    uci5_ndx         = get_rxt_ndx( 'uci5' )  ! HO2 + HO2
    uci6_ndx         = get_rxt_ndx( 'uci6' )  ! OH + HNO3
    uci7_ndx         = get_rxt_ndx( 'uci7' )  ! HNO4r thermal decomp
    uci8_ndx         = get_rxt_ndx( 'uci8' )  ! N2O5r thermal
    uci9_ndx         = get_rxt_ndx( 'uci9' )  ! PANr thermal
    HNO4f_ndx         = get_rxt_ndx( 'HNO4f' )  ! HNO4 formation
    N2O5f_ndx         = get_rxt_ndx( 'N2O5f' )  ! N2O5 formation
    PANf_ndx          = get_rxt_ndx( 'PANf' )   ! PAN formation
    ucih1_ndx        = get_rxt_ndx( 'ucih1' ) ! N2O5 + ASAD (het chem)
    ucih2_ndx        = get_rxt_ndx( 'ucih2' ) ! NO3 + ASAD
    ucih3_ndx        = get_rxt_ndx( 'ucih3' ) ! HO2 + ASAD



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
    use mo_setinv,    only : n2_ndx, o2_ndx, h2o_ndx, o3_ndx  !PJC + MJP

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
    real(r8), parameter :: dg = 0.1_r8            ! mole diffusion =0.1 cm2/s (Dentener, 1993)

    integer  ::  i, k
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  ko(ncol)
    real(r8) ::  kinf(ncol)
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
    real(r8) ::  fc(ncol)
    real(r8) ::  prod_O1D(ncol)
    real(r8) ::  c_n2o5, c_ho2, c_no3
    real(r8) ::  gamma_n2o5, gamma_ho2, gamma_no3

    real(r8) :: n2_rate(ncol,pver)
    real(r8) :: o2_rate(ncol,pver)
    real(r8) :: h2o_rate(ncol,pver)


             sfc    => sfc_array(i,k,:)
             dm_aer => dm_array(i,k,:)
             rxt(i,k,ucih2_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3 , gamma_no3  )

! maybe need these is calling aero_model_surfarea()
    real(r8), pointer :: sfc(:), dm_aer(:)
    integer :: ntot_amode
    real(r8), pointer :: sfc_array(:,:,:), dm_array(:,:,:)


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
   if( uci1_ndx > 0 ) then

!  all this is to get reactive aerosol surface area for std trop het-chem
      if( ucih1_ndx > 0 .or. ucih2_ndx > 0 .or. ucih3_ndx > 0 ) then
         call rad_cnst_get_info(0, nmodes=ntot_amode)
         if (ntot_amode>0) then
            allocate(sfc_array(pcols,pver,ntot_amode), dm_array(pcols,pver,ntot_amode) )
         else
            allocate(sfc_array(pcols,pver,5), dm_array(pcols,pver,5) )
         endif
         sfc_array(:,:,:) = 0._r8
         dm_array(:,:,:) = 0._r8
         sad_total(:,:) = 0._r8
         if( carma_do_hetchem ) then
            sad_total(:ncol,:pver)=strato_sad(:ncol,:pver)
         else
            call aero_model_surfarea( &
               mmr, rm1, relhum, pmid, temp, strato_sad, &
               sulfate, m, ltrop, het1_ndx, pbuf, ncol, sfc_array, dm_array, sad_total )
         endif
      endif


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

! thermal decomposition rates for HNO4(uci7), N2O5(uci8), and PAN(uci9)
       if( uci7_ndx > 0 ) then
          if ( HNO4_ndx > 0 ) then
             rxt(:,k,uci7_ndx) = rxt(:,k,HNO4f_ndx)/rxt(:,k,uci7_ndx)
          else
             rxt(:,k,uci7_ndx) = 0._r8
          endif
       endif

       if( uci8_ndx > 0 ) then
          if ( N2O5_ndx > 0 ) then
             rxt(:,k,uci8_ndx) = rxt(:,k,N2O5f_ndx)/rxt(:,k,uci8_ndx)
          else
             rxt(:,k,uci8_ndx) = 0._r8
          endif
       endif

       if( uci9_ndx > 0 ) then
          if ( PAN_ndx > 0 ) then
             rxt(:,k,uci9_ndx) = rxt(:,k,PANf_ndx)/rxt(:,k,uci9_ndx)
          else
             rxt(:,k,uci9_ndx) = 0._r8
          endif
       endif

!  hetrxtrate sums over all the aerosol types: hydrolysis reactions on wetted aerosols
!    NEEDS surface area(sfc) & diameter(dm_aer) for as many aerosols as there are.
!    SUM  rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
!         rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))

       if( ucih1_ndx > 0 .or. ucih2_ndx > 0 .or. ucih3_ndx > 0 ) then
          uci_loop : do i = 1,ncol
             sfc    => sfc_array(i,k,:)
             dm_aer => dm_array(i,k,:)
             c_n2o5 = 1.40e3_r8 * sqrt_t(i)         ! mean molecular speed of n2o5
             c_no3  = 1.85e3_r8 * sqrt_t(i)         ! mean molecular speed of no3
             c_ho2  = 2.53e3_r8 * sqrt_t(i)         ! mean molecular speed of ho2
             gamma_n2o5 = rxt(i,k,ucih1_ndx)
             gamma_no3  = rxt(i,k,ucih2_ndx)
             gamma_ho2  = rxt(i,k,ucih3_ndx)
             if( ucih1_ndx > 0 ) then
                rxt(i,k,ucih1_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( ucih2_ndx > 0 ) then
                rxt(i,k,ucih2_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3 , gamma_no3  )
             end if
             if( ucih3_ndx > 0 ) then
                rxt(i,k,ucih3_ndx) = hetrxtrate( sfc, dm_aer, dg, c_ho2 , gamma_ho2  )
             end if
          enddo uci_loop
       endif

    end do level_loop
!
!-chemUCI----------------------------------------------------------------

! LLNL ensure that this is not executed if chemUCI because it might destroy the rates set above.
   else if ( jo1d_ndx > 0 ) then

    n2_rate(:,:)  = x1 * Exp ( y1 / temp(:ncol,:)) * invariants(:,:,n2_ndx)
    o2_rate(:,:)  = x2 * Exp ( y2 / temp(:ncol,:)) * invariants(:,:,o2_ndx)
    h2o_rate(:,:) = x3 * Exp ( y3 / temp(:ncol,:)) * invariants(:,:,h2o_ndx)

    rxt(:,:,jo1d_ndx) = rxt(:,:,jo1d_ndx) *   &
                          (h2o_rate(:,:) / (h2o_rate(:,:) + n2_rate(:,:) + o2_rate(:,:)))

   endif


  end subroutine O1D_to_2OH_adj

  !-------------------------------------------------------------------------
  subroutine comp_exp( x, y, n )

      implicit none

      real(r8), intent(out) :: x(:)
      real(r8), intent(in)  :: y(:)
      integer,  intent(in)  :: n

!  #ifdef IBM
!      call vexp( x, y, n )
!  #else
      x(:n) = exp( y(:n) )
!  #endif

  end subroutine comp_exp

    !-------------------------------------------------------------------------
    !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
    !-------------------------------------------------------------------------
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

      real(r8), intent(in) :: sfc(:)
      real(r8), intent(in) :: dm_aer(:)
      real(r8), intent(in) :: dg_gas
      real(r8), intent(in) :: c_gas
      real(r8), intent(in) :: gamma_gas
      real(r8) :: rate

      real(r8),allocatable :: rxt(:)
      integer :: n, i

      n = size(sfc)

      allocate(rxt(n))
      do i=1,n
         rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
      enddo

      rate = sum(rxt)

      deallocate(rxt)

  endfunction hetrxtrate


end module llnl_O1D_to_2OH_adj

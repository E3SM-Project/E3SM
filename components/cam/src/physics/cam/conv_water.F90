
module conv_water

   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       ! 
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !         Andrew Gettelman October 2010  Separate module                !
   !---------------------------------------------------------------------- !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice
  use cam_abortutils,    only: endrun

  use perf_mod
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: conv_water_register, conv_water_4rad, conv_water_init

! pbuf indices

  integer :: icwmrsh_idx, icwmrdp_idx, fice_idx, sh_frac_idx, dp_frac_idx, &
             ast_idx, sh_cldliq1_idx, sh_cldice1_idx, rei_idx

  integer :: ixcldice, ixcldliq

  contains

  !============================================================================ !

  subroutine conv_water_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the fields in the physics buffer.                            !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

  !-----------------------------------------------------------------------

    ! these calls were already done in convect_shallow...so here I add the same fields to the physics buffer with a "1" at the end
! shallow gbm cloud liquid water (kg/kg)
    call pbuf_add_field('SH_CLDLIQ1','physpkg',dtype_r8,(/pcols,pver/),sh_cldliq1_idx)  
! shallow gbm cloud ice water (kg/kg)
    call pbuf_add_field('SH_CLDICE1','physpkg',dtype_r8,(/pcols,pver/),sh_cldice1_idx)  

  end subroutine conv_water_register


  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine conv_water_init()
   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   !   Initializes the pbuf indices required by conv_water
   ! --------------------------------------------------------------------- ! 

   
   use physics_buffer, only : pbuf_get_index
   use cam_history,    only :  addfld

   use constituents,  only: cnst_get_ind

   implicit none

   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)
 
   icwmrsh_idx  = pbuf_get_index('ICWMRSH')
   icwmrdp_idx  = pbuf_get_index('ICWMRDP')
   fice_idx     = pbuf_get_index('FICE')
   sh_frac_idx  = pbuf_get_index('SH_FRAC')
   dp_frac_idx  = pbuf_get_index('DP_FRAC')
   ast_idx      = pbuf_get_index('AST')
   rei_idx      = pbuf_get_index('REI')

   ! Convective cloud water variables.
   call addfld ('ICIMRCU', (/ 'lev' /), 'A', 'kg/kg', 'Convection in-cloud ice mixing ratio '                   )
   call addfld ('ICLMRCU', (/ 'lev' /), 'A', 'kg/kg', 'Convection in-cloud liquid mixing ratio '                )
   call addfld ('ICIMRTOT', (/ 'lev' /), 'A', 'kg/kg', 'Total in-cloud ice mixing ratio '                        )
   call addfld ('ICLMRTOT', (/ 'lev' /), 'A', 'kg/kg', 'Total in-cloud liquid mixing ratio '                     )

   end subroutine conv_water_init

   subroutine conv_water_4rad( state, pbuf,  conv_water_mode, totg_liq, totg_ice )

   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       ! 
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !                                                                       !
   !---------------------------------------------------------------------- !

   
   use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

   use physics_types,   only: physics_state
   use cam_history,     only: outfld
   use phys_control,    only: phys_getopts
   
   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !

   
   type(physics_state), intent(in)    :: state        ! state variables
   type(physics_buffer_desc), pointer :: pbuf(:)


   integer,  intent(in) :: conv_water_mode

   real(r8), intent(out):: totg_ice(pcols,pver)   ! Total GBA in-cloud ice
   real(r8), intent(out):: totg_liq(pcols,pver)   ! Total GBA in-cloud liquid

   ! --------------- !
   ! Local Workspace !
   ! --------------- !

   ! Physics buffer fields

   real(r8), pointer, dimension(:,:) ::  ast      ! Physical liquid+ice stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  rei      ! Ice effective drop size (microns)

   real(r8), pointer, dimension(:,:) ::  dp_icwmr ! Deep conv. cloud water
   real(r8), pointer, dimension(:,:) ::  sh_icwmr ! Shallow conv. cloud water
   real(r8), pointer, dimension(:,:) ::  fice     ! Ice partitioning ratio
   real(r8), pointer, dimension(:,:) ::  sh_cldliq ! shallow convection gbx liq cld mixing ratio for COSP
   real(r8), pointer, dimension(:,:) ::  sh_cldice ! shallow convection gbx ice cld mixing ratio for COSP

   ! Local Variables

   real(r8) :: conv_ice(pcols,pver)               ! Convective contributions to IC cloud ice
   real(r8) :: conv_liq(pcols,pver)               ! Convective contributions to IC cloud liquid
   real(r8) :: tot_ice(pcols,pver)                ! Total IC ice
   real(r8) :: tot_liq(pcols,pver)                ! Total IC liquid

   integer  :: i,k,itim_old                       ! Lon, lev indices buff stuff.
   real(r8) :: cu_icwmr                           ! Convective  water for this grid-box.   
   real(r8) :: ls_icwmr                           ! Large-scale water for this grid-box. 
   real(r8) :: tot_icwmr                          ! Large-scale water for this grid-box.  
   real(r8) :: ls_frac                            ! Large-scale cloud frac for this grid-box. 
   real(r8) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac 
   real(r8) :: kabs, kabsi, kabsl, alpha, dp0, sh0, ic_limit, frac_limit  
   real(r8) :: wrk1         

   integer :: lchnk
   integer :: ncol

   ! --------- !
   ! Parameter !
   ! --------- !

   parameter( kabsl = 0.090361_r8, frac_limit = 0.01_r8, ic_limit = 1.e-12_r8 )
   character(len=16) :: microp_scheme 

   ncol  = state%ncol
   lchnk = state%lchnk

 ! Get microphysics option
   call phys_getopts( microp_scheme_out = microp_scheme )

 ! Get convective in-cloud water and ice/water temperature partitioning.

   call pbuf_get_field(pbuf, icwmrsh_idx, sh_icwmr )
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr )
   call pbuf_get_field(pbuf, fice_idx,    fice )

 ! Get convective in-cloud fraction    

   call pbuf_get_field(pbuf, sh_frac_idx,  sh_frac )
   call pbuf_get_field(pbuf, dp_frac_idx,  dp_frac )
   call pbuf_get_field(pbuf, rei_idx,      rei )

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,  ast,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) ) 

   ! --------------------------------------------------------------- !
   ! Loop through grid-boxes and determine:                          !
   ! 1. Effective mean in-cloud convective ice/liquid (deep+shallow) !
   ! 2. Effective mean in-cloud total ice/liquid (ls+convective)     !
   ! --------------------------------------------------------------- !

   do k = 1, pver
   do i = 1, ncol

      if( sh_frac(i,k) <= frac_limit .or. sh_icwmr(i,k) <= ic_limit ) then
          sh0_frac = 0._r8
      else
          sh0_frac = sh_frac(i,k)
      endif
      if( dp_frac(i,k) <= frac_limit .or. dp_icwmr(i,k) <= ic_limit ) then
          dp0_frac = 0._r8
      else
          dp0_frac = dp_frac(i,k)
      endif
      cu0_frac = sh0_frac + dp0_frac

    ! For the moment calculate the emissivity based upon the ls clouds ice fraction

      wrk1 = min(1._r8,max(0._r8, state%q(i,k,ixcldice)/(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)+1.e-36_r8)))

      if( ( cu0_frac < frac_limit ) .or. ( ( sh_icwmr(i,k) + dp_icwmr(i,k) ) < ic_limit ) ) then

            cu0_frac = 0._r8
            cu_icwmr = 0._r8
         
            ls_frac = ast(i,k)
            if( ls_frac < frac_limit ) then
                ls_frac  = 0._r8
                ls_icwmr = 0._r8
            else
                ls_icwmr = ( state%q(i,k,ixcldliq) + state%q(i,k,ixcldice) )/max(frac_limit,ls_frac) ! Convert to IC value.
            end if

            tot0_frac = ls_frac
            tot_icwmr = ls_icwmr
           
      else

          ! Select radiation constants (effective radii) for emissivity averaging.
            
            if( microp_scheme == 'RK' ) then
               kabsi = 0.005_r8 + 1._r8/rei(i,k)
            else
               kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(i,k)),130._r8)
            endif
            kabs  = kabsl * ( 1._r8 - wrk1 ) + kabsi * wrk1
            alpha = -1.66_r8*kabs*state%pdel(i,k)/gravit*1000.0_r8

          ! Selecting cumulus in-cloud water.            

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted arithmetic average
               cu_icwmr = ( sh0_frac * sh_icwmr(i,k) + dp0_frac*dp_icwmr(i,k))/max(frac_limit,cu0_frac)
            case (2)
               sh0 = exp(alpha*sh_icwmr(i,k))
               dp0 = exp(alpha*dp_icwmr(i,k))               
               cu_icwmr = log((sh0_frac*sh0+dp0_frac*dp0)/max(frac_limit,cu0_frac))
               cu_icwmr = cu_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

          ! Selecting total in-cloud water. 
          ! Attribute large-scale/convective area fraction differently from default.

            ls_frac   = ast(i,k) 
            ls_icwmr  = (state%q(i,k,ixcldliq) + state%q(i,k,ixcldice))/max(frac_limit,ls_frac) ! Convert to IC value.
            tot0_frac = (ls_frac + cu0_frac) 

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted 'arithmetic in emissivity' average
               tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit,tot0_frac)
            case (2)
               tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr)+cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit,tot0_frac))
               tot_icwmr = tot_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

      end if

    ! Repartition convective cloud water into liquid and ice phase.
    ! Currently, this partition is made using the ice fraction of stratus condensate.
    ! In future, we should use ice fraction explicitly computed from the convection scheme.

      conv_ice(i,k) = cu_icwmr * wrk1
      conv_liq(i,k) = cu_icwmr * (1._r8-wrk1)

      tot_ice(i,k)  = tot_icwmr * wrk1
      tot_liq(i,k)  = tot_icwmr * (1._r8-wrk1)

      totg_ice(i,k) = tot0_frac * tot_icwmr * wrk1
      totg_liq(i,k) = tot0_frac * tot_icwmr * (1._r8-wrk1)

   end do
   end do

!add pbuff calls for COSP
   call pbuf_get_field(pbuf, sh_cldliq1_idx, sh_cldliq  )
   call pbuf_get_field(pbuf, sh_cldice1_idx, sh_cldice  )

   sh_cldliq(:ncol,:pver)=sh_icwmr(:ncol,:pver)*(1-fice(:ncol,:pver))*sh_frac(:ncol,:pver)
   sh_cldice(:ncol,:pver)=sh_icwmr(:ncol,:pver)*fice(:ncol,:pver)*sh_frac(:ncol,:pver)

  ! Output convective IC WMRs
   
   call outfld( 'ICLMRCU ', conv_liq  , pcols, lchnk )
   call outfld( 'ICIMRCU ', conv_ice  , pcols, lchnk )
   call outfld( 'ICLMRTOT', tot_liq   , pcols, lchnk )
   call outfld( 'ICIMRTOT', tot_ice   , pcols, lchnk )

  end subroutine conv_water_4rad

end module conv_water

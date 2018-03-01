
subroutine qqflx_fixer (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                        q, rpdel,shflx   ,lhflx   ,qflx    )

!!........................................................................
!! Water conservation fixer 
!! 
!! If QFLX is too negative, the condensation or deposition water vapor at 
!! the surface will take all the available moisture in surface layer. This 
!! will cause problems in the vertical diffusion calculation, where QFLX 
!! is applied. In the original CESM model, QNEG4 is called to correct QFLX 
!! so that it won't take out all the available moisture from the surface layer.
!! 
!! The new fixer, named as qqflx fixer, borrows water vapor from the whole 
!! column above the surface layer proportionally and add moisture into the 
!! surface layer, so that it can compensate the downward (negative) QFLX. 
!!
!! The excess downward (negative) q flux is compared to a theoretical
!! maximum downward q flux.  The theoretical max is based upon the
!! given moisture content of lowest level of the model atmosphere.
!!
!! Author: Kai Zhang (kai.zhang@pnnl.gov) and Phil Rasch 
!!........................................................................

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, pver
   use phys_grid,       only: get_lat_p, get_lon_p, get_wght_all_p 
   use physconst,       only: gravit, latvap
   use constituents,    only: qmin, pcnst
   use cam_logfile,     only: iulog
   use spmd_utils,      only: masterproc
   use cam_abortutils,  only: endrun
   use phys_control,    only: print_fixer_message

   implicit none

   integer, parameter :: i_wv = 1 
   real(r8),parameter :: factor = 1.5_r8 

!! Input arguments
!!........................................................................

   character*8, intent(in) :: subnam         ! name of calling routine

   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns

   real(r8), intent(in) :: ztodt             ! time step
   real(r8), intent(in) :: rpdel(pcols,pver) ! 1./(pint(k+1)-pint(k))

!! Input/Output arguments
!!........................................................................

   real(r8), intent(inout) :: shflx(pcols)        ! surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)        ! surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst)  ! surface water flux (kg/m^2/s)
   real(r8), intent(inout) :: q(pcols,pver,pcnst) ! moisture (kg/kg)

!! Local workspace
!!........................................................................


   real(r8):: excess(pcols)     ! excess downward surface latent heat flux 

   integer :: i, k 
   real(r8) :: wst(pver)        ! local
   real(r8) :: wsp              ! local
   real(r8) :: wpnet            ! local 
   real(r8) :: qori             ! local 
   real(r8) :: qbot             ! local 
   real(r8) :: dqbot            ! local 
   real(r8) :: ratio            ! local 
   real(r8) :: dwpbot           ! local 
   real(r8) :: qflx_org         ! local 

!! begin
!!........................................................................


   !! loop over cols
   !!....................................................................... 

   do i = 1, ncol

      !! check if downward water flux is too large
      !!....................................................................... 

      excess(i) = qflx(i,i_wv) - (qmin(i_wv) - q(i,pver,i_wv))/(ztodt*gravit*rpdel(i,pver))

      if (excess(i) < 0.0_r8 ) then

         qori = q(i,pver,i_wv) 

         !! new q at the bottom to balance the negative qflx 
         !!....................................................................... 

         qbot = qmin(i_wv) - qflx(i,i_wv) * ztodt*gravit*rpdel(i,pver) 

         !! due to process splitting, it could happen that the estimated qbot is 
         !! not sufficiently large to compensate the negative qflx in the vertical 
         !! diffusion scheme, so an adjust factor can be applied. 
         !!....................................................................... 

         qbot = qbot * factor

         !! change of q at the bottom
         !!....................................................................... 

         dqbot = qbot - q(i,pver,i_wv) 

         !! change of q * pdel /g  
         !!....................................................................... 

         dwpbot = dqbot / (gravit*rpdel(i,pver))

         if(dwpbot.lt.0._r8) then 
            call endrun('qflx_fixer: dwpbot < 0 ')  
         end if 

         wsp = 0._r8 

         do k = 1, pver-1 
            wst(k) = q(i,k,i_wv) / (gravit*rpdel(i,k))
            wsp    = wsp + wst(k) 
         end do

         wpnet = wsp - dwpbot

         if(wpnet.gt.0._r8 .and. wsp.gt.0._r8) then 
          
            !! if there is sufficient water vapor in the column, scale q at each level
            !! to compensate the potential water vapor sink in the surface layer. 
            !!....................................................................... 
 
            ratio = wpnet / wsp 

            q(i,1:pver-1,i_wv) = q(i,1:pver-1,i_wv) * ratio  
             
            q(i,pver,i_wv) = qbot 

            write(iulog,*) ' ### qflx_fixer ### ', &
                           ' chunk', lchnk, &
                           ' col', i, &
                           ' lat = ', get_lat_p(lchnk,i), &
                           ' lon = ', get_lon_p(lchnk,i), & 
                           ' qflx = ', qflx(i,i_wv), &
                           ' qori = ', qori, &
                           ' qnew = ', qbot, &
                           ' column wp = ', wsp, & 
                           ' dwp by qflx = ', dwpbot, & 
                           ' q scaled to ', 100.*ratio, "%"

         else 

            write(iulog,*) ' ### qqflx_fixer ### ', &
                           ' chunk', lchnk, &
                           ' col', i, &
                           ' lat = ', get_lat_p(lchnk,i), &
                           ' lon = ', get_lon_p(lchnk,i), & 
                           ' original qflx = ', qflx_org, &  
                           ' qsurf = ', qori, & 
                           ' column wp = ', wsp, & 
                           ' dwp by qflx = ', dwpbot
            write(iulog,*) ' column does not have enough water to compensate water vapor sink due to negative qflx !!! '

            call endrun('qflx_fixer: spurious negative qflx ') 


!!            !! if the column above the surface layer doens't have enough water
!!            !! use qneg4 to fix qflx and printout the conservation error message
!!            !!....................................................................... 
!!          
!!            qflx_org   = qflx (i,i_wv) 
!! 
!!            qflx (i,i_wv) = qflx (i,i_wv) - excess(i)
!!            lhflx(i)   = lhflx(i) - excess(i)*latvap
!!            shflx(i)   = shflx(i) + excess(i)*latvap
!!
!!            write(iulog,*) ' ### qqflx_fixer ### ', &
!!                           ' chunk', lchnk, &
!!                           ' col', i, &
!!                           ' lat = ', get_lat_p(lchnk,i), &
!!                           ' lon = ', get_lon_p(lchnk,i), & 
!!                           ' original qflx = ', qflx_org, &  
!!                           ' modified qflx = ', qflx(i,i_wv), & 
!!                           ' qsurf = ', qori, & 
!!                           ' column wp = ', wsp, & 
!!                           ' dwp by qflx = ', dwpbot

         end if !! if(wpnet.gt.0._r8 .and. wsp.gt.0._r8) 
         
     end if !! if (excess(i) < 0.0_r8 ) 
 
   end do !! i 

!! end
!!........................................................................

   return
end subroutine qqflx_fixer

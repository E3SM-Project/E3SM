module cldfrc2m

! cloud fraction calculations

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols
use physconst,        only: rair
use wv_saturation,    only: qsat_water, svp_water, svp_ice
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public ::           &
   cldfrc2m_readnl, &
   cldfrc2m_init,   &
   astG_PDF_single, &
   astG_PDF,        &
   astG_RHU_single, &
   astG_RHU,        &
   aist_single,     &
   aist_vector,     &
   CAMstfrac,       &
   rhmini_const,    &
   rhmaxi_const

! Namelist variables
real(r8) :: cldfrc2m_rhmini            ! Minimum rh for ice cloud fraction > 0.
real(r8) :: cldfrc2m_rhmaxi

! -------------------------- !
! Parameters for Ice Stratus !
! -------------------------- !
real(r8),  protected :: rhmini_const                 ! Minimum rh for ice cloud fraction > 0.
real(r8),  protected :: rhmaxi_const

real(r8),  parameter :: qist_min     = 1.e-7_r8      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
real(r8),  parameter :: qist_max     = 5.e-3_r8      ! Maximum in-stratus ice IWC constraint [ kg/kg ]

! ----------------------------- !
! Parameters for Liquid Stratus !
! ----------------------------- !

logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
                                                  ! use Slingo (triangular PDF-based) liquid stratus fraction
logical,  parameter  :: freeze_dry   = .false.    ! If .true., use 'freeze dry' in liquid stratus fraction formula
real(r8)             :: rhminl_const              ! Critical RH for low-level  liquid stratus clouds
real(r8)             :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
real(r8)             :: rhminh_const              ! Critical RH for high-level liquid stratus clouds
real(r8)             :: premit                    ! Top    height for mid-level liquid stratus fraction
real(r8)             :: premib                    ! Bottom height for mid-level liquid stratus fraction
integer              :: iceopt                    ! option for ice cloud closure 
                                                  ! 1=wang & sassen 2=schiller (iciwc)  
                                                  ! 3=wood & field, 4=Wilson (based on smith)
                                                  ! 5=modified slingo (ssat & empyt cloud)        
real(r8)             :: icecrit                   ! Critical RH for ice clouds in Wilson & Ballard closure
                                                  ! ( smaller = more ice clouds )

!================================================================================================
contains
!================================================================================================

subroutine cldfrc2m_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldfrc2m_readnl'

   namelist /cldfrc2m_nl/ cldfrc2m_rhmini, cldfrc2m_rhmaxi
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldfrc2m_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldfrc2m_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      rhmini_const = cldfrc2m_rhmini
      rhmaxi_const = cldfrc2m_rhmaxi

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(rhmini_const,      1, mpir8,  0, mpicom)
   call mpibcast(rhmaxi_const,      1, mpir8,  0, mpicom)
#endif

end subroutine cldfrc2m_readnl

!================================================================================================

subroutine cldfrc2m_init()

   use cloud_fraction, only: cldfrc_getparams

   call cldfrc_getparams(rhminl_out=rhminl_const, rhminl_adj_land_out=rhminl_adj_land_const,  &
                         rhminh_out=rhminh_const, premit_out=premit, premib_out=premib, &
                         iceopt_out=iceopt, icecrit_out=icecrit)

   if( masterproc ) then
      write(iulog,*) 'cldfrc2m parameters:'
      write(iulog,*) '  rhminl          = ', rhminl_const
      write(iulog,*) '  rhminl_adj_land = ', rhminl_adj_land_const
      write(iulog,*) '  rhminh          = ', rhminh_const
      write(iulog,*) '  premit          = ', premit
      write(iulog,*) '  premib          = ', premib
      write(iulog,*) '  iceopt          = ', iceopt
      write(iulog,*) '  icecrit         = ', icecrit
      write(iulog,*) '  rhmini          = ', rhmini_const
      write(iulog,*) '  rhmaxi          = ', rhmaxi_const
   end if

end subroutine cldfrc2m_init

!================================================================================================


subroutine astG_PDF_single(U, p, qv, landfrac, snowh, a, Ga, orhmin, &
                           rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: U                     ! Relative humidity
   real(r8), intent(in)  :: p                     ! Pressure [Pa]
   real(r8), intent(in)  :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac              ! Land fraction
   real(r8), intent(in)  :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a                     ! Stratus fraction
   real(r8), intent(out) :: Ga                    ! dU/da
   real(r8), optional, intent(out) :: orhmin      ! Critical RH

   real(r8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   real(r8) :: rhminl
   real(r8) :: rhminl_adj_land
   real(r8) :: rhminh

   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8

   rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

end subroutine astG_PDF_single

!================================================================================================

subroutine astG_PDF(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                    rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: U_in(pcols)           ! Relative humidity
   real(r8), intent(in)  :: p_in(pcols)           ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(pcols)          ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols)    ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)       ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out(pcols)          ! Stratus fraction
   real(r8), intent(out) :: Ga_out(pcols)         ! dU/da
   integer,  intent(in)  :: ncol

   real(r8), optional, intent(in)  :: rhminl_in(pcols)                ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in(pcols)       ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in(pcols)                ! Critical relative humidity for high-level liquid stratus

   real(r8)              :: rhminl                ! Critical relative humidity for low-level  liquid stratus
   real(r8)              :: rhminl_adj_land       ! Adjustment drop of rhminl over the land
   real(r8)              :: rhminh                ! Critical relative humidity for high-level liquid stratus

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8

   rhminl          = rhminl_const
   rhminl_adj_land = rhminl_adj_land_const
   rhminh          = rhminh_const

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out(:)  = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if (present(rhminl_in))          rhminl          = rhminl_in(i)      
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in(i)
   if (present(rhminh_in))          rhminh          = rhminh_in(i)

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

end subroutine astG_PDF
!================================================================================================

subroutine astG_RHU_single(U, p, qv, landfrac, snowh, a, Ga, orhmin, &
                              rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: U               ! Relative humidity
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a               ! Stratus fraction
   real(r8), intent(out) :: Ga              ! dU/da
   real(r8), optional, intent(out) :: orhmin ! Critical RH

   real(r8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhdif                                 ! Factor for stratus fraction
   real(r8) rhwght

   real(r8) :: rhminl
   real(r8) :: rhminl_adj_land
   real(r8) :: rhminh

   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e10_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

end subroutine astG_RHU_single

!================================================================================================

subroutine astG_RHU(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                    rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: U_in(pcols)           ! Relative humidity
   real(r8), intent(in)  :: p_in(pcols)           ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(pcols)          ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols)    ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)       ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out(pcols)          ! Stratus fraction
   real(r8), intent(out) :: Ga_out(pcols)         ! dU/da
   integer,  intent(in)  :: ncol

   real(r8), optional, intent(in)  :: rhminl_in(pcols)          ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in(pcols) ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in(pcols)          ! Critical relative humidity for high-level liquid stratus

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: rhminl                ! Critical relative humidity for low-level  liquid stratus
   real(r8)              :: rhminl_adj_land       ! Adjustment drop of rhminl over the land
   real(r8)              :: rhminh                ! Critical relative humidity for high-level liquid stratus

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer  i
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhdif                                 ! Factor for stratus fraction
   real(r8) rhwght

   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   rhminl          = rhminl_const
   rhminl_adj_land = rhminl_adj_land_const
   rhminh          = rhminh_const

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out(:) = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if (present(rhminl_in))          rhminl          = rhminl_in(i)      
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in(i)
   if (present(rhminh_in))          rhminh          = rhminh_in(i)

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e10_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

end subroutine astG_RHU

!================================================================================================

subroutine aist_single(qv, T, p, qi, landfrac, snowh, aist, &
                       rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in)

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T               ! Temperature
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8), optional, intent(in)  :: rhmaxi_in
   real(r8), optional, intent(in)  :: rhmini_in          ! Critical relative humidity for               ice stratus
   real(r8), optional, intent(in)  :: rhminl_in          ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in          ! Critical relative humidity for high-level liquid stratus

   ! Local variables
   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) es, qs

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
   real(r8) rhdif                           ! working variable for slingo scheme

   real(r8) :: rhmaxi
   real(r8) :: rhmini
   real(r8) :: rhminl
   real(r8) :: rhminl_adj_land
   real(r8) :: rhminh

   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8

   rhmaxi = rhmaxi_const
   if (present(rhmaxi_in)) rhmaxi = rhmaxi_in
   rhmini = rhmini_const
   if (present(rhmini_in)) rhmini = rhmini_in
   rhminl = rhminl_const
   if (present(rhminl_in)) rhminl = rhminl_in
   rhminl_adj_land = rhminl_adj_land_const
   if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in
   rhminh = rhminh_const
   if (present(rhminh_in)) rhminh = rhminh_in

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     call qsat_water(T, p, es, qs)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
        ! set rh ice cloud fraction
        rhi= (qv+qi)/qs * (esl/esi)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        aist = min(1.0_r8, max(rhdif,0._r8)**2)


        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists

        if (qi.lt.minice) then
           aist=0._r8
        else
           aist=max(mincld,aist)
        endif

        ! enforce limits on icimr
        if (qi.ge.minice) then
           icimr=qi/aist

           !minimum
           if (icimr.lt.qist_min) then
              aist = max(0._r8,min(1._r8,qi/qist_min))
           endif
           !maximum
           if (icimr.gt.qist_max) then
              aist = max(0._r8,min(1._r8,qi/qist_max))
           endif

        endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

end subroutine aist_single

!================================================================================================

subroutine aist_vector(qv_in, T_in, p_in, qi_in, ni_in, &
#ifdef CLUBB_SGS
                       do_icesuper, &
#endif
                       landfrac_in, snowh_in, aist_out, ncol, &
                       rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   real(r8), intent(in)  :: qv_in(pcols)       ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T_in(pcols)        ! Temperature
   real(r8), intent(in)  :: p_in(pcols)        ! Pressure [Pa]
   real(r8), intent(in)  :: qi_in(pcols)       ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: ni_in(pcols)       ! Grid-mean ice water number concentration [#/kg]
   real(r8), intent(in)  :: landfrac_in(pcols) ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)    ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: aist_out(pcols)    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )
   integer,  intent(in)  :: ncol 

   real(r8), optional, intent(in)  :: rhmaxi_in
   real(r8), optional, intent(in)  :: rhmini_in(pcols)          ! Critical relative humidity for               ice stratus
   real(r8), optional, intent(in)  :: rhminl_in(pcols)          ! Critical relative humidity for low-level  liquid stratus
   real(r8), optional, intent(in)  :: rhminl_adj_land_in(pcols) ! Adjustment drop of rhminl over the land
   real(r8), optional, intent(in)  :: rhminh_in(pcols)          ! Critical relative humidity for high-level liquid stratus

#ifdef CLUBB_SGS
   logical, intent(in)   :: do_icesuper
#endif

   ! Local variables

   real(r8) qv                              ! Grid-mean water vapor[kg/kg]
   real(r8) T                               ! Temperature
   real(r8) p                               ! Pressure [Pa]
   real(r8) qi                              ! Grid-mean ice water content [kg/kg]
   real(r8) ni
   real(r8) landfrac                        ! Land fraction
   real(r8) snowh                           ! Snow depth (liquid water equivalent)

   real(r8) rhmaxi                          ! Critical relative humidity for               ice stratus
   real(r8) rhmini                          ! Critical relative humidity for               ice stratus
   real(r8) rhminl                          ! Critical relative humidity for low-level  liquid stratus
   real(r8) rhminl_adj_land                 ! Adjustment drop of rhminl over the land
   real(r8) rhminh                          ! Critical relative humidity for high-level liquid stratus

   real(r8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs,ah,bh,ch         ! Fit parameters
   real(r8) nil
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) qs
   real(r8) esat_in(pcols)
   real(r8) qsat_in(pcols)

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
   real(r8) rhdif                           ! working variable for slingo scheme

   integer i

#ifdef CLUBB_SGS
   real(r8) aist2
#endif

   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8

     rhmaxi          = rhmaxi_const
     if (present(rhmaxi_in))          rhmaxi          = rhmaxi_in

     rhmini          = rhmini_const
     rhminl          = rhminl_const
     rhminl_adj_land = rhminl_adj_land_const
     rhminh          = rhminh_const

#ifdef CLUBB_SGS
     if (do_icesuper) then
       iceopt = 7
     endif
#endif

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out(:) = 0._r8
     esat_in(:)  = 0._r8
     qsat_in(:)  = 0._r8

     call qsat_water(T_in(1:ncol), p_in(1:ncol), &
          esat_in(1:ncol), qsat_in(1:ncol))
     
     do i = 1, ncol

     landfrac = landfrac_in(i)     
     snowh = snowh_in(i)   
     T = T_in(i)
     qv = qv_in(i)
     p = p_in(i)
     qi = qi_in(i)
     ni = ni_in(i)
     qs = qsat_in(i)
     esl = svp_water(T)
     esi = svp_ice(T)

     if (present(rhmini_in))          rhmini          = rhmini_in(i)      
     if (present(rhminl_in))          rhminl          = rhminl_in(i)      
     if (present(rhminl_adj_land_in)) rhminl_adj_land = rhminl_adj_land_in(i)
     if (present(rhminh_in))          rhminh          = rhminh_in(i)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land(i) .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
        ! set rh ice cloud fraction
        rhi= (qv+qi)/qs * (esl/esi)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        aist = min(1.0_r8, max(rhdif,0._r8)**2)

#ifdef CLUBB_SGS
     elseif (iceopt.eq.6 .or. iceopt.eq.7) then
#else
     elseif (iceopt.eq.6) then
#endif
        !----- ICE CLOUD OPTION 6: fit based on T and Number (Gettelman: based on Heymsfield obs)
        ! Use observations from Heymsfield et al 2012 of IWC and Ni v. Temp
        ! Multivariate fit follows form of Boudala 2002: ICIWC = a * exp(b*T) * N^c
        ! a=6.73e-8, b=0.05, c=0.349
        ! N is #/L, so need to convert Ni_L=N*rhoa/1000.

#ifdef CLUBB_SGS
     if (iceopt.eq.7) then
        ah= 6.73834e-12_r8
     else
#endif
        ah= 6.73834e-08_r8
#ifdef CLUBB_SGS
     endif
#endif
        bh= 0.0533110_r8
        ch= 0.3493813_r8
        rho=p/(rair*T)
        nil=ni*rho/1000._r8
        icicval = ah * exp(bh*T) * nil**ch
        !result is in g m-3, convert to kg H2O / kg air (icimr...)
        icicval = icicval / rho / 1000._r8
        aist =  max(0._r8,min(qi/icicval,1._r8))
        aist =  min(aist,1._r8)

#ifdef CLUBB_SGS
        rhi= (qv+qi)/qs * (esl/esi)
        rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
        aist2 = min(1.0_r8, max(rhdif,0._r8)**2)

        !aist = max(aist,aist2)
        aist = min(aist,aist2)
#endif
     endif     

     if (iceopt.eq.5 .or. iceopt.eq.6) then


        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists

        if (qi.lt.minice) then
           aist=0._r8
        else
           aist=max(mincld,aist)
        endif

        ! enforce limits on icimr
        if (qi.ge.minice) then
           icimr=qi/aist

           !minimum
           if (icimr.lt.qist_min) then
              aist = max(0._r8,min(1._r8,qi/qist_min))
           endif
           !maximum
           if (icimr.gt.qist_max) then
              aist = max(0._r8,min(1._r8,qi/qist_max))
           endif

        endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

     aist_out(i) = aist

     enddo

end subroutine aist_vector

!================================================================================================

end module cldfrc2m

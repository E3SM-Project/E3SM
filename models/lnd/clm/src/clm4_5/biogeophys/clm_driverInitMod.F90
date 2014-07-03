module clm_driverInitMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialization of clm driver variables needed from previous timestep
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_driverInit
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_driverInit(bounds, num_nolakec, filter_nolakec)
    !
    ! !DESCRIPTION:
    ! Initialization of clm driver variables needed from previous timestep
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevsno
    use clm_varcon   , only : h2osno_max, istice_mec
    use decompMod    , only : bounds_type
    use shr_infnan_mod,only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)   ! column filter for non-lake points
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, p, f, j         ! indices

    !-----------------------------------------------------------------------

   associate(& 
   snl                => cps%snl                 , & ! Input:  [integer (:)]  number of snow layers                    
   h2osno             => cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   h2osoi_ice         => cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                      
   h2osoi_liq         => cws%h2osoi_liq          , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                  
   frac_veg_nosno_alb => pps%frac_veg_nosno_alb  , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1) [-]
   frac_veg_nosno     => pps%frac_veg_nosno      , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1 now) [-] (pft-level)
   h2osno_old         => cws%h2osno_old          , & ! Output: [real(r8) (:)]  snow water (mm H2O) at previous time step
   do_capsnow         => cps%do_capsnow          , & ! Output: [logical (:)]  true => do snow capping                  
   frac_iceold        => cps%frac_iceold         , & ! Output: [real(r8) (:,:)]  fraction of ice relative to the tot water
   qflx_glcice        => cwf%qflx_glcice         , & ! Output: [real(r8) (:)]  flux of new glacier ice (mm H2O/s) [+ = ice grows]
   eflx_bot           => cef%eflx_bot              & ! Output: [real(r8) (:)]  heat flux from beneath soil/ice column (W/m**2)
   )

    do c = bounds%begc,bounds%endc

      l = col%landunit(c)

      ! Save snow mass at previous time step
      h2osno_old(c) = h2osno(c)

      ! Decide whether to cap snow
      if (h2osno(c) > h2osno_max) then
         do_capsnow(c) = .true.
      else
         do_capsnow(c) = .false.
      end if
      eflx_bot(c)    = 0._r8
      
      ! Initialize qflx_glcice everywhere, to zero.
      qflx_glcice(c) = 0._r8     
      
    end do

    ! Initialize fraction of vegetation not covered by snow (pft-level)

    do p = bounds%begp,bounds%endp
       if (pft%active(p)) then
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

    end associate 
   end subroutine clm_driverInit


end module clm_driverInitMod

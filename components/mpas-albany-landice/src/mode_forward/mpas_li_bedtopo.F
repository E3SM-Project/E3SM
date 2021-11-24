! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_bedtopo
!
!> \MPAS land-ice bedtopo driver
!> \author Matt Hoffman
!> \date   20 June 2019
!> \details
!>  This module contains the routines for bed topography for solid earth changes
!>
!
!-----------------------------------------------------------------------

module li_bedtopo

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use li_mask
   use li_setup

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------
   public :: li_bedtopo_init, &
             li_bedtopo_finalize, &
             li_bedtopo_block_init, &
             li_bedtopo_solve

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------



!***********************************************************************

contains

!***********************************************************************
!
!  routine li_bedtopo_init
!
!> \brief   Initializes bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine initializes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_init(domain, err)
#ifdef USE_SEALEVELMODEL
		use spharmt 
      use planets_mod
		use io_mod
		!use user_specs_mod
#endif
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------


      ! No init is needed.
      err = 0

#ifdef USE_SEALEVELMODEL
      call mpas_log_write("<<< 1.Successfully called slmodel_init()! >>>") ! Replace with actual call
		call slmodel_init()
      call mpas_log_write("<<< 2.Successfully called slmodel_init()! >>>") ! Replace with actual call
		
#endif

   !--------------------------------------------------------------------

   end subroutine li_bedtopo_init



!***********************************************************************
!
!  routine li_bedtopo_block_init
!
!> \brief   Initializes blocks for bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine initializes each block of the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_block_init(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: &
         block          !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! No block init needed.
      err = 0

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_block_init



!***********************************************************************
!
!  subroutine li_bedtopo_solve
!
!> \brief   Updates bed topography
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine updates the bed topography.  Currently the only option
!>  is a data field passed in as input.
!
!-----------------------------------------------------------------------
subroutine li_bedtopo_solve(domain, err)

      use li_mask
      use li_advection

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block
      character (len=StrKIND), pointer :: config_uplift_method
      type (mpas_pool_type), pointer :: meshPool          !< mesh information
      type (mpas_pool_type), pointer :: geometryPool      !< geometry information
      type (mpas_pool_type), pointer :: velocityPool      !< velocity information

      real (kind=RKIND), dimension(:), pointer :: bedTopography, upliftRate
      real (kind=RKIND), pointer :: deltat
      integer :: err_tmp

		!HH: variables for the sea-level solver module --------------------
		integer :: i, itersl_sh, iter_sh, dtime_sh
		real :: starttime_sh                ! start time of the simulation 
		integer :: iargc, nargs             ! Arguments read in from a bash script                           
		character(16) :: carg(20)           ! Arguments from a bash script

		! Reading in arguments from a bash script
	!	nargs = iargc()
	!	do i=1,nargs
	!	   call getarg(i, carg(i))
	!	enddo
	!	read (carg(1),*) itersl_sh
	!	read (carg(2),*) iter_sh      ! the coupling time step we are on (in years)
	!	read (carg(3),*) dtime_sh     ! coupling time (in years)
	!	read (carg(4),*) starttime_sh ! start time of the simulation (in years)
		!HH: variables for the sea-level solver module --------------------
		
      err = 0
      err_tmp = 0

      ! Set needed variables and pointers
      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='none') then
         ! do nothing
      elseif (trim(config_uplift_method)=='data') then

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

            call mpas_pool_get_array(meshPool, 'deltat', deltat)
            call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
            call mpas_pool_get_array(geometryPool, 'upliftRate', upliftRate)

            bedTopography(:) = bedTopography(:) + upliftRate(:) * deltat

            call li_update_geometry(geometryPool)
            call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)

            block => block % next
         end do

      else
         call mpas_log_write("Unknown option selected for 'config_uplift_method'", MPAS_LOG_ERR)
      endif


#ifdef USE_SEALEVELMODEL
      ! add calls to code in the sea level model inside ifdefs only.
      call mpas_log_write("<<< Compiled successfully with Sea Level Model! >>>") ! Replace with actual call
		
		!HH
		!-------------------------------------------------------------------------------------------	
     ! elseif (trim(config_uplift_method)=='slmodel') then !HH
			
		!	call sl_solver_checkpoint(itersl_sh, dtime_sh)
		!	call set_planet 
		!	call sl_timewindow(iter_sh)
			
	    !  elseif (iter_sh .gt. 0) then 
	!		   call sl_solver(itersl_sh, iter_sh, dtime_sh, starttime_sh) !dtime_sh should come from nml.
	!		endif
		!-------------------------------------------------------------------------------------------	
#endif


      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_bedtopo_solve.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------

   end subroutine li_bedtopo_solve




!***********************************************************************
!
!  routine li_bedtopo_finalize
!
!> \brief   finalizes bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine finalizes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_finalize(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      err = 0

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_finalize



   ! private subroutines

	subroutine slmodel_init
	   use sl_model_mod !HH
		! declare variables
		integer :: iter_sh, itersl_sh, dtime_sh
		real    :: starttime_sh
		! check if we are on the head processor
		! get initial topography from MALI 
		     ! gather topography to the headnode
		! interpolate init topo from MALI grid to SLM grid
		! pass interpolated topo from MALI to SLM
		iter_sh = 0 
		itersl_sh = 1	   ! * note: hardcode itersl_sh == 1. 
		starttime_sh = 0 ! get itersl_sh,, starttime_sh (from MALI simulation, needs manipulation) 
	   dtime_sh = 100
		
		call sl_solver_checkpoint(itersl_sh, dtime_sh)
		call set_planet 
		call sl_timewindow(iter_sh)
			 ! * note: read in namelist values. SLM should return dtime from its own namelist file and MALI stores the values
	   call sl_solver_init(itersl_sh, starttime_sh)
	end subroutine slmodel_init
	
	subroutine slmodel_solve
		! check if we are on the head node
		! check if SLM should be called. If SLM should be called, then:
		   ! get new iceload from MALI
			    !gather iceload to the head node 
	    	! interpolate iceload from MALI to SLM grid
		   ! pass interpolated iceload from MALI to SLM 
		   ! get itersl_sh, iter_sh, dtime_sh, starttime_sh
			! call set_planet & sl_timewindow
         ! call sl_solver
			! read in and interpolate the new bedrock topography from SLM to MALI grid
			   ! option 1 (initial choice): set bed topography relative to geoid 
				! option 2: keep track of both surfaces or one? 
			! scatter topography to all processors. 
	end subroutine slmodel_solve
	
	! last step (refer to the MPAS Ocean model)
	subroutine interp_li_to_slm
	end subroutine interp_li_to_slm
	
	subroutine interp_slm_to_li
	end subroutine interp_slm_to_li
	




!***********************************************************************

end module li_bedtopo

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection
!
!> \brief MPAS ocean tracer advection driver
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This module contains initialization and driver routines for
!>  computing tracer advection tendencies. It primarily calls submodule
!>  routines based on user choices for advection options.
!
!-----------------------------------------------------------------------

module ocn_tracer_advection

   ! module includes
   use mpas_kind_types
   use mpas_timer

   use ocn_tracer_advection_shared
   use ocn_tracer_advection_std
   use ocn_tracer_advection_mono
   use ocn_tracer_advection_vert

   use ocn_config
   use ocn_mesh

   implicit none
   private
   save

   ! public module method interfaces
   public :: ocn_tracer_advection_init,         &
             ocn_tracer_advection_tend

   ! privat module variables
   logical :: tracerAdvOff  !< flag to turn off tracer advection
   logical :: monotonicOn   !< flag to choose a monotone advection scheme

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine is the driver routine for computing tracer advection
!>  tendencies. It simply calls submodule tendency routines based on
!>  choice of algorithm.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_tend(tracers, normalThicknessFlux, &
                                        w, layerThickness, dt, tend,  &
                                        computeBudgets)!{{{

      !*** Input/Output parameters

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend            !< [in,out] Tracer tendency to which advection added

      !*** Input parameters

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers               !< [in] Current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux, &!< [in] Thickness weighted horiz velocity
         w,                   &!< [in] Vertical velocity
         layerThickness        !< [in] Thickness field

      real (kind=RKIND), intent(in) :: &
         dt                    !< [in] Time step

      logical, intent(in) :: &
         computeBudgets      ! flag to compute active tracer budget

      ! end of preamble
      !----------------
      ! begin code

      ! immediate return if tracer advection not selected
      if (tracerAdvOff) return

      call mpas_timer_start("tracer adv")

      ! call specific advection routine based on choice of monotonicity
      if (monotonicOn) then
         call ocn_tracer_advection_mono_tend(tend, tracers, layerThickness,    &
                                             normalThicknessFlux, w, dt,       &
                                             computeBudgets)

      else
         call ocn_tracer_advection_std_tend(tracers, normalThicknessFlux, &
                                            w, layerThickness, dt, tend)
      endif

      call mpas_timer_stop("tracer adv")

   end subroutine ocn_tracer_advection_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_advection_init
!
!> \brief MPAS ocean tracer advection tendency
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine is the driver routine for initializing various tracer
!>  advection choices and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_advection_init(err)!{{{

      !*** output parameters

      integer, intent(out) :: err !< [out] Error flag

      !*** local variables

      integer :: err1, err2, err3, err4 ! local error flags

      ! end preamble
      !-------------
      ! begin code

      err = 0 ! initialize error code to success

      ! set some basic flags for options
      tracerAdvOff = config_disable_tr_adv
      if (config_flux_limiter == 'monotonic') then
         monotonicOn = .true.
      elseif (config_flux_limiter == 'none') then
         monotonicOn = .false.
      else
         CALL mpas_log_write('Advection limiter not supported', &
            MPAS_LOG_CRIT)
      endif

      ! set all other options from submodule initialization routines
      call ocn_tracer_advection_shared_init(err1)
      call ocn_tracer_advection_std_init   (err2)
      call ocn_tracer_advection_mono_init  (err3)
      call ocn_tracer_advection_vert_init  (err4)

      ! if an error is returned from init routines, write an error
      ! message and return a non-zero error code
      if (err1 /= 0 .or. err2 /= 0 .or. err3 /=0 .or. err4 /= 0) then
         err = 1
         call mpas_log_write(                                 &
            'Error encountered during tracer advection init', &
            MPAS_LOG_ERR, masterOnly=.true.)
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_advection_init!}}}

!***********************************************************************

end module ocn_tracer_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

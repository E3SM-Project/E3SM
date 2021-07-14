! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_test_sht
!
!> \brief MPAS ocean initialize case -- Cosine Bell
!> \author Luke Van Roekel
!> \date   09/01/2020
!> \details
!>  This module contains the routines for initializing the
!>  the cosine bell test case
!>  Reference: Section 2a in Skamarock, W.C. and A. Gassmann, 2011: Conservative
!>  Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for
!>  ODE-Based Time Integration. Mon. Wea. Rev., 139, 2962–2975,
!>  https://doi.org/10.1175/MWR-D-10-05056.1
!
!-----------------------------------------------------------------------

module ocn_init_test_sht

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
   use mpas_dmpar

   use ocn_init_cell_markers
   use ocn_init_vertical_grids
   use ocn_vel_self_attraction_loading

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_test_sht, &
             ocn_init_validate_test_sht

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_test_sht
!
!> \brief   Set-up for test case
!> \author  Luke Van Roekel
!> \date    09/01/2020
!> \details
!>  Reference: Section 2a in Skamarock, W.C. and A. Gassmann, 2011: Conservative
!>  Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for
!>  ODE-Based Time Integration. Mon. Wea. Rev., 139, 2962–2975,
!>  https://doi.org/10.1175/MWR-D-10-05056.1
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_test_sht(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool
      type (mpas_pool_type), pointer :: forcingPool, diagnosticsPool

      integer, pointer :: nCellsSolve
      real (kind=RKIND), dimension(:), pointer :: lonCell, latCell
      real (kind=RKIND), dimension(:), pointer :: ssh, ssh_sal
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), pointer :: config_test_sht_cosine_bell_lat_center, &
                                    config_test_sht_cosine_bell_lon_center, &
                                    config_test_sht_cosine_bell_psi0, &
                                    config_test_sht_cosine_bell_radius
      integer, pointer :: config_test_sht_function_option
      integer, pointer :: config_test_sht_n_iterations
      integer, pointer :: config_parallel_self_attraction_loading_order 
      integer, pointer :: config_parallel_self_attraction_loading_n_cells_per_block
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), pointer :: cellWidthEq, cellWidthMidLat 
      real (kind=RKIND), pointer :: latTransitionStart, latTransitionWidth 

      integer :: i                                                  
      integer :: iCell
      integer :: err
      integer :: nCells

      real (kind=RKIND) :: temp
      real (kind=RKIND) :: cellWidth
      real (kind=RKIND) :: lat_center
      real (kind=RKIND) :: error_local, error_global
      real (kind=RKIND) :: area_local, area_global
      real (kind=RKIND) :: sinlat_highres,sinlat_lowres,sinlon_highres,sinlon_lowres
      real (kind=RKIND) :: trans_width,trans_start,x,w1,w2,rad


      character (len=StrKIND), pointer :: config_init_configuration


      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs,'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('test_sht')) return

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs,'config_test_sht_cosine_bell_lat_center', &
                                config_test_sht_cosine_bell_lat_center)
      call mpas_pool_get_config(domain % configs,'config_test_sht_cosine_bell_lon_center', &
                                config_test_sht_cosine_bell_lon_center)
      call mpas_pool_get_config(domain % configs,'config_test_sht_cosine_bell_psi0', &
                                config_test_sht_cosine_bell_psi0)
      call mpas_pool_get_config(domain % configs,'config_test_sht_cosine_bell_radius', &
                                config_test_sht_cosine_bell_radius)
      call mpas_pool_get_config(domain % configs,'config_test_sht_function_option', &
                                config_test_sht_function_option)
      call mpas_pool_get_config(domain % configs,'config_parallel_self_attraction_loading_order', &
                                config_parallel_self_attraction_loading_order)
      call mpas_pool_get_config(domain % configs,'config_parallel_self_attraction_loading_n_cells_per_block', &
                                config_parallel_self_attraction_loading_n_cells_per_block)
      call mpas_pool_get_config(domain % configs,'config_test_sht_n_iterations', &
                                config_test_sht_n_iterations)
      call mpas_pool_get_config(domain % configs,'config_test_sht_function3_cell_width_equator', &
                                cellWidthEq)
      call mpas_pool_get_config(domain % configs,'config_test_sht_function3_cell_width_pole', &
                                cellWidthMidLat)
      call mpas_pool_get_config(domain % configs,'config_test_sht_function3_lat_transition_start', &
                                latTransitionStart)
      call mpas_pool_get_config(domain % configs,'config_test_sht_function3_lat_transition_width', &
                                latTransitionWidth)

      !load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh',meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state',statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing',forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics',diagnosticsPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve',nCellsSolve)
        call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)

        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

        call mpas_pool_get_array(statePool, 'ssh', ssh)
        call mpas_pool_get_array(diagnosticsPool, 'ssh_sal', ssh_sal)

        
        ! Setup function to approximate
        if (config_test_sht_function_option == 1) then
            do iCell = 1, nCellsSolve
              temp = sphere_radius*acos(sin(config_test_sht_cosine_bell_lat_center)*sin(latCell(iCell)) + &
                      cos(config_test_sht_cosine_bell_lat_center)*cos(latCell(iCell))*cos(lonCell(iCell) - &
                      config_test_sht_cosine_bell_lon_center))
    
                ssh(iCell) = config_test_sht_cosine_bell_psi0/2.0_RKIND * ( 1.0_RKIND + &
                  cos(3.1415926_RKIND*temp/config_test_sht_cosine_bell_radius))
            enddo
    
            lat_center = config_test_sht_cosine_bell_lat_center + pii*0.5_RKIND
            do iCell = 1, nCellsSolve
              temp = sphere_radius*acos(sin(lat_center)*sin(latCell(iCell)) + &
                      cos(lat_center)*cos(latCell(iCell))*cos(lonCell(iCell) - &
                      config_test_sht_cosine_bell_lon_center))
    
                ssh(iCell) = ssh(iCell) + config_test_sht_cosine_bell_psi0/2.0_RKIND * ( 1.0_RKIND + &
                  cos(3.1415926_RKIND*temp/config_test_sht_cosine_bell_radius))
            enddo
        elseif (config_test_sht_function_option == 2) then
               do iCell = 1, nCellsSolve
                 temp = sphere_radius*acos(sin(config_test_sht_cosine_bell_lat_center)*sin(latCell(iCell)) + &
                         cos(config_test_sht_cosine_bell_lat_center)*cos(latCell(iCell))*cos(lonCell(iCell) - &
                         config_test_sht_cosine_bell_lon_center))
       
                 if( temp < config_test_sht_cosine_bell_radius ) then
                  ssh(iCell) = config_test_sht_cosine_bell_psi0/2.0_RKIND * ( 1.0_RKIND + &
                    cos(3.1415926_RKIND*temp/config_test_sht_cosine_bell_radius))
                 else
                   ssh(iCell) = 0.0_RKIND
                 endif
               enddo
        elseif (config_test_sht_function_option == 3) then

          
          rad = sphere_radius/1000.0_RKIND

          do iCell = 1,nCellsSolve
            sinlat_highres = sin(floor(2.0_RKIND*pii*rad/cellWidthEq/30.0_RKIND)*latCell(iCell))
            sinlat_lowres =  sin(floor(2.0_RKIND*pii*rad/cellWidthMidLat/30.0_RKIND)*latCell(iCell))
            sinlon_highres = sin(floor(2.0_RKIND*pii*rad/cellWidthEq/30.0_RKIND)*lonCell(iCell))
            sinlon_lowres =  sin(floor(2.0_RKIND*pii*rad/cellWidthMidLat/30.0_RKIND)*lonCell(iCell))

            x = abs(latCell(iCell))
            trans_start = pii/180.0_RKIND*latTransitionStart
            trans_width = pii/180.0_RKIND*latTransitionWidth
            w1 = 0.5_RKIND*(tanh((x - trans_start - .5 * trans_width) / (.2 *trans_width)) + 1.0)
            w2 = 1.0_RKIND - w1

            ssh(iCell) =  sinlat_highres*w2 + sinlat_lowres*w1 + sinlon_highres*w2 + sinlon_lowres*w1

            trans_start = pii/180.0_RKIND*40.0_RKIND
            trans_width = pii/180.0_RKIND*latTransitionWidth
            w1 = 0.5_RKIND*(tanh((x - trans_start - .5 * trans_width) / (.2 *trans_width)) + 1.0)
            w2 = 1.0_RKIND - w1
            ssh(iCell) = w1*sinlat_lowres + w2*ssh(iCell)

          enddo

        else
          call mpas_log_write('config_test_sht_function_option not supported')
          iErr = 1
        endif


     call ocn_vel_self_attraction_loading_init(domain,err)
     do i = 1,config_test_sht_n_iterations
       call mpas_log_write('Iteration = $i',intArgs=(/i/))
       call ocn_compute_self_attraction_loading(domain, forcingPool, domain % dminfo, ssh, err)
     enddo

     error_local = 0.0_RKIND
     do iCell = 1,nCellsSolve
       error_local = error_local + (ssh(iCell)-ssh_sal(iCell))**2
     enddo 

     call mpas_dmpar_sum_real(domain % dminfo,error_local,error_global)
     call mpas_dmpar_sum_int(domain % dminfo,nCellsSolve,nCells)
     error_global = sqrt(error_global/real(nCells,RKIND))


     call mpas_log_write('nCells = $i',intArgs=(/nCells/))
     call mpas_log_write('function = $i',intArgs=(/config_test_sht_function_option/))
     call mpas_log_write('nOrder = $i',intArgs=(/config_parallel_self_attraction_loading_order/))
     call mpas_log_write('cells per block = $i',intArgs=(/config_parallel_self_attraction_loading_n_cells_per_block/))
     call mpas_log_write('error = $r',realArgs=(/error_global/))


     area_local = 0.0_RKIND
     do iCell = 1,nCellsSolve
      area_local = area_local + areaCell(iCell)
     enddo 

     call mpas_dmpar_sum_real(domain % dminfo, area_local, area_global)
     call mpas_log_write('area error = $r',realArgs=(/4.0_RKIND*pii-area_global/sphere_radius**2/))
     call mpas_log_write('sphere radius = $r',realArgs=(/sphere_radius/))

     block_ptr => block_ptr % next
   end do

  !--------------------------------------------------------------------

   end subroutine ocn_init_setup_test_sht!}}}

!***********************************************************************
!
!  routine ocn_init_validate_test_sht
!
!> \brief   Validation for test case
!> \author  Luke Van Roekel
!> \date    09/01/2020
!> \details
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_test_sht(configPool, packagePool,iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration',config_init_configuration)

      if(config_init_configuration .ne. trim('test_sht')) return

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_test_sht!}}}

!***********************************************************************

end module ocn_init_test_sht

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

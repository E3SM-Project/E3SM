! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_cases

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_abort, only : mpas_dmpar_global_abort
   use mpas_sort
   use atm_advection
   use mpas_RBF_interpolation
   use mpas_vector_reconstruction
   use mpas_timer
   use mpas_init_atm_static
   use mpas_init_atm_surface
   use mpas_atmphys_constants, only: svpt0,svp1,svp2,svp3
   use mpas_atmphys_functions
   use mpas_atmphys_initialize_real

   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping !, only: MPAS_Time_type, MPAS_TimeInterval_type, MPAS_Clock_type, &
                        !        mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti
   


   contains


   subroutine init_atm_setup_case(domain, stream_manager)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the hydrostatic test case
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_stream_manager
      use mpas_init_atm_gwd, only : compute_gwd_fields

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager


      integer :: i
      integer :: ierr
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: fg
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: diag_physics

      integer, pointer :: config_init_case
      logical, pointer :: config_static_interp
      logical, pointer :: config_native_gwd_static
      logical, pointer :: config_met_interp

      character(len=StrKIND), pointer :: mminlu

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: nVertLevels


      call mpas_pool_get_config(domain % blocklist % configs, 'config_init_case', config_init_case)

      !
      ! Do some quick checks to make sure compile options are compatible with the chosen test case
      !
      if (config_init_case == 6) then
#ifndef ROTATED_GRID
         call mpas_dmpar_global_abort('*****************************************************************',  deferredAbort=.true.)
         call mpas_dmpar_global_abort('To initialize and run the mountain wave test case (case 6),',        deferredAbort=.true.)
         call mpas_dmpar_global_abort('   please clean and re-compile init_atmosphere with -DROTATED_GRID', deferredAbort=.true.)
         call mpas_dmpar_global_abort('   added to the specification of MODEL_FORMULATION',                 deferredAbort=.true.)
         call mpas_dmpar_global_abort('   at the top of the Makefile.',                                     deferredAbort=.true.)
         call mpas_dmpar_global_abort('*****************************************************************')
#endif
      else
#ifdef ROTATED_GRID
         call mpas_dmpar_global_abort('*****************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Only test case 6 should use code compiled with -DROTATED_GRID',     deferredAbort=.true.)
         call mpas_dmpar_global_abort('   specified in the Makefile.',                                     deferredAbort=.true.)
         call mpas_dmpar_global_abort('*****************************************************************')
#endif
      end if



      if ((config_init_case == 1) .or. (config_init_case == 2) .or. (config_init_case == 3)) then

         write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
         if (config_init_case == 1) write(0,*) ' no initial perturbation '
         if (config_init_case == 2) write(0,*) ' initial perturbation included '
         if (config_init_case == 3) write(0,*) ' normal-mode perturbation included '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call mpas_pool_get_config(block_ptr % configs, 'config_static_interp', config_static_interp)
            call mpas_pool_get_config(block_ptr % configs, 'config_met_interp', config_met_interp)

            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            write(0,*) ' calling test case setup '
            call init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, block_ptr % configs, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if ((config_init_case == 4) .or. (config_init_case == 5)) then

         write(0,*) ' squall line - super cell test case '
         if (config_init_case == 4) write(0,*) ' squall line test case' 
         if (config_init_case == 5) write(0,*) ' supercell test case'
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            write(0,*) ' calling test case setup '
            call init_atm_case_squall_line(domain % dminfo, mesh, nCells, nVertLevels, state, diag, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 6 ) then

         write(0,*) ' mountain wave test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))

            call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            write(0,*) ' calling test case setup '
            call init_atm_case_mtn_wave(domain % dminfo, mesh, nCells, nVertLevels, state, diag, block_ptr % configs, config_init_case)
            call decouple_variables(mesh, nCells, nVertLevels, state, diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 7 ) then


         write(0,*) ' real-data GFS test case '
         block_ptr => domain % blocklist

         do while (associated(block_ptr))

            call mpas_pool_get_config(block_ptr % configs, 'config_static_interp', config_static_interp)
            call mpas_pool_get_config(block_ptr % configs, 'config_native_gwd_static', config_native_gwd_static)
            call mpas_pool_get_config(block_ptr % configs, 'config_met_interp', config_met_interp)

            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'fg', fg)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)

            call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

            if (config_static_interp) then

               !  
               !  Without a convex mesh partition file, interpolating static fields in parallel 
               !     will give incorrect results. Since it is very unlikely that typical users
               !     will have convex partitions, it's safer to just stop if multiple MPI tasks are
               !     detected when performing the static_interp step.
               !  
               if (domain % dminfo % nprocs > 1) then
                 call mpas_dmpar_global_abort('****************************************************************', deferredAbort=.true.)
                 call mpas_dmpar_global_abort('Error: Interpolation of static fields does not work in parallel.', deferredAbort=.true.)
                 call mpas_dmpar_global_abort('Please run the static_interp step using only a single MPI task.',  deferredAbort=.true.)
                 call mpas_dmpar_global_abort('****************************************************************')
               end if

               call init_atm_static(mesh, block_ptr % dimensions, block_ptr % configs)
               call init_atm_static_orogwd(mesh, block_ptr % dimensions, block_ptr % configs)
            end if

            if (config_native_gwd_static) then
               write(0,*) ' '
               write(0,*) 'Computing GWDO static fields on the native MPAS mesh'
               write(0,*) ' '
               ierr = compute_gwd_fields(domain)
               if (ierr /= 0) then
                  call mpas_dmpar_global_abort('****************************************************************', deferredAbort=.true.)
                  call mpas_dmpar_global_abort('Error while trying to compute sub-grid-scale orography',           deferredAbort=.true.)
                  call mpas_dmpar_global_abort('statistics for use with the GWDO scheme.',                         deferredAbort=.true.)
                  call mpas_dmpar_global_abort('****************************************************************')
               end if
            end if

            !
            ! If at this point the mminlu variable is blank, we assume that the static interp step was
            !   not run, and that we are working with a static file created before there was a choice
            !   of land use datasets; in this case, the dataset was almost necessarily USGS
            !
            call mpas_pool_get_array(mesh, 'mminlu', mminlu)
            if (len_trim(mminlu) == 0) then
                  write(0,*) '****************************************************************'
                  write(0,*) 'No information on land use dataset is available.'
                  write(0,*) 'Assume that we are using ''USGS''.'
                  write(0,*) '****************************************************************'
                  write(mminlu,'(a)') 'USGS'
            end if

            call init_atm_case_gfs(block_ptr, mesh, nCells, nEdges, nVertLevels, fg, state, &
                                   diag, diag_physics, config_init_case, block_ptr % dimensions, block_ptr % configs)
            
            if (config_met_interp) then
               call physics_initialize_real(mesh, fg, domain % dminfo, block_ptr % dimensions, block_ptr % configs)
            end if

            block_ptr => block_ptr % next
         end do

      else if (config_init_case == 8 ) then

         write(0,*) 'real-data surface (SST) update test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'fg', fg)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)

            ! Defined in mpas_init_atm_surface.F
            call init_atm_case_sfc(domain, domain % dminfo, stream_manager, mesh, fg, state, block_ptr % dimensions, block_ptr % configs)
            block_ptr => block_ptr % next
         end do

      else

         call mpas_dmpar_global_abort(' ****************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort(' Only test cases 1 through 8 are currently supported.', deferredAbort=.true.)
         call mpas_dmpar_global_abort(' ****************************************************')

      end if


      !initialization of surface input variables technically not needed to run our current set of
      !idealized test cases:
      if (config_init_case < 7)  then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'fg', fg)

            call physics_idealized_init(mesh, fg)

            block_ptr => block_ptr % next
         end do
      end if

   end subroutine init_atm_setup_case

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, configs, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

!      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:), pointer :: u, ru, w, rw, v
      real (kind=RKIND), dimension(:,:), pointer :: rho, theta
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      
!.. initialization of moisture:
      integer, pointer :: index_qv
      real (kind=RKIND), parameter :: rh_max = 0.40 ! Maximum relative humidity
!      real (kind=RKIND), parameter :: rh_max = 0.70 ! Maximum relative humidity
      real (kind=RKIND), dimension(nVertLevels, nCells) :: qsat, relhum
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
!.. end initialization of moisture.

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
      integer, pointer :: nz1, nCellsSolve, nEdges, maxEdges, nVertices

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, verticesOnEdge, cellsOnCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND) :: es, qvs, xnutr, znut, ptemp 
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(nVertLevels ) :: eta, etav, teta, ppi, tt, temperature_1d

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf
      real (kind=RKIND), pointer :: cf1, cf2, cf3

      !  storage for (lat,z) arrays for zonal velocity calculation

      logical, parameter :: rebalance = .true.
      integer, parameter :: nlat=721
      real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nVertLevels + 1, nlat) :: zgrid_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: u_2d, pp_2d, rho_2d, qv_2d, etavs_2d, zz_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: p_2d, pb_2d, ppb_2d, rr_2d, rb_2d, tb_2d, rtb_2d
      real (kind=RKIND), dimension(nVertLevels, nlat-1) :: zx_2d 
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat, hx_1d
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

!      logical, parameter :: moisture = .true.
      logical, parameter :: moisture = .false.

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:), pointer :: latCell, latVertex, lonVertex, latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: fEdge, fVertex

      real (kind=RKIND), pointer :: sphere_radius
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: config_init_case


      call mpas_pool_get_config(configs, 'config_init_case', config_init_case)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)


      !
      ! Scale all distances and areas from a unit sphere to one with radius sphere_radius
      !
      call mpas_pool_get_array(mesh, 'xCell', xCell)
      call mpas_pool_get_array(mesh, 'yCell', yCell)
      call mpas_pool_get_array(mesh, 'zCell', zCell)
      call mpas_pool_get_array(mesh, 'xEdge', xEdge)
      call mpas_pool_get_array(mesh, 'yEdge', yEdge)
      call mpas_pool_get_array(mesh, 'zEdge', zEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      xCell(:) = xCell(:) * sphere_radius
      yCell(:) = yCell(:) * sphere_radius
      zCell(:) = zCell(:) * sphere_radius
      xVertex(:) = xVertex(:) * sphere_radius
      yVertex(:) = yVertex(:) * sphere_radius
      zVertex(:) = zVertex(:) * sphere_radius
      xEdge(:) = xEdge(:) * sphere_radius
      yEdge(:) = yEdge(:) * sphere_radius
      zEdge(:) = zEdge(:) * sphere_radius
      dvEdge(:) = dvEdge(:) * sphere_radius
      dcEdge(:) = dcEdge(:) * sphere_radius
      areaCell(:) = areaCell(:) * sphere_radius**2.0
      areaTriangle(:) = areaTriangle(:) * sphere_radius**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * sphere_radius**2.0

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nz1)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz = nz1 + 1

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'dss', dss)

      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'v', v)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(state, 'theta_m', t)
      call mpas_pool_get_array(state, 'scalars', scalars)
      call mpas_pool_get_array(state, 'u', u)
      call mpas_pool_get_array(state, 'w', w)

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

!.. initialization of moisture:
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0
!.. end initialization of moisture.

      surface_pressure(:) = 0.0

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv) 

      xnutr = 0.0
      zd = 12000.0
      znut = eta_t

      etavs = (1.0 - 0.252) * pii/2.
      r_earth = sphere_radius
      omega_e = omega
      p0 = 1.0e+05

      write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,nCells
        do k=1,nz
          phi = latCell(iCell)
          hx(k,iCell) = u0/gravity*cos(etavs)**1.5                                   &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        end do
      end do

      !     Metrics for hybrid coordinate and vertical stretching

      str = 1.5
      zt = 45000.
      dz = zt/float(nz1)

      write(0,*) ' hx computation complete '

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            sh(k) = (real(k-1)*dz/zt)**str 

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) and define sh(k) = zc(k)/zt
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

            zw(k) = float(k-1)*dz
!            zw(k) = sh(k)*zt
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
!            ah(k) = 0.
            write(0,*) ' k, sh, zw, ah ',k,sh(k),zw(k),ah(k)
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      write(0,*) ' cf1, cf2, cf3 = ',cf1,cf2,cf3

      do iCell=1,nCells
        do k=1,nz
          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
                         + ah(k) * sh(k)* zt
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

      !do k=1,nz1
      !  write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
      !end do

      !do k=1,nz1
      !  write(0,*) ' k, zx(k,1) ',k,zx(k,1)
      !end do

      write(0,*) ' grid metrics setup complete '

!**************  section for 2d (z,lat) calc for zonal velocity

      dlat = 0.5*pii/float(nlat-1)
      do i = 1,nlat

        lat_2d(i) = float(i-1)*dlat
        phi = lat_2d(i)
        hx_1d    = u0/gravity*cos(etavs)**1.5                           &
                   *((-2.*sin(phi)**6                                   &
                         *(cos(phi)**2+1./3.)+10./63.)                  &
                         *(u0)*cos(etavs)**1.5                          &
                    +(1.6*cos(phi)**3                                   &
                         *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)

        do k=1,nz
          zgrid_2d(k,i) = (1.-ah(k))*(sh(k)*(zt-hx_1d)+hx_1d)  &
                         + ah(k) * sh(k)* zt
        end do
        do k=1,nz1
          zz_2d(k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
          ppb_2d(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb_2d(k,i) = (ppb_2d(k,i)/p0)**(rgas/cp)
          rb_2d(k,i) = ppb_2d(k,i)/(rgas*t0b*zz_2d(k,i))
          tb_2d(k,i) = t0b/pb_2d(k,i)
          rtb_2d(k,i) = rb_2d(k,i)*tb_2d(k,i)
          p_2d(k,i) = pb_2d(k,i)
          pp_2d(k,i) = 0.0
          rr_2d(k,i) = 0.0
        end do


        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb_2d(k,i)+pp_2d(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k) >= znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do

          phi = lat_2d(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))


            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
            ptemp   = ppb_2d(k,i) + pp_2d(k,i)

            !get moisture 
            if (moisture) then
              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
          end do

          do itrp = 1,25
            do k=1,nz1
              rr_2d(k,i)  = (pp_2d(k,i)/(rgas*zz_2d(k,i)) - rb_2d(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                            &
                          *(1.25*(rr_2d(1,i)+rb_2d(1,i))*(1.+qv_2d(1,i))   &
                            -.25*(rr_2d(2,i)+rb_2d(2,i))*(1.+qv_2d(2,i)))

            ppi(1) = ppi(1)-ppb_2d(1,i)
            do k=1,nz1-1

              ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                       &
                            ( (rr_2d(k  ,i)+(rr_2d(k  ,i)+rb_2d(k  ,i))*qv_2d(k  ,i))*fzp(k+1)   &
                            + (rr_2d(k+1,i)+(rr_2d(k+1,i)+rb_2d(k+1,i))*qv_2d(k+1,i))*fzm(k+1) )
            end do

            do k=1,nz1
              pp_2d(k,i) = .2*ppi(k)+.8*pp_2d(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          rho_2d(k,i) = rr_2d(k,i)+rb_2d(k,i)
          etavs_2d(k,i) = ((ppb_2d(k,i)+pp_2d(k,i))/p0 - 0.252)*pii/2.
          u_2d(k,i) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(k,i))**1.5)
        end do

      end do  ! end loop over latitudes for 2D zonal wind field calc

      !SHP-balance:: in case of rebalacing for geostrophic wind component
      if (rebalance) then

        do i=1,nlat-1
          do k=1,nz1
            zx_2d(k,i) = (zgrid_2d(k,i+1)-zgrid_2d(k,i))/(dlat*r_earth)
          end do
        end do

        call init_atm_recompute_geostrophic_wind(u_2d, rho_2d, pp_2d, qv_2d, lat_2d, zz_2d, zx_2d,     &
                                        cf1, cf2, cf3, fzm, fzp, rdzw, nz1, nlat, dlat, sphere_radius)

      end if

!******************************************************************      

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, nCells
        do k=1,nz1
          ztemp    = .5*(zgrid(k+1,i)+zgrid(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do

!       if(i == 1) then
!         do k=1,nz1
!           write(0,*) ' k, ppb, pb, rb, tb (k,1) ',k,ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)
!         end do
!       end if

      200 format(4i6,8(1x,e15.8))
      201 format(3i6,8(1x,e15.8))
      202 format(2i6,10(1x,e15.8))
      203 format(i6,10(1x,e15.8))

!     iterations to converge temperature as a function of pressure
!
        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb(k,i)+pp(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k).ge.znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do
          phi = latCell(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*scalars(index_qv,k,i))

            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture 
            if (moisture) then
 
                !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )

               if(ptemp < 50000.) then
                  relhum(k,i) = 0.0
               elseif(ptemp > p0) then
                  relhum(k,i) = 1.0
               else
                  relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
               end if
               relhum(k,i) = min(rh_max,relhum(k,i))

               !.. calculation of water vapor mixing ratio:
               if (temperature_1d(k) > 273.15) then
                   es  = 1000.*0.6112*exp(17.67*(temperature_1d(k)-273.15)/(temperature_1d(k)-29.65))
               else
                   es  = 1000.*0.6112*exp(21.8745584*(temperature_1d(k)-273.15)/(temperature_1d(k)-7.66))
               end if
               qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
               if(relhum(k,i) .eq. 0.0) qsat(k,i) = 0.0
               scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*scalars(index_qv,k,i))

          end do

          do itrp = 1,25
            do k=1,nz1
              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1

!              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)   &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))

               ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                                  &
                             ( (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))*fzp(k+1)   &
                             + (rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))*fzm(k+1) )

            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          p (k,i) = ((ppb(k,i)+pp(k,i))/p0)**(rgas/cp)
          t (k,i) = tt(k)/p(k,i)
          rt (k,i) = t(k,i)*rr(k,i)+rb(k,i)*(t(k,i)-tb(k,i))
          rho_zz (k,i) = rb(k,i) + rr(k,i)
        end do

        !calculation of surface pressure:
        surface_pressure(i) = 0.5*dzw(1)*gravity                                    &
                        * (1.25*(rr(1,i) + rb(1,i)) * (1. + scalars(index_qv,1,i))  &
                        -  0.25*(rr(2,i) + rb(2,i)) * (1. + scalars(index_qv,2,i)))
        surface_pressure(i) = surface_pressure(i) + pp(1,i) + ppb(1,i)

      end do  ! end loop over cells

      !write(0,*)
      !write(0,*) '--- initialization of water vapor:'
      !do iCell = 1, nCells
      !   if(iCell == 1 .or. iCell == nCells) then
      !      do k = nz1, 1, -1
      !         write(0,202) iCell,k,t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)
      !      end do
      !      write(0,*)
      !   end if
      !end do

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,nEdges

         vtx1 = verticesOnEdge(1,iEdge)
         vtx2 = verticesOnEdge(2,iEdge)
         lat1 = latVertex(vtx1)
         lat2 = latVertex(vtx2)
         iCell1 = cellsOnEdge(1,iEdge)
         iCell2 = cellsOnEdge(2,iEdge)
         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)

         if (config_init_case == 2) then
            r_pert = sphere_distance( latEdge(iEdge), lonEdge(iEdge), &
                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1) * sphere_radius / dvEdge(iEdge)

         else if (config_init_case == 3) then
            lon_Edge = lonEdge(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
         else
            u_pert = 0.0
         end if

         if (rebalance) then

           call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,dvEdge(iEdge),sphere_radius,u0,nz1,nlat)
           do k=1,nVertLevels
             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
             u(k,iEdge) = fluxk + u_pert
           end do

         else 

           do k=1,nVertLevels
             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
             fluxk = u0*flux*(cos(etavs)**1.5)
             u(k,iEdge) = fluxk + u_pert
           end do

         end if

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nz1
            ru(k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,iEdge)
         end do

      !
      ! Generate rotated Coriolis field
      !

         fEdge(iEdge) = 2.0 * omega_e * &
                                       ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha_grid) + &
                                         sin(latEdge(iEdge)) * cos(alpha_grid) &
                                       )
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 2.0 * omega_e * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha_grid) + &
                                          sin(latVertex(iVtx)) * cos(alpha_grid) &
                                         )
      end do

      !
      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
      !

      !
      !     pre-calculation z-metric terms in omega eqn.
      !
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 1, nVertLevels

            if (config_theta_adv_order == 2) then

               z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

            else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

               d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)

!  WCS fix 20120711

                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do             

               z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                             - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               if (config_theta_adv_order == 3) then
                  z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
               else
                  z_edge3 = 0.
               end if

            end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

         end do

      end do

      ! for including terrain
      rw = 0.0
      w = 0.0
      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               rw(k,cell2) = rw(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do


      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
           end do
         end do
      end do

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

        write(0,*) ' i, psurf, lat ',i,psurf,latCell(i)*180./3.1415828
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_jw


   subroutine init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1_in,lat2_in,dvEdge,a,u0,nz1,nlat)

      implicit none
   
      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: u_2d,etavs_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(out) :: flux_zonal
      real (kind=RKIND), intent(in) :: lat1_in, lat2_in, dvEdge, a, u0
   
      integer :: k,i
      real (kind=RKIND) :: lat1, lat2, w1, w2
      real (kind=RKIND) :: dlat,da,db
   
      lat1 = abs(lat1_in)
      lat2 = abs(lat2_in)
      if(lat2 <= lat1) then
        lat1 = abs(lat2_in)
        lat2 = abs(lat1_in)
      end if
   
      do k=1,nz1
        flux_zonal(k) = 0.
      end do
   
      do i=1,nlat-1
        if( (lat1 <= lat_2d(i+1)) .and. (lat2 >= lat_2d(i)) ) then
   
        dlat = lat_2d(i+1)-lat_2d(i)
        da = (max(lat1,lat_2d(i))-lat_2d(i))/dlat
        db = (min(lat2,lat_2d(i+1))-lat_2d(i))/dlat
        w1 = (db-da) -0.5*(db-da)**2
        w2 = 0.5*(db-da)**2
   
        do k=1,nz1
          flux_zonal(k) = flux_zonal(k) + w1*u_2d(k,i) + w2*u_2d(k,i+1)
        end do
   
        end if
   
      end do
   
   !  renormalize for setting cell-face fluxes
   
      do k=1,nz1
        flux_zonal(k) = sign(1.0_RKIND,lat2_in-lat1_in)*flux_zonal(k)*dlat*a/dvEdge/u0
      end do
        
   end subroutine init_atm_calc_flux_zonal


   !SHP-balance
   subroutine init_atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                         cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat,rad)

      implicit none
   
      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(inout) :: u_2d
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: rho_2d, pp_2d, qv_2d, zz_2d
      real (kind=RKIND), dimension(nz1,nlat-1), intent(in) :: zx_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(in) :: fzm, fzp, rdzw
      real (kind=RKIND), intent(in) :: cf1, cf2, cf3, dlat, rad
   
      !local variable
      real (kind=RKIND), dimension(nz1,nlat-1) :: pgrad, ru, u
      real (kind=RKIND), dimension(nlat-1) :: f
      real (kind=RKIND), dimension(nz1+1)  :: dpzx
   
   !   real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e
   
      real (kind=RKIND) :: rdx, qtot, r_earth, phi
      integer :: k,i, itr
   
      r_earth  = rad
      omega_e = omega
      rdx = 1./(dlat*r_earth)
   
      do i=1,nlat-1
        do k=1,nz1
          pgrad(k,i) = rdx*(pp_2d(k,i+1)/zz_2d(k,i+1)-pp_2d(k,i)/zz_2d(k,i))
        end do
   
        dpzx(:) = 0.
   
        k=1
        dpzx(k) = .5*zx_2d(k,i)*(cf1*(pp_2d(k  ,i+1)+pp_2d(k  ,i))        &
                                +cf2*(pp_2d(k+1,i+1)+pp_2d(k+1,i))        &
                                +cf3*(pp_2d(k+2,i+1)+pp_2d(k+2,i)))
        do k=2,nz1
           dpzx(k) = .5*zx_2d(k,i)*(fzm(k)*(pp_2d(k  ,i+1)+pp_2d(k  ,i))   &
                                   +fzp(k)*(pp_2d(k-1,i+1)+pp_2d(k-1,i)))
        end do
   
        do k=1,nz1
           pgrad(k,i) = pgrad(k,i) - rdzw(k)*(dpzx(k+1)-dpzx(k))
        end do
      end do
   
   
      !initial value of v and rv -> that is from analytic sln. 
      do i=1,nlat-1
         do k=1,nz1
            u(k,i) = .5*(u_2d(k,i)+u_2d(k,i+1))
            ru(k,i) = u(k,i)*(rho_2d(k,i)+rho_2d(k,i+1))*.5
         end do
      end do
   
      write(0,*) "MAX U wind before REBALANCING ---->", maxval(abs(u))
   
      !re-calculate geostrophic wind using iteration 
      do itr=1,50
      do i=1,nlat-1
         phi = (lat_2d(i)+lat_2d(i+1))/2.
         f(i) = 2.*omega_e*sin(phi)
         do k=1,nz1
            if (f(i).eq.0.) then
              ru(k,i) = 0.
            else
              qtot = .5*(qv_2d(k,i)+qv_2d(k,i+1))
              ru(k,i) = - ( 1./(1.+qtot)*pgrad(k,i) + tan(phi)/r_earth*u(k,i)*ru(k,i) )/f(i)
            end if
              u(k,i) = ru(k,i)*2./(rho_2d(k,i)+rho_2d(k,i+1))
         end do
      end do
      end do
   
      write(0,*) "MAX U wind after REBALANCING ---->", maxval(abs(u))
   
      !update 2d ru
      do i=2,nlat-1
        do k=1,nz1
          u_2d(k,i) = (ru(k,i-1)+ru(k,i))*.5
        end do
      end do
   
      i=1
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i+1)-u_2d(k,i+2))*.5
      end do
      i=nlat
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i-1)-u_2d(k,i-2))*.5
      end do

   end subroutine init_atm_recompute_geostrophic_wind


   subroutine init_atm_case_squall_line(dminfo, mesh, nCells, nVertLevels, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup squall line and supercell test case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, cell1, cell2
      integer, pointer :: nEdges, nVertices, maxEdges, nCellsSolve
      integer, pointer :: index_qv

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znu, znw, znwc, znwv
      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znuc, znuv

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND), dimension(nVertLevels, nCells) :: relhum, thi, tbi, cqwb

      real (kind=RKIND) ::  r, xnutr
      real (kind=RKIND) ::  ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(nVertLevels ) :: qvb
      real (kind=RKIND), dimension(nVertLevels ) :: t_init_1d

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2
      real (kind=RKIND), pointer :: cf1, cf2, cf3
      real (kind=RKIND) :: ztr, thetar, ttr, thetas, um, us, zts, pitop, pibtop, ptopb, ptop, rcp, rcv, p0
      real (kind=RKIND) :: radx, radz, zcent, xmid, delt, xloc, rad, yloc, ymid, a_scale
      real (kind=RKIND) :: pres, temp, es, qvs

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius

      real (kind=RKIND), dimension(:,:), pointer :: t_init, w, rw, v, rho, theta
      real (kind=RKIND), dimension(:), pointer :: u_init, qv_init, angleEdge, fEdge, fVertex

      call mpas_pool_get_array(mesh, 'xCell', xCell)
      call mpas_pool_get_array(mesh, 'yCell', yCell)
      call mpas_pool_get_array(mesh, 'zCell', zCell)
      call mpas_pool_get_array(mesh, 'xEdge', xEdge)
      call mpas_pool_get_array(mesh, 'yEdge', yEdge)
      call mpas_pool_get_array(mesh, 'zEdge', zEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

      call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      !
      ! Scale all distances
      !

      a_scale = 1.0

      xCell(:) = xCell(:) * a_scale
      yCell(:) = yCell(:) * a_scale
      zCell(:) = zCell(:) * a_scale
      xVertex(:) = xVertex(:) * a_scale
      yVertex(:) = yVertex(:) * a_scale
      zVertex(:) = zVertex(:) * a_scale
      xEdge(:) = xEdge(:) * a_scale
      yEdge(:) = yEdge(:) * a_scale
      zEdge(:) = zEdge(:) * a_scale
      dvEdge(:) = dvEdge(:) * a_scale
      dcEdge(:) = dcEdge(:) * a_scale
      areaCell(:) = areaCell(:) * a_scale**2.0
      areaTriangle(:) = areaTriangle(:) * a_scale**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * a_scale**2.0

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz1 = nVertLevels
      nz = nz1 + 1

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'dss', dss)
      call mpas_pool_get_array(mesh, 't_init', t_init)
      call mpas_pool_get_array(mesh, 'u_init', u_init)
      call mpas_pool_get_array(mesh, 'qv_init', qv_init)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'cqw', cqw)

      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(state, 'theta_m', t)
      call mpas_pool_get_array(state, 'u', u)
      call mpas_pool_get_array(state, 'w', w)
      call mpas_pool_get_array(state, 'scalars', scalars)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)
      
      xnutr = 0.
      zd = 12000.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

     write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,nCells
        do k=1,nz
          hx(k,iCell) = 0.  ! squall line or supercell on flat plane
        end do
      end do

      !     metrics for hybrid coordinate and vertical stretching

      str = 1.0
      zt = 20000.
      dz = zt/float(nz1)

!      write(0,*) ' dz = ',dz
      write(0,*) ' hx computation complete '

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str 

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) 
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!            write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      do iCell=1,nCells
        do k=1,nz
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

!
! convective initialization
!
         ztr    = 12000.
         thetar = 343.
         ttr    = 213.
         thetas = 300.5

!         write(0,*) ' rgas, cp, gravity ',rgas,cp, gravity

      if ( test_case == 4) then ! squall line parameters
         um = 12.
         us = 10.
         zts = 2500.
      else if (test_case == 5) then !supercell parameters
         um = 30.
         us = 15.
         zts = 5000.
      end if

         do i=1,nCells
            do k=1,nz1
               ztemp = .5*(zgrid(k,i)+zgrid(k+1,i))
               if(ztemp .gt. ztr) then
                  t (k,i) = thetar*exp(9.8*(ztemp-ztr)/(1003.*ttr))
                  relhum(k,i) = 0.25
               else
                  t (k,i) = 300.+43.*(ztemp/ztr)**1.25
                  relhum(k,i) = (1.-0.75*(ztemp/ztr)**1.25)
                  if(t(k,i).lt.thetas) t(k,i) = thetas
               end if
               tb(k,i) = t(k,i)
               thi(k,i) = t(k,i)
               tbi(k,i) = t(k,i)
               cqw(k,i) = 1.
               cqwb(k,i) = 1.
            end do
         end do

!         relhum(:,:) = 0.

!  set the velocity field - we are on a plane here.

         do i=1, nEdges
            cell1 = cellsOnEdge(1,i)
            cell2 = cellsOnEdge(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               if(ztemp.lt.zts)  then
                  u(k,i) = um*ztemp/zts
               else
                  u(k,i) = um
               end if
               if(i == 1 ) u_init(k) = u(k,i) - us
               u(k,i) = cos(angleEdge(i)) * (u(k,i) - us)
            end do
            end if
         end do

         call mpas_dmpar_bcast_reals(dminfo, nz1, u_init)

!
!    for reference sounding 
!
     do itr=1,30

      pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
      pibtop = 1.-.5*dzw(1)*gravity*(1.+qvb(1))/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t(k,1)+t(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))
         pibtop = pibtop-dzu(k)*gravity/(cp*cqwb(k,1)*.5*(tb(k,1)+tb(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))

         !write(0,*) k,pitop,tb(k,1),dzu(k),tb(k,1)
      end do
      pitop = pitop-.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
      pibtop = pibtop-.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,1)*zz(nz1,1))

      call mpas_dmpar_bcast_real(dminfo, pitop)
      call mpas_dmpar_bcast_real(dminfo, pibtop)

      ptopb = p0*pibtop**(1./rcp)
      write(6,*) 'ptopb = ',.01*ptopb

      do i=1, nCells
         pb(nz1,i) = pibtop+.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,i))/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*cqwb(k+1,i)*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*cqw(k+1,i)*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            ppb(k,i) = p0*(zz(k,i)*rgas*rtb(k,i)/p0)**(cp/cv)
         end do
      end do

     !
     ! update water vapor mixing ratio from humidity profile
     !
      do i= 1,nCells
         do k=1,nz1
            temp     = p(k,i)*thi(k,i)
            pres     = p0*p(k,i)**(1./rcp)
            qvs      = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
            scalars(index_qv,k,i) = min(0.014_RKIND,relhum(k,i)*qvs)
         end do
      end do

      do k=1,nz1
!*********************************************************************
!           QVB = QV INCLUDES MOISTURE IN REFERENCE STATE
!            qvb(k) = scalars(index_qv,k,1)
                                        
!           QVB = 0 PRODUCES DRY REFERENCE STATE
            qvb(k) = 0.
!*********************************************************************
      end do

      do i= 1,nCells
         do k=1,nz1
            t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
            tb(k,i) = tbi(k,i)*(1.+1.61*qvb(k))
         end do
         do k=2,nz1
            cqw (k,i) = 1./(1.+.5*(scalars(index_qv,k,i)+scalars(index_qv,k-1,i)))
            cqwb(k,i) = 1./(1.+.5*(qvb(k)+qvb(k-1)))
         end do
      end do

      end do !end of iteration loop

      write(0,*) ' base state sounding '
      write(0,*) ' k,     pb,     rb,     tb,     rtb,     t,     rr,      p,    qvb'
      do k=1,nVertLevels
         write (0,'(i2,8(2x,f19.15))') k,pb(k,1),rb(k,1),tb(k,1),rtb(k,1),t(k,1),rr(k,1),p(k,1),qvb(k)
      end do

!
!     potential temperature perturbation
!
!      delt = -10.
!      delt = -0.01
      delt = 3.
      radx  = 10000.
      radz  = 1500.
      zcent = 1500.

      if (test_case == 4) then          ! squall line prameters
         call mpas_dmpar_max_real(dminfo, maxval(xCell(:)), xmid)
         xmid = xmid * 0.5
         ymid = 0.0          ! Not used for squall line
      else if (test_case == 5) then     ! supercell parameters
         call mpas_dmpar_max_real(dminfo, maxval(xCell(:)), xmid)
         call mpas_dmpar_max_real(dminfo, maxval(yCell(:)), ymid)
         xmid = xmid * 0.5
         ymid = ymid * 0.5
      end if

      do i=1, nCells
        xloc = xCell(i) - xmid
        if (test_case == 4) then 
           yloc = 0.                            !squall line setting
        else if (test_case == 5) then
           yloc = yCell(i) - ymid !supercell setting
        end if

        do k = 1,nz1
          ztemp     = .5*(zgrid(k+1,i)+zgrid(k,i))
          rad =sqrt((xloc/radx)**2+(yloc/radx)**2+((ztemp-zcent)/radz)**2)
          if(rad.lt.1)  then
            thi(k,i) = thi(k,i) + delt*cos(.5*pii*rad)**2
          end if
           t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
        end do
      end do

      do itr=1,30

        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t (k,1)+t (k-1,1)) &
                                                  *.5*(zz(k,1)+zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)
        write(0,*) 'ptop  = ',.01*ptop, .01*ptopb

        call mpas_dmpar_bcast_real(dminfo, ptop)

        do i = 1, nCells

          pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
          do k=nz1-1,1,-1
!             pp(k,i) = pp(k+1,i)+.5*dzu(k+1)*gravity*                   &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)  &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))
               pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*(    &
                            fzm(k+1)*(rb(k+1,i)*(scalars(index_qv,k+1,i)-qvb(k+1))    &
                                     +rr(k+1,i)*(1.+scalars(index_qv,k+1,i)))         &
                           +fzp(k+1)*(rb(k  ,i)*(scalars(index_qv,k  ,i)-qvb(k))      &
                                     +rr(k  ,i)*(1.+scalars(index_qv,k  ,i))))
          end do
          if (itr==1.and.i==1) then
          do k=1,nz1
          write(0,*) "pp-check", pp(k,i) 
          end do
          end if
          do k=1,nz1
             rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                     -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)       
             p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
             rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
          end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      do k=1,nz1
         qv_init(k) = scalars(index_qv,k,1)
      end do

      t_init_1d(:) = t(:,1)
      call mpas_dmpar_bcast_reals(dminfo, nz1, t_init_1d)
      call mpas_dmpar_bcast_reals(dminfo, nz1, qv_init)

      do i=1,nCells
         do k=1,nz1
            t_init(k,i) = t_init_1d(k)
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do
      end do

      do i=1,nEdges
        cell1 = cellsOnEdge(1,i)
        cell2 = cellsOnEdge(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)    
          end do
        end if
      end do


      !
      !  we are assuming w and rw are zero for this initialization
      !  i.e., no terrain
      !
       rw = 0.0
       w = 0.0

       zb = 0.0
       zb3 = 0.0

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,nEdges
         fEdge(iEdge) = 0.0
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 0.0
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, nVertLevels
                 v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
              end do
            end if
         end do
      end do

     ! write(0,*) ' k,u_init, t_init, qv_init '
     ! do k=1,nVertLevels
     !   write(0,'(i2,3(2x,f14.10)') k,u_init(k),t_initk),qv_init(k)
     ! end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_squall_line


!----------------------------------------------------------------------------------------------------------


   subroutine init_atm_case_mtn_wave(dminfo, mesh, nCells, nVertLevels, state, diag, configs, init_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(in) :: init_case

      real (kind=RKIND), parameter :: t0=288., hm=250.

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru 
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, deriv_two, zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, cellsOnCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2, nz1
      integer, pointer :: nEdges, maxEdges, nCellsSolve, nVertices
      integer, pointer :: index_qv

      real (kind=RKIND) :: ptop, pitop, ptopb, p0, flux, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND) :: ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(nVertLevels, nCells) :: relhum
      real (kind=RKIND) :: es, qvs, xnutr, ptemp
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2
      real (kind=RKIND) :: um, us,  rcp, rcv
      real (kind=RKIND) :: xmid, temp, pres, a_scale

      real (kind=RKIND) :: xi, xa, xc, xla, zinv, xn2, xn2m, xn2l, sm, dzh, dzht, dzmin, z_edge, z_edge3 

      integer, dimension(nCells, 2) :: next_cell
      real (kind=RKIND),  dimension(nCells) :: hxzt
      logical, parameter :: terrain_smooth = .false.

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_theta_adv_order

      real (kind=RKIND), pointer :: cf1, cf2, cf3

      real (kind=RKIND), dimension(:,:), pointer :: t_init, w, rw, v, rho, theta
      real (kind=RKIND), dimension(:), pointer :: u_init, angleEdge, fEdge, fVertex


      call mpas_pool_get_array(mesh, 'xCell', xCell)
      call mpas_pool_get_array(mesh, 'yCell', yCell)
      call mpas_pool_get_array(mesh, 'zCell', zCell)
      call mpas_pool_get_array(mesh, 'xEdge', xEdge)
      call mpas_pool_get_array(mesh, 'yEdge', yEdge)
      call mpas_pool_get_array(mesh, 'zEdge', zEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

      call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 't_init', t_init)
      call mpas_pool_get_array(mesh, 'u_init', u_init)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)

      !
      ! Scale all distances
      !
      a_scale = 1.0

      xCell(:) = xCell(:) * a_scale
      yCell(:) = yCell(:) * a_scale
      zCell(:) = zCell(:) * a_scale
      xVertex(:) = xVertex(:) * a_scale
      yVertex(:) = yVertex(:) * a_scale
      zVertex(:) = zVertex(:) * a_scale
      xEdge(:) = xEdge(:) * a_scale
      yEdge(:) = yEdge(:) * a_scale
      zEdge(:) = zEdge(:) * a_scale
      dvEdge(:) = dvEdge(:) * a_scale
      dcEdge(:) = dcEdge(:) * a_scale
      areaCell(:) = areaCell(:) * a_scale**2.0
      areaTriangle(:) = areaTriangle(:) * a_scale**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * a_scale**2.0

      
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz1 = nVertLevels
      nz = nz1 + 1

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'dss', dss)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'cqw', cqw)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(state, 'theta_m', t)
      call mpas_pool_get_array(state, 'u', u)
      call mpas_pool_get_array(state, 'w', w)
      call mpas_pool_get_array(state, 'scalars', scalars)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      xnutr = 0.1
      zd = 10500.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

      ! for hx computation
      xa = 5000. !SHP - should be changed based on grid distance 
      xla = 4000.
      xc = maxval (xCell(:))/2. 

      !     metrics for hybrid coordinate and vertical stretching
      str = 1.0
      zt = 21000.
      dz = zt/float(nz1)
!      write(0,*) ' dz = ',dz

      do k=1,nz

!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str 

!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) 
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!            write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      d1  = .5*dzw(1)
      d2  = dzw(1)+.5*dzw(2)
      d3  = dzw(1)+dzw(2)+.5*dzw(3)
      !cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      cof1 = (2.*dzu(2)+dzu(3))/(dzu(2)+dzu(3))*dzw(1)/dzu(2)
      cof2 =     dzu(2)        /(dzu(2)+dzu(3))*dzw(1)/dzu(3)
      cf1  = fzp(2) + cof1
      cf2  = fzm(2) - cof1 - cof2
      cf3  = cof2

! setting for terrain
      do iCell=1,nCells
         xi = xCell(iCell)
         !====1. for pure cosine mountain
         ! if(abs(xi-xc).ge.2.*xa)  then
         !    hx(1,iCell) = 0.
         ! else
         !    hx(1,iCell) = hm*cos(.5*pii*(xi-xc)/(2.*xa))**2.
         ! end if

         !====2. for cosine mountain
         !if(abs(xi-xc).lt.xa)  THEN
         !     hx(1,iCell) = hm*cos(pii*(xi-xc)/xla)**2. *cos(.5*pii*(xi-xc)/xa )**2.
         ! else
         !    hx(1,iCell) = 0.
         ! end if

         !====3. for shock mountain 
         hx(1,iCell) = hm*exp(-((xi-xc)/xa)**2)*cos(pii*(xi-xc)/xla)**2.

         hx(nz,iCell) = zt

!***** SHP -> get the temporary point information for the neighbor cell ->> should be changed!!!!! 
         do i=1,nCells 
            !option 1
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)-sqrt(3.)*dcEdge(1)) next_cell(iCell,1) = i 
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)+sqrt(3.)*dcEdge(1)) next_cell(iCell,2) = i 
            !option 2
            next_cell(iCell,1) = iCell - 8 ! note ny=4
            next_cell(iCell,2) = iCell + 8 ! note ny=4

            if (xCell(iCell).le. 3.*dcEdge(1)) then
                next_cell(iCell,1) = 1
            else if (xCell(iCell).ge. maxval(xCell(:))-3.*dcEdge(1)) then
                next_cell(iCell,2) = 1
            end if

         end do
      end do
      
      write(0,*) ' hx computation complete '

      if (terrain_smooth) then
         call mpas_dmpar_global_abort('***************************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Please contact the MPAS-A developers for up-to-date terrain-smoothing code.', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Otherwise, set terrain_smooth=.false. in the mountain wave test case',        deferredAbort=.true.)
         call mpas_dmpar_global_abort('   initialization routine and re-compile.',                                   deferredAbort=.true.)
         call mpas_dmpar_global_abort('***************************************************************************')
      end if

      do iCell=1,nCells
        do k=1,nz
            if (terrain_smooth) then
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)
            else
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(1,iCell)/zt)+hx(1,iCell)) &
                           + (1.-ah(k)) * zc(k)
            end if
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz1
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

      write(0,*) ' grid metrics setup complete '

!
! mountain wave initialization
!
         !SHP-original
         !zinv = 1000.
         !SHP-schar case
         zinv = 3000.

         xn2  = 0.0001
         xn2m = 0.0000
         xn2l = 0.0001

         um = 10.
         us = 0.

         do i=1,nCells
            do k=1,nz1
               ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
               tb(k,i) =  t0*(1. + xn2m/gravity*ztemp) 
               if(ztemp .le. zinv) then
                  t (k,i) = t0*(1.+xn2l/gravity*ztemp)
               else
                  t (k,i) = t0*(1.+xn2l/gravity*zinv+xn2/gravity*(ztemp-zinv)) 
               end if
                  relhum(k,i) = 0.
            end do
         end do

!  set the velocity field - we are on a plane here.

         do i=1, nEdges
            cell1 = cellsOnEdge(1,i)
            cell2 = cellsOnEdge(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               u(k,i) = um
               if(i == 1 ) u_init(k) = u(k,i) - us
#ifdef ROTATED_GRID
               u(k,i) = sin(angleEdge(i)) * (u(k,i) - us)
#else
               u(k,i) = cos(angleEdge(i)) * (u(k,i) - us)
#endif
            end do
            end if
         end do

!
!     reference sounding based on dry atmosphere
!
      pitop = 1.-.5*dzw(1)*gravity/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*(fzm(k)*tb(k,1)+fzp(k)*tb(k-1,1))   &
                                         *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
      end do
      pitop = pitop-.5*dzw(nz1)*gravity/(cp*tb(nz1,1)*zz(nz1,1))
      ptopb = p0*pitop**(1./rcp)
                
      do i=1, nCells
         pb(nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            cqw(k,i) = 1.
         end do
      end do

       write(0,*) ' ***** base state sounding ***** '
       write(0,*) 'k       pb        p         rb         rtb         rr          tb          t'
       do k=1,nVertLevels
          write(0,'(i2,7(2x,f14.9))') k,pb(k,1),p(k,1),rb(k,1),rtb(k,1),rr(k,1),tb(k,1),t(k,1)
       end do
 
       scalars(index_qv,:,:) = 0.

!-------------------------------------------------------------------
!     ITERATIONS TO CONVERGE MOIST SOUNDING
      do itr=1,30
        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))

        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*(fzm(k)*t (k,1)+fzp(k)*t (k-1,1)) &
                                                   *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)

        do i = 1, nCells

           pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
           do k=nz1-1,1,-1
              pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*                   &
                            (fzm(k)*(rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))  &
                            +fzp(k)*(rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i)))
           end do
           do k=1,nz1
              rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                      -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)
              p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
              rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
           end do
!
!     update water vapor mixing ratio from humitidty profile
!
           do k=1,nz1
              temp   = p(k,i)*t(k,i)
              pres   = p0*p(k,i)**(1./rcp)
              qvs    = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
              scalars(index_qv,k,i) = min(0.014_RKIND,relhum(k,i)*qvs)
           end do
                         
           do k=1,nz1
              t (k,i) = t(k,i)*(1.+1.61*scalars(index_qv,k,i))
           end do
           do k=2,nz1
              cqw(k,i) = 1./(1.+.5*( scalars(index_qv,k-1,i)  &
                                    +scalars(index_qv,k  ,i)))
           end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      write(0,*) ' *** sounding for the simulation ***'
      write(0,*) '    z       theta       pres         qv       rho_m        u        rr'
      do k=1,nz1
         write(0,'(8(f14.9,2x))') .5*(zgrid(k,1)+zgrid(k+1,1))/1000.,   &
                       t(k,1)/(1.+1.61*scalars(index_qv,k,1)),        &
                       .01*p0*p(k,1)**(1./rcp),                       &
                       1000.*scalars(index_qv,k,1),                   &
                       (rb(k,1)+rr(k,1))*(1.+scalars(index_qv,k,1)),  &
                       u_init(k), rr(k,1)
      end do

      do i=1,ncells
         do k=1,nz1
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do

        do k=1,nz1
            t_init(k,i) = t(k,i)
        end do
      end do

      do i=1,nEdges
        cell1 = cellsOnEdge(1,i)
        cell2 = cellsOnEdge(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)    
          end do
        end if
      end do

!
!     pre-calculation z-metric terms in omega eqn.
!
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do             
             
                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. 

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.   
                  else 
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1) 
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2) 
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1) 
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2) 
  
            end do

         end if
       end do

!     for including terrain
      w(:,:) = 0.0
      rw(:,:) = 0.0

!
!     calculation of omega, rw = zx * ru + zz * rw
!

      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))  
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux 
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux 

            if (config_theta_adv_order ==3) then
               rw(k,cell2) = rw(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do


      do iEdge=1,nEdges
         fEdge(iEdge) = 0.
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, nVertLevels
                 v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
              end do
            end if
         end do
      end do

!      do k=1,nVertLevels
!        write(0,*) ' k,u_init, t_init, qv_init ',k,u_init(k),t_init(k),qv_init(k)
!      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_mtn_wave


   subroutine init_atm_case_gfs(block, mesh, nCells, nEdges, nVertLevels, fg, state, diag, diag_physics, init_case, dims, configs)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Real-data test case using GFS data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_dmpar
      use init_atm_read_met
      use init_atm_llxy
      use init_atm_hinterp
      use mpas_hash
      use mpas_atmphys_constants, only : svpt0, svp1, svp2, svp3

      implicit none

      type (block_type), intent(inout), target :: block
      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: fg
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout):: diag_physics
      integer, intent(in) :: init_case
      type (mpas_pool_type), intent(inout):: dims
      type (mpas_pool_type), intent(inout):: configs

      type (parallel_info), pointer :: parinfo
      type (dm_info), pointer :: dminfo

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

!      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      type (met_data) :: field
      type (proj_info) :: proj

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: vert_level, latPoints, lonPoints, ter
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:), pointer :: destField1d
      real (kind=RKIND), dimension(:,:), pointer :: destField2d
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      real (kind=RKIND) :: target_z
      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
      integer, pointer :: nCellsSolve, nz1
      integer :: nInterpPoints, ndims

      integer :: nfglevels_actual
      integer, pointer :: index_qv

      integer, dimension(5) :: interp_list
      real (kind=RKIND) :: maskval
      real (kind=RKIND) :: msgval
      real (kind=RKIND) :: fillval
      integer :: masked

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, edgesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell 
      real (kind=RKIND), dimension(:,:), pointer :: v
      real (kind=RKIND), dimension(:,:), pointer :: sorted_arr

      type (field1DReal), pointer :: tempField
      type (field1DReal), pointer :: ter_field
      type (field1DReal), target :: tempFieldTarget

      real(kind=RKIND), dimension(:), pointer :: hs, hs1, sm0
      real(kind=RKIND) :: hm, hm_global, zh, dzmin, dzmina, dzmina_global, dzminf, dzminf_global, sm
      real(kind=RKIND) :: dcsum
      integer :: nsmterrain, kz, sfc_k
      logical :: hybrid, smooth

      integer :: it
      real (kind=RKIND) :: p_check

      ! For interpolating terrain and land use
      integer :: nx, ny
      integer :: istatus

      real (kind=RKIND), allocatable, dimension(:,:) :: rslab, maskslab
      integer, dimension(:), pointer :: mask_array
      integer, dimension(nEdges), target :: edge_mask
      character (len=StrKIND) :: fname
      logical :: is_sfc_field

      real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r
      real (kind=RKIND) :: lat, lon, x, y

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND), dimension(nVertLevels, nCells) :: rel_hum, temperature, qv
      real (kind=RKIND) :: ptmp, es, rs, rgas_moist, qvs, xnutr, znut, ptemp, rcv
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(nVertLevels + 1) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(nVertLevels + 1) :: sh, zw, ah
      real (kind=RKIND), dimension(nVertLevels) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(nVertLevels) :: eta, etav, teta, ppi, tt

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      integer, parameter :: nlat=361
      real (kind=RKIND), dimension(nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
      real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nlat, nVertLevels) :: u_2d, etavs_2d
      real (kind=RKIND), dimension(nVertLevels + 1) :: fsum
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND) :: alt, als, zetal, zl

      !  calculation of the water vapor mixing ratio:
      real (kind=RKIND) :: sh_max,sh_min,global_sh_max,global_sh_min

      character (len=StrKIND), pointer :: config_met_prefix
      character (len=StrKIND), pointer :: config_start_time
      logical, pointer :: config_met_interp
      logical, pointer :: config_vertical_grid
      integer, pointer :: config_nsmterrain
      integer, pointer :: config_nsm
      real (kind=RKIND), pointer :: config_dzmin
      real (kind=RKIND), pointer :: config_ztop
      logical, pointer :: config_tc_vertical_grid
      logical, pointer :: config_use_spechumd
      integer, pointer :: config_nfglevels
      integer, pointer :: config_nfgsoillevels
      logical, pointer :: config_smooth_surfaces
      integer, pointer :: config_theta_adv_order
      real (kind=RKIND), pointer :: config_coef_3rd_order

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), pointer :: cf1, cf2, cf3
      integer, dimension(:), pointer :: landmask

      real (kind=RKIND), dimension(:,:), pointer :: dzs_fg
      real (kind=RKIND), dimension(:,:), pointer :: zs_fg

      real (kind=RKIND), dimension(:), pointer :: sst
      real (kind=RKIND), dimension(:), pointer :: seaice
      real (kind=RKIND), dimension(:), pointer :: xice
      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: w
      real (kind=RKIND), dimension(:,:), pointer :: theta
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:), pointer :: relhum
      real (kind=RKIND), dimension(:,:), pointer :: spechum
      real (kind=RKIND), dimension(:,:), pointer :: ru
      real (kind=RKIND), dimension(:,:), pointer :: rw
      real (kind=RKIND), dimension(:), pointer :: precipw
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructY
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructZ
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructMeridional

      real (kind=RKIND), dimension(:), pointer :: psfc
      real (kind=RKIND), dimension(:), pointer :: skintemp
      real (kind=RKIND), dimension(:), pointer :: snow
      real (kind=RKIND), dimension(:), pointer :: snowc
      real (kind=RKIND), dimension(:,:), pointer :: u_fg
      real (kind=RKIND), dimension(:,:), pointer :: v_fg
      real (kind=RKIND), dimension(:,:), pointer :: z_fg
      real (kind=RKIND), dimension(:,:), pointer :: t_fg
      real (kind=RKIND), dimension(:,:), pointer :: rh_fg
      real (kind=RKIND), dimension(:,:), pointer :: sh_fg
      real (kind=RKIND), dimension(:,:), pointer :: gfs_z
      real (kind=RKIND), dimension(:,:), pointer :: p_fg
      real (kind=RKIND), dimension(:,:), pointer :: st_fg
      real (kind=RKIND), dimension(:,:), pointer :: sm_fg
      real (kind=RKIND), dimension(:), pointer :: soilz

      type (hashtable) :: level_hash
      logical :: too_many_fg_levs
      integer :: level_value

      ! For outputting surface fields u10, v10, q2, rh2, and t2m from first-guess data
      real (kind=RKIND), dimension(:), pointer :: u10
      real (kind=RKIND), dimension(:), pointer :: v10
      real (kind=RKIND), dimension(:), pointer :: q2
      real (kind=RKIND), dimension(:), pointer :: rh2
      real (kind=RKIND), dimension(:), pointer :: t2m

      character (len=StrKIND) :: errstring

      call mpas_pool_get_config(configs, 'config_met_prefix', config_met_prefix)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_met_interp', config_met_interp)
      call mpas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(configs, 'config_nsmterrain', config_nsmterrain)
      call mpas_pool_get_config(configs, 'config_nsm', config_nsm)
      call mpas_pool_get_config(configs, 'config_dzmin', config_dzmin)
      call mpas_pool_get_config(configs, 'config_ztop', config_ztop)
      call mpas_pool_get_config(configs, 'config_tc_vertical_grid', config_tc_vertical_grid)
      call mpas_pool_get_config(configs, 'config_use_spechumd', config_use_spechumd)
      call mpas_pool_get_config(configs, 'config_nfglevels', config_nfglevels)
      call mpas_pool_get_config(configs, 'config_nfgsoillevels', config_nfgsoillevels)
      call mpas_pool_get_config(configs, 'config_smooth_surfaces', config_smooth_surfaces)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)

      parinfo => block % parinfo
      dminfo => block % domain % dminfo

      call mpas_pool_get_field(mesh, 'ter', ter_field)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)

      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'ter', ter)
      call mpas_pool_get_array(mesh, 'dss', dss)

      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure', pressure)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'relhum', relhum)
      call mpas_pool_get_array(diag, 'spechum', spechum)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(state, 'theta_m', t)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(state, 'scalars', scalars)
      call mpas_pool_get_array(state, 'u', u)
      call mpas_pool_get_array(state, 'w', w)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag_physics, 'precipw', precipw)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)
      call mpas_pool_get_array(mesh, 'landmask', landmask)

      call mpas_pool_get_array(fg, 'dzs_fg', dzs_fg)
      call mpas_pool_get_array(fg, 'zs_fg', zs_fg)
      call mpas_pool_get_array(fg, 'sst', sst)
      call mpas_pool_get_array(fg, 'xice', xice)
      call mpas_pool_get_array(fg, 'seaice', seaice)
      call mpas_pool_get_array(fg, 'st_fg', st_fg)
      call mpas_pool_get_array(fg, 'sm_fg', sm_fg)
      call mpas_pool_get_array(fg, 'psfc', psfc)
      call mpas_pool_get_array(fg, 'skintemp', skintemp)
      call mpas_pool_get_array(fg, 'snow', snow)
      call mpas_pool_get_array(fg, 'snowc', snowc)
      call mpas_pool_get_array(fg, 'u', u_fg)
      call mpas_pool_get_array(fg, 'v', v_fg)
      call mpas_pool_get_array(fg, 'z', z_fg)
      call mpas_pool_get_array(fg, 't', t_fg)
      call mpas_pool_get_array(fg, 'rh', rh_fg)
      call mpas_pool_get_array(fg, 'sh', sh_fg)
      call mpas_pool_get_array(fg, 'gfs_z', gfs_z)
      call mpas_pool_get_array(fg, 'p', p_fg)

      call mpas_pool_get_dimension(dims, 'nVertLevels', nz1)
      call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
      nz = nz1 + 1

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      rcv = rgas/(cp-rgas)
      omega_e = omega
      p0 = 1.e+05

      scalars(:,:,:) = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN ADOPT GFS TERRAIN HEIGHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call mpas_dmpar_global_abort('********************************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Error opening initial meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                            deferredAbort=.true.)
         call mpas_dmpar_global_abort('********************************************************************************')
      end if

      call read_next_met_field(field, istatus)
      do while (istatus == 0)
         if (trim(field % field) == 'SOILHGT') then


write(0,*) 'USING ECMWF TERRAIN...'

            interp_list(1) = FOUR_POINT
            interp_list(2) = SEARCH
            interp_list(3) = 0

            !
            ! Set up projection
            !
            call map_init(proj)
          
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat,RKIND), &
                            loninc = real(field % deltalon,RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            end if


            if (trim(field % field) == 'SOILHGT') then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               destField1d => ter
               ndims = 1
            end if

            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               lat = latPoints(i)*DEG_PER_RAD
               lon = lonPoints(i)*DEG_PER_RAD
               call latlon_to_ij(proj, lat, lon, x, y)
               if (x < 0.5) then
                  lon = lon + 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               else if (x >= real(field%nx)+0.5) then
                  lon = lon - 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               end if
               if (y < 0.5) then
                  y = 1.0
               else if (y >= real(field%ny)+0.5) then
                  y = real(field%ny)
               end if
               if (ndims == 1) then
                  destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, -1.e30_RKIND, interp_list, 1)
               else if (ndims == 2) then
                  destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, -1.e30_RKIND, interp_list, 1)
               end if
            end do

            deallocate(rslab)
         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END ADOPT GFS TERRAIN HEIGHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (config_vertical_grid) then

      !
      ! Vertical grid setup
      !
      allocate(hs (nCells+1))
      allocate(hs1(nCells+1))

!     Fourth order smoother for terrain

      nsmterrain = config_nsmterrain

      do i=1,nsmterrain

         do iCell=1,nCells
            hs(iCell) = 0.
            if(ter(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)
                  hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                        / dcEdge(edgesOnCell(j,iCell))    &
                                        * (ter(cellsOnCell(j,iCell))-ter(iCell))
               end do
            end if
            hs(iCell) = ter(iCell) + 0.216*hs(iCell)
         end do

         do iCell=1,nCells
            ter(iCell) = 0.
            if(hs(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)
                  ter(iCell) = ter(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                          / dcEdge(edgesOnCell(j,iCell))    &
                                          * (hs(cellsOnCell(j,iCell))-hs(iCell))
               end do
            end if
            ter(iCell) = hs(iCell) - 0.216*ter(iCell)
         end do

         ! note that ther variable ter used throughout this section is a pointer to grid % ter % array, here we are passing ter's parent field
         call mpas_dmpar_exch_halo_field(ter_field)

      end do

      do iCell=1,nCells
         hx(:,iCell) = ter(iCell)
      end do

      hm = maxval(ter(1:nCellsSolve))
      call mpas_dmpar_max_real(dminfo, hm, hm_global)
      hm = hm_global
      write(0,*) "max ter = ", hm

!     Metrics for hybrid coordinate and vertical stretching

      if (config_tc_vertical_grid) then

         write(0,*) 'Setting up vertical levels as in 2014 TC experiments'

         zt = config_ztop
         dz = zt/float(nz1)

!... Laura D. Fowler: change the values als,alt,and zetal valid for nVertLevels=41 to values
!    needed for nVertLevels equal or greater than 55.
         if (nVertLevels >= 55) then
            als   =  .075  ! 55 levels.
            alt   = 1.70   ! 55 levels.
            zetal =  .75   ! 55 levels.
         else
            als   =  .075  ! 41 levels.
            alt   = 1.23   ! 41 levels.
            zetal =  .31   ! 41 levels.
         endif
         write(0,*)
         write(0,*) '--- config_tc_vertical_grid = ', config_tc_vertical_grid
         write(0,*) '--- als   = ', als
         write(0,*) '--- alt   = ', alt
         write(0,*) '--- zetal = ',zetal
         if (nVertLevels /= 55 .and. nVertLevels /= 41) then
            write(0,*) '********************************************************************************'
            write(0,*) '* Note: als, alt, and zetal have not been tested for nVertLevels equal ', nVertLevels
            write(0,*) '*       The values of these parameters have only been tested with 41 and 55 layers'
            write(0,*) '********************************************************************************'
         end if
         write(0,*) ' '
!... end Laura D. Fowler / 2016-04-12.

         do k=1,nz
            zl = 1.-alt*(1.-zetal)
            if ((real(k-1)/real(nz1)).LT.zetal) then
                 zw(k) =              (als*real(K-1)/real(nz1)      &
                     +(3.*(1.-alt)+2.*(alt-als)*zetal)              &
                                     *((K-1)*dz/(zt*zetal))**2      &
                     -(2.*(1.-alt) +  (alt-als)*zetal)              &
                                     *(real(k-1)*dz/(zt*zetal))**3) * zt
            else
                 zw(K) = (zl+alt*(real(K-1)/real(nz1)-zetal))*zt
            end if
            if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
         end do

      else

         write(0,*) 'Setting up vertical levels as in MPAS 2.0 and earlier'

         str = 1.5
         zt = config_ztop
         dz = zt/float(nz1)

         do k=1,nz
            zw(k) = (real(k-1)/real(nz1))**str*zt
            if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
         end do
      end if

!     ah(k) governs the transition between terrain-following 
!        and pure height coordinates
!           ah(k) = 1           is a smoothed terrain-following coordinate
!           ah(k) = 1.-zw(k)/zt is the basic terrain-following coordinate
!           ah(k) = 0           is a height coordinate
 
      hybrid = .true.
!      hybrid = .false.

      kz = nz
      if (hybrid) then
      
         zh = 30000.0
!         zh = 0.5*zt

         do k=1,nz
            if (zw(k) < zh) then
               ah(k) = cos(.5*pii*zw(k)/zh)**6

!!!               ah(k) = ah(k)*(1.-zw(k)/zt)

            else
               ah(k) = 0.
               kz = min(kz,k)
            end if
         end do

      else

         do k=1,nz
            ah(k) = 1.-zw(k)/zt
         end do

      end if

      do k=1,nz
         write(0,*) k,zw(k), ah(k)
      end do

      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))


      write(0,*) ' cf1, cf2, cf3 = ',cf1,cf2,cf3

!     Smoothing algorithm for coordinate surfaces 

      smooth = config_smooth_surfaces

      if (smooth) then

         dzmin = config_dzmin

         allocate(sm0(nCells+1))

         do iCell=1,nCells
            dcsum = 0.0
            do j=1,nEdgesOnCell(iCell)
               dcsum = dcsum + dcEdge(edgesOnCell(j,iCell))
            end do
            dcsum = dcsum / real(nEdgesOnCell(iCell))
            sm0(iCell) = max(0.01_RKIND, 0.125 * min(1.0_RKIND, 3000.0_RKIND/dcsum))
         end do

         do k=2,kz-1
            hx(k,:) = hx(k-1,:)
            dzminf = zw(k)-zw(k-1)

            do i=1,config_nsm + k
               do iCell=1,nCells

                  sm = sm0(iCell) * min((3.0*zw(k)/zt)**2.0, 1.0_RKIND)

                  hs1(iCell) = 0.
                  do j = 1,nEdgesOnCell(iCell)

                     hs1(iCell) = hs1(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                           / dcEdge(edgesOnCell(j,iCell))    &
                                           *  (hx(k,cellsOnCell(j,iCell))-hx(k,iCell))
                  end do
                  hs1(iCell) = hx(k,iCell) + sm*hs1(iCell)

                  hs(iCell) = 0.
              !    do j = 1,nEdgesOnCell(iCell)
              !       hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
              !                             / dcEdge(edgesOnCell(j,iCell))    &
              !                             *  (hs1(cellsOnCell(j,iCell))-hs1(iCell))
              !    end do
                  hs(iCell) = hs1(iCell) - 0.*hs(iCell)

               end do

               tempField => tempFieldTarget
               tempField % block => block
               tempField % dimSizes(1) = nCells
               tempField % sendList => parinfo % cellsToSend
               tempField % recvList => parinfo % cellsToRecv
               tempField % copyList => parinfo % cellsToCopy
               tempField % array => hs
               tempField % isActive = .true.
               tempField % prev => null()
               tempField % next => null()

               call mpas_dmpar_exch_halo_field(tempField)

               do iCell=1,nCells
                  dzmina = (zw(k) + ah(k)*hs(iCell)) - (zw(k-1) + ah(k-1)*hx(k-1,iCell))
                  if (dzmina > dzmin*(zw(k)-zw(k-1))) then
                     hx(k,iCell) = hs(iCell)
                     if (dzmina < dzminf) dzminf = dzmina
                  end if 
               end do

            end do
            call mpas_dmpar_min_real(dminfo, dzminf, dzminf_global)
            write(0,*) k,i,sm,dzminf_global/(zw(k)-zw(k-1))
         end do

         deallocate(sm0)

         do k=kz,nz
               hx(k,:) = 0.
         end do

      else

         do k=2,nz1
            dzmina = minval(zw(k)+ah(k)*hx(k,:)-zw(k-1)-ah(k-1)*hx(k-1,:))
            write(0,*) k,dzmina/(zw(k)-zw(k-1))
         end do

      end if

      deallocate(hs )
      deallocate(hs1)

!     Height of coordinate levels (calculation of zgrid)

      do iCell=1,nCells
         do k=1,nz
            zgrid(k,iCell) = zw(k) + ah(k)*hx(k,iCell)
         end do
         do k=1,nz1
            zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
         end do
      end do

      do i=1, nEdges
         iCell1 = cellsOnEdge(1,i)
         iCell2 = cellsOnEdge(2,i)
         do k=1,nz1
            zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
         end do
      end do
      do i=1, nCells
         do k=1,nz1
           ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
           dss(k,i) = 0.
           ztemp = zgrid(k,i)
           if (ztemp.gt.zd+.1)  then
               dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
           end if
         end do
      end do


      ! For z-metric term in omega equation
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                        d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                        d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do             
             
                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. 

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.   
                  else 
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1) 
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2) 
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1) 
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2) 
  
            end do

         end if
      end do

      write(0,*) ' grid metrics setup complete '

      end if    ! config_vertical_grid


      if (config_met_interp) then

      !ldf (2011-11-19): added initialization of the sea-surface temperature, seaice fraction, and
      !seaice flag:
       sst = 0.0
       xice = 0.0
       seaice = 0.0
      !ldf end.

      !
      ! First, try to locate the LANDSEA field for use as an interpolation mask
      !
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call mpas_dmpar_global_abort('********************************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Error opening initial meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                            deferredAbort=.true.)
         call mpas_dmpar_global_abort('********************************************************************************')
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)
         if (trim(field % field) == 'LANDSEA') then

            allocate(maskslab(-2:field % nx+3, field % ny))
            maskslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            maskslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            maskslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            maskslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            maskslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            maskslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            maskslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
write(0,*) 'minval, maxval of LANDSEA = ', minval(maskslab), maxval(maskslab)

         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()

      if (.not. allocated(maskslab)) then
         call mpas_dmpar_global_abort('********************************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('LANDSEA field not found in meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                            deferredAbort=.true.)
         call mpas_dmpar_global_abort('********************************************************************************')
      end if

      edge_mask(:) = 1


      !
      ! Horizontally interpolate meteorological data
      !
      allocate(vert_level(config_nfglevels))
      vert_level(:) = -1.0

      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         call mpas_dmpar_global_abort('********************************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Error opening initial meteorological data file ' &
                                      //trim(config_met_prefix)//':'//config_start_time(1:13),                            deferredAbort=.true.)
         call mpas_dmpar_global_abort('********************************************************************************')
      end if

      call mpas_hash_init(level_hash)
      too_many_fg_levs = .false.

      call read_next_met_field(field, istatus)

      do while (istatus == 0)

         interp_list(1) = FOUR_POINT
         interp_list(2) = SEARCH
         interp_list(3) = 0

         maskval = -1.0
         masked = -1
         fillval = 0.0
         msgval = -1.e30

         mask_array => landmask

         if (trim(field % field) == 'UU' .or. &
             trim(field % field) == 'VV' .or. &
             trim(field % field) == 'TT' .or. &
             trim(field % field) == 'RH' .or. &
             trim(field % field) == 'SPECHUMD' .or. &
             trim(field % field) == 'GHT' .or. &
             trim(field % field) == 'PMSL' .or. &
             trim(field % field) == 'PSFC' .or. &
             trim(field % field) == 'SOILHGT' .or. &
             trim(field % field) == 'SM000010' .or. &
             trim(field % field) == 'SM010040' .or. &
             trim(field % field) == 'SM040100' .or. &
             trim(field % field) == 'SM100200' .or. &
             trim(field % field) == 'SM010200' .or. &
             trim(field % field) == 'SM000007' .or. &
             trim(field % field) == 'SM007028' .or. &
             trim(field % field) == 'SM028100' .or. &
             trim(field % field) == 'SM100255' .or. &
             trim(field % field) == 'ST000010' .or. &
             trim(field % field) == 'ST010040' .or. &
             trim(field % field) == 'ST040100' .or. &
             trim(field % field) == 'ST100200' .or. &
             trim(field % field) == 'ST010200' .or. &
             trim(field % field) == 'ST000007' .or. &
             trim(field % field) == 'ST007028' .or. &
             trim(field % field) == 'ST028100' .or. &
             trim(field % field) == 'ST100255' .or. &
             trim(field % field) == 'PRES' .or. &
             trim(field % field) == 'PRESSURE' .or. &
             trim(field % field) == 'SNOW' .or. &
             trim(field % field) == 'SEAICE' .or. &
             trim(field % field) == 'SKINTEMP') then

            if (trim(field % field) == 'SM000010' .or. &
                trim(field % field) == 'SM010040' .or. &
                trim(field % field) == 'SM040100' .or. &
                trim(field % field) == 'SM100200' .or. &
                trim(field % field) == 'SM010200' .or. &
                trim(field % field) == 'SM000007' .or. &
                trim(field % field) == 'SM007028' .or. &
                trim(field % field) == 'SM028100' .or. &
                trim(field % field) == 'SM100255' .or. &
                trim(field % field) == 'ST000010' .or. &
                trim(field % field) == 'ST010040' .or. &
                trim(field % field) == 'ST040100' .or. &
                trim(field % field) == 'ST100200' .or. &
                trim(field % field) == 'ST010200' .or. &
                trim(field % field) == 'ST000007' .or. &
                trim(field % field) == 'ST007028' .or. &
                trim(field % field) == 'ST028100' .or. &
                trim(field % field) == 'ST100255' .or. &
                trim(field % field) == 'SNOW' .or. &
                trim(field % field) == 'SEAICE' .or. &
                trim(field % field) == 'SKINTEMP') then
               k = 1
            else if (trim(field % field) /= 'PMSL' .and. &
                     trim(field % field) /= 'PSFC' .and. &
                     trim(field % field) /= 'SOILHGT') then

               ! Since the hash table can only store integers, transfer the bit pattern from 
               ! the real-valued xlvl into an integer; that the result is not an integer version
               ! of the level is not important, since we only want to test uniqueness of levels
               level_value = transfer(field % xlvl, level_value)
               if (.not. mpas_hash_search(level_hash, level_value)) then
                  call mpas_hash_insert(level_hash, level_value)
                  if (mpas_hash_size(level_hash) > config_nfglevels) then
                     too_many_fg_levs = .true.
                  end if
               end if

               !
               ! In case we have more than config_nfglevels levels, just keep cycling through
               ! the remaining fields in the intermediate file for the purpose of counting how
               ! many unique levels are found using the code above
               !
               if (too_many_fg_levs) then
                  call read_next_met_field(field, istatus)
                  cycle
               end if

               do k=1,config_nfglevels
                  if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
               end do
               if (vert_level(k) == -1.0) vert_level(k) = field % xlvl
            else
               k = 1
            end if

            !
            ! Set up projection
            !
            call map_init(proj)
          
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat,RKIND), &
                            loninc = real(field % deltalon,RKIND), &
                            knowni = 1.0_RKIND, &
                            knownj = 1.0_RKIND, &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = 360.0_RKIND / real(field % nx,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
!                            nxmax = nint(360.0 / field % deltalon), &
            end if


            !
            ! Horizontally interpolate the field at level k
            !
            if (trim(field % field) == 'UU') then
write(0,*) 'Interpolating U at ', k, vert_level(k)

               ! For U10, interpolate to cell centers
               if (vert_level(k) == 200100.0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call mpas_pool_get_array(fg, 'u10', destField1d)
                  ndims = 1

               ! otherwise to edges
               else
                  mask_array => edge_mask

                  nInterpPoints = nEdges
                  latPoints => latEdge
                  lonPoints => lonEdge
                  call mpas_pool_get_array(fg, 'u', destField2d)
                  ndims = 2
               end if

            else if (trim(field % field) == 'VV') then
write(0,*) 'Interpolating V at ', k, vert_level(k)

               ! For V10, interpolate to cell centers
               if (vert_level(k) == 200100.0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call mpas_pool_get_array(fg, 'v10', destField1d)
                  ndims = 1

               ! otherwise to edges
               else
                  mask_array => edge_mask

                  nInterpPoints = nEdges
                  latPoints => latEdge
                  lonPoints => lonEdge
                  call mpas_pool_get_array(fg, 'v', destField2d)
                  ndims = 2
               end if

            else if (trim(field % field) == 'TT') then
write(0,*) 'Interpolating T at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 't', destField2d)
               ndims = 2
            else if (trim(field % field) == 'RH') then
write(0,*) 'Interpolating RH at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'rh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'SPECHUMD') then
write(0,*) 'Interpolating SPECHUMD at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sh', destField2d)
               ndims = 2
            else if (trim(field % field) == 'GHT') then
write(0,*) 'Interpolating GHT at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'z', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRES') then
write(0,*) 'Interpolating PRES at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PRESSURE') then
write(0,*) 'Interpolating PRESSURE at ', k, vert_level(k)
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'p', destField2d)
               ndims = 2
            else if (trim(field % field) == 'PMSL') then
write(0,*) 'Interpolating PMSL'
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'pmsl', destField1d)
               ndims = 1
            else if (trim(field % field) == 'PSFC') then
write(0,*) 'Interpolating PSFC'
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'psfc', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SOILHGT') then
write(0,*) 'Interpolating SOILHGT'
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'soilz', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SM000010') then
write(0,*) 'Interpolating SM000010'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 10.
               zs_fg(k,:) = 10.
            else if (trim(field % field) == 'SM010200') then
write(0,*) 'Interpolating SM010200'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 200.-10.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'SM010040') then
write(0,*) 'Interpolating SM010040'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 40.-10.
               zs_fg(k,:) = 40.
            else if (trim(field % field) == 'SM040100') then
write(0,*) 'Interpolating SM040100'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-40.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'SM100200') then
write(0,*) 'Interpolating SM100200'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 200.-100.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'SM000007') then
write(0,*) 'Interpolating SM000007'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 7.
               zs_fg(k,:) = 7.
            else if (trim(field % field) == 'SM007028') then
write(0,*) 'Interpolating SM007028'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 28.-7.
               zs_fg(k,:) = 28.
            else if (trim(field % field) == 'SM028100') then
write(0,*) 'Interpolating SM028100'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-28.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'SM100255') then
write(0,*) 'Interpolating SM100255'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'sm_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 255.-100.
               zs_fg(k,:) = 255.
            else if (trim(field % field) == 'ST000010') then
write(0,*) 'Interpolating ST000010'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 10.
               zs_fg(k,:) = 10.
            else if (trim(field % field) == 'ST010200') then
write(0,*) 'Interpolating ST010200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 200.-10.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'ST010040') then
write(0,*) 'Interpolating ST010040'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 40.-10.               
               zs_fg(k,:) = 40.
            else if (trim(field % field) == 'ST040100') then
write(0,*) 'Interpolating ST040100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-40.             
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'ST100200') then
write(0,*) 'Interpolating ST100200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 200.-100.
               zs_fg(k,:) = 200.
            else if (trim(field % field) == 'ST000007') then
write(0,*) 'Interpolating ST000007'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 1
               ndims = 2
               dzs_fg(k,:) = 7.
               zs_fg(k,:) = 7.
            else if (trim(field % field) == 'ST007028') then
write(0,*) 'Interpolating ST007028'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 2
               ndims = 2
               dzs_fg(k,:) = 28.-7.
               zs_fg(k,:) = 28.
            else if (trim(field % field) == 'ST028100') then
write(0,*) 'Interpolating ST028100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 3
               ndims = 2
               dzs_fg(k,:) = 100.-28.
               zs_fg(k,:) = 100.
            else if (trim(field % field) == 'ST100255') then
write(0,*) 'Interpolating ST100255'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'st_fg', destField2d)
               k = 4
               ndims = 2
               dzs_fg(k,:) = 255.-100.
               zs_fg(k,:) = 255.
            else if (trim(field % field) == 'SNOW') then
write(0,*) 'Interpolating SNOW'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0

               masked = 0
               fillval = 0.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'snow', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SEAICE') then
write(0,*) 'Interpolating SEAICE'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0

               maskval = 1.0
               masked = 1
               fillval = 0.0

               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'xice', destField1d)
               ndims = 1
            else if (trim(field % field) == 'SKINTEMP') then
write(0,*) 'Interpolating SKINTEMP'
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'skintemp', destField1d)
               ndims = 1
            end if

            allocate(rslab(-2:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               if (mask_array(i) /= masked) then
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field%ny)
                  end if
                  if (ndims == 1) then
                     if (maskval /= -1.0) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  else if (ndims == 2) then
                     if (maskval /= -1.0) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                     else
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1)
                     end if
                  end if
               else
                  if (ndims == 1) then
                     destField1d(i) = fillval
                  else if (ndims == 2) then
                     destField2d(k,i) = fillval
                  end if
               end if
            end do

            !
            ! In addition to interpolating wind fields to cell edges, we should
            ! also intperolate to cell centers at the surface in order to
            ! produce U10 and V10 fields
            !
            is_sfc_field = .false.
            mask_array => landmask
            if (index(field % field, 'UU') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'u10', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'VV') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'v10', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'TT') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 't2m', destField1d)
               ndims = 1
               is_sfc_field = .true.
            else if (index(field % field, 'RH') /= 0 .and. vert_level(k) == 200100.0) then
               nInterpPoints = nCells
               latPoints => latCell
               lonPoints => lonCell
               call mpas_pool_get_array(fg, 'rh2', destField1d)
               ndims = 1
               is_sfc_field = .true.
            end if

            if (is_sfc_field) then
               do i=1,nInterpPoints
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, &
                                                   msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
               end do
            end if

            deallocate(rslab)
     
         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()
      level_value = mpas_hash_size(level_hash)
      call mpas_hash_destroy(level_hash)

      if (too_many_fg_levs) then
         write(errstring,'(a,i4)') '       Please increase config_nfglevels to at least ', level_value
         call mpas_dmpar_global_abort('*******************************************************************', deferredAbort=.true.)
         call mpas_dmpar_global_abort('Error: The meteorological data file has more than config_nfglevels.', deferredAbort=.true.)
         call mpas_dmpar_global_abort(trim(errstring),                                                       deferredAbort=.true.)
         call mpas_dmpar_global_abort('       in the namelist and re-run.',                                  deferredAbort=.true.)
         call mpas_dmpar_global_abort('*******************************************************************')
      end if


      !
      ! Check how many distinct levels we actually found in the meteorological data
      !
      do k=1,config_nfglevels
         if (vert_level(k) == -1.0) exit 
      end do
      nfglevels_actual = k-1
      write(0,*) '*************************************************'
      write(0,*) 'Found ', nfglevels_actual, ' levels in the first-guess data'
      write(0,*) '*************************************************'

   
      !
      ! Extract surface fields from first-guess
      !
!      call mpas_pool_get_array(fg, 'u10', u10)
!      call mpas_pool_get_array(fg, 'v10', v10)
      call mpas_pool_get_array(fg, 'q2', q2)
      call mpas_pool_get_array(fg, 'rh2', rh2)
      call mpas_pool_get_array(fg, 't2m', t2m)
      call mpas_pool_get_array(fg, 'u', u_fg)
      call mpas_pool_get_array(fg, 'v', v_fg)
      call mpas_pool_get_array(fg, 't', t_fg)
      call mpas_pool_get_array(fg, 'rh', rh_fg)
!      u10(:) = 0.0
!      v10(:) = 0.0
      q2(:) = 0.0
      rh2(:) = 0.0
      t2m(:) = 0.0

      do k=1,config_nfglevels
         if (vert_level(k) == 200100.0) then
!            u10(:) = u_fg(k,:)
!            v10(:) = v_fg(k,:)
            t2m(:) = t_fg(k,:)
            rh2(:) = rh_fg(k,:)
            do iCell = 1, nCells
               es = svp1 * 10.0_RKIND * exp(svp2 * (t2m(iCell)-svpt0) / (t2m(iCell)-svp3))
               es = min(es, 0.99_RKIND * 0.01_RKIND * psfc(iCell))
               rs = 0.622_RKIND * es * 100.0_RKIND / (psfc(iCell) - es * 100.0_RKIND)
               q2(iCell) = 0.01_RKIND * rs * rh2(iCell)
             end do
         end if
      end do



      ! 
      ! For isobaric data, fill in the 3-d pressure field; otherwise, ensure
      ! that the surface pressure and height fields are filled in
      ! 
      if (minval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0 .and. &
          maxval(p_fg(1:nfglevels_actual,1:nCellsSolve)) == 0.0) then
         write(0,*) 'Setting pressure field for isobaric data'
         do k=1,config_nfglevels
            if (vert_level(k) /= 200100.0) then
               p_fg(k,:) = vert_level(k)
            else
               p_fg(k,:) = psfc(:)
            end if
         end do
      else
         call mpas_pool_get_array(fg, 'z', z_fg)
         call mpas_pool_get_array(fg, 'soilz', soilz)
         write(0,*) 'Assuming model-level input data'
         do k=1,config_nfglevels
            if (vert_level(k) == 200100.0) then
               p_fg(k,:) = psfc(:)
               z_fg(k,:) = soilz(:)
            end if
         end do
      end if

      ! Set SST based on SKINTEMP field if it wasn't found in input data
      if (minval(sst) == 0.0 .and. maxval(sst) == 0.0) then
         write(0,*) 'Setting SST from SKINTEMP'
         !where (landmask == 0) sst = skintemp
         sst = skintemp
      end if

      ! Set SNOWC (snow-cover flag) based on SNOW
      snowc(:) = 0.0
      where (snow > 0.0) snowc = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MGD CHECK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iCell=1,nCells
   if (landmask(iCell) == 1) then

      do k = 1, config_nfgsoillevels
         if (st_fg(k,iCell) <= 0.0) write(0,*) 'Bad st_fg ', k, iCell
      end do

      do k = 1, config_nfgsoillevels
         if (sm_fg(k,iCell) <= 0.0) then
            write(0,*) 'Bad sm_fg ', sm_fg(k,iCell), k, iCell
            sm_fg(k,iCell) = 0.001
         end if
      end do
      !LDF end.

   end if
end do
write(0,*) 'Done with soil consistency check'


      !
      ! Get SEAICE from a separate file
      !
      call read_met_init('SEAICE_FRACTIONAL', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'SEAICE_FRACTIONAL file not found...'
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (trim(field % field) == 'SEAICE') then

write(0,*) 'PROCESSING SEAICE'

               !
               ! Set up projection
               !
               call map_init(proj)
          
               if (field % iproj == PROJ_PS) then
                  call map_set(PROJ_PS, proj, &
                               dx = real(field % dx,RKIND), &
                               truelat1 = real(field % truelat1,RKIND), &
                               stdlon = real(field % xlonc,RKIND), &
                               knowni = real(field % nx / 2.0,RKIND), &
                               knownj = real(field % ny / 2.0,RKIND), &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               end if

               if (trim(field % field) == 'SEAICE') then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call mpas_pool_get_array(fg, 'xice', destField1d)
                  ndims = 1
               end if
   
               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0
   
               masked = 1
               fillval = 0.0
               msgval = 1.01
               mask_array => landmask


               allocate(rslab(field % nx, field % ny))
               rslab(:,:) = field % slab(:,:)
               do i=1,nInterpPoints
                  if (mask_array(i) /= masked) then
                     lat = latPoints(i)*DEG_PER_RAD
                     lon = lonPoints(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     else if (x >= real(field%nx)+0.5) then
                        lon = lon - 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (y < 0.5) then
                        y = 1.0
                     else if (y >= real(field%ny)+0.5) then
                        y = real(field%ny)
                     end if
                     if (ndims == 1) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField1d(i) == msgval) destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField2d(k,i) == msgval) destField2d(k,i) = fillval
                     end if
                  else
                     if (ndims == 1) then
                        destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = fillval
                     end if
                  end if
               end do
               deallocate(rslab)

            end if
      
            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()


      !
      ! Get OMLD climatology from a separate file
      !
      call read_met_init('OMLD', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'OMLD file not found...'
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (index(field % field, 'OMLD') /= 0) then

write(0,*) 'PROCESSING OMLD'

               !
               ! Set up projection
               !
               call map_init(proj)
          
               if (field % iproj == PROJ_LATLON) then
                  call map_set(PROJ_LATLON, proj, &
                               latinc = real(field % deltalat,RKIND), &
                               loninc = real(field % deltalon,RKIND), &
                               knowni = 1.0_RKIND, &
                               knownj = 1.0_RKIND, &
                               lat1 = real(field % startlat,RKIND), &
                               lon1 = real(field % startlon,RKIND))
               else
                  write(stderrUnit,*) 'ERROR: We were expecting OMLD field to be on a lat-lon projection...'
               end if

               if (index(field % field, 'OMLD') /= 0) then
                  nInterpPoints = nCells
                  latPoints => latCell
                  lonPoints => lonCell
                  call mpas_pool_get_array(state, 'h_oml_initial', destField1d)
                  ndims = 1
               end if
   
               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = SEARCH
               interp_list(4) = 0
   
               masked = 1
               fillval = 0.0
               msgval = 0.0
               mask_array => landmask


               allocate(rslab(field % nx, field % ny))
               rslab(:,:) = field % slab(:,:)
               do i=1,nInterpPoints
                  if (mask_array(i) /= masked) then
                     lat = latPoints(i)*DEG_PER_RAD
                     lon = lonPoints(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     else if (x >= real(field%nx)+0.5) then
                        lon = lon - 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (y < 0.5) then
                        y = 1.0
                     else if (y >= real(field%ny)+0.5) then
                        y = real(field%ny)
                     end if
                     if (ndims == 1) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField1d(i) == msgval) destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField2d(k,i) == msgval) destField2d(k,i) = fillval
                     end if
                  else
                     if (ndims == 1) then
                        destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = fillval
                     end if
                  end if
               end do
               deallocate(rslab)

            end if
      
            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()


      if (allocated(maskslab)) deallocate(maskslab)

      ! Freeze really cold ocean
      where (sst < 271.0 .and. landmask == 0) xice = 1.0

      ! Limit XICE to values between 0 and 1. Although the input meteorological field is between 0.
      ! and 1., interpolation to the MPAS grid can yield values of XiCE less than 0. and greater
      ! than 1.:
      where (xice < 0._RKIND) xice = 0._RKIND
      where (xice > 1._RKIND) xice = 1._RKIND

      ! Set SEAICE (0/1 flag) based on XICE (fractional ice coverage)
      seaice(:) = 0.0
      where (xice >= 0.5) seaice = 1.0


      !  
      ! Compute normal wind component and store in fg % u
      !  
      do iEdge=1,nEdges
         do k=1,nfglevels_actual
            u_fg(k,iEdge) = cos(angleEdge(iEdge)) * u_fg(k,iEdge) &
                          + sin(angleEdge(iEdge)) * v_fg(k,iEdge)
         end do
      end do

      !  
      ! Vertically interpolate meteorological data
      !  
      allocate(sorted_arr(2,nfglevels_actual))

      do iCell=1,nCells

         ! T
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = t_fg(k,iCell)
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            t(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                  sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
         end do


         ! RH
         sorted_arr(:,:) = -999.0
         relhum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = rh_fg(k,iCell)
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            relhum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                       sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
            if (target_z < z_fg(1,iCell) .and. k < nVertLevels) relhum(k,iCell) = relhum(k+1,iCell)
         end do


         ! SPECHUM: if first-guess values are negative, set those values to zero before
         ! vertical interpolation.
         sorted_arr(:,:) = -999.0
         spechum(:,iCell) = 0._RKIND
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = max(0._RKIND,sh_fg(k,iCell))
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k = nVertLevels, 1, -1
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            spechum(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                        sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=0)
            if (target_z < z_fg(1,iCell) .and. k < nVertLevels) spechum(k,iCell) = spechum(k+1,iCell)
         end do


         ! GHT
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = z_fg(k,iCell)
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            gfs_z(k,iCell) = vertical_interp(target_z, nfglevels_actual-1, &
                                      sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
         end do


         ! PRESSURE
         sorted_arr(:,:) = -999.0
         do k = 1, nfglevels_actual
            sorted_arr(1,k) = z_fg(k,iCell)
            if (vert_level(k) == 200100.0) then 
               sorted_arr(1,k) = 99999.0
               sfc_k = k
            end if
            sorted_arr(2,k) = log(p_fg(k,iCell))
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k = 1, nVertLevels
            target_z = 0.5 * (zgrid(k,iCell) + zgrid(k+1,iCell))
            pressure(k,iCell) = exp(vertical_interp(target_z, nfglevels_actual-1, &
                                    sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1))
         end do

      end do


      do iEdge=1,nEdges

         ! U
         sorted_arr(:,:) = -999.0
         do k=1,nfglevels_actual
            sorted_arr(1,k) = 0.5 * (z_fg(k,cellsOnEdge(1,iEdge)) + z_fg(k,cellsOnEdge(2,iEdge)))
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 0.5 * (fg % soilz % array(cellsOnEdge(1,iEdge)) + fg % soilz % array(cellsOnEdge(2,iEdge)))
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = u_fg(k,iEdge)
         end do
         call mpas_quicksort(nfglevels_actual, sorted_arr)
         do k=1,nVertLevels
            target_z = 0.25 * (zgrid(k,cellsOnEdge(1,iEdge)) + zgrid(k+1,cellsOnEdge(1,iEdge)) + zgrid(k,cellsOnEdge(2,iEdge)) + zgrid(k+1,cellsOnEdge(2,iEdge)))
!           u(k,iEdge) = vertical_interp(target_z, nfglevels_actual, sorted_arr, order=1, extrap=0)
            u(k,iEdge) = vertical_interp(target_z, nfglevels_actual-1, sorted_arr(:,1:nfglevels_actual-1), order=1, extrap=1)
         end do

      end do


      !
      ! Reconstruct zonal and meridional winds for diagnostic puposes:
      !
      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)
      call mpas_reconstruct(mesh, u,                 &
                            uReconstructX,           &
                            uReconstructY,           &
                            uReconstructZ,           &
                            uReconstructZonal,       &
                            uReconstructMeridional   &
                           )
   

      !
      ! Adjust surface pressure for difference in topography
      !
      do sfc_k=1,nfglevels_actual
         if (vert_level(sfc_k) == 200100.) exit
      end do 
      do iCell=1,nCells

         ! We need to extrapolate
            sorted_arr(:,:) = -999.0
            do k=1,nfglevels_actual
               sorted_arr(1,k) = z_fg(k,iCell)
               if (vert_level(k) == 200100.0) then 
!NOSFC                  sorted_arr(1,k) = fg % soilz % array(iCell)
                  sorted_arr(1,k) = 99999.0
               end if
               sorted_arr(2,k) = log(p_fg(k,iCell))
            end do
            call mpas_quicksort(nfglevels_actual, sorted_arr)
            target_z = zgrid(1,iCell)
            psfc(iCell) = exp(vertical_interp(target_z, nfglevels_actual, sorted_arr, order=1, extrap=1))

      end do

      deallocate(sorted_arr)


      ! Diagnose the water vapor mixing ratios:
      global_sh_min = 0._RKIND
      global_sh_max = 0._RKIND
      if(config_use_spechumd) then
         sh_min = minval(spechum(:,1:nCellsSolve))
         sh_max = maxval(spechum(:,1:nCellsSolve))
         call mpas_dmpar_min_real(dminfo,sh_min,global_sh_min)
         call mpas_dmpar_max_real(dminfo,sh_max,global_sh_max)
      endif
      write(0,*)
      write(0,*) '--- global_sh_min = ', global_sh_min
      write(0,*) '--- global_sh_max = ', global_sh_max
      write(0,*)

      write(0,*) '--- config_use_spechumd = ', config_use_spechumd 
      if(.not. config_use_spechumd .or. (global_sh_min==0._RKIND .and. global_sh_max==0._RKIND)) then
         !--- calculate the saturation mixing ratio and interpolated first-guess relative humidity:
         if (config_use_spechumd) then
            write(0,*) 'WARNING: config_use_spechumd=T, but specific humidity was not found in '//trim(config_met_prefix)//':'//config_start_time(1:13)
         end if
         write(0,*) ' *** initializing water vapor mixing ratio using first-guess relative humidity' 
         write(0,*)

         do k = 1, nVertLevels
            do iCell = 1, nCells
!              es = svp1*10.*exp(svp2*(t(k,iCell)-svpt0)/(t(k,iCell)-svp3))
!              es = min(es,0.99*0.01*pressure(k,iCell))
!              rs = 0.622*es*100. / (pressure(k,iCell)-es*100.)
               rs = rslf(pressure(k,iCell),t(k,iCell))
               if(t(k,iCell) .lt. svpt0) rs = rsif(pressure(k,iCell),t(k,iCell))
               scalars(index_qv,k,iCell) = 0.01_RKIND*rs*relhum(k,iCell)
             enddo
         enddo
      else
         !--- use the interpolated first-guess specific humidity:
         write(0,*) ' *** initializing water vapor mixing ratio using first-guess specific humidity' 
         write(0,*)
         do k = 1, nVertLevels
            do iCell = 1, nCells
               scalars(index_qv,k,iCell) = spechum(k,iCell)/(1._RKIND-spechum(k,iCell))
            enddo
         enddo
      endif

      !
      ! Diagnose fields needed in initial conditions file (u, w, rho, theta)
      ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
      !
      do iCell=1,nCells

         ! Q2
         es = 6.112 * exp((17.27*(t2m(iCell) - 273.16))/(t2m(iCell) - 35.86))
         rs = 0.622 * es * 100. / (psfc(iCell) - es * 100.)
         q2(iCell) = 0.01 * rs * rh2(iCell)
    

         do k=1,nVertLevels
            ! PI
            p(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)

            ! THETA - can compute this using PI instead
!            t(k,iCell) = t(k,iCell) / p(k,iCell)
            t(k,iCell) = t(k,iCell) * (p0 / pressure(k,iCell)) ** (rgas / cp)

            ! RHO_ZZ
            rho_zz(k,iCell) = pressure(k,iCell) / rgas / (p(k,iCell) * t(k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / (1.0 + scalars(index_qv,k,iCell))
         end do
      end do


      !
      ! Calculation of the initial precipitable water:
      ! 
      do iCell = 1,nCells
         precipw(iCell) = 0.0
         do k = 1,nVertLevels
            precipw(iCell) = precipw(iCell) + rho_zz(k,iCell)*scalars(index_qv,k,iCell)*(zgrid(k+1,iCell)-zgrid(k,iCell))
         end do
      end do

      !
      ! Reference state based on a dry isothermal atmosphere
      !
      do iCell=1,nCells
         do k=1,nz1
            ztemp    = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
            ppb(k,iCell) = p0*exp(-gravity*ztemp/(rgas*t0b))      ! pressure_base
            pb (k,iCell) = (ppb(k,iCell)/p0)**(rgas/cp)           ! exner_base
!            rb (k,iCell) = ppb(k,iCell)/(rgas*t0b*zz(k,iCell))    ! rho_base
            rb (k,iCell) = ppb(k,iCell)/(rgas*t0b)                ! rho_base
            tb (k,iCell) = t0b/pb(k,iCell)                        ! theta_base
            rtb(k,iCell) = rb(k,iCell)*tb(k,iCell)                ! rtheta_base
            p  (k,iCell) = pb(k,iCell)                            ! exner
            pp (k,iCell) = 0.                                     ! pressure_p
            rr (k,iCell) = 0.                                     ! rho_p
         end do
      end do

      do iCell=1,nCells
         do k=1,nVertLevels

!  WCS 20130821 - couple with vertical metric

            rb(k,iCell) = rb(k,iCell) / zz(k,iCell)
            rho_zz(k,iCell) = rho_zz(k,iCell) / zz(k,iCell)

            pp(k,iCell) = pressure(k,iCell) - ppb(k,iCell) 
            rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         end do
      end do

      do iCell=1,nCells
         k = 1
!  WCS 20130821 - couple with vertical metric, note: rr is coupled here
         rho_zz(k,iCell) = ((pressure(k,iCell) / p0)**(cv / cp)) * (p0 / rgas) &
                            / (t(k,iCell)*(1.0 + 1.61*scalars(index_qv,k,iCell))) / zz(k,iCell)
         rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         do k=2,nVertLevels
            it = 0
            p_check = 2.0 * 0.0001
            do while ( (it < 30) .and. (p_check > 0.0001) )

               p_check = pp(k,iCell)
!  WCS 20130821 - MPAS hydrostatic relation
               pp(k,iCell) = pp(k-1,iCell) - (fzm(k)*rr(k,iCell) + fzp(k)*rr(k-1,iCell))*gravity*dzu(k) &
                                           - (fzm(k)*rho_zz(k,iCell)*scalars(index_qv,k,iCell) &
                                                  + fzp(k)*rho_zz(k-1,iCell)*scalars(index_qv,k-1,iCell))*gravity*dzu(k)
               pressure(k,iCell) = pp(k,iCell) + ppb(k,iCell)
               p(k,iCell) = (pressure(k,iCell) / p0) ** (rgas / cp)
!  WCS 20130821 - couple with vertical metric
               rho_zz(k,iCell) = pressure(k,iCell) / rgas &
                     / (p(k,iCell)*t(k,iCell)*(1.0 + 1.61*scalars(index_qv,k,iCell)))/zz(k,iCell)
               rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

               p_check = abs(p_check - pp(k,iCell))
                
               it = it + 1
            end do
         end do
      end do

      ! Compute theta_m and rho-tilde
      do iCell=1,nCells
         do k=1,nVertLevels
            t(k,iCell) = t(k,iCell) * (1.0 + 1.61*scalars(index_qv,k,iCell))
!!  WCS 20130821 - coupling with vertical metric already accomplished...
!!            rho_zz(k,iCell) = rho_zz(k,iCell) / zz(k,iCell)
!!            rb(k,iCell) = rb(k,iCell) / zz(k,iCell)
!  WCS 20130821 - decouple rr from vertical metric
            rr(k,iCell) = rr(k,iCell)*zz(k,iCell)
         end do
      end do

      do iEdge=1,nEdges
         do k=1,nVertLevels
            ru(k,iEdge) = u(k,iEdge) * 0.5*(rho_zz(k,cellsOnEdge(1,iEdge)) + rho_zz(k,cellsOnEdge(2,iEdge)))
         end do
      end do


      rw= 0.0
      w = 0.0
      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               rw(k,cell2) = rw(k,cell2)    &
                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do
   
      deallocate(vert_level)

     
      ! Calculate surface pressure (This is an ad-hoc calculation. The actual surface pressure is actually re-calculated at
      !the top of the subroutine MPAS_to_physics in ../core_atmos_physics/mpas_atmphys_interface_nhyd.F
      do iCell=1,nCells
         surface_pressure(iCell) = 0.5*gravity/rdzw(1)                                              &
                                 * (1.25* rho_zz(1,iCell) * (1. + scalars(index_qv, 1, iCell))  &
                                 -  0.25* rho_zz(2,iCell) * (1. + scalars(index_qv, 2, iCell)))
         surface_pressure(iCell) = surface_pressure(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do


      end if    ! config_met_interp

   end subroutine init_atm_case_gfs


   integer function nearest_edge(target_lat, target_lon, &
                                 start_edge, &
                                 nCells, nEdges, maxEdges, nEdgesOnCell, edgesOnCell, cellsOnEdge, latCell, lonCell, latEdge, lonEdge)

      implicit none

      real (kind=RKIND), intent(in) :: target_lat, target_lon
      integer, intent(in) :: start_edge
      integer, intent(in) :: nCells, nEdges, maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
      real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
      real (kind=RKIND), dimension(nEdges), intent(in) :: latEdge, lonEdge

      integer :: i, cell1, cell2, iCell
      integer :: iEdge
      integer :: current_edge
      real (kind=RKIND) :: cell1_dist, cell2_dist
      real (kind=RKIND) :: current_distance, d
      real (kind=RKIND) :: nearest_distance

      nearest_edge = start_edge
      current_edge = -1

      do while (nearest_edge /= current_edge)
         current_edge = nearest_edge
         current_distance = sphere_distance(latEdge(current_edge), lonEdge(current_edge), target_lat, target_lon, 1.0_RKIND)
         nearest_edge = current_edge
         nearest_distance = current_distance
         cell1 = cellsOnEdge(1,current_edge)
         cell2 = cellsOnEdge(2,current_edge)
         cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), target_lat, target_lon, 1.0_RKIND)
         cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), target_lat, target_lon, 1.0_RKIND)
         if (cell1_dist < cell2_dist) then
            iCell = cell1
         else
            iCell = cell2
         end if
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            if (iEdge <= nEdges) then
               d = sphere_distance(latEdge(iEdge), lonEdge(iEdge), target_lat, target_lon, 1.0_RKIND)
               if (d < nearest_distance) then
                  nearest_edge = iEdge
                  nearest_distance = d
               end if
            end if
         end do
      end do

   end function nearest_edge


   real (kind=RKIND) function vertical_interp(target_z, nz, zf, order, extrap, surface_val, sealev_val)

      implicit none

      real (kind=RKIND), intent(in) :: target_z
      integer, intent(in) :: nz 
      real (kind=RKIND), dimension(2,nz), intent(in) :: zf      ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
      integer, intent(in), optional :: order
      integer, intent(in), optional :: extrap
      real (kind=RKIND), intent(in), optional :: surface_val
      real (kind=RKIND), intent(in), optional :: sealev_val

      integer :: k, lm, lp
      real (kind=RKIND) :: wm, wp
      real (kind=RKIND) :: slope

      integer :: interp_order, extrap_type
      real (kind=RKIND) :: surface, sealevel


      if (present(order)) then
         interp_order = order
      else
         interp_order = 2
      end if

      if (present(extrap)) then
         extrap_type = extrap
      else
         extrap_type = 1
      end if

      if (present(surface_val)) then
         surface = surface_val
      else
         surface = 200100.0
      end if

      if (present(sealev_val)) then
         sealevel = sealev_val
      else
         sealevel = 201300.0
      end if

      !
      ! Extrapolation required
      !
      if (target_z < zf(1,1)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,1)
         else if (extrap_type == 1) then
            slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
            vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
         else if (extrap_type == 2) then
            vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
         end if
         return
      end if
      if (target_z >= zf(1,nz)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,nz)
         else if (extrap_type == 1) then
            slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
            vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
         end if
         return
      end if


      !
      ! No extrapolation required
      !
      do k=1,nz-1
         if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
            lm = k
            lp = k+1
            wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
            wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
            exit
         end if
      end do

      vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

      return

   end function vertical_interp


!----------------------------------------------------------------------------------------------------------

   real (kind=RKIND) function env_qv( z, temperature, pressure, rh_max )

      implicit none
      real (kind=RKIND) :: z, temperature, pressure, ztr, es, qvs, p0, rh_max

      p0 = 100000.

!      ztr = 5000.
!
!      if(z .gt. ztr) then
!         env_qv = 0.
!      else
!         if(z.lt.2000.) then
!            env_qv = .5
!         else
!            env_qv = .5*(1.-(z-2000.)/(ztr-2000.))
!         end if
!      end if

       if (pressure .lt. 50000. ) then
           env_qv = 0.0
       else
           env_qv = (1.-((p0-pressure)/50000.)**1.25)
       end if

       env_qv = min(rh_max,env_qv)

! env_qv is the relative humidity, turn it into mixing ratio
       if (temperature .gt. 273.15) then
           es  = 1000.*0.6112*exp(17.67*(temperature-273.15)/(temperature-29.65))
       else
           es  = 1000.*0.6112*exp(21.8745584*(temperature-273.16)/(temperature-7.66))
       end if
       qvs = (287.04/461.6)*es/(pressure-es)

       ! qvs =  380.*exp(17.27*(temperature-273.)/(temperature-36.))/pressure

        env_qv = env_qv*qvs

   end function env_qv


   subroutine physics_idealized_init(mesh, fg)
   
      implicit none
      
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(inout) :: fg
      
      !local variables:
      integer :: iCell, iMonth, iSoil
      integer, pointer :: nCells, nSoilLevels, nMonths
      integer, dimension(:), pointer :: landmask, lu_index, soilcat_top
      real (kind=RKIND), dimension(:), pointer :: ter, xice, shdmin, shdmax, vegfra, sfc_albbck, xland, seaice
      real (kind=RKIND), dimension(:), pointer :: snow, snowc, snoalb, snowh, skintemp, sst, tmn
      real (kind=RKIND), dimension(:,:), pointer :: tslb, smcrel, sh2o, smois, dzs, albedo12m, greenfrac
      
      !---------------------------------------------------------------------------------------------

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nSoilLevels', nSoilLevels)
      call mpas_pool_get_dimension(mesh, 'nMonths', nMonths)

      call mpas_pool_get_array(mesh, 'ter', ter)
      call mpas_pool_get_array(mesh, 'landmask', landmask)
      call mpas_pool_get_array(mesh, 'lu_index', lu_index)
      call mpas_pool_get_array(mesh, 'soilcat_top', soilcat_top)
      call mpas_pool_get_array(mesh, 'shdmin', shdmin)
      call mpas_pool_get_array(mesh, 'shdmax', shdmax)
      call mpas_pool_get_array(mesh, 'snoalb', snoalb)
      call mpas_pool_get_array(mesh, 'albedo12m', albedo12m)
      call mpas_pool_get_array(mesh, 'greenfrac', greenfrac)

      call mpas_pool_get_array(fg, 'xice', xice)
      call mpas_pool_get_array(fg, 'vegfra', vegfra)
      call mpas_pool_get_array(fg, 'sfc_albbck', sfc_albbck)
      call mpas_pool_get_array(fg, 'xland', xland)
      call mpas_pool_get_array(fg, 'seaice', seaice)
      call mpas_pool_get_array(fg, 'snow', snow)
      call mpas_pool_get_array(fg, 'snowc', snowc)
      call mpas_pool_get_array(fg, 'snowh', snowh)
      call mpas_pool_get_array(fg, 'skintemp', skintemp)
      call mpas_pool_get_array(fg, 'sst', sst)
      call mpas_pool_get_array(fg, 'tmn', tmn)
      call mpas_pool_get_array(fg, 'tslb', tslb)
      call mpas_pool_get_array(fg, 'smcrel', smcrel)
      call mpas_pool_get_array(fg, 'sh2o', sh2o)
      call mpas_pool_get_array(fg, 'smois', smois)
      call mpas_pool_get_array(fg, 'dzs', dzs)
      
      !initialization of surface input variables that are not needed if we run the current set of
      !idealized test cases:
      
      
      do iCell = 1, nCells
      
         !terrain,soil type, and vegetation:
         ter(iCell) = 0.0
         xice(iCell) = 0.0
         landmask(iCell) = 0
         lu_index(iCell) = 0
         soilcat_top(iCell) = 0
         shdmin(iCell) = 0.0
         shdmax(iCell) = 0.0
         vegfra(iCell) = 0.0
         sfc_albbck(iCell) = 0.0
         xland(iCell) = 0.0
         seaice(iCell) = 0.0
      
         !snow coverage:
         snow(iCell) = 0.0
         snowc(iCell) = 0.0
         snoalb(iCell) = 0.08
         snowh(iCell) = 0.0
      
         !surface and sea-surface temperatures:
         skintemp(iCell) = 288.0
         sst(iCell) = 288.0
      
         !soil layers:
         tmn(iCell) = 288.0
         do iSoil = 1, nSoilLevels
            tslb(iSoil,iCell)   = 288.0
            smcrel(iSoil,iCell) =   0.0
            sh2o(iSoil,iCell) =   0.0
            smois(iSoil,iCell) =   0.0
            dzs(iSoil,iCell) =   0.0
         end do
      
         !monthly climatological surface albedo and greeness fraction:
         do iMonth = 1, nMonths
            albedo12m(iMonth,iCell) = 0.08
            greenfrac(iMonth,iCell) = 0.0
         end do
      
      end do
   
   end subroutine physics_idealized_init
   
   
   subroutine decouple_variables(mesh, nCells, nVertLevels, state, diag)

      implicit none

      type (mpas_pool_type), intent(in) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag

      integer :: iCell, iEdge, k

      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: rdzw
      real (kind=RKIND), dimension(:,:), pointer :: zz, pp, ppb, rho, rho_zz, theta, theta_m
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer, pointer :: index_qv

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz)
      call mpas_pool_get_array(state, 'theta_m', theta_m)
      call mpas_pool_get_array(state, 'scalars', scalars)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
     
      ! Compute surface pressure
      do iCell=1,nCells
         surface_pressure(iCell) = 0.5*gravity/rdzw(1)                                        &
                                   * (1.25* rho_zz(1,iCell) * (1. + scalars(index_qv, 1, iCell))  &
                                      -  0.25* rho_zz(2,iCell) * (1. + scalars(index_qv, 2, iCell)))
         surface_pressure(iCell) = surface_pressure(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do


      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = theta_m(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine decouple_variables

end module init_atm_cases

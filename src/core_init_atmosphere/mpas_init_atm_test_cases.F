! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_test_cases

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use atm_advection
   use mpas_atmphys_initialize_real
   use mpas_RBF_interpolation
   use mpas_vector_reconstruction
   use mpas_timer
   use mpas_init_atm_static
   use mpas_init_atm_surface

   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping !, only: MPAS_Time_type, MPAS_TimeInterval_type, MPAS_Clock_type, &
                        !        mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti
   


   contains


   subroutine init_atm_setup_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the hydrostatic test case
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr


      !
      ! Do some quick checks to make sure compile options are compatible with the chosen test case
      !
      if (config_test_case == 6) then
#ifndef ROTATED_GRID
         write(0,*) '*****************************************************************'
         write(0,*) 'To initialize and run the mountain wave test case (case 6),'
         write(0,*) '   please clean and re-compile init_atmosphere with -DROTATED_GRID'
         write(0,*) '   added to the specification of MODEL_FORMULATION'
         write(0,*) '   at the top of the Makefile.'
         write(0,*) '*****************************************************************'
         call mpas_dmpar_abort(domain % dminfo)
#endif
      else
#ifdef ROTATED_GRID
         write(0,*) '*****************************************************************'
         write(0,*) 'Only test case 6 should use code compiled with -DROTATED_GRID'
         write(0,*) '   specified in the Makefile.'
         write(0,*) '*****************************************************************'
         call mpas_dmpar_abort(domain % dminfo)
#endif
      end if



      if ((config_test_case == 1) .or. (config_test_case == 2) .or. (config_test_case == 3)) then

         write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
         if (config_test_case == 1) write(0,*) ' no initial perturbation '
         if (config_test_case == 2) write(0,*) ' initial perturbation included '
         if (config_test_case == 3) write(0,*) ' normal-mode perturbation included '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call init_atm_test_case_jw(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            call decouple_variables(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if ((config_test_case == 4) .or. (config_test_case ==5)) then

         write(0,*) ' squall line - super cell test case '
         if (config_test_case == 4) write(0,*) ' squall line test case' 
         if (config_test_case == 5) write(0,*) ' supercell test case'
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call init_atm_test_case_squall_line(domain % dminfo, block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            call decouple_variables(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6 ) then

         write(0,*) ' mountain wave test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call init_atm_test_case_mtn_wave(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            call decouple_variables(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag)
            write(0,*) ' returned from test case setup '
            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 7 ) then

         write(0,*) ' real-data GFS test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            if (config_static_interp) then

               !  
               !  Without a convex mesh partition file, interpolating static fields in parallel 
               !     will give incorrect results. Since it is very unlikely that typical users
               !     will have convex partitions, it's safer to just stop if multiple MPI tasks are
               !     detected when performing the static_interp step.
               !  
               if (domain % dminfo % nprocs > 1) then
                  write(0,*) ' '
                  write(0,*) '****************************************************************'
                  write(0,*) 'Error: Interpolation of static fields does not work in parallel.'
                  write(0,*) 'Please run the static_interp step using only a single MPI task.'
                  write(0,*) '****************************************************************'
                  write(0,*) ' '
                  call mpas_dmpar_abort(domain % dminfo)
               end if

               call init_atm_static(block_ptr % mesh)
               call init_atm_static_orogwd(block_ptr % mesh)
            endif
            call init_atm_test_case_gfs(block_ptr % mesh, block_ptr % fg, & 
                                        block_ptr % state % time_levs(1) % state, block_ptr % diag, &
                                        block_ptr % diag_physics, config_test_case)
            if (config_met_interp) call physics_initialize_real(block_ptr % mesh, block_ptr % fg, domain % dminfo)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 8 ) then

         write(0,*) 'real-data surface (SST) update test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            ! Defined in mpas_init_atm_surface.F
            call init_atm_test_case_sfc(domain, domain % dminfo, block_ptr % mesh,block_ptr % fg, block_ptr % state % time_levs(1) % state)
            block_ptr => block_ptr % next
         end do

      else

         write(0,*) ' '
         write(0,*) ' ****************************************************'
         write(0,*) ' Only test cases 1 through 8 are currently supported.'
         write(0,*) ' ****************************************************'
         write(0,*) ' '
         call mpas_dmpar_abort(domain % dminfo)

      end if


      ! Copy initialized state to all time levels
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         do i=2,nTimeLevs
            call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
         end do
         block_ptr => block_ptr % next
      end do

      !initialization of surface input variables technically not needed to run our current set of
      !idealized test cases:
      if (config_test_case < 7)  then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call physics_idealized_init(block_ptr % mesh, block_ptr % fg)
            block_ptr => block_ptr % next
         end do
      endif

   end subroutine init_atm_setup_test_case

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_test_case_jw(grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      !type (diag_physics_type), intent(inout) :: diag_physics
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
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      
!.. initialization of moisture:
      integer:: index_qv
      real (kind=RKIND),parameter :: rh_max = 0.40 ! Maximum relative humidity
!      real (kind=RKIND),parameter :: rh_max = 0.70 ! Maximum relative humidity
      real (kind=RKIND),dimension(grid % nVertLevels, grid % nCells) :: qsat, relhum
      real (kind=RKIND),dimension(:,:,:),pointer:: scalars
!.. end initialization of moisture.

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, AreaCell 
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: u, v, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND) :: es, qvs, xnutr, znut, ptemp 
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(grid % nVertLevels ) :: eta, etav, teta, ppi, tt, temperature_1d

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      logical, parameter :: rebalance = .true.
      integer, parameter :: nlat=721
      real (kind=RKIND), dimension(grid % nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(grid % nVertLevels + 1, nlat) :: zgrid_2d
      real (kind=RKIND), dimension(grid % nVertLevels, nlat) :: u_2d, pp_2d, rho_2d, qv_2d, etavs_2d, zz_2d
      real (kind=RKIND), dimension(grid % nVertLevels, nlat-1) :: zx_2d 
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat, hx_1d
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

      logical, parameter :: moisture = .true.
!      logical, parameter :: moisture = .false.

      !
      ! Scale all distances and areas from a unit sphere to one with radius sphere_radius
      !
      grid % xCell % array = grid % xCell % array * grid % sphere_radius
      grid % yCell % array = grid % yCell % array * grid % sphere_radius
      grid % zCell % array = grid % zCell % array * grid % sphere_radius
      grid % xVertex % array = grid % xVertex % array * grid % sphere_radius
      grid % yVertex % array = grid % yVertex % array * grid % sphere_radius
      grid % zVertex % array = grid % zVertex % array * grid % sphere_radius
      grid % xEdge % array = grid % xEdge % array * grid % sphere_radius
      grid % yEdge % array = grid % yEdge % array * grid % sphere_radius
      grid % zEdge % array = grid % zEdge % array * grid % sphere_radius
      grid % dvEdge % array = grid % dvEdge % array * grid % sphere_radius
      grid % dcEdge % array = grid % dcEdge % array * grid % sphere_radius
      grid % areaCell % array = grid % areaCell % array * grid % sphere_radius**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * grid % sphere_radius**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * grid % sphere_radius**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      dvEdge            => grid % dvEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array

      deriv_two  => grid % deriv_two % array
      zb  => grid % zb % array
      zb3 => grid % zb3% array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array

      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array

      ppb => diag % pressure_base % array
      pp  => diag % pressure_p % array

      rho_zz => state % rho_zz % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array

      surface_pressure => diag % surface_pressure % array

!.. initialization of moisture:
      scalars => state % scalars % array
      !qsat    => diag_physics % qsat % array
      !relhum  => diag_physics % relhum % array
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0
!.. end initialization of moisture.

      surface_pressure(:) = 0.0

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      index_qv = state % index_qv

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      r_earth = grid % sphere_radius
      omega_e = omega
      p0 = 1.e+05

      write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,grid % nCells
        do k=1,nz
          phi = grid % latCell % array (iCell)
          hx(k,iCell) = u0/gravity*cos(etavs)**1.5                                   &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        enddo
      enddo

      !     Metrics for hybrid coordinate and vertical stretching

      str = 1.8
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

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

      do iCell=1,grid % nCells
        do k=1,nz	
          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

      !do k=1,nz1
      !  write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
      !enddo

      !do k=1,nz1
      !  write(0,*) ' k, zx(k,1) ',k,zx(k,1)
      !enddo

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
          zz_2d (k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz_2d(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do


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

          phi = lat_2d (i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))


            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture 
            if (moisture) then
              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
          end do

          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz_2d(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                            &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+qv_2d(1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+qv_2d(2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1

!              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                        &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*qv_2d(k  ,i)   &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*qv_2d(k+1,i))

              ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                       &
                            ( (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*qv_2d(k  ,i))*fzp(k+1)   &
                            + (rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*qv_2d(k+1,i))*fzm(k+1) )
            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          rho_2d(k,i) = rr(k,i)+rb(k,i)
          pp_2d(k,i) = pp(k,i)
          etavs_2d(k,i) = ((ppb(k,i)+pp(k,i))/p0 - 0.252)*pii/2.
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

        call init_atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                        cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat,grid%sphere_radius)

      end if

!******************************************************************      

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, grid % nCells
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
!         enddo
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
          phi = grid % latCell % array (i)
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
               endif
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
      !do iCell = 1, grid % nCells
      !   if(iCell == 1 .or. iCell == grid % nCells) then
      !      do k = nz1, 1, -1
      !         write(0,202) iCell,k,t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)
      !      enddo
      !      write(0,*)
      !   endif
      !enddo

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,grid % nEdges

         vtx1 = grid % VerticesOnEdge % array (1,iEdge)
         vtx2 = grid % VerticesOnEdge % array (2,iEdge)
         lat1 = grid%latVertex%array(vtx1)
         lat2 = grid%latVertex%array(vtx2)
         iCell1 = grid % cellsOnEdge % array(1,iEdge)
         iCell2 = grid % cellsOnEdge % array(2,iEdge)
         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * grid % sphere_radius / grid % dvEdge % array(iEdge)

         if (config_test_case == 2) then
            r_pert = sphere_distance( grid % latEdge % array (iEdge), grid % lonEdge % array (iEdge), &
                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1) * grid % sphere_radius / grid % dvEdge % array(iEdge)

         else if (config_test_case == 3) then
            lon_Edge = grid % lonEdge % array(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * grid % sphere_radius / grid % dvEdge % array(iEdge)
         else
            u_pert = 0.0
         end if

         if (rebalance) then

           call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,grid % dvEdge % array(iEdge),grid%sphere_radius,u0,nz1,nlat)
           do k=1,grid % nVertLevels
             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
             state % u % array(k,iEdge) = fluxk + u_pert
           end do

         else 

           do k=1,grid % nVertLevels
             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
             fluxk = u0*flux*(cos(etavs)**1.5)
             state % u % array(k,iEdge) = fluxk + u_pert
           end do

         end if

         cell1 = grid % CellsOnEdge % array(1,iEdge)
         cell2 = grid % CellsOnEdge % array(2,iEdge)
         do k=1,nz1
            diag % ru % array (k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*state % u % array (k,iEdge)
         end do

      !
      ! Generate rotated Coriolis field
      !

         grid % fEdge % array(iEdge) = 2.0 * omega_e * &
                                       ( -cos(grid%lonEdge%array(iEdge)) * cos(grid%latEdge%array(iEdge)) * sin(alpha_grid) + &
                                         sin(grid%latEdge%array(iEdge)) * cos(alpha_grid) &
                                       )
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 2.0 * omega_e * &
                                         (-cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha_grid) + &
                                          sin(grid%latVertex%array(iVtx)) * cos(alpha_grid) &
                                         )
      end do

      !
      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
      !

      !
      !     pre-calculation z-metric terms in omega eqn.
      !
      do iEdge = 1,grid % nEdges
         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         do k = 1, grid%nVertLevels

            if (config_theta_adv_order == 2) then

               z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

            else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

               d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)

!  WCS fix 20120711

                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell2))
                  end do             

               z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                             - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               if (config_theta_adv_order == 3) then
                  z_edge3 =  - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
               else
                  z_edge3 = 0.
               end if

            end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/AreaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/AreaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell2)

         end do

      end do

      ! for including terrain
      diag % rw % array = 0.
      state % w % array = 0.
      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*diag % ru % array(k,iEdge)+fzp(k)*diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,grid%nCells
         do k=2,grid%nVertLevels
            state % w % array(k,iCell) = diag % rw % array(k,iCell)     &
                                       / (fzp(k) * state % rho_zz % array(k-1,iCell) + fzm(k) * state % rho_zz % array(k,iCell))
         end do
      end do


      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, grid%nVertLevels
               diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
           end do
         end do
      end do

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

        write(0,*) ' i, psurf, lat ',i,psurf,grid%latCell%array(i)*180./3.1415828
      enddo

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_jw


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


   subroutine init_atm_test_case_squall_line(dminfo, grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup squall line and supercell test case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, cell1, cell2, nCellsSolve
      integer :: index_qv

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znu, znw, znwc, znwv
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rh, thi, tbi, cqwb

      real (kind=RKIND) ::  r, xnutr
      real (kind=RKIND) ::  ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(grid % nVertLevels ) :: qvb
      real (kind=RKIND), dimension(grid % nVertLevels ) :: t_init_1d

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2
      real (kind=RKIND) :: ztr, thetar, ttr, thetas, um, us, zts, pitop, pibtop, ptopb, ptop, rcp, rcv, p0
      real (kind=RKIND) :: radx, radz, zcent, xmid, delt, xloc, rad, yloc, ymid, a_scale
      real (kind=RKIND) :: pres, temp, es, qvs

      !
      ! Scale all distances
      !

      a_scale = 1.0

      grid % xCell % array = grid % xCell % array * a_scale
      grid % yCell % array = grid % yCell % array * a_scale
      grid % zCell % array = grid % zCell % array * a_scale
      grid % xVertex % array = grid % xVertex % array * a_scale
      grid % yVertex % array = grid % yVertex % array * a_scale
      grid % zVertex % array = grid % zVertex % array * a_scale
      grid % xEdge % array = grid % xEdge % array * a_scale
      grid % yEdge % array = grid % yEdge % array * a_scale
      grid % zEdge % array = grid % zEdge % array * a_scale
      grid % dvEdge % array = grid % dvEdge % array * a_scale
      grid % dcEdge % array = grid % dcEdge % array * a_scale
      grid % areaCell % array = grid % areaCell % array * a_scale**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * a_scale**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * a_scale**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array

      ppb => diag % pressure_base % array
      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array
      cqw => diag % cqw % array

      rho_zz => state % rho_zz % array

      pp => diag % pressure_p % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array
      u => state % u % array
      ru => diag % ru % array

      scalars => state % scalars % array

      index_qv = state % index_qv

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      xnutr = 0.
      zd = 12000.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

     write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,grid % nCells
        do k=1,nz
          hx(k,iCell) = 0.  ! squall line or supercell on flat plane
        enddo
      enddo

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
!	    write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)			
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

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

      do iCell=1,grid % nCells
        do k=1,nz	
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)	
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

!
! convective initialization
!
         ztr    = 12000.
         thetar = 343.
         ttr    = 213.
         thetas = 300.5

!         write(0,*) ' rgas, cp, gravity ',rgas,cp, gravity

      if ( config_test_case == 4) then ! squall line parameters
         um = 12.
         us = 10.
         zts = 2500.
      else if (config_test_case == 5) then !supercell parameters
         um = 30.
         us = 15.
         zts = 5000.
      end if

         do i=1,grid % nCells
            do k=1,nz1
               ztemp = .5*(zgrid(k,i)+zgrid(k+1,i))
               if(ztemp .gt. ztr) then
                  t (k,i) = thetar*exp(9.8*(ztemp-ztr)/(1003.*ttr))
                  rh(k,i) = 0.25
               else
                  t (k,i) = 300.+43.*(ztemp/ztr)**1.25
                  rh(k,i) = (1.-0.75*(ztemp/ztr)**1.25)
                  if(t(k,i).lt.thetas) t(k,i) = thetas
               end if
               tb(k,i) = t(k,i)
               thi(k,i) = t(k,i)
               tbi(k,i) = t(k,i)
               cqw(k,i) = 1.
               cqwb(k,i) = 1.
            end do
         end do

!         rh(:,:) = 0.

!  set the velocity field - we are on a plane here.

         do i=1, grid % nEdges
            cell1 = grid % CellsOnEdge % array(1,i)
            cell2 = grid % CellsOnEdge % array(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               if(ztemp.lt.zts)  then
                  u(k,i) = um*ztemp/zts
               else
                  u(k,i) = um
               end if
               if(i == 1 ) grid % u_init % array(k) = u(k,i) - us
               u(k,i) = cos(grid % angleEdge % array(i)) * (u(k,i) - us)
            end do
            end if
         end do

         call mpas_dmpar_bcast_reals(dminfo, nz1, grid % u_init % array)

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

      do i=1, grid % nCells
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
      do i= 1,grid%nCells
         do k=1,nz1
            temp     = p(k,i)*thi(k,i)
            pres     = p0*p(k,i)**(1./rcp)
            qvs      = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
            scalars(index_qv,k,i) = min(0.014_RKIND,rh(k,i)*qvs)
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

      do i= 1,grid%nCells
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
      do k=1,grid%nVertLevels
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

      if (config_test_case == 4) then          ! squall line prameters
         call mpas_dmpar_max_real(dminfo, maxval(grid % xCell % array(:)), xmid)
         xmid = xmid * 0.5
         ymid = 0.0          ! Not used for squall line
      else if (config_test_case == 5) then     ! supercell parameters
         call mpas_dmpar_max_real(dminfo, maxval(grid % xCell % array(:)), xmid)
         call mpas_dmpar_max_real(dminfo, maxval(grid % yCell % array(:)), ymid)
         xmid = xmid * 0.5
         ymid = ymid * 0.5
      end if

      do i=1, grid % nCells
        xloc = grid % xCell % array(i) - xmid
        if (config_test_case == 4) then 
           yloc = 0.                            !squall line setting
        else if (config_test_case == 5) then
           yloc = grid % yCell % array(i) - ymid !supercell setting
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

        do i = 1, grid % nCells

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
        grid % qv_init % array(k) = scalars(index_qv,k,1)
      end do

      t_init_1d(:) = t(:,1)
      call mpas_dmpar_bcast_reals(dminfo, nz1, t_init_1d)
      call mpas_dmpar_bcast_reals(dminfo, nz1, grid % qv_init % array)

      do i=1,grid % ncells
         do k=1,nz1
            grid % t_init % array(k,i) = t_init_1d(k)
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do
      end do

      do i=1,grid % nEdges
        cell1 = grid % CellsOnEdge % array(1,i)
        cell2 = grid % CellsOnEdge % array(2,i)
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
       diag % rw % array = 0.
       state % w % array = 0.

       grid % zb % array = 0.
       grid % zb3% array = 0.

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 0.
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, grid%nVertLevels
                 diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
              end do
            end if
         end do
      end do

     ! write(0,*) ' k,u_init, t_init, qv_init '
     ! do k=1,grid%nVertLevels
     !   write(0,'(i2,3(2x,f14.10)') k,grid % u_init % array(k),grid % t_init% array(k),grid % qv_init % array(k)
     ! end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_squall_line


!----------------------------------------------------------------------------------------------------------


   subroutine init_atm_test_case_mtn_wave(grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: t0=288., hm=250.

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru 
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, deriv_two, zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, AreaCell, xCell, yCell 
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve
      integer :: index_qv

      real (kind=RKIND) :: ptop, pitop, ptopb, p0, flux, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND) :: ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rh
      real (kind=RKIND) :: es, qvs, xnutr, ptemp
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, cf1, cf2, cf3
      real (kind=RKIND) :: um, us,  rcp, rcv
      real (kind=RKIND) :: xmid, temp, pres, a_scale

      real (kind=RKIND) :: xi, xa, xc, xla, zinv, xn2, xn2m, xn2l, sm, dzh, dzht, dzmin, z_edge, z_edge3 

      integer, dimension(grid % nCells, 2) :: next_cell
      real (kind=RKIND),  dimension(grid % nCells) :: hxzt
      logical, parameter :: terrain_smooth = .false.

      !
      ! Scale all distances
      !

      a_scale = 1.0

      grid % xCell % array = grid % xCell % array * a_scale
      grid % yCell % array = grid % yCell % array * a_scale
      grid % zCell % array = grid % zCell % array * a_scale
      grid % xVertex % array = grid % xVertex % array * a_scale
      grid % yVertex % array = grid % yVertex % array * a_scale
      grid % zVertex % array = grid % zVertex % array * a_scale
      grid % xEdge % array = grid % xEdge % array * a_scale
      grid % yEdge % array = grid % yEdge % array * a_scale
      grid % zEdge % array = grid % zEdge % array * a_scale
      grid % dvEdge % array = grid % dvEdge % array * a_scale
      grid % dcEdge % array = grid % dcEdge % array * a_scale
      grid % areaCell % array = grid % areaCell % array * a_scale**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * a_scale**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * a_scale**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array  
      dvEdge            => grid % dvEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array
      deriv_two         => grid % deriv_two % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      zb => grid % zb % array
      zb3 => grid % zb3 % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array
 
      xCell => grid % xCell % array
      yCell => grid % yCell % array

      ppb => diag % pressure_base % array
      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array
      cqw => diag % cqw % array

      rho_zz => state % rho_zz % array

      pp => diag % pressure_p % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array
      u => state % u % array
      ru => diag % ru % array

      scalars => state % scalars % array

      index_qv = state % index_qv

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      xnutr = 0.1
      zd = 10500.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

      ! for hx computation
      xa = 5000. !SHP - should be changed based on grid distance 
      xla = 4000.
      xc = maxval (grid % xCell % array(:))/2. 

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
!	    write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)			
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

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

! setting for terrain
      do iCell=1,grid % nCells
         xi = grid % xCell % array(iCell)
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
         do i=1,grid % nCells 
            !option 1
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)-sqrt(3.)*grid % dcEdge % array(1)) next_cell(iCell,1) = i 
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)+sqrt(3.)*grid % dcEdge % array(1)) next_cell(iCell,2) = i 
            !option 2
            next_cell(iCell,1) = iCell - 8 ! note ny=4
            next_cell(iCell,2) = iCell + 8 ! note ny=4

            if (xCell(iCell).le. 3.*grid % dcEdge % array(1)) then
                next_cell(iCell,1) = 1
            else if (xCell(iCell).ge. maxval(xCell(:))-3.*grid % dcEdge % array(1)) then
                next_cell(iCell,2) = 1
            end if

         end do
      enddo
      
      write(0,*) ' hx computation complete '

      if (terrain_smooth) then
         write(0,*) '***************************************************************************'
         write(0,*) 'Please contact the MPAS-A developers for up-to-date terrain-smoothing code.'
         write(0,*) 'Otherwise, set terrain_smooth=.false. in the mountain wave test case'
         write(0,*) '   initialization routine and re-compile.'
         write(0,*) '***************************************************************************'
         call mpas_dmpar_abort(grid % block % domain % dminfo)
      end if

      do iCell=1,grid % nCells
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

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

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

         do i=1,grid % nCells
            do k=1,nz1
               ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
               tb(k,i) =  t0*(1. + xn2m/gravity*ztemp) 
               if(ztemp .le. zinv) then
                  t (k,i) = t0*(1.+xn2l/gravity*ztemp)
               else
                  t (k,i) = t0*(1.+xn2l/gravity*zinv+xn2/gravity*(ztemp-zinv)) 
               end if
                  rh(k,i) = 0. 
            end do
         end do

!  set the velocity field - we are on a plane here.

         do i=1, grid % nEdges
            cell1 = grid % CellsOnEdge % array(1,i)
            cell2 = grid % CellsOnEdge % array(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               u(k,i) = um
               if(i == 1 ) grid % u_init % array(k) = u(k,i) - us
#ifdef ROTATED_GRID
               u(k,i) = sin(grid % angleEdge % array(i)) * (u(k,i) - us)
#else
               u(k,i) = cos(grid % angleEdge % array(i)) * (u(k,i) - us)
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
                
      do i=1, grid % nCells
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
       do k=1,grid%nVertLevels
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

        do i = 1, grid % nCells

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
              scalars(index_qv,k,i) = min(0.014_RKIND,rh(k,i)*qvs)
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
                       grid % u_init % array(k), rr(k,1)
      end do

      do i=1,grid % ncells
         do k=1,nz1
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do

        do k=1,nz1
            grid % t_init % array(k,i) = t(k,i)
        end do
      end do

      do i=1,grid % nEdges
        cell1 = grid % CellsOnEdge % array(1,i)
        cell2 = grid % CellsOnEdge % array(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)    
          end do
        end if
      end do

!
!     pre-calculation z-metric terms in omega eqn.
!
      do iEdge = 1,grid % nEdges
         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, grid%nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell2))
                  end do             
             
                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. 

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.   
                  else 
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/AreaCell(cell1) 
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/AreaCell(cell2) 
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell1) 
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell2) 
  
            end do

         end if
       end do

!     for including terrain
      state % w % array(:,:) = 0.0
      diag % rw % array(:,:) = 0.0

!
!     calculation of omega, rw = zx * ru + zz * rw
!

      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))  
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux 
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux 

            if (config_theta_adv_order ==3) then
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,grid%nCells
         do k=2,grid%nVertLevels
            state % w % array(k,iCell) = diag % rw % array(k,iCell)     & 
                                       / (fzp(k) * state % rho_zz % array(k-1,iCell) + fzm(k) * state % rho_zz % array(k,iCell))
         end do
      end do


      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 0.
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, grid%nVertLevels
                 diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
              end do
            end if
         end do
      end do

!      do k=1,grid%nVertLevels
!        write(0,*) ' k,u_init, t_init, qv_init ',k,grid % u_init % array(k),grid % t_init% array(k),grid % qv_init % array(k)
!      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_mtn_wave


   subroutine init_atm_test_case_gfs(grid, fg, state, diag, diag_physics, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Real-data test case using GFS data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_dmpar
      use init_atm_read_met
      use init_atm_llxy
      use init_atm_hinterp

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (fg_type), intent(inout) :: fg
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (diag_physics_type), intent(inout):: diag_physics
      integer, intent(in) :: test_case

      type (block_type), pointer :: block
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
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:), pointer :: destField1d
      real (kind=RKIND), dimension(:,:), pointer :: destField2d
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      real (kind=RKIND) :: target_z
      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve
      integer :: nInterpPoints, ndims

      integer, dimension(5) :: interp_list
      real (kind=RKIND) :: maskval
      real (kind=RKIND) :: msgval
      real (kind=RKIND) :: fillval
      integer :: masked

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, edgesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, AreaCell 
      real (kind=RKIND), dimension(:,:), pointer :: v
      real (kind=RKIND), dimension(:,:), pointer :: sorted_arr

      type (field1DReal), pointer :: tempField
      type (field1DReal), target :: tempFieldTarget

      real(kind=RKIND), dimension(:), pointer :: hs, hs1
      real(kind=RKIND) :: hm, hm_global, zh, dzmin, dzmina, dzmina_global, dzminf, sm
      integer :: nsmterrain, kz, sfc_k
      logical :: hybrid, smooth

      integer :: it
      real (kind=RKIND) :: p_check

      ! For interpolating terrain and land use
      integer :: nx, ny
      integer :: istatus

      real (kind=RKIND), allocatable, dimension(:,:) :: rslab, maskslab
      integer, dimension(:), pointer :: mask_array
      integer, dimension(grid % nEdges), target :: edge_mask
      character (len=StrKIND) :: fname

      real (kind=RKIND) :: u, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r
      real (kind=RKIND) :: lat, lon, x, y

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rel_hum, temperature, qv
      real (kind=RKIND) :: ptmp, es, rs, rgas_moist, qvs, xnutr, znut, ptemp, rcv
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: sh, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(grid % nVertLevels) :: eta, etav, teta, ppi, tt

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      integer, parameter :: nlat=361
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
      real (kind=RKIND), dimension(grid % nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nlat, grid % nVertLevels) :: u_2d, etavs_2d
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: fsum
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2


      block => grid % block
      parinfo => block % parinfo
      dminfo => block % domain % dminfo

      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      dvEdge            => grid % dvEdge % array
      dcEdge            => grid % dcEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array
      cellsOnCell       => grid % cellsOnCell % array

      deriv_two  => grid % deriv_two % array
      zb  => grid % zb % array
      zb3 => grid % zb3% array

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      ter => grid % ter % array
      dss => grid % dss % array

      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array

      ppb => diag % pressure_base % array
      pp  => diag % pressure_p % array

      rho_zz => state % rho_zz % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array

      scalars => state % scalars % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      rcv = rgas/(cp-rgas)
      r_earth = grid % sphere_radius
      omega_e = omega
      p0 = 1.e+05

      scalars(:,:,:) = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN ADOPT GFS TERRAIN HEIGHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)
      do while (istatus == 0)
         if (index(field % field, 'SOILHGT') /= 0) then

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


            if (index(field % field, 'SOILHGT') /= 0) then
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => ter
               ndims = 1
            end if

            do i=1,nInterpPoints
               lat = latPoints(i)*DEG_PER_RAD
               lon = lonPoints(i)*DEG_PER_RAD
               call latlon_to_ij(proj, lat, lon, x, y)
               if (x < 0.5) then
                  lon = lon + 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               end if
               if (ndims == 1) then
                  destField1d(i) = interp_sequence(x, y, 1, field % slab, 1, field % nx, 1, field % ny, 1, 1, -1.e30_RKIND, interp_list, 1)
               else if (ndims == 2) then
                  destField2d(k,i) = interp_sequence(x, y, 1, field % slab, 1, field % nx, 1, field % ny, 1, 1, -1.e30_RKIND, interp_list, 1)
               end if
            end do
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
      allocate(hs (grid%nCells+1))
      allocate(hs1(grid%nCells+1))

!     Fourth order smoother for terrain

      nsmterrain = config_nsmterrain

      do i=1,nsmterrain

         do iCell=1,grid%nCells
            hs(iCell) = 0.
            if(ter(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)
                  hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                        / dcEdge(edgesOnCell(j,iCell))    &
                                        * (ter(cellsOnCell(j,iCell))-ter(iCell))
               end do
            endif
            hs(iCell) = ter(iCell) + 0.125*hs(iCell)
         end do

         do iCell=1,grid %nCells
            ter(iCell) = 0.
            if(hs(iCell) .ne. 0.) then
               do j = 1,nEdgesOnCell(iCell)
                  ter(iCell) = ter(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                          / dcEdge(edgesOnCell(j,iCell))    &
                                          * (hs(cellsOnCell(j,iCell))-hs(iCell))
               end do
            endif
!           ter(iCell) = hs(iCell) - 0.25*ter(iCell)
            ter(iCell) = hs(iCell) - 0.125*ter(iCell)
         end do

         ! note that ther variable ter used throughout this section is a pointer to grid % ter % array, here we are passing ter's parent field
         call mpas_dmpar_exch_halo_field(grid % ter)

      end do

      do iCell=1,grid % nCells
         hx(:,iCell) = ter(iCell)
      end do

      hm = maxval(ter(1:nCellsSolve))
      call mpas_dmpar_max_real(dminfo, hm, hm_global)
      hm = hm_global
      write(0,*) "max ter = ", hm

!     Metrics for hybrid coordinate and vertical stretching

      str = 1.5
!      str = 1.
      zt = config_ztop
      dz = zt/float(nz1)

      do k=1,nz
         zw(k) = (real(k-1)/real(nz1))**str*zt
         if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
      end do

!     ah(k) governs the transition between terrain-following 
!        and pure height coordinates
!           ah(k) = 1           is a smoothed terrain-following coordinate
!           ah(k) = 1.-zw(k)/zt is the basic terrain-following coordinate
!           ah(k) = 0           is a height coordinate
 
      hybrid = .true.
!      hybrid = .false.

      kz = nz
      if (hybrid) then
      
         zh = zt
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

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

!     Smoothing algorithm for coordinate surfaces 

      smooth = config_smooth_surfaces
!      smooth = .false.

      if (smooth) then

         dzmin = 0.5

         do k=2,kz-1
            hx(k,:) = hx(k-1,:)
            dzminf = zw(k)-zw(k-1)

!            dzmin = max(0.5_RKIND,1.-.5*zw(k)/hm)

            sm = .02*min(0.5_RKIND*zw(k)/hm,1.0_RKIND)
          
            do i=1,30
               do iCell=1,grid % nCells
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
               tempField % dimSizes(1) = grid % nCells
               tempField % sendList => parinfo % cellsToSend
               tempField % recvList => parinfo % cellsToRecv
               tempField % copyList => parinfo % cellsToCopy
               tempField % array => hs
               tempField % prev => null()
               tempField % next => null()

               call mpas_dmpar_exch_halo_field(tempField)

             !  dzmina = minval(hs(:)-hx(k-1,:))
               dzmina = minval(zw(k)+ah(k)*hs(1:nCellsSolve)-zw(k-1)-ah(k-1)*hx(k-1,1:nCellsSolve))
               call mpas_dmpar_min_real(dminfo, dzmina, dzmina_global)
             !  write(0,*) ' k,i, dzmina, dzmin, zw(k)-zw(k-1) ', k,i, dzmina, dzmin, zw(k)-zw(k-1)
               if (dzmina_global >= dzmin*(zw(k)-zw(k-1))) then
                  hx(k,:)=hs(:)
                  dzminf = dzmina_global
               else
                  exit
               end if
            end do
            write(0,*) k,i,sm,dzminf/(zw(k)-zw(k-1)),dzmina/(zw(k)-zw(k-1))
         end do

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

      do iCell=1,grid % nCells
         do k=1,nz	
            zgrid(k,iCell) = zw(k) + ah(k)*hx(k,iCell)
         end do
         do k=1,nz1
            zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
         end do
      end do

      do i=1, grid % nEdges
         iCell1 = grid % CellsOnEdge % array(1,i)
         iCell2 = grid % CellsOnEdge % array(2,i)
         do k=1,nz
            zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
         end do
      end do
      do i=1, grid % nCells
         do k=1,nz1
           ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
           dss(k,i) = 0.
           ztemp = zgrid(k,i)
           if (ztemp.gt.zd+.1)  then
               dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
           end if
         end do
      enddo

!      do k=1,nz1
!         write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
!      enddo

!      do k=1,nz1
!         write(0,*) ' k, zx(k,1) ',k,zx(k,1)
!      enddo


      ! For z-metric term in omega equation
      do iEdge = 1,grid % nEdges
         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, grid%nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell2))
                  end do             
             
                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. 

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.   
                  else 
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/AreaCell(cell1) 
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/AreaCell(cell2) 
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell1) 
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell2) 
  
            end do

         end if
      end do

      write(0,*) ' grid metrics setup complete '

      end if    ! config_vertical_grid


      if (config_met_interp) then

      !ldf (2011-11-19): added initialization of the sea-surface temperature, seaice fraction, and
      !seaice flag:
       fg % sst % array    = 0.
       fg % xice % array   = 0.
       fg % seaice % array = 0.
      !ldf end.

      !
      ! First, try to locate the LANDSEA field for use as an interpolation mask
      !
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)
         if (index(field % field, 'LANDSEA') /= 0) then

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

      edge_mask(:) = 1


      !
      ! Horizontally interpolate meteorological data
      !
      allocate(vert_level(config_nfglevels))
      vert_level(:) = -1.0

      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)

         interp_list(1) = FOUR_POINT
         interp_list(2) = SEARCH
         interp_list(3) = 0

         maskval = -1.0
         masked = -1
         fillval = 0.0
         msgval = -1.e30

         mask_array => grid % landmask % array

         if (index(field % field, 'UU') /= 0 .or. &
             index(field % field, 'VV') /= 0 .or. &
             index(field % field, 'TT') /= 0 .or. &
             index(field % field, 'RH') /= 0 .or. &
             index(field % field, 'GHT') /= 0 .or. &
             index(field % field, 'PMSL') /= 0 .or. &
             index(field % field, 'PSFC') /= 0 .or. &
             index(field % field, 'SOILHGT') /= 0 .or. &
             index(field % field, 'SM000010') /= 0 .or. &
             index(field % field, 'SM010040') /= 0 .or. &
             index(field % field, 'SM040100') /= 0 .or. &
             index(field % field, 'SM100200') /= 0 .or. &
             index(field % field, 'SM010200') /= 0 .or. &
             index(field % field, 'ST000010') /= 0 .or. &
             index(field % field, 'ST010040') /= 0 .or. &
             index(field % field, 'ST040100') /= 0 .or. &
             index(field % field, 'ST100200') /= 0 .or. &
             index(field % field, 'ST010200') /= 0 .or. &
             index(field % field, 'PRES') /= 0 .or. &
             index(field % field, 'SNOW') /= 0 .or. &
             index(field % field, 'SEAICE') /= 0 .or. &
             index(field % field, 'SKINTEMP') /= 0) then

            if (index(field % field, 'SM000010') /= 0 .or. &
                index(field % field, 'SM010040') /= 0 .or. &
                index(field % field, 'SM040100') /= 0 .or. &
                index(field % field, 'SM100200') /= 0 .or. &
                index(field % field, 'SM010200') /= 0 .or. &
                index(field % field, 'ST000010') /= 0 .or. &
                index(field % field, 'ST010040') /= 0 .or. &
                index(field % field, 'ST040100') /= 0 .or. &
                index(field % field, 'ST100200') /= 0 .or. &
                index(field % field, 'ST010200') /= 0 .or. &
                index(field % field, 'SNOW') /= 0 .or. &
                index(field % field, 'SEAICE') /= 0 .or. &
                index(field % field, 'SKINTEMP') /= 0) then
               k = 1
            else if (index(field % field, 'PMSL') == 0 .and. &
                     index(field % field, 'PSFC') == 0 .and. &
                     index(field % field, 'SOILHGT') == 0) then
               do k=1,config_nfglevels
                  if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
               end do
               if (k > config_nfglevels) write(0,*) 'ERROR: We seem to have more levels than we thought we should!'
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
                            loninc = real(field % deltalon,RKIND), &
                            lat1 = real(field % startlat,RKIND), &
                            lon1 = real(field % startlon,RKIND))
!                            nxmax = nint(360.0 / field % deltalon), &
            end if


            !
            ! Horizontally interpolate the field at level k
            !
            if (index(field % field, 'UU') /= 0) then
write(0,*) 'Interpolating U at ', k, vert_level(k)

               mask_array => edge_mask

               nInterpPoints = grid % nEdges
               latPoints => grid % latEdge % array
               lonPoints => grid % lonEdge % array
               destField2d => fg % u % array
               ndims = 2
            else if (index(field % field, 'VV') /= 0) then
write(0,*) 'Interpolating V at ', k, vert_level(k)

               mask_array => edge_mask

               nInterpPoints = grid % nEdges
               latPoints => grid % latEdge % array
               lonPoints => grid % lonEdge % array
               destField2d => fg % v % array
               ndims = 2
            else if (index(field % field, 'TT') /= 0) then
write(0,*) 'Interpolating T at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % t % array
               ndims = 2
            else if (index(field % field, 'RH') /= 0) then
write(0,*) 'Interpolating RH at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % rh % array
               ndims = 2
            else if (index(field % field, 'GHT') /= 0) then
write(0,*) 'Interpolating GHT at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % z % array
               ndims = 2
            else if (index(field % field, 'PRES') /= 0) then
write(0,*) 'Interpolating PRES at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % p % array
               ndims = 2
            else if (index(field % field, 'PMSL') /= 0) then
write(0,*) 'Interpolating PMSL'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % pmsl % array
               ndims = 1
            else if (index(field % field, 'PSFC') /= 0) then
write(0,*) 'Interpolating PSFC'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % psfc % array
               ndims = 1
            else if (index(field % field, 'SOILHGT') /= 0) then
write(0,*) 'Interpolating SOILHGT'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % soilz % array
               ndims = 1
            else if (index(field % field, 'SM000010') /= 0) then
write(0,*) 'Interpolating SM000010'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 1
               ndims = 2
               fg % dzs_fg % array(k,:) = 10.
               fg % zs_fg  % array(k,:) = 10.
            else if (index(field % field, 'SM010200') /= 0) then
write(0,*) 'Interpolating SM010200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-10.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'SM010040') /= 0) then
write(0,*) 'Interpolating SM010040'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 40.-10.
               fg % zs_fg  % array(k,:) = 40.
            else if (index(field % field, 'SM040100') /= 0) then
write(0,*) 'Interpolating SM040100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 3
               ndims = 2
               fg % dzs_fg % array(k,:) = 100.-40.
               fg % zs_fg  % array(k,:) = 100.
            else if (index(field % field, 'SM100200') /= 0) then
write(0,*) 'Interpolating SM100200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 4
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-100.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'ST000010') /= 0) then
write(0,*) 'Interpolating ST000010'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 1
               ndims = 2
               fg % dzs_fg % array(k,:) = 10.
               fg % zs_fg  % array(k,:) = 10.
            else if (index(field % field, 'ST010200') /= 0) then
write(0,*) 'Interpolating ST010200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-10.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'ST010040') /= 0) then
write(0,*) 'Interpolating ST010040'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 40.-10.               
               fg % zs_fg  % array(k,:) = 40.
            else if (index(field % field, 'ST040100') /= 0) then
write(0,*) 'Interpolating ST040100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 3
               ndims = 2
               fg % dzs_fg % array(k,:) = 100.-40.             
               fg % zs_fg  % array(k,:) = 100.
            else if (index(field % field, 'ST100200') /= 0) then
write(0,*) 'Interpolating ST100200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 4
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-100.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'SNOW') /= 0) then
write(0,*) 'Interpolating SNOW'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0

               masked = 0
               fillval = 0.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % snow % array
               ndims = 1
            else if (index(field % field, 'SEAICE') /= 0) then
write(0,*) 'Interpolating SEAICE'

               !interp_list(1) = SIXTEEN_POINT
               interp_list(1) = FOUR_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 1.0
               masked = 1
               fillval = 0.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % xice % array
               ndims = 1
            else if (index(field % field, 'SKINTEMP') /= 0) then
write(0,*) 'Interpolating SKINTEMP'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % skintemp % array
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
                  end if
                  if (ndims == 1) then
                     destField1d(i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                  else if (ndims == 2) then
                     destField2d(k,i) = interp_sequence(x, y, 1, rslab, -2, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
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

      call read_met_close()

      ! Fix for isobaric data
      if (minval(fg % p % array(:,:)) == 0.0) then
         write(0,*) 'Setting pressure field for isobaric data'
         do k=1,config_nfglevels
            if (vert_level(k) /= 200100.0) then
               fg % p % array(k,:) = vert_level(k)
            else
               fg % p % array(k,:) = fg % psfc % array(:)
            end if
         end do
      end if

      ! Set SST based on SKINTEMP field if it wasn't found in input data
      if (minval(fg % sst % array) == 0.0 .and. maxval(fg % sst % array) == 0.0) then
         write(0,*) 'Setting SST from SKINTEMP'
         !where (grid % landmask % array == 0) fg % sst % array = fg % skintemp % array
         fg % sst % array = fg % skintemp % array
      end if

      ! Set SNOWC (snow-cover flag) based on SNOW
      fg % snowc % array(:) = 0.0
      where (fg % snow % array > 0.0) fg % snowc % array = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MGD CHECK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iCell=1,grid%nCells
   if (grid % landmask % array(iCell) == 1) then

      do k = 1, config_nfgsoillevels
         if (fg % st_fg % array(k,iCell) <= 0.0) write(0,*) 'Bad st_fg ', k, iCell
      enddo

      do k = 1, config_nfgsoillevels
         if (fg % sm_fg % array(k,iCell) <= 0.0) write(0,*) 'Bad sm_fg ', k, iCell
      enddo
      !LDF end.

   end if
end do
write(0,*) 'Done with soil consistency check'


      !
      ! Get SEAICE from a separate file
      !
      call read_met_init('SEAICE_FRACTIONAL', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading SEAICE_FRACTIONAL data'
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (index(field % field, 'SEAICE') /= 0) then

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

               if (index(field % field, 'SEAICE') /= 0) then
                  nInterpPoints = grid % nCells
                  latPoints => grid % latCell % array
                  lonPoints => grid % lonCell % array
                  destField1d => fg % xice % array
                  ndims = 1
               end if
   
               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0
   
               masked = 1
               fillval = 0.0
               msgval = 1.01
               mask_array => grid % landmask % array


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
      where (fg % sst % array < 271.0 .and. grid % landmask % array == 0) fg % xice % array = 1.0

      ! Limit XICE to values between 0 and 1. Although the input meteorological field is between 0.
      ! and 1., interpolation to the MPAS grid can yield values of XiCE less than 0. and greater
      ! than 1.:
      where (fg % xice % array < 0._RKIND) fg % xice % array = 0._RKIND
      where (fg % xice % array > 1._RKIND) fg % xice % array = 1._RKIND

      ! Set SEAICE (0/1 flag) based on XICE (fractional ice coverage)
      fg % seaice % array(:) = 0.0
      where (fg % xice % array >= 0.5) fg % seaice % array = 1.0


      !  
      ! Compute normal wind component and store in fg%u
      !  
      do iEdge=1,grid%nEdges
         do k=1,config_nfglevels
            fg % u % array(k,iEdge) = cos(grid % angleEdge % array(iEdge)) * fg % u % array(k,iEdge) &
                                    + sin(grid % angleEdge % array(iEdge)) * fg % v % array(k,iEdge)
         end do
      end do


      !  
      ! Vertically interpolate meteorological data
      !  
      allocate(sorted_arr(2,config_nfglevels))

      do iCell=1,grid%nCells

         ! T
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % t % array(k,iCell)
         end do
         call mpas_quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
!           state % theta_m % array(k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1)
            state % theta_m % array(k,iCell) = vertical_interp(target_z, config_nfglevels-1, &
                                      sorted_arr(:,1:config_nfglevels-1), order=1, extrap=1)
         end do


         ! RH
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % rh % array(k,iCell)
         end do
         call mpas_quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
!           state % scalars % array(state % index_qv,k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=0)
            state % scalars % array(state % index_qv,k,iCell) = vertical_interp(target_z, config_nfglevels-1, &
                                                       sorted_arr(:,1:config_nfglevels-1), order=1, extrap=1)
            diag % rh % array(k,iCell) = state % scalars % array(state % index_qv,k,iCell)
         end do


         ! GHT
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % z % array(k,iCell)
         end do
         call mpas_quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
!           fg % gfs_z % array(k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1)
            fg % gfs_z % array(k,iCell) = vertical_interp(target_z, config_nfglevels-1, &
                                 sorted_arr(:,1:config_nfglevels-1), order=1, extrap=1)
         end do


         ! PRESSURE
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
            if (vert_level(k) == 200100.0) then 
!NOSFC               sorted_arr(1,k) = fg % soilz % array(iCell)
               sorted_arr(1,k) = 99999.0
               sfc_k = k
            end if
            sorted_arr(2,k) = log(fg % p % array(k,iCell))
         end do
         call mpas_quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
!           diag % pressure % array(k,iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))
            diag % pressure % array(k,iCell) = exp(vertical_interp(target_z, config_nfglevels-1, &
                                         sorted_arr(:,1:config_nfglevels-1), order=1, extrap=1))
         end do


         ! PRESSURE
!         sorted_arr(:,:) = -999.0
!         do k=1,config_nfglevels
!            sorted_arr(1,k) = fg % z % array(k,iCell)
!            if (vert_level(k) == 200100.0) then 
!!NOSFC               sorted_arr(1,k) = fg % soilz % array(iCell)
!               sorted_arr(1,k) = 99999.0
!               sfc_k = k
!            end if
!            sorted_arr(2,k) = log(fg % p % array(k,iCell))
!         end do
!         call mpas_quicksort(config_nfglevels, sorted_arr)
!         do k=1,grid%nVertLevels+1
!            target_z = grid % zgrid % array(k,iCell)
!            fg % gfs_p % array(k,iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))
!         end do

      end do


      do iEdge=1,grid%nEdges

         ! U
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = 0.5 * (fg % z % array(k,cellsOnEdge(1,iEdge)) + fg % z % array(k,cellsOnEdge(2,iEdge)))
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 0.5 * (fg % soilz % array(cellsOnEdge(1,iEdge)) + fg % soilz % array(cellsOnEdge(2,iEdge)))
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % u % array(k,iEdge)
         end do
         call mpas_quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.25 * (grid % zgrid % array(k,cellsOnEdge(1,iEdge)) + grid % zgrid % array(k+1,cellsOnEdge(1,iEdge)) + grid % zgrid % array(k,cellsOnEdge(2,iEdge)) + grid % zgrid % array(k+1,cellsOnEdge(2,iEdge)))
!           state % u % array(k,iEdge) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=0)
            state % u % array(k,iEdge) = vertical_interp(target_z, config_nfglevels-1, & 
                                sorted_arr(:,1:config_nfglevels-1), order=1, extrap=1)
         end do

      end do


      !
      ! Reconstruct zonal and meridional winds for diagnostic puposes:
      !
      call mpas_rbf_interp_initialize(grid)
      call mpas_init_reconstruct(grid)
      call mpas_reconstruct(grid, state % u % array,                        &
                            diag % uReconstructX % array,                   &
                            diag % uReconstructY % array,                   &
                            diag % uReconstructZ % array,                   &
                            diag % uReconstructZonal % array,               &
                            diag % uReconstructMeridional % array           &
                           )
   

      !
      ! Adjust surface pressure for difference in topography
      !
      do sfc_k=1,config_nfglevels
         if (vert_level(sfc_k) == 200100.) exit
      end do 
      do iCell=1,grid%nCells

         ! We need to extrapolate
            sorted_arr(:,:) = -999.0
            do k=1,config_nfglevels
               sorted_arr(1,k) = fg % z % array(k,iCell)
               if (vert_level(k) == 200100.0) then 
!NOSFC                  sorted_arr(1,k) = fg % soilz % array(iCell)
                  sorted_arr(1,k) = 99999.0
               end if
               sorted_arr(2,k) = log(fg % p % array(k,iCell))
            end do
            call mpas_quicksort(config_nfglevels, sorted_arr)
            target_z = grid % zgrid % array(1,iCell)
            fg % psfc % array(iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))

      end do

      deallocate(sorted_arr)



      !
      ! Diagnose fields needed in initial conditions file (u, w, rho, theta, scalars)
      ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
      !
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels

            ! QV
            es = 6.112 * exp((17.27*(state % theta_m % array(k,iCell) - 273.16))/(state % theta_m % array(k,iCell) - 35.86))
            rs = 0.622 * es * 100. / (diag % pressure % array(k,iCell) - es * 100.)
            scalars(state % index_qv,k,iCell) = 0.01 * rs * scalars(state % index_qv,k,iCell)

            ! PI
            p(k,iCell) = (diag % pressure % array(k,iCell) / p0) ** (rgas / cp)

            ! THETA - can compute this using PI instead
!            t(k,iCell) = t(k,iCell) / p(k,iCell)
            t(k,iCell) = t(k,iCell) * (p0 / diag % pressure % array(k,iCell)) ** (rgas / cp)

            ! RHO_ZZ
            rho_zz(k,iCell) = diag % pressure % array(k,iCell) / rgas / (p(k,iCell) * t(k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / (1.0 + scalars(state % index_qv,k,iCell))
         end do
      end do


      !
      ! Calculation of the initial precipitable water:
      ! 
      do iCell = 1,grid%nCells
         diag_physics%precipw%array(iCell) = 0.0
         do k = 1,grid%nVertLevels
            diag_physics%precipw%array(iCell) = diag_physics%precipw%array(iCell) &
                         + rho_zz(k,iCell)*scalars(state%index_qv,k,iCell)*(zgrid(k+1,iCell)-zgrid(k,iCell))
         enddo
      enddo

      !
      ! Reference state based on a dry isothermal atmosphere
      !
      do iCell=1,grid % nCells
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

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            pp(k,iCell) = diag % pressure % array(k,iCell) - ppb(k,iCell) 
            rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)
         end do
      end do

      do iCell=1,grid%nCells
         k = 1
         rho_zz(k,iCell) = ((diag % pressure % array(k,iCell) / p0)**(cv / cp)) * (p0 / rgas) / (t(k,iCell)*(1.0 + 1.61*scalars(state % index_qv,k,iCell)))
         rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         do k=2,grid % nVertLevels
            it = 0
            p_check = 2.0 * 0.0001
            do while ( (it < 30) .and. (p_check > 0.0001) )

               p_check = pp(k,iCell)
               dz = (zgrid(k,iCell) - zgrid(k-1,iCell))
               pp(k,iCell) = pp(k-1,iCell) - 0.5 * (rr(k,iCell) + rr(k-1,iCell))*gravity*dz &
                                           - 0.5 * (rho_zz(k,iCell)*scalars(state % index_qv,k,iCell) + rho_zz(k-1,iCell)*scalars(state % index_qv,k-1,iCell))*gravity*dz
               diag % pressure % array(k,iCell) = pp(k,iCell) + ppb(k,iCell)
               p(k,iCell) = (diag % pressure % array(k,iCell) / p0) ** (rgas / cp)
               rho_zz(k,iCell) = diag % pressure % array(k,iCell) / rgas / (p(k,iCell)*t(k,iCell)*(1.0 + 1.61*scalars(state % index_qv,k,iCell)))
               rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

               p_check = abs(p_check - pp(k,iCell))
                
               it = it + 1
            end do
         end do
      end do

      ! Compute theta_m and rho-tilde
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            t(k,iCell) = t(k,iCell) * (1.0 + 1.61*scalars(state % index_qv,k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / zz(k,iCell)
            rb(k,iCell) = rb(k,iCell) / zz(k,iCell)
         end do
      end do

      do iEdge=1,grid%nEdges
         do k=1,grid%nVertLevels
            diag % ru % array(k,iEdge) = state % u % array(k,iEdge) * 0.5*(state % rho_zz % array(k,cellsOnEdge(1,iEdge)) + state % rho_zz % array(k,cellsOnEdge(2,iEdge)))
         end do
      end do


      diag % rw % array = 0.
      state % w % array = 0.
      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*diag % ru % array(k,iEdge)+fzp(k)*diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,grid%nCells
         do k=2,grid%nVertLevels
            state % w % array(k,iCell) = diag % rw % array(k,iCell)     &
                                       / (fzp(k) * state % rho_zz % array(k-1,iCell) + fzm(k) * state % rho_zz % array(k,iCell))
         end do
      end do
   
      deallocate(vert_level)

     
      ! Calculate surface pressure (This is an ad-hoc calculation. The actual surface pressure is actually re-calculated at
      !the top of the subroutine MPAS_to_physics in ../core_atmos_physics/mpas_atmphys_interface_nhyd.F
      do iCell=1,grid%nCells
         diag % surface_pressure % array(iCell) = 0.5*gravity/rdzw(1)                                              &
                                                * (1.25* rho_zz(1,iCell) * (1. + scalars(state % index_qv, 1, iCell))  &
                                                -  0.25* rho_zz(2,iCell) * (1. + scalars(state % index_qv, 2, iCell)))
         diag % surface_pressure % array(iCell) = diag % surface_pressure % array(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do


      end if    ! config_met_interp

   end subroutine init_atm_test_case_gfs


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
      
      !input and output arguments:
      type(mesh_type),intent(inout):: mesh
      type (fg_type), intent(inout) :: fg
      
      !local variables:
      integer:: iCell,iMonth,iSoil
      
      !---------------------------------------------------------------------------------------------
      
      !initialization of surface input variables that are not needed if we run the current set of
      !idealized test cases:
      
      
      do iCell = 1, mesh % nCells
      
         !terrain,soil type, and vegetation:
         mesh % ter      % array(iCell) = 0.
         fg % xice       % array(iCell) = 0.
         mesh % landmask % array(iCell) = 0
         mesh % lu_index  % array(iCell) = 0
         mesh % soilcat_top % array(iCell) = 0
         mesh % shdmin   % array(iCell) = 0.
         mesh % shdmax   % array(iCell) = 0.
         fg % vegfra   % array(iCell) = 0.
         fg % sfc_albbck % array(iCell) = 0.
         fg % xland % array(iCell) = 0.
         fg % seaice % array(iCell) = 0.
      
         !snow coverage:
         fg % snow     % array(iCell) = 0.
         fg % snowc    % array(iCell) = 0.
         mesh % snoalb % array(iCell) = 0.08
         fg % snowh % array(iCell) = 0.
      
         !surface and sea-surface temperatures:
         fg % skintemp % array(iCell) = 288.0
         fg % sst      % array(iCell) = 288.0
      
         !soil layers:
         fg % tmn % array(iCell) = 288.0
         do iSoil = 1, mesh % nSoilLevels
            fg % tslb % array(iSoil,iCell)   = 288.0
            fg % smcrel % array(iSoil,iCell) =   0.0
            fg % sh2o   % array(iSoil,iCell) =   0.0
            fg % smois  % array(iSoil,iCell) =   0.0
            fg % dzs    % array(iSoil,iCell) =   0.0
         enddo
      
         !monthly climatological surface albedo and greeness fraction:
         do iMonth = 1, mesh % nMonths
            mesh % albedo12m % array(iMonth,iCell) = 0.08
            mesh % greenfrac % array(iMonth,iCell) = 0.
         enddo
      
      enddo
   
   end subroutine physics_idealized_init
   
   
   subroutine decouple_variables(grid, state, diag)

      implicit none

      type (mesh_type), intent(in) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag

      integer :: iCell, iEdge, k

      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, rdzw
      real (kind=RKIND), dimension(:,:), pointer :: zz, pp, ppb
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      cellsOnEdge => grid % CellsOnEdge % array
      fzp => grid % fzm % array
      fzp => grid % fzp % array
      rdzw => grid % rdzw % array
      zz => grid % zz % array

      pp  => diag % pressure_p % array
      ppb => diag % pressure_base % array

      scalars => state % scalars % array
 
     
      ! Compute surface pressure
      do iCell=1,grid%nCells
         diag % surface_pressure % array(iCell) = 0.5*gravity/rdzw(1)                                        &
                                                * (1.25* state % rho_zz % array(1,iCell) * (1. + scalars(state % index_qv, 1, iCell))  &
                                                -  0.25* state % rho_zz % array(2,iCell) * (1. + scalars(state % index_qv, 2, iCell)))
         diag % surface_pressure % array(iCell) = diag % surface_pressure % array(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do


      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do

   end subroutine decouple_variables


end module init_atm_test_cases

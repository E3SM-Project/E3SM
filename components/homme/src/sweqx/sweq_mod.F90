#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module sweq_mod
contains
  subroutine sweq(elem,fvm,edge1,edge2,edge3,red,par,ithr,nets,nete)
    !-----------------
    use kinds, only : real_kind, longdouble_kind
    !-----------------
    use parallel_mod, only : parallel_t, syncmp, abortmp
    !-----------------
    use thread_mod, only : nthreads
    !-----------------
    use hybrid_mod, only : hybrid_t, hybrid_create
    !-----------------
    use time_mod, only : timelevel_t , tstep, secpday, time_at, nmax, timelevel_update, timelevel_init
    !-----------------
    use derivative_mod, only : derivative_t, derivinit, allocate_subcell_integration_matrix
    !-----------------
    use dimensions_mod, only : np, nlev, npsq, npsq, nelemd, nvar, nc, ntrac
    !-----------------
    use shallow_water_mod, only : tc1_init_state, tc2_init_state, tc5_init_state, tc6_init_state, tc5_invariants, &
         tc8_init_state, vortex_init_state, vortex_errors, sj1_init_state, tc6_errors, &
         tc1_errors, tc2_errors, tc5_errors, sweq_invariants, swirl_init_state, swirl_errors, sj1_errors, tc1_velocity,&
         toy_chemistry_forcing
    !-----------------
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
    !-----------------
    use global_norms_mod, only : test_global_integral, print_cfl
    !-----------------
    use quadrature_mod, only : quadrature_t, gausslobatto
    !-----------------
    use edgetype_mod, only : EdgeBuffer_t
    ! ----------------
    use reduction_mod, only : ReductionBuffer_ordered_1d_t
    !-----------------
    use element_mod, only : element_t
    !-----------------
    use state_mod, only : printstate
    !-----------------
    use solver_mod, only : solver_test
    !-----------------
    use cg_mod, only : cg_t, cg_create
    !-----------------
    use restart_io_mod, only : readrestart, writerestart
    !-----------------
    use advance_mod, only : advance_nonstag
    !-----------------
#ifdef TRILINOS
    use implicit_mod, only : advance_imp_nonstag
    !-----------------
    use, intrinsic :: iso_c_binding 
    !-----------------
    use derived_type_mod ,only : derived_type, initialize


    use precon_type_mod ,only : precon_type, init_precon
    use solver_mod, only : solver_test, solver_test_ml
#endif
    !-----------------
    use control_mod, only : integration, filter_mu, filter_type, transfer_type, debug_level,  &
         restartfreq, statefreq, runtype, topology,   &
         test_case, sub_case, qsplit, nu, nu_s, limiter_option, hypervis_subcycle, test_cfldep, g_sw_output, &
         tstep_type, toy_chemistry
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    use bndry_mod, only : compute_ghost_corner_orientation, &
         sort_neighbor_buffer_mapping
    use checksum_mod, only : test_ghost, test_bilin_phys2gll

    use fvm_control_volume_mod, only : fvm_struct
    
    use fvm_mod, only : fvm_init2,fvm_init3    
    use fvm_bsp_mod, only: fvm_init_tracer

    use reduction_mod, only : parallelmax
    use mesh_mod, only : MeshUseMeshFile
    use unit_tests_mod, only : test_ibyp, test_edge_flux, &
                              test_sub_integration,       &
                              test_subcell_dss_fluxes,    &
                              test_subcell_div_fluxes,    &
                              test_subcell_Laplace_fluxes, &
                              test_subcell_div_fluxes_again, &
                              test_subcell_dss_fluxes_again, &
                              test_subcell_Laplace_fluxes_again  ! dont remove

    
    implicit none

    integer, parameter :: facs = 4            ! starting face number to print
    integer, parameter :: face = 4            ! ending  face number to print
    type (element_t), intent(inout) :: elem(:)
    type (fvm_struct), intent(inout) :: fvm(:)
    
    type (EdgeBuffer_t), intent(in)             :: edge1 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(in)             :: edge2 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(inout)             :: edge3 ! edge buffer entity             (shared)
    type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer               (shared)
    type (parallel_t), intent(in)               :: par   ! distributed parallel structure (shared)
    integer, intent(in)                         :: ithr  ! thread number                  (private)
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)

    ! ==================================
    ! Local thread (private) memory
    ! ==================================

    real (kind=real_kind)       :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn          ! dynamics timestep
    real (kind=real_kind) :: dt_tracers      ! tracer timestep

    real (kind=real_kind)       :: pmean           ! mean geopotential
    type (derivative_t)         :: deriv           ! derivative struct
    type (TimeLevel_t)          :: tl              ! time level struct
    type (hybrid_t)             :: hybrid
    type (quadrature_t)         :: gll,gs          ! gauss-lobatto and gauss wts and pts

    type (quadrature_t)   :: gp           ! quadratures on velocity and pressure grids

#ifdef TRILINOS
    integer :: lenx
    real (c_double) ,allocatable ,dimension(:) :: xstate
! state_object is a derived data type passed thru noxinit as a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    type(derived_type) ,target          :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    type(derived_type) ,target          :: jac_object
    type(derived_type) ,pointer        :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac

    type (element_t)  :: pc_elem(size(elem))
    type (element_t)  :: jac_elem(size(elem))


    real (kind=real_kind), dimension(np,np)  :: utemp1,utemp2




#endif

    integer :: simday

    integer :: point
    integer :: i,j,iptr
    integer :: it,ie,k
    integer :: ntmp
    integer :: nm1,n0,np1
    integer :: nstep

    real*8  :: tot_iter
    logical, parameter :: Debug = .FALSE.

  real (kind=longdouble_kind)                    :: fvm_corners(nc+1)
  real(kind=longdouble_kind)                     :: fvm_points(nc)     ! fvm cell centers on reference element
  
  real (kind=real_kind)                          :: xtmp
  real (kind=real_kind)                          :: maxcflx, maxcfly  


#ifdef TRILINOS
  interface 
     subroutine noxinit(vectorSize, vector, comm, v_container, p_container) &
          bind(C,name='noxinit')
       use ,intrinsic :: iso_c_binding
       integer(c_int)                :: vectorSize,comm
       real(c_double)  ,dimension(*) :: vector
       type(c_ptr)                   :: v_container
       type(c_ptr)                   :: p_container  !precon ptr
     end subroutine noxinit

    subroutine noxfinish() bind(C,name='noxfinish')
      use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
    end subroutine noxfinish

  end interface
#endif

#if 0
     call allocate_subcell_integration_matrix(np,  6)
#endif




    if(Debug) print *,'homme: point #1'

    ! ==========================
    ! begin executable code
    ! ==========================
    call t_startf('sweq')
    hybrid = hybrid_create(par,ithr,NThreads)
    simday=0
    if (topology == "cube") then
       call test_global_integral(elem,hybrid,nets,nete)

       dtnu = 2.0d0*tstep*max(nu,nu_s)/hypervis_subcycle
       call print_cfl(elem,hybrid,nets,nete,dtnu)

       if (.not. MeshUseMeshFile) then
          ! orientation code assumes only one corner element neighbor
          ! orientation algorithm only works for cubed-sphere meshes
          call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
          call test_ghost(hybrid,elem,nets,nete)
       endif
       call sort_neighbor_buffer_mapping(hybrid%par,elem,nets,nete)
    end if

    if(Debug) print *,'homme: point #2'
    ! ==================================
    ! Initialize derivative structure
    ! ==================================

    ! Initialize derivative structure
    ! fvm nodes are equally spaced in alpha/beta
    ! HOMME with equ-angular gnomonic projection maps alpha/beta space
    ! to the reference element via simple scale + translation
    ! thus, fvm nodes in reference element [-1,1] are a tensor product of
    ! array 'fvm_nodes(:)' computed below:
    xtmp=nc 
    do i=1,nc+1
      fvm_corners(i)= 2*(i-1)/xtmp - 1
    end do
    do i=1,nc
       fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    end do
    call derivinit(deriv,fvm_corners,fvm_points)
    if (ntrac>0) then
       call fvm_init2(elem,fvm,hybrid,nets,nete,tl)
       call test_bilin_phys2gll(elem,fvm,hybrid,nets,nete)
    endif

    ! ========================================
    ! Initialize velocity and pressure grid
    ! quadrature points...
    ! ========================================

    
    gp =gausslobatto(np)

!   some test code
#if 0
    if (hybrid%masterthread) print *,'running CG solver test'
    call solver_test(elem,edge1,red,hybrid,deriv,nets,nete)
#ifdef TRILINOS
    call solver_test_ml(elem,edge1,red,hybrid,deriv,nets,nete)
#endif
    if (hybrid%masterthread) print *,'running global integration-by-parts checks'
    call test_ibyp(elem,hybrid,nets,nete)
    if (hybrid%masterthread) print *,'running element divergence/edge flux checks'
    call test_edge_flux(elem,deriv,nets,nete)
    call test_sub_integration(elem,deriv,nets,nete)
    call test_subcell_dss_fluxes(elem,deriv,nets,nete)
    call test_subcell_div_fluxes(elem,deriv,nets,nete)
    call test_subcell_Laplace_fluxes(elem,deriv,nets,nete)
    call test_subcell_div_fluxes_again(elem,deriv,nets,nete)
!   call test_subcell_dss_fluxes_again(elem,deriv,nets,nete)
    call test_subcell_Laplace_fluxes_again(elem,deriv,nets,nete)
    stop
#endif



    ! =================================
    ! Slow start leapfrog...
    ! =================================
    ! mt 4/2007:  added better leapfrog bootstrap procedure (due to J. Tribbia)
    ! original algorithm looses 2 digits in the Energy because of first 2 timesteps.
    !
    !   original algorithm                                        better version
    !
    !  nm = u(0)                                                  nm = u(0)
    !  n0 = u(0)                                                  n0 = u(0)           
    !  np = undefined                                             np = undefined
    ! 
    !  call advance (dt/2)    np = nm + 2*(dt/2)*n0               call advance(dt/4)   np = nm + 2*dt/4*n0
    !
    !  nm = u(0)                                                  nm=u(0)
    !  n0 = u(0)                                                  n0=u(0)             
    !  np = u(dt)                                                 np = u(dt/2)
    !
    !  call timelevel_update(forward)                             call timelevel_updatate(forward)
    ! 
    !  nm = u(0)                                                  nm=u(0)             
    !  n0 = u(dt)                                                 n0=u(dt/2)          
    !  np = undefined                                             np = undefiend
    !
    !  call advance (dt)                                          call advance (dt/2)
    !
    !  nm = u(0)                                                  nm=u(0)
    !  n0 = u(dt)                                                 n0=u(dt/2)          
    !  np = u(2dt)                                                np=u(dt)
    !
    !  call timelevel_updatate(leapfrog)                          call timelevel_updatate(forward)
    !  
    !  nm = u(dt)                                                 nm=u(0)
    !  n0 = u(2dt)                                                n0=u(dt)            
    !  np = undefined                                             np=undefined
    !
    !                                                             call advance (dt) 
    !
    !                                                             nm=u(0)
    !                                                             n0=u(dt)            
    !                                                             np=u(2dt)
    !
    !                                                             call timelevel_updatate(leapfrog)
    !
    !                                                             nm=u(dt)
    !                                                             n0=u(2dt)           
    !                                                             np=undefined
    !

    call TimeLevel_init(tl) !note: hard wired in tc1_init_state so if chg here then chg there too

    if(Debug) print *,'homme: point #5'
    ! =================================================================
    ! Initialize geopotential and velocity for different test cases...
    ! =================================================================

    if (topology == "cube") then
       if (runtype .eq. 1) then 
          if (hybrid%masterthread) then
             print *,'runtype: RESTART of Shallow Water equations'
          end if
          if (test_case(1:5) == "swtc1") then
             if (hybrid%masterthread) print *,"Restarting swtc1..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc1_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc2") then
             if (hybrid%masterthread) print *,"Restarting swtc2..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc2_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc5") then
             if (hybrid%masterthread) print *,"Restarting swtc5..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc5_init_state(elem, nets,nete,pmean,deriv)
             call tc5_invariants(elem,90,tl,pmean,edge2,deriv,hybrid,nets,nete)
             call tc5_errors(elem,7,tl,pmean,"ref_tc5_imp",simday,hybrid,nets,nete,par)
          else if (test_case(1:5) == "swtc6") then
             if (hybrid%masterthread) print *,"Restarting swtc6..."
             call tc6_init_state(elem,nets,nete,pmean)
             simday=0
             call tc6_errors(elem,7,tl,pmean,"ref_tc6_imp",simday,hybrid,nets,nete,par)
          else if (test_case(1:5) == "swtc8") then
             if (hybrid%masterthread) print *,"Restarting swtc8..."
             call tc8_init_state(elem,nets,nete,hybrid,pmean)
          else if (test_case(1:6) == "vortex") then
             if (hybrid%masterthread) print *,"Restarting vortex..."
             call vortex_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swirl") then
             if (hybrid%masterthread) print *,"Restarting swirl..."
             call swirl_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swsj1") then
             if (hybrid%masterthread) print *,"Restarting swsj1..."
             call sj1_init_state(elem,nets,nete,hybrid,pmean,deriv)
             simday=0
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if
          !============================
          ! Read in the restarted state 
          !============================
          call ReadRestart(elem,ithr,nete,nets,tl)
          !================================================
          ! Print out the state variables 
          !================================================
          !DBG print *,'homme: right after ReadRestart pmean is: ',pmean

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,-1)
          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       else 
          if (hybrid%masterthread) then
             print *,'runtype: INITIAL of Shallow Water equations'
          end if
          dt = tstep/4
          if (test_case(1:5) == "swtc1") then
             if (hybrid%masterthread) print *,"initializing swtc1..."
             call tc1_init_state(elem,nets,nete,pmean)
             call tc1_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          else if (test_case(1:5) == "swtc2") then
             if (hybrid%masterthread) print *,"initializing swtc2..."
             call tc2_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc5") then
             if (hybrid%masterthread) print *,"initializing swtc5..."
             call tc5_init_state(elem,nets,nete,pmean,deriv)
             call tc5_invariants(elem,90,tl,pmean,edge2,deriv,hybrid,nets,nete)
             call tc5_errors(elem,7,tl,pmean,"ref_tc5_imp",simday,hybrid,nets,nete,par)
             
             if (ntrac>0) then
             do ie=nets,nete
               call fvm_init_tracer(fvm(ie),tl)
             end do
             call fvm_init3(elem,fvm,hybrid,nets,nete,tl%n0)
             if (hybrid%masterthread) print *,"initializing fvm tracers for swtc5..."
             endif

             
          else if (test_case(1:5) == "swtc6") then
             if (hybrid%masterthread)  print *,"initializing swtc6..."
             call tc6_init_state(elem,nets,nete,pmean)
             simday=0
             call tc6_errors(elem,7,tl,pmean,"ref_tc6_imp",simday,hybrid,nets,nete,par)
          else if (test_case(1:5) == "swtc8") then
             if (hybrid%masterthread) print *,"initializing swtc8..."
             call tc8_init_state(elem,nets,nete,hybrid,pmean)
          else if (test_case(1:6) == "vortex") then
             if (hybrid%masterthread) print *,"initializing vortex..."
             call vortex_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swirl") then
             if (hybrid%masterthread) print *,"initializing swirl..."
             call swirl_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swsj1") then
             if (hybrid%masterthread) print *,"initializing swsj1..."
             call sj1_init_state(elem,nets,nete,hybrid,pmean,deriv)
             simday=0
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if  ! test case init calls

          ! ==============================================
          ! Output initial picture of geopotential...
          ! ============================================== 
          ! I/O routines are not thread safe
          ! we need to stop threads and restart after I/O
          if (NThreads>1) then
             call abortmp('Shallow water model I/O can not be called with threads')
          endif
#ifdef PIO_INTERP
	  call interp_movie_init(elem,par,tl=tl)
          call interp_movie_output(elem,tl, par, pmean, fvm)     
#else
          call shal_movie_init(elem,hybrid,fvm)
          call shal_movie_output(elem,tl, hybrid, pmean, nets, nete,deriv,fvm)
#endif
          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,-1)
          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
          if(Debug) print *,'homme: point #6'

          if (toy_chemistry==1) call toy_chemistry_forcing(elem,nets,nete,tl,dt)
          if(Debug) print *,'homme: point #7'


          ! =================================
          ! Call advance
          ! =================================
          call advance_nonstag(elem, edge2, edge3,    deriv,   hybrid, &
                  dt,    pmean,     tl,  nets,   nete)
          call TimeLevel_update(tl,"forward")
          if(Debug) print *,'homme: point #8'

          ! ============================================================
          ! Shallow Water Test Case 1 outputs
          ! L1,L2,Linf error norms every timestep
          ! ============================================================

          if (test_case(1:5) == "swtc1") then
             call tc1_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          else if (test_case(1:5) == "swtc2") then
             call tc2_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #9'

          ! ===============================================================
          ! Print Min/Max/Sum of State vars after first timestep
          ! ===============================================================
          if (integration /= "full_imp") then

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,-1)

          endif  ! if time step taken
          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       endif  ! if initial run 
    end if ! if topology == "cube"

    ! reset timestep counter.  New more accurate leapfrog bootstrap routine takes
    ! one extra timestep to get started.  dont count that timestep, otherwise
    ! times will all be off by tstep. Also, full_imp is just starting.

    tl%nstep=0 

    ! =========================================
    ! Set up for leapfrog time integration...
    ! =========================================

    dt = tstep/2

    if (integration == "full_imp") then
      if (hybrid%masterthread) print *,'initializing Trilinos solver info'

#ifdef TRILINOS
      lenx=np*np*nlev*nvar*(nete-nets+1)
      allocate(xstate(lenx))
      xstate(:) = 0
       call initialize(state_object, lenx, elem, pmean,edge1,edge2, edge3, &
        hybrid, deriv, dt, tl, nets, nete)

    
       pc_elem=elem
       jac_elem=elem


       call initialize(pre_object, lenx, pc_elem, pmean, edge1,edge2,edge3, &
        hybrid, deriv, dt, tl, nets, nete)

       call initialize(jac_object, lenx, jac_elem, pmean, edge1,edge2,edge3, &
        hybrid, deriv, dt, tl, nets, nete)


       fptr => state_object
       c_ptr_to_object =  c_loc(fptr)
       pptr => pre_object
       c_ptr_to_pre =  c_loc(pptr)
       jptr => jac_object
       c_ptr_to_jac =  c_loc(jptr)





        call noxinit(size(xstate), xstate, 1, c_ptr_to_object, c_ptr_to_pre)
#endif

    end if
    
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (ithr==0) then
       call syncmp(par)
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

    ! ===================================
    ! Main timestepping loop
    ! ===================================
    tot_iter=0.0       
    do while(tl%nstep<nmax)

       ! =================================
       ! Call advance
       ! =================================
       if(Debug) print *,'homme: point #12'

       if      (integration == "explicit") then
          call advance_nonstag( elem,  edge2, edge3, deriv, hybrid, &
               dt   , pmean, tl   , nets, nete)
       else if (integration == "full_imp") then
            dt=tstep
#ifdef TRILINOS

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1


        do ie=nets,nete
         do k=1,nlev
          do i=1,np
           do j=1,np
           utemp1(i,j)= elem(ie)%D(i,j,1,1)*elem(ie)%state%v(i,j,1,k,np1) + elem(ie)%D(i,j,1,2)*elem(ie)%state%v(i,j,2,k,np1)
           utemp2(i,j)= elem(ie)%D(i,j,2,1)*elem(ie)%state%v(i,j,1,k,np1) + elem(ie)%D(i,j,2,2)*elem(ie)%state%v(i,j,2,k,np1)
           elem(ie)%state%v(i,j,1,k,np1)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,np1)=utemp2(i,j)

           utemp1(i,j)= elem(ie)%D(i,j,1,1)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%D(i,j,1,2)*elem(ie)%state%v(i,j,2,k,n0)
           utemp2(i,j)= elem(ie)%D(i,j,2,1)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%D(i,j,2,2)*elem(ie)%state%v(i,j,2,k,n0) 
           elem(ie)%state%v(i,j,1,k,n0)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,n0)=utemp2(i,j)

           utemp1(i,j)= elem(ie)%D(i,j,1,1)*elem(ie)%state%v(i,j,1,k,nm1) + elem(ie)%D(i,j,1,2)*elem(ie)%state%v(i,j,2,k,nm1)
           utemp2(i,j)= elem(ie)%D(i,j,2,1)*elem(ie)%state%v(i,j,1,k,nm1) + elem(ie)%D(i,j,2,2)*elem(ie)%state%v(i,j,2,k,nm1)

           elem(ie)%state%v(i,j,1,k,nm1)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,nm1)=utemp2(i,j)
! need  to  convert  xstate  to  latlon

           end do !nv
          end do !nv 
         end do !nlev 
        end do !ie 




            call advance_imp_nonstag(elem, edge1, edge2, edge3, red, deriv,  &
               cg, hybrid, dt, pmean, tl, nets, nete, xstate)

            ! TODO update with vortex and swirl possibly using set_prescribed_velocity

        do ie=nets,nete

         if (topology == "cube" .and. test_case=="swtc1") then
             do k=1,nlev
                elem(ie)%state%v(:,:,:,k,np1)=tc1_velocity(elem(ie)%spherep,elem(ie)%Dinv)
                elem(ie)%state%v(:,:,:,k,n0)=elem(ie)%state%v(:,:,:,k,np1)
                elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,nm1)
             end do 
         else
         do k=1,nlev
          do i=1,np
           do j=1,np
           utemp1(i,j)= elem(ie)%Dinv(i,j,1,1)*elem(ie)%state%v(i,j,1,k,np1) + elem(ie)%Dinv(i,j,1,2)*elem(ie)%state%v(i,j,2,k,np1)
           utemp2(i,j)= elem(ie)%Dinv(i,j,2,1)*elem(ie)%state%v(i,j,1,k,np1) + elem(ie)%Dinv(i,j,2,2)*elem(ie)%state%v(i,j,2,k,np1)

           elem(ie)%state%v(i,j,1,k,np1)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,np1)=utemp2(i,j)

           utemp1(i,j)= elem(ie)%Dinv(i,j,1,1)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%Dinv(i,j,1,2)*elem(ie)%state%v(i,j,2,k,n0)
           utemp2(i,j)= elem(ie)%Dinv(i,j,2,1)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%Dinv(i,j,2,2)*elem(ie)%state%v(i,j,2,k,n0)

           elem(ie)%state%v(i,j,1,k,n0)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,n0)=utemp2(i,j)

           utemp1(i,j)= elem(ie)%Dinv(i,j,1,1)*elem(ie)%state%v(i,j,1,k,nm1) + elem(ie)%Dinv(i,j,1,2)*elem(ie)%state%v(i,j,2,k,nm1)
           utemp2(i,j)= elem(ie)%Dinv(i,j,2,1)*elem(ie)%state%v(i,j,1,k,nm1) + elem(ie)%Dinv(i,j,2,2)*elem(ie)%state%v(i,j,2,k,nm1)

           elem(ie)%state%v(i,j,1,k,nm1)=utemp1(i,j)
           elem(ie)%state%v(i,j,2,k,nm1)=utemp2(i,j)
! need ! to ! convert ! xstate ! to ! contravariant 
           end do !np 
          end do !np 
         end do !nlev 
        end if
        end do !ie

#else
!           Check /utils/trilinos/README for more details
            call abortmp('Need to include -DTRILINOS at compile time to execute FI solver')
#endif
       end if
       if(Debug) print *,'homme: point #13'
!      if (hybrid%masterthread) print *,'post solve same ts'

       ! =================================
       ! update time level pointers
       ! =================================
        if (integration == "full_imp") then
           if (tstep_type == 13) then
              call TimeLevel_update(tl,"leapfrog") ! for BDF2
           else     
              call TimeLevel_update(tl,"forward") ! second order Crank Nicolson
           end if
        else
           if (tl%nstep==0) then
              call TimeLevel_update(tl,"forward")
              dt=dt*2    
           else
              call TimeLevel_update(tl,"leapfrog")
           endif
        end if 

       ! ============================================================
       ! Instrumentation alley:
       !
       ! Shallow Water Test Case output files
       ! ============================================================
#ifdef PIO_INTERP
        call interp_movie_output(elem,tl, par, pmean, fvm)
#else     
        call shal_movie_output(elem,tl, hybrid, pmean, nets, nete,deriv,fvm)
#endif
       ! ==================================================
       ! Shallow Water Test Cases:
       ! ==================================================
       if(Debug) print *,'homme: point #14'

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       if (test_case(1:5) == "swtc1") then

          ! ==================================================
          ! Shallow Water Test Case 1:  cosine bell
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call tc1_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          end if

       else if (test_case(1:5) == "swtc2") then
          if(Debug) print *,'homme: point #15'

          ! ==================================================
          ! Shallow Water Test Case 2:  cosine bell
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call tc2_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16'

       else if (test_case(1:5) == "swtc5") then

          ! ===============================================================
          ! Shallow Water Test Case 5: Rossby Haurwitz Waves
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ===============================================================

          if(Debug) print *,'homme: point #16.1'
          if (MODULO(tl%nstep,statefreq)==0) then
             call tc5_invariants(elem, 90, tl, pmean, edge2, deriv, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16.2'

          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
!          if (MODULO(tl%nstep,statefreq)==0) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call tc5_errors(elem,7,tl,pmean,"ref_tc5_imp",simday,hybrid,nets,nete,par)
          end if
          if(Debug) print *,'homme: point #16.3'

       else if (test_case(1:5) == "swtc6") then

          ! ===============================================================
          ! Shallow Water Test Case 6: Rossby Haurwitz Waves
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ===============================================================

          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call tc6_errors(elem,7,tl,pmean,"ref_tc6_imp",simday,hybrid,nets,nete,par)
          end if

       else if (test_case(1:5) == "swtc8") then
!!! RDL            call tc8_invariants(90, tl, pmean, edge2, deriv, hybrid, nets, nete)
       else if (test_case(1:6) == "vortex") then

          ! ===============================================================
          ! Shallow Water Test Case 9: vortex
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! detect day rollover
          ! ===============================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call vortex_errors(elem, 7, tl,  hybrid, nets, nete)
          end if
       else if (test_case(1:5) == "swirl") then

          if (MODULO(tl%nstep,statefreq)==0) then
             call swirl_errors(elem, 7, tl,  hybrid, nets, nete)
          end if
       else if (test_case(1:5) == "swirl") then

          if (MODULO(tl%nstep,statefreq)==0) then
             call swirl_errors(elem, 7, tl,  hybrid, nets, nete)
          end if
       else if (test_case(1:5) == "swsj1") then

          ! ===============================================================
          ! Shallow Water Test Case SJ:
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ===============================================================

          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if

       end if

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       if(Debug) print *,'homme: point #17'
       if (MODULO(tl%nstep,statefreq)==0 ) then 
          if(hybrid%masterthread) then
             print *,tl%nstep,"time=",Time_at(tl%nstep)/secpday," days"
          end if

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,-1)

          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       end if

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif

       ! ============================================================
       ! Write restart files if required
       ! ============================================================

       if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then
          call WriteRestart(elem, ithr,nets,nete,tl)
          if (par%masterproc) print *, "after restart"
       endif

    end do

     if (integration == "full_imp") then ! closeout trilinos assignments
#ifdef TRILINOS
       call noxfinish()
       deallocate(xstate)
#endif
     end if

    ! ======================================================
    ! compute and report times...
    ! ======================================================

    ! ==========================
    ! end of the hybrid program
    ! ==========================
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
    call t_stopf('sweq')
  end subroutine sweq

  subroutine sweq_rk(elem, edge1,edge2,edge3,red,par,ithr,nets,nete)
    !-----------------
    use kinds, only : real_kind
    !-----------------
    use physical_constants, only : dd_pi
    !-----------------
    use parallel_mod, only : parallel_t, syncmp
    !-----------------
    use thread_mod, only : nthreads
    !-----------------
    use hybrid_mod, only : hybrid_t, hybrid_create
    !-----------------
    use time_mod, only : timelevel_t , ndays, tstep, secpday, time_at, nmax, timelevel_update, timelevel_init
    !-----------------
    use derivative_mod, only : derivative_t, derivinit
    !-----------------
    use dimensions_mod, only :  np, nlev, npsq
    !-----------------
    use shallow_water_mod, only : tc1_init_state, tc2_init_state, tc5_init_state, &
         tc6_init_state, tc5_invariants, tc8_init_state, vortex_init_state, &
         vortex_errors, sj1_init_state, tc6_errors, tc1_errors, tc2_errors, &
         tc5_errors, sweq_invariants, swirl_init_state, swirl_errors, sj1_errors,&
         kmass_swirl, toy_chemistry_forcing
    !-----------------
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish

#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
    


#endif
    !-----------------
    use global_norms_mod, only : test_global_integral, print_cfl
    !-----------------
    use quadrature_mod, only : quadrature_t, gausslobatto
    !-----------------
    use edgetype_mod, only : EdgeBuffer_t
    ! ----------------
    use reduction_mod, only : ReductionBuffer_ordered_1d_t
    !-----------------
    use element_mod, only : element_t
    !-----------------
    use state_mod, only : printstate
    !-----------------
    use restart_io_mod, only : readrestart, writerestart
    !-----------------
    use advance_mod, only : advance_nonstag_rk
    !-----------------
    use control_mod, only : integration, debug_level,  &
         restartfreq, statefreq, runtype, topology, &
         rk_stage_user, test_case, sub_case, kmass, qsplit, nu, nu_s, limiter_option, &
         hypervis_subcycle, g_sw_output, toy_chemistry
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL

    use rk_mod, only     : RkInit
    use types_mod, only : rk_t
    
    implicit none

    integer, parameter :: facs = 4            ! starting face number to print
    integer, parameter :: face = 4            ! ending  face number to print
    type (element_t), intent(inout) :: elem(:)
    type (EdgeBuffer_t), intent(in)             :: edge1 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(in)             :: edge2 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(inout)             :: edge3 ! edge buffer entity             (shared)
    type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer               (shared)
    type (parallel_t), intent(in)               :: par   ! distributed parallel structure (shared)
    integer, intent(in)                         :: ithr  ! thread number                  (private)
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)

    ! ==================================
    ! Local thread (private) memory
    ! ==================================

    real (kind=real_kind)       :: dt,dt_rk        ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn          ! dynamics timestep
    real (kind=real_kind) :: dt_tracers      ! tracer timestep

    real (kind=real_kind)       :: pmean           ! mean geopotential
    type (derivative_t)         :: deriv           ! derivative struct
    type (TimeLevel_t)          :: tl              ! time level struct
    type (hybrid_t)             :: hybrid
    type (quadrature_t)         :: gll,gs          ! gauss-lobatto and gauss wts and pts


    type (quadrature_t)   :: gp           ! quadratures on velocity and pressure grids
    
    real (kind=real_kind) :: mindx,dt_gv

    integer  :: simday

    integer :: point
    integer :: i,j,iptr
    integer :: it,ie,k
    integer :: ntmp
    integer :: nm1,n0,np1
    integer :: nstep
    integer :: cfl

    type (rk_t) :: RungeKutta

    real*8  :: tot_iter
    logical, parameter :: Debug = .FALSE.


    if(Debug) print *,'e: point #1'
    ! ==========================
    ! begin executable code
    ! ==========================
    call t_startf('sweq')

    hybrid = hybrid_create(par,ithr,NThreads)

    if (topology == "cube") then
       call test_global_integral(elem,hybrid,nets,nete,mindx)
       dtnu = (tstep/rk_stage_user)*max(nu,nu_s)/hypervis_subcycle
       call print_cfl(elem,hybrid,nets,nete,dtnu)
    end if

    ! Find time-step to gravity wave speed
    ! 2012: broken because mindx=0.  also, should be updated to use true eigenvalue,
    ! not mindx.  
    dt_gv = (mindx/dd_pi)/300.0D0    
    dt = tstep ! this is the "user" time-step
    if (hybrid%masterthread) then
       print *,"dt grv = ", dt_gv
       print *,"dt user= ", dt
    endif
  
    if(Debug) print *,'homme: point #2'
    ! ==================================
    ! Initialize derivative structure
    ! ==================================

    call derivinit(deriv)

    ! ========================================
    ! Initialize velocity and pressure grid
    ! quadrature points...
    ! ========================================

    gp=gausslobatto(np)

    if(Debug) print *,'homme: point #4'
    call TimeLevel_init(tl)

    ! =================================================================
    ! Initialize geopotential and velocity for different test cases...
    ! =================================================================

    if (topology == "cube") then
       if (runtype .eq. 1) then 
          if (hybrid%masterthread) then
             print *,'runtype: RESTART of Shallow Water equations'
          end if
          if (test_case(1:5) == "swtc1") then
             if (hybrid%masterthread) print *,"Restarting swtc1..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc1_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc2") then
             if (hybrid%masterthread) print *,"Restarting swtc2..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc2_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc5") then
             if (hybrid%masterthread) print *,"Restarting swtc5..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call tc5_init_state(elem, nets,nete,pmean,deriv)
             call tc5_invariants(elem,90,tl,pmean,edge2,deriv,hybrid,nets,nete)
             call tc5_errors(elem,7, tl, pmean, "ref_tc5_imp", simday, hybrid, nets, nete,par)
          else if (test_case(1:5) == "swtc6") then
             if (hybrid%masterthread) print *,"Restarting swtc6..."
             call tc6_init_state(elem,nets,nete,pmean)
             simday=0
             call tc6_errors(elem,7, tl, pmean, "ref_tc6_imp", simday, hybrid, nets, nete,par)
          else if (test_case(1:5) == "swtc8") then
             if (hybrid%masterthread) print *,"Restarting swtc8..."
             call tc8_init_state(elem,nets,nete,hybrid,pmean)
          else if (test_case(1:6) == "vortex") then
             if (hybrid%masterthread) print *,"Restarting vortex..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call vortex_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swirl") then
             if (hybrid%masterthread) print *,"Restarting swirl..."
             !==================================================
             ! Recover the initial state for diagnostic purposes
             !==================================================
             call swirl_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swsj1") then
             if (hybrid%masterthread) print *,"Restarting swsj1..."
             call sj1_init_state(elem,nets,nete,hybrid,pmean,deriv)
             simday=0
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if
          !============================
          ! Read in the restarted state 
          !============================
          call ReadRestart(elem,ithr,nete,nets,tl)
          !================================================
          ! Print out the state variables 
          !================================================
          !DBG print *,'homme: right after ReadRestart pmean is: ',pmean

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,kmass)

          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       else 
          if (hybrid%masterthread) then
             print *,'runtype: INITIAL of Shallow Water equations'
          end if
          
          if (test_case(1:5) == "swtc1") then
             if (hybrid%masterthread) print *,"initializing swtc1..."
             call tc1_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc2") then
             if (hybrid%masterthread) print *,"initializing swtc2..."
             call tc2_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swtc5") then
             if (hybrid%masterthread) print *,"initializing swtc5..."
             call tc5_init_state(elem,nets,nete,pmean,deriv)
             call tc5_invariants(elem,90,tl,pmean,edge2,deriv,hybrid,nets,nete)
             call tc5_errors(elem,7, tl, pmean, "ref_tc5_imp", simday, hybrid, nets, nete,par)
          else if (test_case(1:5) == "swtc6") then
             if (hybrid%masterthread)  print *,"initializing swtc6..."
             call tc6_init_state(elem,nets,nete,pmean)
             simday=0
             call tc6_errors(elem,7, tl, pmean, "ref_tc6_imp", simday, hybrid, nets, nete,par)
          else if (test_case(1:5) == "swtc8") then
             if (hybrid%masterthread) print *,"initializing swtc8..."
             call tc8_init_state(elem,nets,nete,hybrid,pmean)
          else if (test_case(1:6) == "vortex") then
             if (hybrid%masterthread) print *,"initializing vortex..."
             call vortex_init_state(elem,nets,nete,pmean)
          else if (test_case(1:5) == "swirl") then
             if (hybrid%masterthread) print *,"initializing swirl..."
             call swirl_init_state(elem,nets,nete,pmean,hybrid,edge3)
          else if (test_case(1:5) == "swsj1") then
             if (hybrid%masterthread) print *,"initializing swsj1..."
             call sj1_init_state(elem,nets,nete,hybrid,pmean,deriv)
             simday=0
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if

          ! ==============================================
          ! Output initial picture of geopotential...
          ! ============================================== 
#ifdef PIO_INTERP
	  call interp_movie_init(elem,par,tl=tl)
          call interp_movie_output(elem,tl, par, pmean )
#else
	  call shal_movie_init(elem,hybrid)
          call shal_movie_output(elem,tl, hybrid, pmean, nets, nete,deriv)
#endif

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,kmass)

          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
          if(Debug) print *,'homme: point #6'

          ! ============================================================
          ! Shallow Water Test Case outputs
          ! L1,L2,Linf error norms every timestep
          ! ============================================================

          if (test_case(1:5) == "swtc1") then
             call tc1_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          else if (test_case(1:5) == "swtc2") then
             call tc2_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          else if (test_case(1:6) == "vortex") then
             call vortex_errors(elem, 7, tl, hybrid, nets, nete)
          else if (test_case(1:5) == "swirl") then
             call swirl_errors(elem, 7, tl, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #9'

          ! ===============================================================
          ! Print Min/Max/Sum of State vars after first timestep
          ! ===============================================================

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,kmass)

          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       endif  ! if initial run 
    end if ! if topology == "cube"

    ! reset timestep counter.  New more accurate leapfrog bootstrap routine takes
    ! one extra timestep to get started.  dont count that timestep, otherwise
    ! times will all be off by tstep. 
    tl%nstep=0 


    if(Debug) print *,'homme: point #11'

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (ithr==0) then
       call syncmp(par)
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

    ! ===================================
    ! Main timestepping loop
    ! ===================================
    
    do while(tl%nstep<nmax)  

       ! =================================
       ! Call advance
       ! =================================
       if (rk_stage_user > 0) then
          ! user specified number of stages.  has to be >= 2
          cfl = max(2,rk_stage_user)
          cfl = cfl - 1  
       else
          cfl = ceiling(dt/dt_gv + 1.0D0)
       endif
       dt_rk = dt

       call RkInit(cfl,RungeKutta)  ! number of stages: cfl+1

       if (toy_chemistry==1) call toy_chemistry_forcing(elem,nets,nete,tl,dt_rk)
       call advance_nonstag_rk(RungeKutta, elem,  edge2, edge3, deriv,  hybrid, &
            dt_rk   , pmean, tl   , nets, nete)
       
       call TimeLevel_update(tl,"leapfrog")       
       ! ============================================================
       ! Instrumentation alley:
       !
       ! Shallow Water Test Case output files
       ! ============================================================
#ifdef PIO_INTERP
       call interp_movie_output(elem, tl, par, pmean)
#else       
          call shal_movie_output(elem,tl, hybrid, pmean, nets, nete,deriv)
#endif
       ! ==================================================
       ! Shallow Water Test Cases:
       ! ==================================================
       if(Debug) print *,'homme: point #14'

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       if (test_case(1:5) == "swtc1") then

          ! ==================================================
          ! Shallow Water Test Case 1:  cosine bell
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call tc1_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          end if

       else if (test_case(1:5) == "swtc2") then
          if(Debug) print *,'homme: point #15'

          ! ==================================================
          ! Shallow Water Test Case 2:  cosine bell
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call tc2_errors(elem, 7, tl, pmean, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16'

       else if (test_case(1:5) == "swtc5") then

          ! ===============================================================
          ! Shallow Water Test Case 5: Rossby Haurwitz Waves
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ===============================================================

          if(Debug) print *,'homme: point #16.1'
          if (MODULO(tl%nstep,statefreq)==0) then
             call tc5_invariants(elem, 90, tl, pmean, edge2, deriv, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16.2'

          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call tc5_errors(elem, 7,tl, pmean, "ref_tc5_imp", simday, hybrid, nets, nete,par)
          end if
          if(Debug) print *,'homme: point #16.3'

       else if (test_case(1:5) == "swtc6") then

          ! ===============================================================
          ! Shallow Water Test Case 6: Rossby Haurwitz Waves
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ===============================================================

          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call tc6_errors(elem, 7,tl, pmean, "ref_tc6_imp", simday, hybrid, nets, nete,par)
          end if

       else if (test_case(1:5) == "swtc8") then
          ! Nothing !
       else if (test_case(1:6) == "vortex") then

          ! ==================================================
          ! Shallow Water Test Case 9:  vortex
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================

          if (MODULO(tl%nstep,statefreq)==0) then
             call vortex_errors(elem, 7, tl, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16'

       else if (test_case(1:5) == "swirl") then

          ! ==================================================
          ! Shallow Water Test Case:  swirl
          ! L1,L2,Linf error norms every statefreq timestep
          ! ==================================================
          if (MODULO(tl%nstep,statefreq)==0) then
             call swirl_errors(elem, 7, tl, hybrid, nets, nete)
          end if
          if(Debug) print *,'homme: point #16'

       else if (test_case(1:5) == "swsj1") then

          ! ==================================================
          ! Shallow Water Test Case:  swsj1
          ! L1,L2,Linf error norms every model day (compared to real soln)
          ! 
          ! detect day rollover
          ! ==================================================
          if (MODULO(Time_at(tl%nstep),secpday) <= 0.5*tstep) then
             simday=NINT(Time_at(tl%nstep)/secpday)
             call sj1_errors(elem,7,tl,pmean,"ref_sj1_imp",simday,hybrid,nets,nete,par)
          end if
          if(Debug) print *,'homme: point #16'

       end if

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       if(Debug) print *,'homme: point #17'
       if (MODULO(tl%nstep,statefreq)==0 ) then 
          if(hybrid%masterthread) then
             print *,tl%nstep,"time=",Time_at(tl%nstep)/secpday," days"
          !   print *,"Integrating at ",dt/dt_gv," times gravity wave restriction"
          end if

          call printstate(elem,pmean,g_sw_output,tl%n0,hybrid,nets,nete,kmass)

          call sweq_invariants(elem,190,tl,pmean,edge3,deriv,hybrid,nets,nete)
       end if
#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif

       ! ============================================================
       ! Write restart files if required
       ! ============================================================

       if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then
          call WriteRestart(elem, ithr,nets,nete,tl)
       endif

    end do

    ! ======================================================
    ! compute and report times...
    ! ======================================================


    ! ==========================
    ! end of the hybrid program
    ! ==========================
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
    call t_stopf('sweq')
  end subroutine sweq_rk

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! fvm driver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine Shal_Advec_Tracers_fvm(elem, fvm, deriv,hybrid,&
                                        dt,tl,nets,nete)
      use element_mod, only : element_t
      use fvm_control_volume_mod, only : fvm_struct
      use derivative_mod, only : derivative_t
      use kinds, only : real_kind
      use hybrid_mod, only : hybrid_t
      use time_mod, only : timelevel_t
      use perf_mod, only : t_startf, t_stopf, t_barrierf            ! _EXTERNAL
      use derivative_mod, only : divergence_sphere, ugradv_sphere
      use fvm_mod, only :  cslam_runairdensity, edgeveloc, fvm_mcgregor,fvm_mcgregordss
      use bndry_mod, only : bndry_exchangev
      use edge_mod, only  : edgevpack, edgevunpack
      use dimensions_mod, only : np, nlev
      
      implicit none
      type (element_t), intent(inout)               :: elem(:)
      type (fvm_struct), intent(inout)              :: fvm(:)
      type (derivative_t), intent(in)               :: deriv
      type (hybrid_t),     intent(in)               :: hybrid
      type (TimeLevel_t), intent(in)                :: tl

      real(kind=real_kind) , intent(in)             :: dt
      integer,intent(in)                            :: nets,nete

      integer :: ie,k

      real (kind=real_kind), dimension(np, np,2) :: vstar, vhat
      real (kind=real_kind), dimension(np, np) :: v1, v2


      call t_barrierf('sync_shal_advec_tracers_fvm', hybrid%par%comm)
      call t_startf('shal_advec_tracers_fvm')

      ! using McGregor AMS 1993 scheme: Economical Determination of Departure Points for
      ! Semi-Lagrangian Models 
!-BEGIN McGregor Without DSS
!       do ie=nets,nete
!         do k=1,nlev
!            ! Convert wind to lat-lon
!           v1     = (elem(ie)%state%v(:,:,1,k,tl%n0) + elem(ie)%state%v(:,:,1,k,tl%np1))/2.0D0  ! contra
!           v2     = (elem(ie)%state%v(:,:,2,k,tl%n0) + elem(ie)%state%v(:,:,2,k,tl%np1))/2.0D0   ! contra 
!           vhat(:,:,1)=elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2   ! contra->latlon
!           vhat(:,:,2)=elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2   ! contra->latlon
!           
!            ! Convert wind to lat-lon
!           v1     = elem(ie)%state%v(:,:,1,k,tl%np1)  ! contra
!           v2     = elem(ie)%state%v(:,:,2,k,tl%np1)   ! contra 
!           vstar(:,:,1)=elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2   ! contra->latlon
!           vstar(:,:,2)=elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2   ! contra->latlon
!           
!           ! calculate high order approximation
!           call fvm_mcgregor(elem(ie), deriv, dt, vhat,vstar, 1)
!           ! apply DSS to make vstar C0
!           elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%spheremp(:,:)*vstar(:,:,1) 
!           elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%spheremp(:,:)*vstar(:,:,2) 
!         enddo 
!         call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,ie)
!         call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,ie)
!       enddo 
!       call bndry_exchangeV(hybrid,edgeveloc)
!       do ie=nets,nete
!          call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,ie)
!          call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,ie)
!          do k=1, nlev  
!            elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
!            elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
!          end do
!       end do
!-END McGregor Without DSS

!-------BEGIN McGregor scheme
    do ie=nets,nete
      do k=1,nlev
         ! Convert wind to lat-lon
        v1     = elem(ie)%state%v(:,:,1,k,tl%n0)  ! contra
        v2     = elem(ie)%state%v(:,:,2,k,tl%n0)  ! contra 
        fvm(ie)%vn0(:,:,1,k)=elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2   ! contra->latlon
        fvm(ie)%vn0(:,:,2,k)=elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2   ! contra->latlon
        
         ! Convert wind to lat-lon
        v1     = elem(ie)%state%v(:,:,1,k,tl%np1)  ! contra
        v2     = elem(ie)%state%v(:,:,2,k,tl%np1)   ! contra 
        elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2   ! contra->latlon
        elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2   ! contra->latlon
        
      enddo  
    end do  
      call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
!-------END new McGregor scheme--------

      ! fvm departure calcluation should use vstar.
      ! from c(n0) compute c(np1): 
      ! call cslam_run(elem,fvm,hybrid,deriv,dt,tl,nets,nete)
      
      call cslam_runairdensity(elem,fvm,hybrid,deriv,dt,tl,nets,nete,0.0_real_kind)

      call t_stopf('shal_advec_tracers_fvm')
    end subroutine shal_advec_tracers_fvm  


end module sweq_mod



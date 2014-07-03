!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module diags_on_lat_aux_grid 

!BOP
! !MODULE: diags_on_lat_aux_grid

! !DESCRIPTION:
!  This module contains routines to compute some diagnostics
!  on an auxilary latitudinal grid. These diagnostics include
!  the meridional overturning circulation, northward T and S
!  transports, and zonal averages.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use io
   use io_tools
   use exit_mod
   use grid
   use global_reductions
   use gather_scatter
   use constants
   use registry
   use timers
   use shr_pio_mod, only : shr_pio_getioroot

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public ::   &
     init_lat_aux_grid,              &
     init_moc_ts_transport_arrays,   &
     compute_moc,                    &
     compute_tracer_transports

! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:),public,allocatable ::  &
      lat_aux_center,     &! cell center latitude values (degrees north) 
      lat_aux_edge         ! cell edge   latitude values (degrees north)

   integer (int_kind),public ::  &
      n_lat_aux_grid,     &! auxilary grid dimension  
      n_transport_reg,    &! number of regions for all transport diagnostics
      n_moc_comp,         &! number of MOC components 
      n_transport_comp,   &! number of T & S transport components
                           !  (see init_moc_ts_transport_arrays for region
                           !   and component details/definitions)
      nreg2_transport      ! number of basins in transports region 2

   integer(int_kind), private ::  &
      timer_moc, timer_tracer_transports

   integer (int_kind), dimension(:), allocatable ::  &
      lat_aux_region_start ! starting latitude indices for 
                           !  regions (not used for "global" region)

   real (r8), dimension(:,:), allocatable ::  &
      TLATD_G              ! latitude of T points in degrees (global array)

   real (r4), dimension(:,:,:,:), public, allocatable ::  &
      TAVG_MOC_G             ! meridional overturning circulation (ioroot only)

   real (r4), dimension(:,:,:), public, allocatable ::  &
      TR_TRANS_G,            &! tracer transports; used to compute both heat & salt
      TAVG_N_HEAT_TRANS_G,   &! northward heat transport (ioroot only)
      TAVG_N_SALT_TRANS_G     ! northward salt transport (ioroot only)

   real (r8),dimension(:,:), allocatable ::  &
      trans_s               ! southern boundary transports 

   integer (int_kind), dimension(:,:,:), allocatable ::  &
      REGION_MASK_LAT_AUX   ! latitude-longitude region mask 
                            !  for these diagnostics (ioroot only)

   logical (log_kind), dimension(:,:,:), allocatable ::  &
      MASK_LAT_DEPTH        ! latitude-depth mask for these diagnostics 
                            !  (ioroot only)

   type (regions), dimension(max_regions),public ::  &
      transport_region_info
 
 
   logical (log_kind), public ::  &
      moc_requested,              &! true if meridional overturning circulation 
                                   !  output is requested
      n_heat_trans_requested,     &! true if northward heat transport output is
                                   !  requested
      n_salt_trans_requested       ! true if northward salt transport output is
                                   !  requested
 
   character (char_len),dimension(max_regions),public ::  &
      transport_reg2_names
    
!EOP
!BOC


!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_lat_aux_grid
! !INTERFACE:

   subroutine init_lat_aux_grid
 

!
! !DESCRIPTION:
! This subroutine initialize the auxilary latitudinal grid. The grid choices are
! (i.e. lat_aux_grid_type =)
!
!   'southern' assumes that the model grid is regular lat-lon in the
!              Southern Hemisphere only and uses it identically in
!              the Southern Hemisphere. it is flipped across the
!              Equator and padded at northern high latitudes if 
!              necessary.
!   'full'     assumes that the entire model grid is regular
!              lat-lon. simply copies this grid into the axilary
!              grid arrays. 
!   'user'     allows the user to specify an equally-spaced grid,
!              starting at lat_aux_begin and ending at lat_aux_end.
!              the grid dimension/resolution is specified with
!              n_lat_aux_grid. lat_aux_begin and lat_aux_end specify
!              the egde coordinates in degrees north.
!
! the default namelist choice is 'southern'.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
 
!-----------------------------------------------------------------------
!
! local variables 
!
!-----------------------------------------------------------------------

   integer (int_kind) ::   &
      ioroot,              &! Task number of PIO root process
      nml_error,           &! namelist i/o error flag
      grid_error,          &! auxilary grid choice error
      j, jj, n,            &! loop indices
      lat_aux_grid,        &! index for the chosen grid type
      i_copy,              &! TLATD_G(i_copy,*) and ULATD_G(i_copy,*) 
                            !  (global arrays) are used to create 
                            !  the auxilary grid
      j_dim_sh              ! number of southern hemisphere TLATD_G(i_copy,*)
                            !  grid points

   real (r8) ::               &
      dlat,                   &! work variable for auxilary grid spacing (degrees)
      southern_edge,          &! latitude of the southern-most edge point 
      np_minus_northern_edge, &! latitudinal range for padding 
      eps_grid = 1.0e-7        ! epsilon difference allowed in regular grid

   integer (int_kind), parameter ::  &
      lat_aux_grid_sh   = 1,         &
      lat_aux_grid_full = 2,         &
      lat_aux_grid_user = 3

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK

   real (r8), dimension(:,:), allocatable ::  &
      ULATD_G             ! latitude of U points in degrees (global array) 

   character (char_len) :: string
!-----------------------------------------------------------------------
!
!  input namelist 
!
!  input values for lat_aux_begin, lat_aux_end, and n_lat_aux_grid
!  are used only when a user defined auxilary grid is requested. 
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      lat_aux_grid_type      ! type of the auxilary latitudinal grid,
                             ! i.e. how it is generated 
 
   real (r8) ::             &
      lat_aux_begin,        &! beginning latitude for the auxilary 
                             !    grid (degrees north)
      lat_aux_end            ! ending latitude for the auxilary
                             !    grid (degrees north)
 
   namelist /transports_nml/ lat_aux_grid_type, lat_aux_begin,  & 
                             lat_aux_end, n_lat_aux_grid,       &
                             moc_requested, n_heat_trans_requested, n_salt_trans_requested,   &
                             transport_reg2_names,              &
                             n_transport_reg

 
!-----------------------------------------------------------------------
!
!  set defaults and then read the namelist 
!
!-----------------------------------------------------------------------

   lat_aux_grid_type      =  'southern'   
   lat_aux_begin          = -90.0_r8
   lat_aux_end            =  90.0_r8
   n_lat_aux_grid         =  180
   moc_requested          = .false.
   n_heat_trans_requested = .false.
   n_salt_trans_requested = .false.
   transport_reg2_names   =  char_blank
   n_transport_reg        = 2 

!-----------------------------------------------------------------------
!
!  read options from namelist input file
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=transports_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading transports_nml')
   endif


   call broadcast_scalar (lat_aux_grid_type,      master_task)
   call broadcast_scalar (lat_aux_begin,          master_task)
   call broadcast_scalar (lat_aux_end,            master_task)
   call broadcast_scalar (n_lat_aux_grid,         master_task)
   call broadcast_scalar (moc_requested,          master_task)
   call broadcast_scalar (n_heat_trans_requested, master_task)
   call broadcast_scalar (n_salt_trans_requested, master_task)
   call broadcast_scalar (n_transport_reg,        master_task)
   do n=1,max_regions
     call broadcast_scalar(transport_reg2_names(n), master_task)
   end do

   if ( nml_error /= 0 ) then
     call exit_POP (SigAbort,'(init_lat_aux_grid) reading transports_nml')
   endif

!-----------------------------------------------------------------------
!
!  document namelist parameters
!
!-----------------------------------------------------------------------
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Transport Diagnostics:'
      write(stdout,blank_fmt)
      write(stdout,*) ' transports_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,transports_nml)
      write(stdout,blank_fmt)
   endif

 
!-----------------------------------------------------------------------
!
!  determine if transport diagnostics need to be computed
!
!-----------------------------------------------------------------------
 
   if (.not. (moc_requested .or. n_heat_trans_requested .or. n_salt_trans_requested) ) return
 
 
!-----------------------------------------------------------------------
!
!  document transport diagnostics selections
!
!-----------------------------------------------------------------------
 
   if ( my_task == master_task ) then
     write (stdout,*)
     write (stdout,'(a)  ') 'Transport Diagnostics Information'
     write (stdout,'(a,/)') '_________________________________'
     write (stdout,'(a)  ') 'Latitudinal Auxiliary Grid Choice is:'

     select case (lat_aux_grid_type(1:3))

     case ('sou')
       lat_aux_grid = lat_aux_grid_sh
       write (stdout,'(a)')  &
         ' ... based on flipped Southern Hemisphere latititude grid'
     case ('ful')
       lat_aux_grid = lat_aux_grid_full
       write (stdout,'(a)') ' ... full model latitudinal grid'
     case ('use')
       lat_aux_grid = lat_aux_grid_user
       write (stdout,'(a)') ' ... user-specified w/ equal spacing '
     case default
       lat_aux_grid = -1000
     end select
   endif
 
!-----------------------------------------------------------------------
!
!  test -- is lat_aux_grid defined properly?
!
!-----------------------------------------------------------------------
 
   call broadcast_scalar (lat_aux_grid, master_task)
 
   if ( lat_aux_grid == -1000 ) then
     call exit_POP (SigAbort, &
       '(init_lat_aux_grid) unknown auxilary latitudinal grid choice')
   endif


!-----------------------------------------------------------------------
!
!  continue documenting transport diagnostics
!
!-----------------------------------------------------------------------
 
   if ( my_task == master_task ) then
     write (stdout,*)
     write (stdout,'(a)  ') 'Transport diagnostics include:'
     if (moc_requested)          write (stdout,'(a)') 'MOC'
     if (n_heat_trans_requested) write (stdout,'(a)') 'N_HEAT'
     if (n_salt_trans_requested) write (stdout,'(a)') 'N_SALT'
     write (stdout,*)
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif
 

 
!-----------------------------------------------------------------------
!
!  more error checking
!
!-----------------------------------------------------------------------

   if (partial_bottom_cells) then
       string ='if (partial_bottom_cells), then 1st modify ' /&
                 &/' subroutines compute_moc and compute_tracer_transports'
       call exit_POP (SigAbort,'(init_lat_aux_grid): '// trim(string))
   endif
 
   ioroot = shr_pio_getioroot(inst_name)
 
!-----------------------------------------------------------------------
!
!  allocate arrays
!
!-----------------------------------------------------------------------
 
   allocate ( TLATD_G(nx_global,ny_global) )

   
   call gather_global (TLATD_G, TLATD, ioroot,distrb_clinic)

   if ( lat_aux_grid == lat_aux_grid_sh  .or.  &
        lat_aux_grid == lat_aux_grid_full ) then

     allocate ( ULATD_G(nx_global,ny_global) )

     WORK = ULAT * radian
     call gather_global (ULATD_G, WORK, ioroot,distrb_clinic)

     i_copy = 1

     if ( my_task == ioroot ) then
       dlat          = c2 * (ULATD_G(i_copy,1)-TLATD_G(i_copy,1)) 
       southern_edge = ULATD_G(i_copy,1) - dlat
     endif

     call broadcast_scalar (dlat,          ioroot)
     call broadcast_scalar (southern_edge, ioroot)

   endif

   if ( lat_aux_grid == lat_aux_grid_sh ) then

     if ( my_task == ioroot )  &
       j_dim_sh = count( TLATD_G(i_copy,:) < c0 )

     call broadcast_scalar (j_dim_sh, ioroot) 

     if ( j_dim_sh == 0 ) then
       call exit_POP (SigAbort,  &
         '(init_lat_aux_grid) there are no Southern Hemisphere grid points')
     endif

     grid_error = 0
     if ( my_task == ioroot ) then
       do j=1,j_dim_sh
         if(any(abs(TLATD_G(:,j)-TLATD_G(i_copy,j)) > eps_grid))then
           grid_error = -1000
         endif
       enddo
     endif

     call broadcast_scalar (grid_error, ioroot)
           
     if ( grid_error /= 0 ) then 
       string = 'SH is not a regular at-lon grid. '  /&
                &/'Use a different choice for lat_aux_grid_type.'
       call exit_POP (SigAbort,'(init_lat_aux_grid): '// trim(string))
     endif

     np_minus_northern_edge = 90.0_r8 + southern_edge        

     if ( np_minus_northern_edge >= p5*dlat ) then

       n_lat_aux_grid = nint( np_minus_northern_edge / dlat )

       if ( n_lat_aux_grid == 0 ) &
         call exit_POP (SigAbort,  &
            '(init_lat_aux_grid) n_lat_aux_grid is zero')
   
       dlat = np_minus_northern_edge / dble(n_lat_aux_grid)

       n_lat_aux_grid = 2 * j_dim_sh + n_lat_aux_grid 

     else

       n_lat_aux_grid = 2 * j_dim_sh

     endif

   endif

   if ( lat_aux_grid == lat_aux_grid_full ) then

     n_lat_aux_grid = ny_global

     grid_error = 0
     if ( my_task == ioroot ) then
       do j=1,n_lat_aux_grid
         if ( any(TLATD_G(:,j) /= TLATD_G(i_copy,j) ) ) grid_error = -1000
       enddo
     endif

     call broadcast_scalar (grid_error, ioroot)

     if ( grid_error /= 0 ) then 
       string = 'The model grid is not a regular lat-lon grid. ' /&
                &/'Use a different choice for lat_aux_grid_type.'
       call exit_POP (SigAbort,'(init_lat_aux_grid): '// trim(string))
     endif

   endif

   if ( lat_aux_grid == lat_aux_grid_user ) then

     if ( lat_aux_end <= lat_aux_begin ) then
       call exit_POP (SigAbort,  &
          '(init_lat_aux_grid) lat_aux_end should be > lat_aux_begin')
     endif

   endif

 
!-----------------------------------------------------------------------
!
!  allocate arrays
!
!-----------------------------------------------------------------------
 
   allocate ( lat_aux_edge  (n_lat_aux_grid+1),  &
              lat_aux_center(n_lat_aux_grid  ) )

   if ( my_task == ioroot ) then

     select case (lat_aux_grid)
 
     case (lat_aux_grid_sh)

       lat_aux_edge(1)            = southern_edge
       lat_aux_edge(2:j_dim_sh+1) = ULATD_G(i_copy,1:j_dim_sh)

       jj = j_dim_sh
       do j=j_dim_sh+2,2*j_dim_sh+1
         lat_aux_edge(j) = - lat_aux_edge(jj)
         jj = jj - 1 
       enddo 

       lat_aux_center(1:j_dim_sh) = TLATD_G(i_copy,1:j_dim_sh)

       jj = j_dim_sh
       do j=j_dim_sh+1,2*j_dim_sh
         lat_aux_center(j) = - lat_aux_center(jj)
         jj = jj - 1
       enddo

       if ( n_lat_aux_grid == 2 * j_dim_sh ) then
            lat_aux_edge(n_lat_aux_grid+1) = 90.0_r8
       else

         do j=2*j_dim_sh+2,n_lat_aux_grid+1
           lat_aux_edge(j) = lat_aux_edge(2*j_dim_sh+1)+(j-2*j_dim_sh-1)*dlat
         enddo

         do j=2*j_dim_sh+1,n_lat_aux_grid
           lat_aux_center(j) = lat_aux_edge(2*j_dim_sh+1)+(j-2*j_dim_sh-p5)*dlat
         enddo

       endif 

     case (lat_aux_grid_full)

       lat_aux_edge(1)               = southern_edge  
       lat_aux_edge(2:n_lat_aux_grid+1) = ULATD_G(i_copy,:)

       lat_aux_center = TLATD_G(i_copy,:)

     case (lat_aux_grid_user)

       dlat = (lat_aux_end - lat_aux_begin) / dble(n_lat_aux_grid)
    
       do j=1,n_lat_aux_grid+1 
         lat_aux_edge(j) = lat_aux_begin + dble(j-1)*dlat
       enddo

       do j=1,n_lat_aux_grid
         lat_aux_center(j) = lat_aux_begin + p5*dlat + dble(j-1)*dlat
       enddo
   
     end select

   endif

   call broadcast_array (lat_aux_edge,   ioroot)
   call broadcast_array (lat_aux_center, ioroot)
          
   if ( lat_aux_grid == lat_aux_grid_sh  .or.   &
        lat_aux_grid == lat_aux_grid_full ) deallocate ( ULATD_G )

!-----------------------------------------------------------------------
!  initialize timers
!-----------------------------------------------------------------------
 
   if (moc_requested)  &
   call get_timer(timer_moc,'MOC',  &
                  nblocks_clinic, distrb_clinic%nprocs)
   if (n_heat_trans_requested .or. n_salt_trans_requested)  &
   call get_timer(timer_tracer_transports,'TRACER_TRANSPORTS',  &
                  nblocks_clinic, distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

 end subroutine init_lat_aux_grid

!***********************************************************************
!BOP
! !IROUTINE: init_moc_ts_transport_arrays
! !INTERFACE:
 subroutine init_moc_ts_transport_arrays 

! !DESCRIPTION:
!  This routine allocates necessary arrays for the meridional overturning
!  circulation and northward T & S transports, define the
!  associated region masks, and find the "Atlantic" region
!  starting latitude.  
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8)    ::  &
     eps_grid = 1.0e-7   ! epsilon difference allowed in regular grid
 
   integer (int_kind) ::   &
      ioroot,            &! Task number of PIO root process
      nml_error,         &! namelist i/o error flag
      i, j, k, n, m,     &! loop indices
      j_southern,        &! loop indices
      nrtr,              &! loop index
      grid_error          ! grid error

   integer (int_kind), dimension(:,:), allocatable ::  &
      WORK_G              ! global work array

   if (.not. (moc_requested .or. n_heat_trans_requested .or. n_salt_trans_requested) ) return
 

!-----------------------------------------------------------------------
!
!  for now we consider only 2 lat-lon regions for these diagnostics:
!
!    n_transport_reg = 1 ----> global minus marginal seas
!    n_transport_reg = 2 ----> REGION_MASK = Atlantic + Mediterranean
!    (optional)               (if not a marginal sea) + Labrador +
!                            + GIN + Arctic + Hudson 
!                             (if not a marginal sea)
!
!-----------------------------------------------------------------------
 
   if (n_transport_reg < 1 ) then
      call exit_POP (SigAbort,'(init_moc_ts_transport_arrays) ' /&
                 &/ ' n_transport_reg must be > 0  -- '         /&
                 &/ ' check namelist transports_nml')
   elseif (n_transport_reg > 2) then
      call exit_POP (SigAbort,'(init_moc_ts_transport_arrays) ' /&
                 &/ ' n_transport_reg must be < 3  -- '            /&
                 &/ ' check namelist transports_nml')
   endif

 
      
!-----------------------------------------------------------------------
!
!  determine the region numbers associated with the selected input
!   region names, transport_reg2_names, which are defined when
!   n_transport_reg = 2
!
!-----------------------------------------------------------------------
 
   ioroot = shr_pio_getioroot(inst_name)
   nreg2_transport = 0
   transport_region_info(:)%name   = char_blank
   transport_region_info(:)%number = 9999

   if (n_transport_reg > 1) then
 
   do n=1,max_regions
     do nrtr = 1, max_regions         
      if (len_trim(transport_reg2_names(nrtr)) > 0 ) then
       if ( trim(transport_reg2_names(nrtr)) == &
            trim(region_info(n)%name)  .and.    &
            .not.region_info(n)%marginal_sea) then
         nreg2_transport = nreg2_transport + 1
         transport_region_info(nreg2_transport)%name    = region_info(n)%name   
         transport_region_info(nreg2_transport)%number  = region_info(n)%number
         transport_region_info(nreg2_transport)%marginal_sea = region_info(n)%  &
                                                     marginal_sea
       endif
      endif
     enddo
   enddo
 
   if (nreg2_transport <= 0) then
      call exit_POP (SigAbort,'(init_moc_ts_transport_arrays) '    /&
                 &/ ' no transport regions have been detected -- ' /&
                 &/ ' check namelist transports_nml')
   endif

   if ( my_task == master_task ) then
     write(stdout,*) 'The following ',nreg2_transport,  &
          ' regions will be included in the n_transport_reg = 2 transports:'
     do nrtr = 1, nreg2_transport
       write(stdout,1000) transport_region_info(nrtr)%name,  &
                          transport_region_info(nrtr)%number 
     enddo
1000    format (2x, a35, '(',i2,')')
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

   endif  ! n_transport_reg > 1
 
!-----------------------------------------------------------------------
!
!  allocate lat_aux_region_start, REGION_MASK_LAT_AUX
!
!-----------------------------------------------------------------------
 
   allocate (lat_aux_region_start(n_transport_reg) )

   allocate (REGION_MASK_LAT_AUX(nx_global,ny_global,n_transport_reg) )

 
   lat_aux_region_start = 0
   call gather_global (REGION_MASK_LAT_AUX(:,:,1),REGION_MASK,  &
                       ioroot,distrb_clinic)


   if ( my_task == ioroot ) then


     if (n_transport_reg > 1) then
        REGION_MASK_LAT_AUX(:,:,2) = REGION_MASK_LAT_AUX(:,:,1)
     endif

     REGION_MASK_LAT_AUX(:,:,1) =   &
                    merge( 1, 0, REGION_MASK_LAT_AUX(:,:,1) > 0 )

    if (n_transport_reg > 1) then
 
     do j=1,ny_global
       do i=1,nx_global
        if (any(abs(REGION_MASK_LAT_AUX(i,j,2)) ==   &
            transport_region_info(:)%number)) then
            REGION_MASK_LAT_AUX(i,j,2) = 1
         else
            REGION_MASK_LAT_AUX(i,j,2) = 0
         endif
       enddo
     enddo
 
     j_southern = ny_global
     do j=1,ny_global
       do i=1,nx_global
         if ( REGION_MASK_LAT_AUX(i,j,2) == 1 ) then
           j_southern = j
           goto 100
         endif
       enddo
     enddo
100  continue
 
     grid_error = 0
     do i=1,nx_global
       if ( any(abs(TLATD_G(:,j_southern)-TLATD_G(i,j_southern)) > eps_grid))then
         grid_error = -1000
       endif
     enddo

     lat_aux_region_start(2) = j_southern - 1 
            
   endif  ! n_transport_reg
   endif  ! ioroot

   if (n_transport_reg > 1) then
     call broadcast_scalar (grid_error,           ioroot)
     call broadcast_array  (lat_aux_region_start, ioroot)

     if ( grid_error /= 0 ) then
       call exit_POP (SigAbort,'(init_moc_ts_transport_arrays) SH is' /&
                   &/ ' not a regular lat-lon grid. The'           /&
                   &/ ' southern boundary for region 2'           /&
                   &/ ' ("Atlantic") cannot be specified.')
     endif
   endif  ! n_transport_reg
 
!-----------------------------------------------------------------------
!
!  determine the latitude-depth mask
!
!-----------------------------------------------------------------------

   allocate ( WORK_G(nx_global,ny_global) ) 

   call gather_global (WORK_G, KMT, ioroot,distrb_clinic)

   if ( my_task == ioroot ) then

     allocate ( MASK_LAT_DEPTH(n_lat_aux_grid,km,n_transport_reg) )

     MASK_LAT_DEPTH = .true.
 
     do k=1,km
       do n=2,n_lat_aux_grid+1
         do j=1,ny_global
           do i=1,nx_global

             if ( TLATD_G(i,j) >= lat_aux_edge(n-1) .and.  &
                  TLATD_G(i,j) <  lat_aux_edge(n  ) .and.  &
                  k <= WORK_G(i,j)                  ) then 
               do m=1,n_transport_reg
                 if ( REGION_MASK_LAT_AUX(i,j,m) == 1 )  &
                   MASK_LAT_DEPTH(n-1,k,m) = .false.
               enddo
             endif

           enddo
         enddo
       enddo
     enddo

   endif

   deallocate ( WORK_G )

   if ( moc_requested ) then
!-----------------------------------------------------------------------
!
!  MOC may have 3 components:
!
!    n_moc_comp = 1 ----> Eulerian-mean
!    n_moc_comp = 2 ----> Eddy-induced (bolus) if diag_gm_bolus is true 
!                          and GM is on
!    n_moc_comp = 3 ----> Submesoscale contribution if GM is on and
!                          submesoscale_mixing is true
!
!-----------------------------------------------------------------------

     n_moc_comp = 1
     if ( registry_match('diag_gm_bolus') )  n_moc_comp = 2
     if ( registry_match('init_submeso') )   n_moc_comp = 3

!-----------------------------------------------------------------------
!
!  allocate TAVG_MOC_G on ioroot only
!
!-----------------------------------------------------------------------
 
     if ( my_task == ioroot ) then
       allocate (TAVG_MOC_G(n_lat_aux_grid+1,km+1,n_moc_comp,n_transport_reg))
     endif

   endif

!-----------------------------------------------------------------------  
!
!  T and S transports may have 5 components: 
!
!(1) n_transport_comp = 1 ----> total [i.e. (2) + (3)]
!
!(2) n_transport_comp = 2 ----> Eulerian-mean advection
!
!(3) n_transport_comp = 3 ----> Eddy-induced advection plus diffusion
!                               if GM is on
!                                              OR
!                               Eddy-induced advection plus submesoscale advection
!                               plus diffusion if GM is on and 
!                               submesoscale_mixing is true
!                                              OR
!                               diffusion if hmix_tracer_choice is not GM
!
!(4) n_transport_comp = 4 ----> Eddy-induced advection if diag_gm_bolus
!                               is true and GM is on (diagnostic computation)
!
!(5) n_transport_comp = 5 ----> Submesoscale advection if GM is on and
!                               submesoscale_mixing is true (diagnostic
!                               computation)
!
!-----------------------------------------------------------------------

   if ( n_heat_trans_requested .or. n_salt_trans_requested ) then
     n_transport_comp = 3
     if ( registry_match('diag_gm_bolus') )  n_transport_comp = 4
     if ( registry_match('init_submeso') )   n_transport_comp = 5
   endif

!-----------------------------------------------------------------------
!
!  allocate TAVG_N_HEAT_TRANS_G, TAVG_N_SALT_TRANS_G, and TR_TRANS_G, on
!  ioroot only
!
!-----------------------------------------------------------------------
 
   if ( n_heat_trans_requested .and. my_task == ioroot ) then
    allocate (  &
     TAVG_N_HEAT_TRANS_G(n_lat_aux_grid+1, n_transport_comp,n_transport_reg))
   endif

   if ( n_salt_trans_requested .and. my_task == ioroot ) then
    allocate (  &
     TAVG_N_SALT_TRANS_G(n_lat_aux_grid+1,n_transport_comp,n_transport_reg))
   endif

   if ((n_heat_trans_requested .or. n_salt_trans_requested) ) then
     allocate (trans_s (n_transport_comp,n_transport_reg) )
     if (my_task == ioroot ) then
       allocate (  &
         TR_TRANS_G (n_lat_aux_grid+1,n_transport_comp,n_transport_reg))
     endif
     call document ('init_moc_ts_transport_arrays','allocate TR_TRANS_G')
   endif

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) 'End of transport regions initialization'
     write(stdout,blank_fmt)
     write(stdout,ndelim_fmt)
     write(stdout,blank_fmt)
     call POP_IOUnitsFlush(POP_stdout);  call POP_IOUnitsFlush(stdout)
   endif
 
!-----------------------------------------------------------------------
!EOC

 end subroutine init_moc_ts_transport_arrays 

!***********************************************************************
!BOP
! !IROUTINE: compute_moc
! !INTERFACE:
 subroutine compute_moc ( W_E, V_E, W_I, V_I, W_SM, V_SM )

! !DESCRIPTION:
! This subroutine computes meridional overturning circulation
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
!  input variables. the present design assumes that these are 
!  time-averaged inputs from tavg.F.
!
   real (rtavg), dimension(:,:,:,:), intent(in) ::  &
      W_E,    &! Eulerian-mean vertical velocity component 
      V_E      ! Eulerian-mean velocity component in the grid-y direction

   real (rtavg), dimension(:,:,:,:), optional, intent(in) ::  &
      W_I,    &! Eddy-induced (bolus) vertical velocity component
      V_I,    &! Eddy-induced (bolus) velocity component in the
               !  grid-y direction
      W_SM,   &! Submeso vertical velocity component
      V_SM     ! Submeso velocity component in the grid-y direction

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      ioroot,             &     ! Task number of PIO root process
      i, j, k, m, n, iblock     ! loop indices

   real (r8), dimension(km,n_moc_comp,n_transport_reg) ::  &
      moc_s                     ! southern boundary values of moc 

   real (r8), dimension(:,:,:), allocatable ::  &
      WORK1, WORK2              ! work arrays

   real (r8), dimension(:,:,:), allocatable ::  &
      WORK1_G, WORK2_G, WORK3_G ! global work arrays

   logical (log_kind) ::  &
      ldiag_gm_bolus,     &     ! local logical for diag_gm_bolus
      lsubmeso                  ! local logical for submesoscale_mixing

   if (.not. moc_requested) return

   if ( (      present(W_I) .and. .not.present(V_I))  .or.  &
        ( .not.present(W_I) .and.      present(V_I)) ) then
     call exit_POP (SigAbort,'(compute_moc) both W_I and V_I'   &
                // ' are necessary fields for the eddy-induced' &
                // ' transport computations, but one is missing.')
   endif

   if ( (      present(W_SM) .and. .not.present(V_SM))  .or.  &
        ( .not.present(W_SM) .and.      present(V_SM)) ) then
     call exit_POP (SigAbort,'(compute_moc) both W_SM and V_SM'   &
                // ' are necessary fields for the submeso' &
                // ' transport computations, but one is missing.')
   endif

   ioroot = shr_pio_getioroot(inst_name)

   ldiag_gm_bolus = .false.
   if ( present(W_I) )  ldiag_gm_bolus = .true.
 
   lsubmeso = .false.
   if ( present(W_SM) )  lsubmeso = .true.
 
   call timer_start (timer_moc)

   allocate ( WORK1(nx_block,ny_block,nblocks_clinic), &
              WORK2(nx_block,ny_block,nblocks_clinic) )

   allocate ( WORK1_G(nx_global,ny_global,km) )

   if ( ldiag_gm_bolus ) allocate ( WORK2_G(nx_global,ny_global,km) )
   if ( lsubmeso )       allocate ( WORK3_G(nx_global,ny_global,km) )

   do k=1,km
    !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
        WORK1(:,:,iblock) = merge(W_E(:,:,k,iblock)*TAREA(:,:,iblock), c0, k <= KMT(:,:,iblock))
     enddo
    !$OMP END PARALLEL DO

     call gather_global (WORK1_G(:,:,k), WORK1, ioroot,distrb_clinic)

     if ( ldiag_gm_bolus ) then
       !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock = 1,nblocks_clinic
          WORK1(:,:,iblock) = merge(W_I(:,:,k,iblock)*TAREA(:,:,iblock), c0, k <= KMT(:,:,iblock))
        enddo
       !$OMP END PARALLEL DO
       call gather_global (WORK2_G(:,:,k), WORK1, ioroot,distrb_clinic)
     endif

     if ( lsubmeso ) then
       !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock = 1,nblocks_clinic
          WORK1(:,:,iblock) = merge(W_SM(:,:,k,iblock)*TAREA(:,:,iblock), c0, k <= KMT(:,:,iblock))
        enddo
       !$OMP END PARALLEL DO
       call gather_global (WORK3_G(:,:,k), WORK1, ioroot,distrb_clinic)
     endif

   enddo

   if ( my_task == ioroot )  then
    TAVG_MOC_G = c0


     do n=2,n_lat_aux_grid+1
       TAVG_MOC_G(n,:,:,:) = TAVG_MOC_G(n-1,:,:,:)
       do j=1,ny_global
        do i=1,nx_global
          if ( TLATD_G(i,j) >= lat_aux_edge(n-1) .and.  &
               TLATD_G(i,j) <  lat_aux_edge(n  ) ) then 
           do m=1,n_transport_reg
             if ( REGION_MASK_LAT_AUX(i,j,m) == 1 ) then 
                do k=1,km
                  TAVG_MOC_G(n,k,1,m)=TAVG_MOC_G(n,k,1,m)+WORK1_G(i,j,k)
                enddo
                if ( ldiag_gm_bolus ) then
                 do k=1,km
                  TAVG_MOC_G(n,k,2,m)=TAVG_MOC_G(n,k,2,m)+WORK2_G(i,j,k)
                 enddo
                endif 
                if ( lsubmeso ) then
                 do k=1,km
                  TAVG_MOC_G(n,k,3,m)=TAVG_MOC_G(n,k,3,m)+WORK3_G(i,j,k)
                 enddo
                endif
             endif
           enddo ! m
          endif ! n
        enddo ! i
       enddo ! j
     enddo ! n
   endif ! ioroot

!-----------------------------------------------------------------------
!
!  determine the southern boundary transports for all regional mocs 
!
!-----------------------------------------------------------------------

   moc_s = c0

   do k=1,km
    !$OMP PARALLEL DO PRIVATE(iblock,i,j)
     do iblock = 1,nblocks_clinic
       WORK1(:,:,iblock) = p5 * V_E(:,:,k,iblock) * DXU(:,:,iblock)
       do j=1,ny_block
       do i=2,nx_block
         WORK2(i,j,iblock)=WORK1(i-1,j,iblock)
       enddo  
       enddo
       WORK2(1,:,iblock)=c0
       WORK1(:,:,iblock) = WORK1(:,:,iblock) + WORK2(:,:,iblock) 
     enddo ! iblock
    !$OMP END PARALLEL DO
     call gather_global (WORK1_G(:,:,k), WORK1, ioroot,distrb_clinic)
   enddo

   if ( ldiag_gm_bolus ) then
     do k=1,km
       !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock = 1,nblocks_clinic
          WORK1(:,:,iblock) = V_I(:,:,k,iblock) * HTN(:,:,iblock) 
        enddo ! iblock
       !$OMP END PARALLEL DO
       call gather_global (WORK2_G(:,:,k), WORK1, ioroot,distrb_clinic)
     enddo
   endif

   if ( lsubmeso ) then
     do k=1,km
       !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock = 1,nblocks_clinic
          WORK1(:,:,iblock) = V_SM(:,:,k,iblock) * HTN(:,:,iblock)
        enddo ! iblock
       !$OMP END PARALLEL DO
       call gather_global (WORK3_G(:,:,k), WORK1, ioroot,distrb_clinic)
     enddo
   endif
       
   if ( my_task == ioroot ) then
     do m=2,n_transport_reg
       j = lat_aux_region_start(m)
       do i=1,nx_global
         if ( REGION_MASK_LAT_AUX(i,j+1,m) == 1 ) then
           do k=1,km
             moc_s(k,1,m) = moc_s(k,1,m) + WORK1_G(i,j,k)
           enddo ! k
           if ( ldiag_gm_bolus ) then
             do k=1,km
               moc_s(k,2,m) = moc_s(k,2,m) + WORK2_G(i,j,k)
             enddo ! k
           endif 
           if ( lsubmeso ) then
             do k=1,km
               moc_s(k,3,m) = moc_s(k,3,m) + WORK3_G(i,j,k)
             enddo ! k
           endif
         endif ! REGION_MASK_LAT_AUX
       enddo ! i
     enddo ! m

  
     moc_s(km,:,:) = - dz(km) * moc_s(km,:,:) 
     do k=km-1,1,-1
       moc_s(k,:,:) = moc_s(k+1,:,:) - dz(k) * moc_s(k,:,:)
     enddo

!-----------------------------------------------------------------------
!
!  add the southern boundary transports for all regional mocs 
!
!-----------------------------------------------------------------------

     do m=2,n_transport_reg
       do k=1,km
         do n=1,n_lat_aux_grid+1
           TAVG_MOC_G(n,k,1,m) = TAVG_MOC_G(n,k,1,m) + moc_s(k,1,m)
         enddo
         if ( ldiag_gm_bolus ) then
           do n=1,n_lat_aux_grid+1
             TAVG_MOC_G(n,k,2,m) = TAVG_MOC_G(n,k,2,m) + moc_s(k,2,m) 
           enddo
         endif
         if ( lsubmeso ) then
           do n=1,n_lat_aux_grid+1
             TAVG_MOC_G(n,k,3,m) = TAVG_MOC_G(n,k,3,m) + moc_s(k,3,m)
           enddo
         endif
       enddo
     enddo

!-----------------------------------------------------------------------
!
!  convert MOC to Sverdrups, prior to masking
!
!-----------------------------------------------------------------------
     TAVG_MOC_G = TAVG_MOC_G*1.0e-12_r8
 
!-----------------------------------------------------------------------
!
! use masks to determine 0 and missing transport values 
! not used in pop version; not translated to pop2
!
!-----------------------------------------------------------------------

!     interior
!
!      do m=1,n_transport_reg
!        do k=1,km-1
!          do n=1,n_lat_aux_grid-1 
!            if ( MASK_LAT_DEPTH(n  ,k  ,m) .or.
!     &           MASK_LAT_DEPTH(n+1,k  ,m) .or.
!     &           MASK_LAT_DEPTH(n  ,k+1,m) .or.
!     &           MASK_LAT_DEPTH(n+1,k+1,m) ) 
!     &        TAVG_MOC_G(n+1,k+1,:,m) = c0 
!            if ( MASK_LAT_DEPTH(n  ,k  ,m) .and.
!     &           MASK_LAT_DEPTH(n+1,k  ,m) .and.
!     &           MASK_LAT_DEPTH(n  ,k+1,m) .and.
!     &           MASK_LAT_DEPTH(n+1,k+1,m) ) 
!     &          TAVG_MOC_G(n+1,k+1,:,m) = undefined_nf 
!            enddo
!          enddo
!        enddo
!
!     top and bottom boundaries
!
!        do m=1,n_transport_reg
!          do n=1,n_lat_aux_grid-1
!            if ( MASK_LAT_DEPTH(n  ,1,m) .or.
!     &           MASK_LAT_DEPTH(n+1,1,m) )
!     &        TAVG_MOC_G(n+1,1,:,m) = c0
!            if ( MASK_LAT_DEPTH(n  ,1,m) .and.
!     &           MASK_LAT_DEPTH(n+1,1,m) )
!     &        TAVG_MOC_G(n+1,1,:,m) = undefined_nf 
!            if ( MASK_LAT_DEPTH(n  ,km,m) .or.
!     &           MASK_LAT_DEPTH(n+1,km,m) )
!     &        TAVG_MOC_G(n+1,km+1,:,m) = c0
!            if ( MASK_LAT_DEPTH(n  ,km,m) .and.
!     &           MASK_LAT_DEPTH(n+1,km,m) )
!     &        TAVG_MOC_G(n+1,km+1,:,m) = undefined_nf 
!          enddo
!        enddo
!
!     southern boundary
!
!        do m=1,n_transport_reg
!          do k=1,km-1
!            if ( MASK_LAT_DEPTH(1,k  ,m) .or.
!     &           MASK_LAT_DEPTH(1,k+1,m) )
!     &        TAVG_MOC_G(1,k+1,:,m) = c0
!            if ( MASK_LAT_DEPTH(1,k  ,m) .and.
!     &           MASK_LAT_DEPTH(1,k+1,m) )
!     &        TAVG_MOC_G(1,k+1,:,m) = undefined_nf
!          enddo
!        enddo
!
!        do m=1,n_transport_reg
!          if ( MASK_LAT_DEPTH(1,1 ,m) )
!     &      TAVG_MOC_G(1,1   ,:,m) = undefined_nf
!          if ( MASK_LAT_DEPTH(1,km,m) )
!     &      TAVG_MOC_G(1,km+1,:,m) = undefined_nf
!        enddo
!
!     northern boundary
!
!        do m=1,n_transport_reg
!          do k=1,km-1
!            if ( MASK_LAT_DEPTH(n_lat_aux_grid,k  ,m) .or.
!     &           MASK_LAT_DEPTH(n_lat_aux_grid,k+1,m) )
!     &        TAVG_MOC_G(n_lat_aux_grid+1,k+1,:,m) = c0
!            if ( MASK_LAT_DEPTH(n_lat_aux_grid,k  ,m) .and.
!     &           MASK_LAT_DEPTH(n_lat_aux_grid,k+1,m) )
!     &        TAVG_MOC_G(n_lat_aux_grid+1,k+1,:,m) = undefined_nf 
!          enddo
!        enddo
!
!        do m=1,n_transport_reg
!          if ( MASK_LAT_DEPTH(n_lat_aux_grid,1 ,m) )
!     &      TAVG_MOC_G(n_lat_aux_grid+1,1   ,:,m) = undefined_nf 
!          if ( MASK_LAT_DEPTH(n_lat_aux_grid,km,m) )
!     &      TAVG_MOC_G(n_lat_aux_grid+1,km+1,:,m) = undefined_nf 
!        enddo
!
   endif ! ioroot

 
   deallocate ( WORK1, WORK2, WORK1_G )

   if ( ldiag_gm_bolus )  deallocate ( WORK2_G )
   if ( lsubmeso )        deallocate ( WORK3_G )

   call timer_stop  (timer_moc)

!-----------------------------------------------------------------------
!EOC
 end subroutine compute_moc

!***********************************************************************
!BOP
! !IROUTINE: compute_tracer_transports
! !INTERFACE:
 subroutine compute_tracer_transports (tracer_index, ADV, HDIF, FN, &
                                       ADV_I, FN_I, ADV_SM, FN_SM ) 
! !DESCRIPTION
!  This subroutine computes northward tracer (T and S) transports 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
!  input variables. the present design assumes that these are 
!  time-averaged inputs from tavg.F.

   real (rtavg), dimension(:,:,:), intent(in) ::  &
      ADV,    &! vertically-integrated tracer Eulerian-mean advection tendency
      HDIF     ! vertically-integrated horz diff tracer tendency (when GM
               !  and submesoscale mixing are on, it includes eddy-induced and
               !  submeso velocity contributions, respectively)

   real (rtavg), dimension(:,:,:,:), intent(in) ::  &
      FN       ! flux of tracer in grid-y direction due to the Eulerian-mean 
               !  transport velocity

   integer (int_kind), optional, intent(in) ::  &
      tracer_index    
 
   real (rtavg), dimension(:,:,:), optional, intent(in) :: &
      ADV_I,  &! vertically-integrated tracer eddy-induced advection
               !  tendency (diagnostic)
      ADV_SM   ! vertically-integrated tracer submeso advection tendency
               !  (diagnostic)

   real (rtavg), dimension(:,:,:,:), optional, intent(in) ::  &
      FN_I,   &! flux of tracer in grid-y direction due to the
               !  eddy-induced velocity (diagnostic)
      FN_SM    ! flux of tracer in grid-y direction due to the 
               !  submeso velocity (diagnostic)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!   local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::      &
      ioroot,                 &! Task number of PIO root process
      i, j, k, m, n, iblock    ! loop indices

   real (r8) :: conversion
 
   real (r8), dimension(:,:,:), allocatable ::  &
      WORK1                    ! work arrays
   real (r8), dimension(:,:), allocatable ::  & 
      WORK1_G, WORK2_G,       &! global work arrays
      WORK3_G, WORK4_G

   logical (log_kind) ::  &
      ldiag_gm_bolus,     &    ! local logical for diag_gm_bolus
      lsubmeso                 ! local logical for submesoscale_mixing

!-----------------------------------------------------------------------
!
!  determine if this subroutine needs to be executed
!
!-----------------------------------------------------------------------
   if (.not. (n_heat_trans_requested .or. n_salt_trans_requested) ) return
 
!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------
   if      (tracer_index == 1  .and. n_heat_trans_requested ) then    
!           ok, continue
   else if (tracer_index == 2  .and. n_salt_trans_requested ) then    
!           ok, continue
   else
     call document('compute_tracer_transports','return upon entry ')
     return
   endif
 
   if ( (     present(ADV_I) .and. .not.present(FN_I))  .or.  &
        (.not.present(ADV_I) .and.      present(FN_I)) ) then
        call exit_POP (SigAbort,'(compute_tracer_transports)'   &
         // ' both ADV_I and FN_I are necessary fields for'     &
         // ' the related transport computations, but one is missing.')
   endif  

   if ( (     present(ADV_SM) .and. .not.present(FN_SM))  .or.  &
        (.not.present(ADV_SM) .and.      present(FN_SM)) ) then
        call exit_POP (SigAbort,'(compute_tracer_transports)'   &
         // ' both ADV_SM and FN_SM are necessary fields for'     &
         // ' the related transport computations, but one is missing.')
   endif

   ioroot = shr_pio_getioroot(inst_name)

   ldiag_gm_bolus = .false.
   if ( present(ADV_I) )  ldiag_gm_bolus = .true.
 
   lsubmeso = .false.
   if ( present(ADV_SM) ) lsubmeso = .true.
 
   call timer_start (timer_tracer_transports)

 
!-----------------------------------------------------------------------
!
!  allocate WORK arrays
!
!-----------------------------------------------------------------------
 
   allocate ( WORK1(nx_block,ny_block,nblocks_clinic) )

   allocate ( WORK1_G(nx_global,ny_global), WORK2_G(nx_global,ny_global) )

   do iblock = 1,nblocks_clinic
     WORK1(:,:,iblock) = - ADV(:,:,iblock) * TAREA(:,:,iblock)
   enddo
   call gather_global (WORK1_G, WORK1, ioroot,distrb_clinic)

   do iblock = 1,nblocks_clinic
     WORK1(:,:,iblock) = - HDIF(:,:,iblock) * TAREA(:,:,iblock)
   enddo
   call gather_global (WORK2_G, WORK1, ioroot,distrb_clinic)

   if ( ldiag_gm_bolus ) then
     allocate ( WORK3_G(nx_global,ny_global) )

     do iblock = 1,nblocks_clinic
       WORK1(:,:,iblock) = - ADV_I(:,:,iblock) * TAREA(:,:,iblock)
     enddo
     call gather_global (WORK3_G, WORK1, ioroot,distrb_clinic)
   endif

   if ( lsubmeso ) then
     allocate ( WORK4_G(nx_global,ny_global) )

     do iblock = 1,nblocks_clinic
       WORK1(:,:,iblock) = - ADV_SM(:,:,iblock) * TAREA(:,:,iblock)
     enddo
     call gather_global (WORK4_G, WORK1, ioroot,distrb_clinic)
   endif
     
   if ( my_task == ioroot ) then

     TR_TRANS_G = c0

     do n=2,n_lat_aux_grid+1

       TR_TRANS_G(n,:,:) = TR_TRANS_G(n-1,:,:)

       do m=1,n_transport_reg
         do j=1,ny_global
           do i=1,nx_global
             if ( TLATD_G(i,j) >= lat_aux_edge(n-1) .and.  &
                  TLATD_G(i,j) <  lat_aux_edge(n  ) ) then 
               if ( REGION_MASK_LAT_AUX(i,j,m) == 1 ) then
                 TR_TRANS_G(n,1,m) = TR_TRANS_G(n,1,m) + WORK1_G(i,j) + WORK2_G(i,j)
                 TR_TRANS_G(n,2,m) = TR_TRANS_G(n,2,m) + WORK1_G(i,j)
                 TR_TRANS_G(n,3,m) = TR_TRANS_G(n,3,m) + WORK2_G(i,j)
                 if ( ldiag_gm_bolus ) then
                   TR_TRANS_G(n,4,m) = TR_TRANS_G(n,4,m) + WORK3_G(i,j)
                 endif
                 if ( lsubmeso ) then
                   TR_TRANS_G(n,5,m) = TR_TRANS_G(n,5,m) + WORK4_G(i,j)
                 endif
               endif
             endif
           enddo
         enddo
       enddo

     enddo
   endif

!-----------------------------------------------------------------------
!
!  determine the southern boundary transports for all regional
!  transports 
!
!-----------------------------------------------------------------------

   trans_s = c0

   do k=1,km

    !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
        WORK1(:,:,iblock) = FN(:,:,k,iblock) * TAREA(:,:,iblock) * dz(k) 
     enddo
    !$OMP END PARALLEL DO
     call gather_global (WORK1_G, WORK1, ioroot,distrb_clinic)

     if ( ldiag_gm_bolus ) then
    !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
       WORK1(:,:,iblock) = FN_I(:,:,k,iblock) * TAREA(:,:,iblock) * dz(k)
     enddo
    !$OMP END PARALLEL DO
       call gather_global (WORK3_G, WORK1, ioroot,distrb_clinic)
     endif

     if ( lsubmeso ) then
    !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
       WORK1(:,:,iblock) = FN_SM(:,:,k,iblock) * TAREA(:,:,iblock) * dz(k)
     enddo
    !$OMP END PARALLEL DO
       call gather_global (WORK4_G, WORK1, ioroot,distrb_clinic)
     endif

     if ( my_task == ioroot ) then

       do m=2,n_transport_reg
         j = lat_aux_region_start(m)
         do i=1,nx_global
           if ( REGION_MASK_LAT_AUX(i,j+1,m) == 1 ) then
             trans_s(2,m) = trans_s(2,m) + WORK1_G(i,j)
             if ( ldiag_gm_bolus ) trans_s(4,m) = trans_s(4,m)  &
                                                 + WORK3_G(i,j)
             if ( lsubmeso )       trans_s(5,m) = trans_s(5,m)  &
                                                 + WORK4_G(i,j)
           endif
         enddo
       enddo

     endif

   enddo

   if ( my_task == ioroot ) then
 
!-----------------------------------------------------------------------
!
!  add the southern boundary transports for all regional transports
!
!-----------------------------------------------------------------------

     do m=2,n_transport_reg
       do n=1,n_lat_aux_grid+1
         TR_TRANS_G(n,2,m) = TR_TRANS_G(n,2,m) + trans_s(2,m) 
         if ( ldiag_gm_bolus )   TR_TRANS_G(n,4,m) = TR_TRANS_G(n,4,m) &
                                                    + trans_s(4,m)
         if ( lsubmeso )         TR_TRANS_G(n,5,m) = TR_TRANS_G(n,5,m) &
                                                    + trans_s(5,m)
       enddo
     enddo

!-----------------------------------------------------------------------
!
!  apply conversion factor to heat transport, prior to masking
!
!-----------------------------------------------------------------------
     conversion = 1.0e-19_r8/hflux_factor
     if (tracer_index == 1) TR_TRANS_G = TR_TRANS_G*conversion
 
!-----------------------------------------------------------------------
!
!  mask the transports
!
!  - any missing component is filled with undefined_nf.  
!  - because southern boundary diffusive transports are not available,
!    the total and diffusive transport components are not computed
!    for regions.  
!
!-----------------------------------------------------------------------

     do m=1,n_transport_reg
       do n=1,n_lat_aux_grid-1
         if ( MASK_LAT_DEPTH(n  ,1,m) .and. MASK_LAT_DEPTH(n+1,1,m) ) then
           TR_TRANS_G(n+1,:,m) = undefined_nf
         endif
       enddo
       if ( MASK_LAT_DEPTH(1,1,m) ) TR_TRANS_G(1,:,m) = undefined_nf
       if ( MASK_LAT_DEPTH(n_lat_aux_grid,1,m) ) then
           TR_TRANS_G(n_lat_aux_grid+1,:,m) = undefined_nf
       endif
     enddo

     do m=2,n_transport_reg
       TR_TRANS_G(:,1,m) = undefined_nf
       TR_TRANS_G(:,3,m) = undefined_nf
     enddo
 
     if (tracer_index == 1) TAVG_N_HEAT_TRANS_G = TR_TRANS_G
     if (tracer_index == 2) TAVG_N_SALT_TRANS_G = TR_TRANS_G

   endif
 
   deallocate ( WORK1, WORK1_G, WORK2_G)
   if ( ldiag_gm_bolus )  deallocate ( WORK3_G )
   if ( lsubmeso )        deallocate ( WORK4_G )
 
   call timer_stop (timer_tracer_transports)

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_tracer_transports

!***********************************************************************

 end module diags_on_lat_aux_grid
      
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

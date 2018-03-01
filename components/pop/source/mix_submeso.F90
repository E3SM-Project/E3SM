!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module mix_submeso 

!BOP
! !MODULE: mix_submeso 

! !DESCRIPTION:
!  This module contains routines for computing submesoscale mixing
!  using the Fox-Kemper, Ferrari, and Hallberg parameterization
!  for restratification by mixed layer eddies.

! !REVISION HISTORY:
!  SVN:$Id: mix_submeso.F90

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks
   use domain
   use constants
   use broadcast
   use grid
   use io
   use vertical_mix
   use vmix_kpp
   use time_management
   use tavg
   use exit_mod
   use registry
   use communicate
   use hmix_gm_submeso_share
   
#ifdef CCSMCOUPLED
   use shr_sys_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_submeso,   &
             submeso_sf,     &
	     submeso_flux 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  variables to save from one call to next
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ktp = 1, kbt = 2      ! refer to the top and bottom halves of a
                               !  grid cell, respectively
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SF_SUBM_X,  &       ! components of the submesoscale 
      SF_SUBM_Y           !  streamfunction
   real (r8), dimension(:,:,:,:), allocatable :: &
         FZTOP_SUBM 

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to submesoscale mixing.
!  Zonal and meridional refer here to logical space only.
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_USUBM,        &   ! zonal      submeso velocity
      tavg_VSUBM,        &   ! meridional submeso velocity 
      tavg_WSUBM,        &   ! vertical   submeso velocity
      tavg_ADVT_SUBM,    &   ! vertically-integrated T submeso 
                             !  advection tendency
      tavg_ADVS_SUBM,    &   ! vertically-integrated S submeso 
                             !  advection tendency
      tavg_VNT_SUBM,     &   ! heat flux tendency in grid-y direction
                             !  due to submeso velocity
      tavg_VNS_SUBM,     &   ! salt flux tendency in grid-y direction
                             !  due to submeso velocity
      tavg_HLS_SUBM          ! horizontal length scale used in horizontal
                             !  buoyancy gradient scaling in submeso

   real (r8), dimension(:,:,:), allocatable :: &
      TIME_SCALE             ! time scale used in horizontal length scale
                             !  calculation

   real (r8) :: &
      max_hor_grid_scale     ! maximum horizontal grid scale allowed

!-----------------------------------------------------------------------
!
!  namelist variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      luse_const_horiz_len_scale     ! if .true., then use a constant
                                     !  horizontal length scale given by
                                     !  hor_length_scale, otherwise the
                                     !  horizontal length scale varies both
                                     !  in space and time

   real (r8) :: &
      efficiency_factor,   &         ! 0.06 <= efficiency factor <= 0.08
      time_scale_constant, &         ! 1 day <= time scale constant <= 1 week
      hor_length_scale               ! constant horizontal length scale used
                                     !  if luse_const_horiz_len_scale is true.
                                     !  if luse_const_horiz_len_scale is false,
                                     !  then hor_length_scale is used as the 
                                     !  lower limit.

!-----------------------------------------------------------------------
!
!  misc module variables
!
!-----------------------------------------------------------------------
   real (r8) :: &
      sqrt_grav              ! sqrt(grav)
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_submeso
! !INTERFACE:

   subroutine init_submeso

! !DESCRIPTION:
!  Initializes various submesoscale mixing options and allocates necessary
!  space.  Also computes some time-independent arrays.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables and namelist
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error,         &   ! error flag for namelist
      iblock                 ! block index

   namelist /mix_submeso_nml/ efficiency_factor,             &
                              time_scale_constant,           &
                              luse_const_horiz_len_scale,    &
                              hor_length_scale

!-----------------------------------------------------------------------
!
!  register init_submeso
!
!-----------------------------------------------------------------------

   call register_string ('init_submeso')

!-----------------------------------------------------------------------
!
!  default namelist values 
!
!-----------------------------------------------------------------------

   efficiency_factor             = 0.07_r8
   time_scale_constant           = 3.456e5_r8       ! 4 days, in seconds
   luse_const_horiz_len_scale    = .false.
   hor_length_scale              = 5.0e5_r8         ! 5 km

   max_hor_grid_scale            = 111.0e5_r8       ! about 1 degree

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old', iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error =  1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=mix_submeso_nml, iostat=nml_error)
     enddo
     if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort, 'ERROR reading mix_submeso_nml')
   endif

   if (my_task == master_task) then

     write(stdout,*) ' '
     write(stdout,*) ' Document Namelist Parameters:'
     write(stdout,*) ' ============================ '
     write(stdout,*) ' '
     write(stdout,  mix_submeso_nml)
     write(stdout,*) ' '
     write(stdout,*) ' Submesoscale mixing options:'
     write(stdout,'(a21,1pe13.6)') ' efficiency factor = ', efficiency_factor
     write(stdout,'(a23,1pe13.6)') ' time scale constant = ', time_scale_constant
     if ( luse_const_horiz_len_scale ) then
       write(stdout,'(a45,1pe13.6)')  &
         ' using a constant horizontal length scale of ', hor_length_scale
     else
       write(stdout,'(a54)') ' horizontal length scale varies both in space and time'
     endif

     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif

   call broadcast_scalar (efficiency_factor,          master_task)
   call broadcast_scalar (time_scale_constant,        master_task)
   call broadcast_scalar (luse_const_horiz_len_scale, master_task)
   call broadcast_scalar (hor_length_scale,           master_task)

!-----------------------------------------------------------------------
!
!  define scalar constants
!
!-----------------------------------------------------------------------

   sqrt_grav  = sqrt(grav)

!-----------------------------------------------------------------------
!
!  allocate necessary arrays
!
!-----------------------------------------------------------------------

   allocate (SF_SUBM_X(nx_block,ny_block,2,2,km,nblocks_clinic),  &
             SF_SUBM_Y(nx_block,ny_block,2,2,km,nblocks_clinic))

   allocate (TIME_SCALE(nx_block,ny_block,nblocks_clinic))

   allocate (FZTOP_SUBM(nx_block,ny_block,nt,nblocks_clinic))

   SF_SUBM_X  = c0
   SF_SUBM_Y  = c0
   TIME_SCALE = c0
   FZTOP_SUBM = c0

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic
     TIME_SCALE(:,:,iblock) = c1 / sqrt( FCORT(:,:,iblock)**2  &
                            + c1 / (time_scale_constant**2) )
   enddo

!-----------------------------------------------------------------------
!
!  define tavg fields related to the submesoscale parameterization 
!
!-----------------------------------------------------------------------

   call define_tavg_field (tavg_USUBM, 'USUBM', 3,                  &
    long_name='Submeso velocity in grid-x direction (diagnostic)',  &
                 units='cm/s', grid_loc='3211',                     &
                 coordinates='ULONG TLAT z_t time')

   call define_tavg_field (tavg_VSUBM, 'VSUBM', 3,                  &
    long_name='Submeso velocity in grid-y direction (diagnostic)',  &
                 units='cm/s', grid_loc='3121',                     &
                 coordinates='TLONG ULAT z_t time')

   call define_tavg_field (tavg_WSUBM, 'WSUBM', 3,                &
    long_name='Vertical submeso velocity (diagnostic)',           &
                 units='cm/s', grid_loc='3112',                   &
                 coordinates='TLONG TLAT z_w time')

   call define_tavg_field (tavg_HLS_SUBM, 'HLS_SUBM', 2,          &
    long_name='Horizontal length scale used in submeso',          &
                 units='cm', grid_loc='2110',                     &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_ADVT_SUBM, 'ADVT_SUBM', 2,                       &
    long_name='Vertically-Integrated T submeso Advection Tendency (diagnostic)', &
                 units='cm degC/s', grid_loc='2110',                             &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_ADVS_SUBM, 'ADVS_SUBM', 2,                       &
    long_name='Vertically-Integrated S submeso Advection Tendency (diagnostic)', &
                 scale_factor=1000.0_rtavg,                                      &
                 units='cm gram/kilogram/s', grid_loc='2110',                    &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_VNT_SUBM, 'VNT_SUBM', 3,                          &
    long_name='Heat Flux Tendency in grid-y Dir due to submeso Vel (diagnostic)', &
                 units='degC/s', grid_loc='3121',                                 &
                 coordinates='TLONG ULAT z_t time')

   call define_tavg_field (tavg_VNS_SUBM, 'VNS_SUBM', 3,                          &
    long_name='Salt Flux Tendency in grid-y Dir due to submeso Vel (diagnostic)', &
                 scale_factor=1000.0_rtavg,                                       &
                 units='gram/kilogram/s', grid_loc='3121',                        &
                 coordinates='TLONG ULAT z_t time')

!-----------------------------------------------------------------------
!EOC

   end subroutine init_submeso

!***********************************************************************
!BOP
! !IROUTINE: submeso_sf 
! !INTERFACE:

   subroutine submeso_sf ( TMIX, this_block )

! !DESCRIPTION:
!  The Fox-Kemper, Ferrari, and Hallberg [2008] submesoscale parameterization
!  for restratification by mixed layer eddies.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block            ! block info for this sub block

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX                  ! tracers at all vertical levels
                            !   at mixing time level
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, k, kk,       &  ! dummy loop counters
      kp1,               &
      bid                   ! local block address for this sub block

   real (r8), dimension(nx_block,ny_block) :: &
      ML_DEPTH,          &  ! mixed layer depth
      HLS,               &  ! horizontal length scale
      WORK1, WORK2,      &  ! work arrays
      WORK3,             &
      USMT, VSMT,        &  ! arrays for submeso velocity computation
      USMB, VSMB,        &
      U_SUBM, V_SUBM,    &
      WBOT_SUBM,         &
      WTOP_SUBM

   real (r8), dimension(2) :: &
      reference_depth

   logical (log_kind), dimension(nx_block,ny_block) :: &
      CONTINUE_INTEGRAL     ! flag

   real (r8), dimension(nx_block,ny_block,2) :: &
      BX_VERT_AVG,       &  ! horizontal buoyancy differences vertically 
      BY_VERT_AVG           !  averaged within the mixed layer   

   real (r8) :: &
      zw_top, factor
      
!-----------------------------------------------------------------------
!
!  initialize various quantities
!
!-----------------------------------------------------------------------
         
   bid = this_block%local_id
   
   

   WORK1     = c0
   WORK2     = c0
   WORK3     = c0
   HLS       = c0
   USMT      = c0
   USMB      = c0
   VSMT      = c0
   VSMB      = c0
   U_SUBM    = c0
   V_SUBM    = c0
   WTOP_SUBM = c0
   WBOT_SUBM = c0

   BX_VERT_AVG = c0
   BY_VERT_AVG = c0

   SF_SUBM_X(:,:,:,:,:,bid) = c0
   SF_SUBM_Y(:,:,:,:,:,bid) = c0

   ML_DEPTH = zw(1)
   if ( vmix_itype == vmix_type_kpp )  &
     ML_DEPTH(:,:) = HMXL(:,:,bid) 
   

   CONTINUE_INTEGRAL = .true.
   where ( KMT(:,:,bid) == 0 ) 
     CONTINUE_INTEGRAL = .false.
   endwhere

!-----------------------------------------------------------------------
!
!  compute vertical averages of horizontal buoyancy differences 
!  within the mixed layer 
!
!-----------------------------------------------------------------------

   do k=1,km
   
     zw_top = c0
     if ( k > 1 )  zw_top = zw(k-1)

     WORK3 = c0
     where ( CONTINUE_INTEGRAL  .and.  ML_DEPTH > zw(k) )
        WORK3 = dz(k)
     endwhere
     where ( CONTINUE_INTEGRAL  .and.  ML_DEPTH <= zw(k)  &
             .and.  ML_DEPTH > zw_top )
       WORK3 = ML_DEPTH - zw_top
     endwhere

     where ( CONTINUE_INTEGRAL )
       BX_VERT_AVG(:,:,1) = BX_VERT_AVG(:,:,1)        &
                           + RX(:,:,1,k,bid) * WORK3
       BX_VERT_AVG(:,:,2) = BX_VERT_AVG(:,:,2)        &
                           + RX(:,:,2,k,bid) * WORK3
       BY_VERT_AVG(:,:,1) = BY_VERT_AVG(:,:,1)        &
                           + RY(:,:,1,k,bid) * WORK3
       BY_VERT_AVG(:,:,2) = BY_VERT_AVG(:,:,2)        &
                           + RY(:,:,2,k,bid) * WORK3
     endwhere

     where ( CONTINUE_INTEGRAL .and.  ML_DEPTH <= zw(k)  &
             .and.  ML_DEPTH > zw_top )
       CONTINUE_INTEGRAL = .false.
     endwhere  

   enddo

#ifdef CCSMCOUPLED
   if ( any(CONTINUE_INTEGRAL) ) then
     call shr_sys_abort ('Incorrect mixed layer depth in submeso subroutine (I)')
   endif
#endif

   where ( KMT(:,:,bid) > 0 )
     BX_VERT_AVG(:,:,1) = - grav * BX_VERT_AVG(:,:,1) / ML_DEPTH
     BX_VERT_AVG(:,:,2) = - grav * BX_VERT_AVG(:,:,2) / ML_DEPTH
     BY_VERT_AVG(:,:,1) = - grav * BY_VERT_AVG(:,:,1) / ML_DEPTH
     BY_VERT_AVG(:,:,2) = - grav * BY_VERT_AVG(:,:,2) / ML_DEPTH
   endwhere

!-----------------------------------------------------------------------
!
!  compute horizontal length scale if necessary
!
!-----------------------------------------------------------------------

   if ( luse_const_horiz_len_scale ) then

     where ( KMT(:,:,bid) > 0 ) 
       HLS = hor_length_scale
     endwhere

   else

     WORK1 = c0

     where ( KMT(:,:,bid) > 0 )
       WORK1 = sqrt( p5 * (                                       &
               ( BX_VERT_AVG(:,:,1)**2 + BX_VERT_AVG(:,:,2)**2 )  &
                 / DXT(:,:,bid)**2                                &
             + ( BY_VERT_AVG(:,:,1)**2 + BY_VERT_AVG(:,:,2)**2 )  &
                 / DYT(:,:,bid)**2 ) )
       WORK1 = WORK1 * ML_DEPTH * (TIME_SCALE(:,:,bid)**2)
     endwhere

     CONTINUE_INTEGRAL = .true.
     where ( KMT(:,:,bid) == 0 ) 
       CONTINUE_INTEGRAL = .false.
     endwhere

     WORK2 = c0

     do k=2,km

       WORK3 = c0
       where ( CONTINUE_INTEGRAL  .and.  ML_DEPTH > zt(k) )
         WORK3 = dzw(k-1) 
       endwhere
       where ( CONTINUE_INTEGRAL  .and.  ML_DEPTH <= zt(k)  &
            .and.  ML_DEPTH >= zt(k-1) )
         WORK3 = ( (ML_DEPTH - zt(k-1))**2 ) * dzwr(k-1)
       endwhere

       where ( CONTINUE_INTEGRAL )
         WORK2 = WORK2 + sqrt(-RZ_SAVE(:,:,k,bid) * WORK3)
       endwhere

       where ( CONTINUE_INTEGRAL .and.  ML_DEPTH <= zt(k)  &
               .and.  ML_DEPTH >= zt(k-1) )
         CONTINUE_INTEGRAL = .false.
       endwhere

     enddo

#ifdef CCSMCOUPLED
     if ( any(CONTINUE_INTEGRAL) ) then
       call shr_sys_abort ('Incorrect mixed layer depth in submeso subroutine (II)')
     endif
#endif

     where ( KMT(:,:,bid) > 0 )

       WORK2 = sqrt_grav * WORK2 * TIME_SCALE(:,:,bid)

       HLS = max ( WORK1, WORK2, hor_length_scale )

     endwhere

   endif

!-----------------------------------------------------------------------
!
!  compute streamfunction due to submesoscale parameterization 
!
!-----------------------------------------------------------------------

   do k=1,km

     reference_depth(ktp) = zt(k) - p25 * dz(k)
     reference_depth(kbt) = zt(k) + p25 * dz(k)

     do kk=ktp,kbt

       where ( reference_depth(kk) < ML_DEPTH  .and.  &
            KMT(:,:,bid) >= k )

         WORK3 = ( c1 - ( c2 * reference_depth(kk) / ML_DEPTH ) )**2
            
         WORK2 = ( c1 - WORK3 )  &
               * ( c1 + ( 5.0_r8 / 21.0_r8 ) * WORK3 )

         WORK1 = efficiency_factor * (ML_DEPTH**2) * WORK2  &
                * TIME_SCALE(:,:,bid) / HLS

!     in the following negative sign is omitted to be consistent with
!     the GM implementation in hmix_gm subroutine. also, DXT and
!     DYT usage is approximate. 

         SF_SUBM_X(:,:,1,kk,k,bid) = WORK1 * BX_VERT_AVG(:,:,1)  &
                                    * min(DXT(:,:,bid),max_hor_grid_scale)
         SF_SUBM_X(:,:,2,kk,k,bid) = WORK1 * BX_VERT_AVG(:,:,2)  &
                                    * min(DXT(:,:,bid),max_hor_grid_scale)
         SF_SUBM_Y(:,:,1,kk,k,bid) = WORK1 * BY_VERT_AVG(:,:,1)  &
                                    * min(DYT(:,:,bid),max_hor_grid_scale)
         SF_SUBM_Y(:,:,2,kk,k,bid) = WORK1 * BY_VERT_AVG(:,:,2)  &
                                    * min(DYT(:,:,bid),max_hor_grid_scale)

       endwhere

     enddo

   enddo

   USMT = c0
   VSMT = c0

   do k=1,km

!-----------------------------------------------------------------------
!
!  diagnostic computation of the submeso velocities
!
!-----------------------------------------------------------------------

     kp1 = k+1
     factor = c1
     if ( k == km ) then
       kp1 = k
       factor = c0
     endif

     do j=1,ny_block-1
       do i=1,nx_block-1

         WORK1(i,j) = (   SF_SUBM_X(i  ,j,1,kbt,k,  bid)    &
               + factor * SF_SUBM_X(i  ,j,1,ktp,kp1,bid)    &
                        + SF_SUBM_X(i+1,j,2,kbt,k,  bid)    &
               + factor * SF_SUBM_X(i+1,j,2,ktp,kp1,bid) )  &
                * p25 * HYX(i,j,bid)

         WORK2(i,j) = (   SF_SUBM_Y(i,j  ,1,kbt,k,  bid)    &
               + factor * SF_SUBM_Y(i,j  ,1,ktp,kp1,bid)    &
                        + SF_SUBM_Y(i,j+1,2,kbt,k,  bid)    &
               + factor * SF_SUBM_Y(i,j+1,2,ktp,kp1,bid) )  &
               * p25 * HXY(i,j,bid)

       enddo
     enddo

     USMB = merge( WORK1, c0, k < KMT(:,:,bid) .and. k < KMTE(:,:,bid) )
     VSMB = merge( WORK2, c0, k < KMT(:,:,bid) .and. k < KMTN(:,:,bid) )

     WORK1 = merge( USMT - USMB, c0, k <= KMT(:,:,bid)  &
                               .and. k <= KMTE(:,:,bid) )
     WORK2 = merge( VSMT - VSMB, c0, k <= KMT(:,:,bid)  &
                               .and. k <= KMTN(:,:,bid) )

     U_SUBM = WORK1 * dzr(k) / HTE(:,:,bid)
     V_SUBM = WORK2 * dzr(k) / HTN(:,:,bid)

     do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie

         if ( k < KMT(i,j,bid)  .and.  ( zw(k) < max( ML_DEPTH(i,j),  &
              ML_DEPTH(i+1,j), ML_DEPTH(i-1,j), ML_DEPTH(i,j+1),      &
              ML_DEPTH(i,j-1) ) ) ) then
           WBOT_SUBM(i,j) = WTOP_SUBM(i,j)             &
                           + TAREA_R(i,j,bid)          &
                       * ( WORK1(i,j) - WORK1(i-1,j)   &
                         + WORK2(i,j) - WORK2(i,j-1) )
         else
           WBOT_SUBM(i,j) = c0
         endif

       enddo
     enddo

!-----------------------------------------------------------------------
!
!  accumulate time average if necessary; checking is internal to routine
!
!-----------------------------------------------------------------------

     if ( mix_pass /= 1 ) then

       if ( k == 1 ) then
         call accumulate_tavg_field (HLS, tavg_HLS_SUBM, bid, 1)  
       endif

       call accumulate_tavg_field (U_SUBM, tavg_USUBM, bid, k)
       call accumulate_tavg_field (V_SUBM, tavg_VSUBM, bid, k)
       call accumulate_tavg_field (WTOP_SUBM, tavg_WSUBM, bid, k) 

       if (accumulate_tavg_now(tavg_ADVT_SUBM) ) then

         WORK1 = p5 * HTE(:,:,bid) * U_SUBM * ( TMIX(:,:,k,1)  &
                   + eoshift(TMIX(:,:,k,1), dim=1, shift=1) )
         WORK2 = eoshift(WORK1, dim=1, shift=-1)
         WORK3 = WORK1 - WORK2

         WORK1 = p5 * HTN(:,:,bid) * V_SUBM * ( TMIX(:,:,k,1)  &
                   + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )
         WORK2 = eoshift(WORK1, dim=2, shift=-1)
         WORK3 = WORK3 + WORK1 - WORK2

         WORK1 = c0
         do j=this_block%jb,this_block%je
           do i=this_block%ib,this_block%ie
             if ( k <= KMT(i,j,bid) ) then
               WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
             endif
           enddo
         enddo

         call accumulate_tavg_field (WORK1, tavg_ADVT_SUBM, bid, k)

       endif

       if (accumulate_tavg_now(tavg_ADVS_SUBM) ) then

         WORK1 = p5 * HTE(:,:,bid) * U_SUBM * ( TMIX(:,:,k,2)  &
                   + eoshift(TMIX(:,:,k,2), dim=1, shift=1) )
         WORK2 = eoshift(WORK1, dim=1, shift=-1)
         WORK3 = WORK1 - WORK2

         WORK1 = p5 * HTN(:,:,bid) * V_SUBM * ( TMIX(:,:,k,2)  &
                   + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
         WORK2 = eoshift(WORK1, dim=2, shift=-1)
         WORK3 = WORK3 + WORK1 - WORK2

         WORK1 = c0
         do j=this_block%jb,this_block%je
           do i=this_block%ib,this_block%ie
             if ( k <= KMT(i,j,bid) ) then
               WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
             endif
           enddo
         enddo

         call accumulate_tavg_field (WORK1, tavg_ADVS_SUBM, bid, k)

       endif

       if ( accumulate_tavg_now(tavg_VNT_SUBM)  .or.  &
            accumulate_tavg_now(tavg_VNS_SUBM) ) then

         WORK1 = p5 * V_SUBM * HTN(:,:,bid) * TAREA_R(:,:,bid)

         if (accumulate_tavg_now(tavg_VNT_SUBM) ) then

           WORK2 = WORK1 * (    TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )

           call accumulate_tavg_field (WORK2, tavg_VNT_SUBM, bid, k)

         endif

         if (accumulate_tavg_now(tavg_VNS_SUBM) ) then

           WORK2 = WORK1 * (    TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )

           call accumulate_tavg_field (WORK2, tavg_VNS_SUBM, bid, k)

         endif

       endif

     endif ! mix_pass ne 1

!-----------------------------------------------------------------------
!
!  update remaining bottom-face fields to top-face fields for next pass
!
!-----------------------------------------------------------------------

     USMT = USMB
     VSMT = VSMB

     WTOP_SUBM = WBOT_SUBM

   enddo

!-----------------------------------------------------------------------
!EOC

   end subroutine submeso_sf 

!***********************************************************************

! !IROUTINE: submeso_flux
! !INTERFACE:

   subroutine submeso_flux (k, GTK, TMIX, tavg_HDIFE_TRACER, &
                     tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Computes the fluxes due to submesoscale mixing of tracers
!
!INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index
   
   type (block), intent(in) :: &
      this_block            ! block info for this sub block
   
   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level   

   integer (int_kind), dimension(nt), intent(in) :: &
      tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
      tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
      tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer
   
! !OUTPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,nt), intent(out) :: &
         GTK     ! submesoscale tracer flux at level k
      
!-----------------------------------------------------------------------
!
!      local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
	 jnorth = 1, jsouth = 2
      integer (int_kind)  :: &
	 i,j,n,                &!dummy loop counters
	 bid,                  &! local block address for this sub block
	 kp1
      
      real (r8) :: &
         fz, factor 
	 
      real (r8), dimension(nx_block,ny_block) :: &
         CX, CY,                  &
	 WORK1, WORK2,            &! local work space
	 KMASK                     ! ocean mask
	 
      real (r8), dimension(nx_block,ny_block,nt)  :: &
         FX, FY                    ! fluxes across east, north faces

     bid = this_block%local_id

     if ( k == 1) FZTOP_SUBM(:,:,:,bid) = c0        ! zero flux B.C. at the surface
     
      CX = merge(HYX(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTE(:,:,bid)))
      CY = merge(HXY(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTN(:,:,bid)))
      
      KMASK = merge(c1, c0, k < KMT(:,:,bid))
            
      kp1 = k + 1
      if ( k == km )  kp1 = k

      if ( k < km ) then
        factor    = c1
      else
        factor    = c0
      endif
      
	do n = 1,nt
          if (.not.registry_match('init_gm')) then
           if (n > 2 .and. k < km)  &
             TZ(:,:,k+1,n,bid) = TMIX(:,:,k  ,n) - TMIX(:,:,k+1,n)
          endif

          do j=1,ny_block
            do i=1,nx_block-1
              FX(i,j,n) = CX(i,j)                          &
               * ( SF_SUBM_X(i  ,j,ieast,ktp,k,bid) * TZ(i,j,k,n,bid)                        &
                 + SF_SUBM_X(i  ,j,ieast,kbt,k,bid) * TZ(i,j,kp1,n,bid)                    &
                 + SF_SUBM_X(i+1,j,iwest,ktp,k,bid) * TZ(i+1,j,k,n,bid)                    &
                 + SF_SUBM_X(i+1,j,iwest,kbt,k,bid) * TZ(i+1,j,kp1,n,bid) )
            enddo
          enddo

        end do
	
       do n = 1,nt

          do j=1,ny_block-1
            do i=1,nx_block
              FY(i,j,n) =  CY(i,j)                          &
               * ( SF_SUBM_Y(i,j  ,jnorth,ktp,k,bid) * TZ(i,j,k,n,bid)                        &
                 + SF_SUBM_Y(i,j  ,jnorth,kbt,k,bid) * TZ(i,j,kp1,n,bid)                    &
                 + SF_SUBM_Y(i,j+1,jsouth,ktp,k,bid) * TZ(i,j+1,k,n,bid)                    &
                 + SF_SUBM_Y(i,j+1,jsouth,kbt,k,bid) * TZ(i,j+1,kp1,n,bid) )
            enddo
          enddo

       end do
	
      ! zero halo regions so accumulate_tavg_field calls do not trap
      WORK1 = c0
      WORK2 = c0
	    
      do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate vertical submesoscale fluxes thru horizontal faces of T-cell
!
!-----------------------------------------------------------------------
 
        GTK(:,:,n) = c0

        if ( k < km ) then

              do j=this_block%jb,this_block%je
                do i=this_block%ib,this_block%ie
               
                  WORK1(i,j) = SF_SUBM_X(i  ,j  ,ieast ,kbt,k  ,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,kbt,k  ,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,kbt,k  ,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,kbt,k  ,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,k  ,n,bid)
                enddo
              enddo

              do j=this_block%jb,this_block%je
                do i=this_block%ib,this_block%ie
                  WORK2(i,j) = factor                                    &
                           * ( SF_SUBM_X(i  ,j  ,ieast ,ktp,kp1,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,ktp,kp1,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,ktp,kp1,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,ktp,kp1,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1,n,bid) ) 
                enddo
              enddo

	    
	  do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                fz = -KMASK(i,j) * p25                                &
                     * (WORK1(i,j) + WORK2(i,j))

                GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                             + FY(i,j,n) - FY(i,j-1,n)  &
                      + FZTOP_SUBM(i,j,n,bid) - fz )*dzr(k)*TAREA_R(i,j,bid)

                FZTOP_SUBM(i,j,n,bid) = fz

              enddo
          enddo
	

	 else                 ! k = km

          do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie

              GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                           + FY(i,j,n) - FY(i,j-1,n)  &
                    + FZTOP_SUBM(i,j,n,bid) )*dzr(k)*TAREA_R(i,j,bid)

              FZTOP_SUBM(i,j,n,bid) = c0 

            enddo
          enddo

        endif

!-----------------------------------------------------------------------
!
!     accumulate time average if necessary
!
!-----------------------------------------------------------------------

      if ( mix_pass /= 1 ) then
          if (accumulate_tavg_now(tavg_HDIFE_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FX(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FY(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FZTOP_SUBM(i,j,n,bid)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            call accumulate_tavg_field(WORK1,tavg_HDIFB_TRACER(n),bid,k)
          endif
      endif   ! mix_pass ne 1

	   
!-----------------------------------------------------------------------
!
!     end of tracer loop
!
!-----------------------------------------------------------------------

      enddo !ends the do n loop 

!-----------------------------------------------------------------------

   end subroutine submeso_flux 

!***********************************************************************

 end module mix_submeso

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

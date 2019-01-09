!#define single_run
! JUNG: comment out for use in CAM branch

MODULE q3d_interface_module
! Subprograms control the VVM channels

USE shr_kind_mod,    only: r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs
USE vGCM_data_types, only: vGCM_state_t,vGCM_tend_t,vGCM_out_t,channel_vGCM_t
USE vvm_data_types,  only: channel_t,allocate_channel_data

USE parmsld, only: nVGCM_seg,nLevel, &
                   nhalo_vGCM    ! JUNG: remove after debugging

IMPLICIT NONE
PRIVATE

public :: crm_init

#ifndef single_run
public :: crm_init_diag
#endif

public :: crm_run
public :: crm_final

CONTAINS

! Local Subroutines:
!------------------------------------------------------------------------------
! Subroutine crm_init         : Initialize data
! Subroutine crm_init_diag    : Initialize diagnostic output
!
! Subroutine crm_run          : Control the CRM run
!   Subroutine crm_prediction : Control the time marching
!   Subroutine crm_out_tend   : Prepare tendency terms (feedback to GCM)
!   Subroutine crm_out        : Prepare CRM data for output (for single_run)
!
! Subroutine crm_final        : For any needed cleanup or output
!
! Subroutine get_xl : Arrange data for the vGCM_out grid structure (x-channel with layers)
! Subroutine get_yl : Arrange data for the vGCM_out grid structure (y-channel with layers)
! Subroutine get_x  : Arrange data for the vGCM_out grid structure (x-channel)
! Subroutine get_y  : Arrange data for the vGCM_out grid structure (y-channel)
! Subroutine vvm_layers_to_vgcm_layers : Vertical interpolation
!------------------------------------------------------------------------------

!==============================================================================
  subroutine crm_init (begchan, endchan, channel, vGCM)
!==============================================================================
    USE parmsld,  only: netsz,nhalo_cal,nhalo_adv,nhalo,nk2,nk2_ext,channel_seg_l
    USE timeinfo, only: rjday0

    USE constld,  only: rad2deg       ! JUNG: remove after debugging

    USE q3d_init_module,    only: init_common
    USE q3d_mapping_module, only: find_channel_info,mapping_channel
    USE time_manager,       only: get_start_date,timemgr_get_calendar_cf  ! share with CAM
    USE shr_cal_mod,        only: shr_cal_ymd2julian                      ! share with CAM

    INTEGER, INTENT(IN) :: begchan    ! first channel on this node
    INTEGER, INTENT(IN) :: endchan    ! last channel on this node

    type(channel_t),      pointer, intent(inout) :: channel(:) ! CRM data for multiple channels
    type(channel_vGCM_t), pointer, intent(inout) :: vGCM(:)    ! vGCM data for multiple channels

!   Local
    LOGICAL :: TEST_WRITE = .FALSE.               ! JUNG: remove after debugging (vGCM)
    LOGICAL :: TEST_WRITE_CRM = .FALSE.           ! JUNG: remove after debugging (CRM)
    INTEGER :: ifort_fix = 90                     ! JUNG: remove after debugging  

    INTEGER :: num_seg, num_gcm, k

    INTEGER :: ifortw, jj                   ! JUNG: remove after debugging

    character (len=cs) :: calendar
    INTEGER :: year0, month0, day0, tod0

!   Find Julian day of the initial timestep (rjday0)
    call get_start_date(year0, month0, day0, tod0)
    calendar = timemgr_get_calendar_cf()
    call shr_cal_ymd2julian(year0, month0, day0, tod0, rjday0, trim(calendar))

!   Initialize common parameters and vertical profiles in CONSTLD.f90
    CALL INIT_common

!--------------------------------------------- JUNG_DEBUG
    IF (TEST_WRITE) THEN
    ! write information for all channels
    
     DO k = begchan, endchan
       ifortw = k
       if (ifortw == 5) ifortw=55
       if (ifortw == 6) ifortw=65
       write(ifortw,*)
       write(ifortw,'(a15,i5,2(a3,i3),a5,i6,a8,f8.3)') &
         'Start date: yr=',year0,' m=',month0,' d=',day0,' sec=',tod0,' rjday0=',rjday0
       write(ifortw,'(a4,i3,2x,a8,i3)') 'nk2=',nk2,'nk2_ext=',nk2_ext
       write(ifortw,*)
     ENDDO
     
    ELSE   ! TEST_WRITE
     ! write information for only one channel
     
     IF (begchan == 1) THEN
       write(ifort_fix,*)
       write(ifort_fix,'(a15,i5,2(a3,i3),a5,i6,a8,f8.3)') &
         'Start date: yr=',year0,' m=',month0,' d=',day0,' sec=',tod0,' rjday0=',rjday0
       write(ifort_fix,'(a4,i3,2x,a8,i3)') 'nk2=',nk2,'nk2_ext=',nk2_ext
       write(ifort_fix,*)
     ENDIF 
     
    ENDIF  ! TEST_WRITE 
!--------------------------------------------- JUNG_DEBUG    

    allocate(channel(begchan:endchan))

!   FIND_CHANNEL_INFO is recalculated at RESTART
    do k = begchan, endchan

      ! Find group-ID of each channel, and size & face number of each channel segment
      CALL FIND_CHANNEL_INFO (vGCM(k)%vGCM_state,vGCM(k)%vGCM_map,channel(k))
      
      channel(k)%num_chn = k
      
    enddo

    CALL allocate_channel_data(begchan, endchan, channel)

!   The below is recalculated at RESTART
    ! Set up segment parameters: channel(:)%seg(:)

    DO k = begchan, endchan
    DO num_seg = 1, 4

    channel(k)%seg(num_seg)%mi1 = channel(k)%seg(num_seg)%nx_size
    channel(k)%seg(num_seg)%mj1 = channel(k)%seg(num_seg)%ny_size

    channel(k)%seg(num_seg)%mim    = 1 - nhalo
    channel(k)%seg(num_seg)%mim_c  = 1 - nhalo_cal
    channel(k)%seg(num_seg)%mim_ce = channel(k)%seg(num_seg)%mim_c - 1
    channel(k)%seg(num_seg)%mim_a  = 1 - nhalo_adv

    channel(k)%seg(num_seg)%mip    = channel(k)%seg(num_seg)%mi1 + nhalo
    channel(k)%seg(num_seg)%mip_c  = channel(k)%seg(num_seg)%mi1 + nhalo_cal
    channel(k)%seg(num_seg)%mip_ce = channel(k)%seg(num_seg)%mip_c + 1
    channel(k)%seg(num_seg)%mip_a  = channel(k)%seg(num_seg)%mi1 + nhalo_adv

    channel(k)%seg(num_seg)%mjm    = 1 - nhalo
    channel(k)%seg(num_seg)%mjm_c  = 1 - nhalo_cal
    channel(k)%seg(num_seg)%mjm_ce = channel(k)%seg(num_seg)%mjm_c - 1
    channel(k)%seg(num_seg)%mjm_a  = 1 - nhalo_adv

    channel(k)%seg(num_seg)%mjp    = channel(k)%seg(num_seg)%mj1 + nhalo
    channel(k)%seg(num_seg)%mjp_c  = channel(k)%seg(num_seg)%mj1 + nhalo_cal
    channel(k)%seg(num_seg)%mjp_ce = channel(k)%seg(num_seg)%mjp_c + 1
    channel(k)%seg(num_seg)%mjp_a  = channel(k)%seg(num_seg)%mj1 + nhalo_adv

    DO num_gcm = 1, nVGCM_seg
     channel(k)%seg(num_seg)%lsta(num_gcm) = 1 + (num_gcm-1)*netsz
     channel(k)%seg(num_seg)%lend(num_gcm) = netsz + (num_gcm-1)*netsz
    ENDDO

    DO num_gcm = 1, nVGCM_seg
     channel(k)%seg(num_seg)%lcen(num_gcm) = &
                0.5_r8*(channel(k)%seg(num_seg)%lsta(num_gcm) &
               +channel(k)%seg(num_seg)%lend(num_gcm))
    ENDDO

!--------------------------------------------- JUNG_DEBUG
    IF (TEST_WRITE) THEN
    ! write information for all channels
     
     ifortw = k
     if (ifortw == 5) ifortw=55
     if (ifortw == 6) ifortw=65

     write(ifortw,'(a8,i5,2x,a8,i2,2x,a5,i2,2x,a6,i2)') &
         'chn_num=',k,'num_seg=',num_seg,'G_id=',channel(k)%num_chg, &
         'nface=',channel(k)%seg(num_seg)%nface
     write(ifortw,*)

     DO jj = -nhalo_vGCM,nhalo_vGCM
     write(ifortw,'(a10,i3)') 'chn_cross=',jj
     DO num_gcm = 1-nhalo_vGCM, nVGCM_seg+nhalo_vGCM
     write(ifortw,'(a7,i4,2x,a4,F8.2,2x,a4,F8.2)') 'n_vGCM=',num_gcm, &
             'lon=',vGCM(k)%vGCM_state(num_seg)%lon(jj,num_gcm)*rad2deg,&
             'lat=',vGCM(k)%vGCM_state(num_seg)%lat(jj,num_gcm)*rad2deg
     ENDDO
     write(ifortw,*)
     ENDDO

     write(ifortw,*)
     write(ifortw,'(a4,i5,2x,a4,i5)') &
         'mi1=',channel(k)%seg(num_seg)%mi1,'mj1=',channel(k)%seg(num_seg)%mj1
     write(ifortw,*)

     DO num_gcm = 1, nVGCM_seg
      write(ifortw,'(a8,i3,2i5,F8.2)') 'num_gcm=',num_gcm, &
                                      channel(k)%seg(num_seg)%lsta(num_gcm), &
                                      channel(k)%seg(num_seg)%lend(num_gcm), &
                                      channel(k)%seg(num_seg)%lcen(num_gcm)
     ENDDO
     write(ifortw,*)

    ELSE   ! TEST_WRITE
    ! write information for only one channel
    
     IF (k == 1) THEN
       write(ifort_fix,'(a8,i5,2x,a8,i2,2x,a5,i2,2x,a6,i2)') &
          'chn_num=',k,'num_seg=',num_seg,'G_id=',channel(k)%num_chg, &
          'nface=',channel(k)%seg(num_seg)%nface
       write(ifort_fix,*)

       DO jj = -nhalo_vGCM,nhalo_vGCM
       write(ifort_fix,'(a10,i3)') 'chn_cross=',jj
       DO num_gcm = 1-nhalo_vGCM, nVGCM_seg+nhalo_vGCM
       write(ifort_fix,'(a7,i4,2x,a4,F8.2,2x,a4,F8.2)') 'n_vGCM=',num_gcm, &
             'lon=',vGCM(k)%vGCM_state(num_seg)%lon(jj,num_gcm)*rad2deg,&
             'lat=',vGCM(k)%vGCM_state(num_seg)%lat(jj,num_gcm)*rad2deg
       ENDDO
       write(ifort_fix,*)
       ENDDO

       write(ifort_fix,*)
       write(ifort_fix,'(a4,i5,2x,a4,i5)') &
         'mi1=',channel(k)%seg(num_seg)%mi1,'mj1=',channel(k)%seg(num_seg)%mj1
       write(ifort_fix,*)

       DO num_gcm = 1, nVGCM_seg
        write(ifort_fix,'(a8,i3,2i5,F8.2)') 'num_gcm=',num_gcm, &
                                      channel(k)%seg(num_seg)%lsta(num_gcm), &
                                      channel(k)%seg(num_seg)%lend(num_gcm), &
                                      channel(k)%seg(num_seg)%lcen(num_gcm)
       ENDDO
       write(ifort_fix,*)     
     ENDIF 
     
    ENDIF  ! TEST_WRITE
!--------------------------------------------- JUNG_DEBUG    

    ENDDO   ! num_seg-loop
    ENDDO   ! k-loop

!------------------------------------------
    do k = begchan, endchan
!------------------------------------------    

      ! Prepare mapping information.
      CALL MAPPING_CHANNEL (vGCM(k)%vGCM_state, channel(k))

!------------------------------------- JUNG: remove after debugging
      IF (TEST_WRITE_CRM) THEN
      ! write information for all channels
      
       ifortw = k
       if (ifortw == 5) ifortw=55
       if (ifortw == 6) ifortw=65

       DO num_seg = 1, 4
       write(ifortw,'(a8,i5,2x,a8,i2,2x,a5,i2,2x,a6,i2)') &
           'chn_num=',k,'num_seg=',num_seg,'G_id=',channel(k)%num_chg, &
           'nface=',channel(k)%seg(num_seg)%nface
       write(ifortw,*)

       write(ifortw,'(4(a4,i5,2x))') 'mim=',channel(k)%seg(num_seg)%mim, &
                                     'mip=',channel(k)%seg(num_seg)%mip, &
                                     'mjm=',channel(k)%seg(num_seg)%mjm, &
                                     'mjp=',channel(k)%seg(num_seg)%mjp

       DO jj = channel(k)%seg(num_seg)%mjm,channel(k)%seg(num_seg)%mjp
       write(ifortw,'(a7,i5)') 'n_crmj=',jj
       DO num_gcm = channel(k)%seg(num_seg)%mim, channel(k)%seg(num_seg)%mip
       write(ifortw,'(a7,i5,2x,a4,F8.2,2x,a4,F8.2)') 'n_crmi=',num_gcm, &
               'lon=',channel(k)%seg(num_seg)%rlon_t(num_gcm,jj)*rad2deg,&
               'lat=',channel(k)%seg(num_seg)%rlat_t(num_gcm,jj)*rad2deg
       ENDDO
       write(ifortw,*)
       ENDDO

       write(ifortw,*)
       write(ifortw,'(a4,i5,2x,a4,i5)') &
           'mi1=',channel(k)%seg(num_seg)%mi1,'mj1=',channel(k)%seg(num_seg)%mj1
       write(ifortw,*)

       ENDDO   ! num_seg-loop

      ELSE   ! TEST_WRITE_CRM
      ! write information for only one channel
      
       IF (k == 1) THEN
       DO num_seg = 1, 4
       write(ifort_fix,'(a8,i5,2x,a8,i2,2x,a5,i2,2x,a6,i2)') &
           'chn_num=',k,'num_seg=',num_seg,'G_id=',channel(k)%num_chg, &
           'nface=',channel(k)%seg(num_seg)%nface
       write(ifort_fix,*)

       write(ifort_fix,'(4(a4,i5,2x))') 'mim=',channel(k)%seg(num_seg)%mim, &
                                     'mip=',channel(k)%seg(num_seg)%mip, &
                                     'mjm=',channel(k)%seg(num_seg)%mjm, &
                                     'mjp=',channel(k)%seg(num_seg)%mjp

       DO jj = channel(k)%seg(num_seg)%mjm,channel(k)%seg(num_seg)%mjp
       write(ifort_fix,'(a7,i5)') 'n_crmj=',jj
       DO num_gcm = channel(k)%seg(num_seg)%mim, channel(k)%seg(num_seg)%mip
       write(ifort_fix,'(a7,i5,2x,a4,F8.2,2x,a4,F8.2)') 'n_crmi=',num_gcm, &
               'lon=',channel(k)%seg(num_seg)%rlon_t(num_gcm,jj)*rad2deg,&
               'lat=',channel(k)%seg(num_seg)%rlat_t(num_gcm,jj)*rad2deg
       ENDDO
       write(ifort_fix,*)
       ENDDO

       write(ifort_fix,*)
       write(ifort_fix,'(a4,i5,2x,a4,i5)') &
           'mi1=',channel(k)%seg(num_seg)%mi1,'mj1=',channel(k)%seg(num_seg)%mj1
       write(ifort_fix,*)  
       
       ENDDO   ! num_seg-loop      
       ENDIF 

      ENDIF  ! TEST_WRITE_CRM

      IF (TEST_WRITE)  THEN
      ! write information for all channels  
      
        ifortw = k
        if (ifortw == 5) ifortw=55
        if (ifortw == 6) ifortw=65

        CALL OUTCON (ifortw)
        
      ELSE   ! TEST_WRITE
      ! write information for only one channel 
      
      IF (k == 1) CALL OUTCON (ifort_fix)
      
      ENDIF  ! TEST_WRITE
!-------------------------------------

      ! Set the horizontal coordinates (lon/lat) of CRM grids,
      ! which fit to the "crm_grid" output.
      DO num_seg = 1, 4
        if (channel(k)%seg(num_seg)%nx_size > 1) then
           vGCM(k)%vGCM_out(num_seg)%clon(1:channel_seg_l) = &
                       channel(k)%seg(num_seg)%RLON_T(1:channel_seg_l,1)
           vGCM(k)%vGCM_out(num_seg)%clat(1:channel_seg_l) = &
                       channel(k)%seg(num_seg)%RLAT_T(1:channel_seg_l,1)
        else
           vGCM(k)%vGCM_out(num_seg)%clon(1:channel_seg_l) = &
                       channel(k)%seg(num_seg)%RLON_T(1,1:channel_seg_l)
           vGCM(k)%vGCM_out(num_seg)%clat(1:channel_seg_l) = &
                       channel(k)%seg(num_seg)%RLAT_T(1,1:channel_seg_l)
        endif
      ENDDO   ! num_seg-loop
!--------------------
    enddo  ! k-loop
!--------------------    

  end subroutine crm_init

#ifndef single_run
!   JUNG: calculation in CAM branch

!==============================================================================
  subroutine crm_init_diag ()
!==============================================================================
! Initialize diagnostic fields
!------------------------------------------------------------------------
     use cam_history_support, only: add_hist_coord
     use cam_history,         only: addfld, add_default, horiz_only
     use parmsld,             only: nhalo_vGCM,nk1,nk2,nk1_ext,nk2_ext
     use constld,             only: zz,zt,pbar,pbarz,rho,rhoz

     integer              :: halo_dim_size
     integer              :: index
     integer,  pointer    :: halo_index(:)
     real(r8), parameter  :: R_UNDEF = -1.0E30 ! Missing value for diags

     ! Setup a dimension for output of variables which include the halo
     halo_dim_size = (2 * nhalo_vGCM) + 1
     allocate(halo_index(halo_dim_size))

     do index = 0, halo_dim_size - 1
        halo_index(index + 1) = index - nhalo_vGCM
     end do

     call add_hist_coord('halo', halo_dim_size, 'vGCM halo indices',               &
          '1', halo_index)
     nullify(halo_index) ! Belongs to coord now

!---------------------------------
!    Add reference profiles
!---------------------------------
     call add_hist_coord('lev_c', nk1, 'CRM level at midpoints',                   &
          'm', zt(2:nk2))

     call add_hist_coord('ilev_c', nk2, 'CRM level at interfaces',                 &
          'm', zz(1:nk2))

     call add_hist_coord('lev_c_ext', nk1_ext, 'CRM level at midpoints (ext)',     &
          'm', zt(2:nk2_ext))

     call add_hist_coord('ilev_c_ext', nk2_ext, 'CRM level at interfaces (ext)',   &
          'm', zz)

     call add_hist_coord('pbar_ext', nk1_ext, 'CRM pressure at midpoints (ext)',   &
          'Pa', pbar(2:nk2_ext))

     call add_hist_coord('pbarz_ext', nk2_ext, 'CRM pressure at interfaces (ext)', &
          'Pa', pbarz)

     call add_hist_coord('rho_ext', nk1_ext, 'CRM density at midpoints (ext)',     &
          'kg/m3', rho(2:nk2_ext))

     call add_hist_coord('rhoz_ext', nk2_ext, 'CRM density at interfaces (ext)',   &
          'kg/m3', rhoz)

!-------------------------------------------------
!    Add CRM diagnostic variables on 'crm_grid'
!-------------------------------------------------
    !------------------------
    ! CRM background fields
    !------------------------
     call addfld ('TH3D_bg', (/'lev_c'/), 'I', 'K',                           &
          'Potential Temperature', gridname='crm_grid',                       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('TH3D_bg', 3, 'I')

     call addfld ('QV3D_bg', (/'lev_c'/), 'I', 'kg/kg',                       &
          'CRM water vapor mixing ratio', gridname='crm_grid',                &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QV3D_bg', 3, 'I')

     call addfld ('U3DX_bg', (/'lev_c'/), 'I', 'm/s',                         &
          'Zonal Wind', gridname='crm_grid',                                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U3DX_bg', 3, 'I')

     call addfld ('U3DY_bg', (/'lev_c'/), 'I', 'm/s',                         &
          'Meridional Wind', gridname='crm_grid',                             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U3DY_bg', 3, 'I')
     
     call addfld ('W3D_bg', (/'ilev_c'/), 'I', 'm/s',                         &
          'Vertical Wind', gridname='crm_grid',                               &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('W3D_bg', 3, 'I')     
     
     call addfld ('Z3DX_bg', (/'ilev_c'/), 'I', '1/s',                         &
          'Zonal Vorticity', gridname='crm_grid',                              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DX_bg', 3, 'I')   
     
     call addfld ('Z3DY_bg', (/'ilev_c'/), 'I', '1/s',                         &
          'Meridional Vorticity', gridname='crm_grid',                         &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DY_bg', 3, 'I')          
     
     call addfld ('Z3DZ_bg_top',horiz_only, 'I', '1/s',                       &
          'Vertical vorticity at top layer', gridname='crm_grid',             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DZ_bg_top', 3, 'I')          

    !------------------------
    ! CRM fields
    !------------------------     
     call addfld ('TH3D', (/'lev_c'/), 'I', 'K',                              &
          'Potential Temperature', gridname='crm_grid',                       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('TH3D', 3, 'I')

     call addfld ('QV3D', (/'lev_c'/), 'I', 'kg/kg',                          &
          'CRM water vapor mixing ratio', gridname='crm_grid',                &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QV3D', 3, 'I')

     call addfld ('U3DX', (/'lev_c'/), 'I', 'm/s',                            &
          'Zonal Wind', gridname='crm_grid',                                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U3DX', 3, 'I')

     call addfld ('U3DY', (/'lev_c'/), 'I', 'm/s',                            &
          'Meridional Wind', gridname='crm_grid',                             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U3DY', 3, 'I')
     
     call addfld ('RXTQ', (/'lev_c'/), 'I', 's',                              &
          'Relaxation time scale', gridname='crm_grid',                       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('RXTQ', 3, 'I')     
     
     call addfld ('W3D', (/'ilev_c'/), 'I', 'm/s',                            &
          'Vertical Wind', gridname='crm_grid',                               &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('W3D', 3, 'I')          

     call addfld ('Z3DX', (/'ilev_c'/), 'I', '1/s',                            &
          'Zonal Vorticity', gridname='crm_grid',                              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DX', 3, 'I')   
     
     call addfld ('Z3DY', (/'ilev_c'/), 'I', '1/s',                            &
          'Meridional Vorticity', gridname='crm_grid',                         &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DY', 3, 'I')          
     
     call addfld ('Z3DZ_top',horiz_only, 'I', '1/s',                           &
          'Vertical vorticity at top layer', gridname='crm_grid',              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Z3DZ_top', 3, 'I')     

     call addfld ('PSI_top',horiz_only, 'I', ' ',                          &
          'Stream function at top layer', gridname='crm_grid',             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('PSI_top', 3, 'I') 
                
!--------------------------------------------------------------
!    Add vGCM tendency diagnostic variables on 'vGCM_grid'
!--------------------------------------------------------------
     call addfld ('dT_vGCM', (/'lev'/), 'I', 'K/s',                           &
          'vGCM temperature tendency', gridname='vGCM_grid',                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dT_vGCM', 2, 'I')

     call addfld ('dQV_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM water vapor tendency', gridname='vGCM_grid',                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQV_vGCM', 2, 'I')
     
     call addfld ('dQC_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM cloud water mixing ratio tendency', gridname='vGCM_grid',     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQC_vGCM', 2, 'I')
     
     call addfld ('dQI_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM cloud ice mixing ratio tendency', gridname='vGCM_grid',       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQI_vGCM', 2, 'I')

     call addfld ('dQR_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM rain mixing ratio tendency', gridname='vGCM_grid',            &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQR_vGCM', 2, 'I')
     
     call addfld ('dQS_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM snow mixing ratio tendency', gridname='vGCM_grid',            &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQS_vGCM', 2, 'I')
     
     call addfld ('dQG_vGCM', (/'lev'/), 'I', 'kg/kg/s',                      &
          'vGCM Graupel mixing ratio tendency', gridname='vGCM_grid',         &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dQG_vGCM', 2, 'I')
     
     call addfld ('dU_vGCM', (/'lev'/), 'I', 'm/s/s',                         &
          'vGCM zonal momentum tendency', gridname='vGCM_grid',               &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dU_vGCM', 2, 'I')
     
     call addfld ('dV_vGCM', (/'lev'/), 'I', 'm/s/s',                         &
          'vGCM meridional momentum tendency', gridname='vGCM_grid',          &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('dV_vGCM', 2, 'I')                         
               
!--------------------------------------------------------------
!    Add vGCM diagnostic variables on 'vGCM_grid'
!--------------------------------------------------------------
     ! 2D Variables
     call addfld ('SPREC', horiz_only, 'I', 'kg/m^2/s',                      &
          'vGCM Surface Precipitation Rate', gridname='vGCM_grid',           &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('SPREC', 2, 'I')

     call addfld ('WTH', horiz_only, 'I', 'kg K/m^2/s',                      &
          'vGCM Surface Sensible Heat Flux', gridname='vGCM_grid',           &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('WTH', 2, 'I')

     call addfld ('WQV ', horiz_only, 'I', 'kg/m^2/s',                       &
          'vGCM Surface Moisture Flux', gridname='vGCM_grid',                &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('WQV', 2, 'I')

     call addfld ('UW', horiz_only, 'I', 'kg/m/s^2',                         &
          'vGCM Surface U-Momentum Flux', gridname='vGCM_grid',              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('UW', 2, 'I')

     call addfld ('WV', horiz_only, 'I', 'kg/m/s^2',                         &
          'vGCM Surface V-Momentum Flux', gridname='vGCM_grid',              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('WV', 2, 'I')                    
     
     ! 3D Variables
     call addfld ('T_vGCMi', (/'lev'/), 'I', 'K',                             &
          'vGCM Potential Temperature', gridname='vGCM_grid',                 &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('T_vGCMi', 2, 'I')

     call addfld ('QV_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM water vapor mixing ratio', gridname='vGCM_grid',              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QV_vGCMi', 2, 'I')

     call addfld ('QC_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM cloud liquid water mixing ratio', gridname='vGCM_grid',       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QC_vGCMi', 2, 'I')
     
     call addfld ('QI_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM cloud ice mixing ratio', gridname='vGCM_grid',                &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QI_vGCMi', 2, 'I')

     call addfld ('QR_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM rain mixing ratio', gridname='vGCM_grid',                     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QR_vGCMi', 2, 'I')
     
     call addfld ('QS_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM snow mixing ratio', gridname='vGCM_grid',                     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QS_vGCMi', 2, 'I')
     
     call addfld ('QG_vGCMi', (/'lev'/), 'I', 'kg/kg',                        &
          'vGCM graupel mixing ratio', gridname='vGCM_grid',                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QG_vGCMi', 2, 'I')                         

     call addfld ('U_vGCMi', (/'lev'/), 'I', 'm/s',                           &
          'vGCM Zonal Wind', gridname='vGCM_grid',                            &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U_vGCMi', 2, 'I')

     call addfld ('V_vGCMi', (/'lev'/), 'I', 'm/s',                           &
          'vGCM Meridional Wind', gridname='vGCM_grid',                       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('V_vGCMi', 2, 'I')

     call addfld ('OMEGA_vGCMi', (/'lev'/), 'I', 'Pa/s',                       &
          'vGCM vertical pressure velocity', gridname='vGCM_grid',             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('OMEGA_vGCMi', 2, 'I')

     call addfld ('Pint_vGCMi', (/'ilev'/), 'I', 'Pa',                         &
          'vGCM dry pressure at interfaces', gridname='vGCM_grid',             &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('Pint_vGCMi', 2, 'I')
     
     call addfld ('ZM_int_vGCMi', (/'ilev'/), 'I', 'm',                        &
          'vGCM geopotential height at interfaces', gridname='vGCM_grid',      &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('ZM_int_vGCMi', 2, 'I')
     
!-----------------------------------------------------------

     call addfld ('T_vGCM', (/'lev'/), 'I', 'K',                              &
          'vGCM Temperature', gridname='vGCM_grid',                           &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('T_vGCM', 2, 'I')

     call addfld ('QV_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM water vapor mixing ratio', gridname='vGCM_grid',              &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QV_vGCM', 2, 'I')

     call addfld ('QC_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM cloud liquid water mixing ratio', gridname='vGCM_grid',       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QC_vGCM', 2, 'I')

     call addfld ('QI_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM cloud ice mixing ratio', gridname='vGCM_grid',                &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QI_vGCM', 2, 'I')

     call addfld ('QR_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM rain mixing ratio', gridname='vGCM_grid',                     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QR_vGCM', 2, 'I')

     call addfld ('QS_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM snow mixing ratio', gridname='vGCM_grid',                     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QS_vGCM', 2, 'I')

     call addfld ('QG_vGCM', (/'lev'/), 'I', 'kg/kg',                         &
          'vGCM graupel mixing ratio', gridname='vGCM_grid',                  &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('QG_vGCM', 2, 'I')

     call addfld ('U_vGCM', (/'lev'/), 'I', 'm/s',                            &
          'vGCM Zonal Wind', gridname='vGCM_grid',                            &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('U_vGCM', 2, 'I')

     call addfld ('V_vGCM', (/'lev'/), 'I', 'm/s',                            &
          'vGCM Meridional Wind', gridname='vGCM_grid',                       &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('V_vGCM', 2, 'I')
     
     call addfld ('W_vGCM', (/'lev'/), 'I', 'm/s',                            &
          'vGCM vertical velocity', gridname='vGCM_grid',                     &
          flag_xyfill=.true., fill_value=R_UNDEF)
     call add_default('W_vGCM', 2, 'I')

#ifdef JUNG_TEST
!--------------------------------------------------------------
!    Add vGCM diagnostic variables on 'vGCM_ext_grid' (with halo)
!    Note: keep the trailing blank of 'lev'
!--------------------------------------------------------------
     call addfld ('T_halo_vGCM', (/'halo', 'lev '/), 'I', 'K',             &
          'vGCM Temperature', gridname='vGCM_ext_grid',                    &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('T_halo_vGCM', 2, 'I')
           
     call addfld ('QV_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM water vapor mixing ratio', gridname='vGCM_ext_grid',       &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QV_halo_vGCM', 2, 'I')      

     call addfld ('QC_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM cloud water mixing ratio', gridname='vGCM_ext_grid',       &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QC_halo_vGCM', 2, 'I')      

     call addfld ('QI_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM cloud ice mixing ratio', gridname='vGCM_ext_grid',         &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QI_halo_vGCM', 2, 'I')      
     
     call addfld ('QR_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM rain mixing ratio', gridname='vGCM_ext_grid',              &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QR_halo_vGCM', 2, 'I')      
     
     call addfld ('QS_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM snow mixing ratio', gridname='vGCM_ext_grid',              &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QS_halo_vGCM', 2, 'I')      
     
     call addfld ('QG_halo_vGCM', (/'halo', 'lev '/), 'I', 'kg/kg',        &
          'vGCM graupel mixing ratio', gridname='vGCM_ext_grid',           &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('QG_halo_vGCM', 2, 'I')                               

     call addfld ('U_halo_vGCM', (/'halo', 'lev '/), 'I', 'm/s',           &
          'vGCM Zonal Wind', gridname='vGCM_ext_grid',                     &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('U_halo_vGCM', 2, 'I')      

     call addfld ('V_halo_vGCM', (/'halo', 'lev '/), 'I', 'm/s',           &
          'vGCM Meridional Wind', gridname='vGCM_ext_grid',                &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('V_halo_vGCM', 2, 'I')      
     
     call addfld ('OMEGA_halo_vGCM', (/'halo', 'lev '/), 'I', 'Pa/s',      &
          'vGCM Meridional Wind', gridname='vGCM_ext_grid',                &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('OMEGA_halo_vGCM', 2, 'I')  
     
     call addfld ('Pint_halo_vGCM', (/'halo', 'ilev'/), 'I', 'Pa',         &
          'Dry pressure at interfaces', gridname='vGCM_ext_grid',          &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('Pint_halo_vGCM', 2, 'I')  

     call addfld ('ZM_int_halo_vGCM', (/'halo', 'ilev'/), 'I', 'm',        &
          'Geopotential height at interfaces', gridname='vGCM_ext_grid',   &
           flag_xyfill=.true.,fill_value=R_UNDEF)
     call add_default('ZM_int_halo_vGCM', 2, 'I')       
#endif                   

  end subroutine crm_init_diag

#endif

!==============================================================================
  subroutine crm_run (isRestart, iyear, rday_beg, rday_end, &
                      begchan, endchan, channel, vGCM)
!==============================================================================
!   main subprogram handling the VVM channels (called by the main program)
!------------------------------------------------------------------------
    USE q3d_bg_module,   only: cal_bg
    USE time_manager,    only: is_first_step

    USE parmsld,         only: nhalo_vGCM,channel_seg_l,nk1,nk2
    USE timeinfo,        only: rjday0
    
    USE q3d_init_module,   only: init_channel
    USE q3d_effect_module, only: cal_out_vgcm     

#ifndef single_run
    USE q3d_runtime, only: crm_chunk
    USE cam_history, only: outfld
#endif

    INTEGER, PARAMETER :: n_Cgrid = 5, i_Cgrid = 3 
    
    LOGICAL, INTENT(IN) :: isRestart        ! .true. iff restart run 
    INTEGER, INTENT(IN) :: iyear            ! Year of current GCM timestep
    REAL (kind=r8), INTENT(IN) :: rday_beg  ! Julian day at beginning of current GCM timestep
    REAL (kind=r8), INTENT(IN) :: rday_end  ! Julian day at end of current GCM timestep

    INTEGER, INTENT(IN) :: begchan  ! number of channels per node
    INTEGER, INTENT(IN) :: endchan  ! number of channels per node

    type(channel_t), pointer, intent(inout) :: channel(:)      ! CRM data for multiple channels
    type(channel_vGCM_t), pointer, intent(inout) :: vGCM(:)    ! vGCM data for multiple channels

    ! Local
    LOGICAL :: FIRST
    INTEGER,  parameter :: order(3) = (/ 2, 1, 3 /)

    INTEGER :: num_seg, k, h
    INTEGER :: VG_SIZE          ! size of vGCM segment plus halo along the channel
    INTEGER :: rshape(3),ishape(3)

    REAL(kind=r8), allocatable :: obufferl(:,:,:)
    REAL(kind=r8), allocatable :: obufferi(:,:,:)
    
    REAL(kind=r8), allocatable :: obuffer2d(:,:)

    FIRST = is_first_step()

    vg_size = nVGCM_seg + (2 * nhalo_vGCM)

    allocate(obufferl(channel_seg_l,nk1,n_Cgrid))
    allocate(obufferi(channel_seg_l,nk2,i_Cgrid))
    
    allocate(obuffer2d(channel_seg_l,2))

    rshape(1) = vg_size
    rshape(2) = (2 * nhalo_vGCM) + 1
    rshape(3) = nLevel

    ishape(1) = vg_size
    ishape(2) = (2 * nhalo_vGCM) + 1
    ishape(3) = nLevel + 1

!=======================
    IF (FIRST) THEN
!   Initialization
!=======================    
       rjday0 = rday_end   ! Fractional Julian day at the start of run: beginning of nstep = 1    
       
       do k = begchan, endchan
        
          channel(k)%isFirst_psi = .TRUE.
          channel(k)%isFirst_chi = .TRUE. 
         
          ! Calculate background fields for initialization
          CALL CAL_BG (vGCM(k)%vGCM_state,vGCM(k)%vGCM_map,channel(k))
         
          ! Initialize the channel data.
          CALL INIT_channel (isRestart,channel(k))
          
          CALL CAL_OUT_vGCM (channel(k))        
          CALL CRM_OUT (channel(k),vGCM(k)%vGCM_out)   
          
          if (channel(k)%num_chn == 1) write(6,*) 'JUNG: Finish FIRST Initialization' 
            
       end do
!============       
    ENDIF
!============   

!**************************************************
    do k = begchan, endchan
!**************************************************

!==================================
      IF (.NOT. FIRST) THEN  
!==================================  
!     For nstep = 0 (CAM-SE initialization step), the CRM components are not called.         
       
      ! Calculate background fields
      CALL CAL_BG (vGCM(k)%vGCM_state,vGCM(k)%vGCM_map,channel(k))
       
      ! Initialize the channel data.
      IF (isRestart) THEN
        
        channel(k)%isFirst_psi = .FALSE.
        channel(k)%isFirst_chi = .FALSE.  
             
        CALL INIT_channel (isRestart,channel(k))
        
        if (channel(k)%num_chn == 1) write(6,*) 'JUNG: Finish RESTART Initialization' 
      ENDIF   
                 
      ! Control the CRM integration
      CALL CRM_PREDICTION (iyear, rday_beg, channel(k))  
      
      ! Prepare the feedback (from channel data to vGCM_tend data)
      CALL CRM_OUT_TEND (channel(k), vGCM(k)%vGCM_tend)
      
      CALL CAL_OUT_vGCM (channel(k))                             
      CALL CRM_OUT (channel(k),vGCM(k)%vGCM_out)        

!============        
      ENDIF 
!============        

#ifndef single_run

      CALL CAL_CRM_WIND_RLL (channel(k))

      do num_seg = 1, 4

!-----------------------------------------------------------------------
! Output CRM diagnostic variables on 'crm_grid'
! (CRM fields need to be buffered due to limitations of outfld call.)
!-----------------------------------------------------------------------

         !------------------------
         ! CRM background fields
         !------------------------

         if (channel(k)%seg(num_seg)%nx_size > 1) then
         
          do h = 1, nk1
           obufferl(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%TH3D_bg(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%QV3D_bg(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%U3DX_ll_bg(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,4) = channel(k)%seg(num_seg)%U3DY_ll_bg(1:channel_seg_l,1,h+1)
          enddo
          do h = 1, nk2
           obufferi(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%W3D_bg(1:channel_seg_l,1,h)
           obufferi(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%Z3DX_bg(1:channel_seg_l,1,h)
           obufferi(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%Z3DY_bg(1:channel_seg_l,1,h)
          enddo 
          
           obuffer2d(1:channel_seg_l,1) = channel(k)%seg(num_seg)%Z3DZ_bg(1:channel_seg_l,1,1)
           
         else
         
          do h = 1, nk1
           obufferl(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%TH3D_bg(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%QV3D_bg(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%U3DX_ll_bg(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,4) = channel(k)%seg(num_seg)%U3DY_ll_bg(1,1:channel_seg_l,h+1)
          enddo
          do h = 1, nk2
           obufferi(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%W3D_bg(1,1:channel_seg_l,h)
           obufferi(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%Z3DX_bg(1,1:channel_seg_l,h)
           obufferi(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%Z3DY_bg(1,1:channel_seg_l,h)
          enddo 
          
          obuffer2d(1:channel_seg_l,1) = channel(k)%seg(num_seg)%Z3DZ_bg(1,1:channel_seg_l,1)
                   
         endif 

         ! values at mid-layers
         call outfld('TH3D_bg',obufferl(:,:,1),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('QV3D_bg',obufferl(:,:,2),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('U3DX_bg',obufferl(:,:,3),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('U3DY_bg',obufferl(:,:,4),channel_seg_l,crm_chunk(k,num_seg))
         
         ! values at interfaces
         call outfld('W3D_bg',obufferi(:,:,1),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('Z3DX_bg',obufferi(:,:,2),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('Z3DY_bg',obufferi(:,:,3),channel_seg_l,crm_chunk(k,num_seg))
         
         ! values at top (active part)
         call outfld('Z3DZ_bg_top',obuffer2d(:,1),channel_seg_l,crm_chunk(k,num_seg))
         
         !------------------------
         ! CRM fields
         !------------------------

         if (channel(k)%seg(num_seg)%nx_size > 1) then
         
          do h = 1, nk1
           obufferl(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%TH3D(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%QV3D(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%U3DX_ll(1:channel_seg_l,1,h+1)
           obufferl(1:channel_seg_l,h,4) = channel(k)%seg(num_seg)%U3DY_ll(1:channel_seg_l,1,h+1)
           
           obufferl(1:channel_seg_l,h,5) = channel(k)%seg(num_seg)%TAU_RX_TQ(1:channel_seg_l,1,h+1)
          enddo
          do h = 1, nk2
           obufferi(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%W3D(1:channel_seg_l,1,h)
           obufferi(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%Z3DX(1:channel_seg_l,1,h)
           obufferi(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%Z3DY(1:channel_seg_l,1,h)
          enddo 
          
          obuffer2d(1:channel_seg_l,1) = channel(k)%seg(num_seg)%Z3DZ(1:channel_seg_l,1,nk2) 
          obuffer2d(1:channel_seg_l,2) = channel(k)%seg(num_seg)%PSI(1:channel_seg_l,1)  
                 
         else
         
          do h = 1, nk1
           obufferl(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%TH3D(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%QV3D(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%U3DX_ll(1,1:channel_seg_l,h+1)
           obufferl(1:channel_seg_l,h,4) = channel(k)%seg(num_seg)%U3DY_ll(1,1:channel_seg_l,h+1)
           
           obufferl(1:channel_seg_l,h,5) = channel(k)%seg(num_seg)%TAU_RX_TQ(1,1:channel_seg_l,h+1)
          enddo
          do h = 1, nk2
           obufferi(1:channel_seg_l,h,1) = channel(k)%seg(num_seg)%W3D(1,1:channel_seg_l,h)
           obufferi(1:channel_seg_l,h,2) = channel(k)%seg(num_seg)%Z3DX(1,1:channel_seg_l,h)
           obufferi(1:channel_seg_l,h,3) = channel(k)%seg(num_seg)%Z3DY(1,1:channel_seg_l,h)
          enddo 
          
          obuffer2d(1:channel_seg_l,1) = channel(k)%seg(num_seg)%Z3DZ(1,1:channel_seg_l,nk2)
          obuffer2d(1:channel_seg_l,2) = channel(k)%seg(num_seg)%PSI(1,1:channel_seg_l) 
                   
         endif 

         ! values at mid-layers
         call outfld('TH3D',obufferl(:,:,1),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('QV3D',obufferl(:,:,2),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('U3DX',obufferl(:,:,3),channel_seg_l,crm_chunk(k,num_seg))
         call outfld('U3DY',obufferl(:,:,4),channel_seg_l,crm_chunk(k,num_seg))
         
         call outfld('RXTQ',obufferl(:,:,5),channel_seg_l,crm_chunk(k,num_seg))
         
         ! values at interfaces
         call outfld('W3D' ,obufferi(:,:,1),channel_seg_l,crm_chunk(k,num_seg)) 
         call outfld('Z3DX',obufferi(:,:,2),channel_seg_l,crm_chunk(k,num_seg)) 
         call outfld('Z3DY',obufferi(:,:,3),channel_seg_l,crm_chunk(k,num_seg))     
         
         ! values at top (active)
         call outfld('Z3DZ_top' ,obuffer2d(:,1),channel_seg_l,crm_chunk(k,num_seg))  
         call outfld('PSI_top' ,obuffer2d(:,2),channel_seg_l,crm_chunk(k,num_seg))         
 
!-----------------------------------------------------------------------
! Output the vGCM tendency fields on 'vGCM_grid'
!-----------------------------------------------------------------------
         call outfld('dT_vGCM',vGCM(k)%vGCM_tend(num_seg)%dT(:,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQV_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQV(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQC_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQC(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQI_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQI(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQR_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQR(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQS_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQS(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dQG_vGCM',vGCM(k)%vGCM_tend(num_seg)%dQG(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         
         call outfld('dU_vGCM',vGCM(k)%vGCM_tend(num_seg)%dU(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('dV_vGCM',vGCM(k)%vGCM_tend(num_seg)%dV(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))                                                                                    

!-----------------------------------------------------------------------
! Output vGCM diagnostic variables on 'vGCM_grid'
!-----------------------------------------------------------------------

         call outfld('SPREC',vGCM(k)%vGCM_out(num_seg)%SPREC(:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('WTH',vGCM(k)%vGCM_out(num_seg)%WTH(:),       &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('WQV',vGCM(k)%vGCM_out(num_seg)%WQV(:),       &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('UW',vGCM(k)%vGCM_out(num_seg)%UW(:),         &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('WV',vGCM(k)%vGCM_out(num_seg)%WV(:),         &
                     nVGCM_seg,crm_chunk(k, num_seg))                                                                                    
                              
         call outfld('T_vGCM',vGCM(k)%vGCM_out(num_seg)%T(:,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QV_vGCM',vGCM(k)%vGCM_out(num_seg)%QV(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QC_vGCM',vGCM(k)%vGCM_out(num_seg)%QC(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QI_vGCM',vGCM(k)%vGCM_out(num_seg)%QI(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QR_vGCM',vGCM(k)%vGCM_out(num_seg)%QR(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QS_vGCM',vGCM(k)%vGCM_out(num_seg)%QS(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QG_vGCM',vGCM(k)%vGCM_out(num_seg)%QG(:,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))                                                            

         call outfld('U_vGCM',vGCM(k)%vGCM_out(num_seg)%U(:,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('V_vGCM',vGCM(k)%vGCM_out(num_seg)%V(:,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('W_vGCM',vGCM(k)%vGCM_out(num_seg)%W(:,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))

         call outfld('T_vGCMi',vGCM(k)%vGCM_state(num_seg)%T(0,1:nVGCM_seg,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QV_vGCMi',vGCM(k)%vGCM_state(num_seg)%QV(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QC_vGCMi',vGCM(k)%vGCM_state(num_seg)%QC(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QI_vGCMi',vGCM(k)%vGCM_state(num_seg)%QI(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QR_vGCMi',vGCM(k)%vGCM_state(num_seg)%QR(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QS_vGCMi',vGCM(k)%vGCM_state(num_seg)%QS(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('QG_vGCMi',vGCM(k)%vGCM_state(num_seg)%QG(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))                                                                                                         
         call outfld('U_vGCMi',vGCM(k)%vGCM_state(num_seg)%U(0,1:nVGCM_seg,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('V_vGCMi',vGCM(k)%vGCM_state(num_seg)%V(0,1:nVGCM_seg,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))
         call outfld('OMEGA_vGCMi',vGCM(k)%vGCM_state(num_seg)%omega(0,1:nVGCM_seg,:), &
                     nVGCM_seg,crm_chunk(k, num_seg))  
         call outfld('ZM_int_vGCMi',vGCM(k)%vGCM_state(num_seg)%zm_int(0,1:nVGCM_seg,:),   &
                     nVGCM_seg,crm_chunk(k, num_seg))     
         call outfld('Pint_vGCMi',vGCM(k)%vGCM_state(num_seg)%pint(0,1:nVGCM_seg,:),      &
                     nVGCM_seg,crm_chunk(k, num_seg))   

#ifdef JUNG_TEST
!-----------------------------------------------------------------------
! Output vGCM diagnostic variables on 'vGCM_ext_grid' (with halo)
! (Note this is inefficient due to differences in array ordering.)
!-----------------------------------------------------------------------

         call outfld('T_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%T,   &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QV_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QV, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QC_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QC, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QI_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QI, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QR_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QR, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QS_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QS, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('QG_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%QG, &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))                                                                                    

         call outfld('U_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%U,   &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('V_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%V,   &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('OMEGA_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%OMEGA,   &
                     rshape, order=order),vg_size,crm_chunk(k, num_seg))
                     
         call outfld('Pint_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%pint,     &
                     ishape, order=order),vg_size,crm_chunk(k, num_seg))

         call outfld('ZM_int_halo_vGCM',reshape(vGCM(k)%vGCM_state(num_seg)%zm_int, &
                     ishape, order=order),vg_size,crm_chunk(k, num_seg)) 
#endif                                                              

      enddo  ! num_seg-loop

#endif

!**************************************************
    enddo ! k-loop (channels)
!**************************************************    

    CONTAINS
     
    SUBROUTINE CAL_CRM_WIND_RLL (channel)
!   Calculate the rll components of wind
    USE constld, only: val_missing
    USE utils,   only: con2rll
    
    type(channel_t), intent(inout) :: channel      ! channel data
   
    ! Local
    REAL(KIND=r8) :: TEMP_U, TEMP_V, VMAP(2)
   
    integer num_seg,mi1,mj1
    integer I,J,K
   
!**************************************************************
    DO num_seg = 1, 4
!**************************************************************

    mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
    mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

    channel%seg(num_seg)%U3DX_ll(:,:,1) = val_missing
    channel%seg(num_seg)%U3DY_ll(:,:,1) = val_missing
      
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         
          TEMP_U = 0.5_r8*(channel%seg(num_seg)%U3DX(i,j,k)   &
                         + channel%seg(num_seg)%U3DX(i-1,j,k))
          TEMP_V = 0.5_r8*(channel%seg(num_seg)%U3DY(i,j,k)   &
                         + channel%seg(num_seg)%U3DY(i,j-1,k))
           
          VMAP = CON2RLL(channel%seg(num_seg)%AM_T(1,i,j),    &
                         channel%seg(num_seg)%AM_T(2,i,j),    &
                         channel%seg(num_seg)%AM_T(3,i,j),    &
                         channel%seg(num_seg)%AM_T(4,i,j),    &
                         TEMP_U,TEMP_V)
          
          channel%seg(num_seg)%U3DX_ll(I,J,K) = VMAP(1) 
          channel%seg(num_seg)%U3DY_ll(I,J,K) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO              

!**************************************************************
    ENDDO   ! num_seg
!*************************************************************
    END SUBROUTINE cal_crm_wind_rll  

  end subroutine crm_run

!==============================================================================
  subroutine crm_final (itt, begchan, endchan, channel, vGCM)
!==============================================================================
    integer, intent(in) :: itt     ! vGCM time step count
    INTEGER, INTENT(IN) :: begchan ! number of channels per node
    INTEGER, INTENT(IN) :: endchan ! number of channels per node

    type(channel_t), pointer, intent(inout) :: channel(:)      ! CRM data for multiple channels
    type(channel_vGCM_t), pointer, intent(inout) :: vGCM(:)    ! vGCM data for multiple channels

    ! This routine is for any needed cleanup or output (e.g., run statistics)
    ! It is called after all run timesteps are completed.

  end subroutine crm_final

!===================================================================================
    SUBROUTINE CRM_PREDICTION (iyear, rday_beg, channel)
!===================================================================================
    USE timeinfo,  only: iyr,rjday
    USE constld,   only: itt_max,a,b,dt,physics,sec_cday,rxtau_q,rxtau_d

    USE q3d_effect_module, only: ini_tendency,ini_crm_eft,cal_eddy_trans, &
                                 cal_tendency,ini_sfc_eft,cal_sfc_eft,cal_feedback, &
                                 eddy_prepare,cal_out_vgcm
    USE q3d_rtime_module,  only: cal_rtime
    USE ab_3d_module,      only: ab_3d

    INTEGER, INTENT(IN) :: iyear               ! Year of current GCM timestep
    REAL (kind=r8), INTENT(IN) :: rday_beg     ! Julian day at beginning of current GCM timestep

    type(channel_t), intent(inout) :: channel  ! CRM data for one channel

    ! Local
    LOGICAL :: DEBUG =.FALSE.
    
    LOGICAL :: MAKE_AVG
    INTEGER :: ITT_CRM, N1, N2
           
    MAKE_AVG = .FALSE.

    iyr = iyear    ! TIMEINFO for radiation calculation
     
!   Initialize the variables related to the calculation of mean eddy effects (for time average)
    CALL INI_CRM_EFT (channel)

    IF (physics) THEN
!    Initialize the variables related to the calculation of mean sfc fluxes (for time average)
     CALL INI_SFC_EFT (channel)
    ENDIF
     
!     Calculate the relaxation timescale (CRM toward GCM)
      IF (RXTAU_q .EQ. 0.0_r8 .AND. RXTAU_d .EQ. 0.0_r8) CALL CAL_RTIME (channel) 

!-----------------------------------------------------------------------
      DO ITT_CRM = 1, itt_max         ! CRM subcycle time marching
!-----------------------------------------------------------------------

      IF(ITT_CRM.eq.1) THEN
        a = dt                         ! forward scheme to initialize the CRM
        b = 0.0_r8                     ! forward scheme to initialize the CRM
      ELSE
        a =   3.0_r8 / 2.0_r8 * dt     ! AB second-order scheme
        b = - 1.0_r8 / 2.0_r8 * dt     ! AB second-order scheme
      ENDIF

      N1 = MOD ( ITT_CRM    , 2 ) + 1
      N2 = MOD ( ITT_CRM - 1, 2 ) + 1

!     Define the fractional Julian day at CRM time-step
      rjday = rday_beg + (ITT_CRM-1) * DT / sec_cday     ! TIMEINFO for radiation calculation
      
      if (channel%num_chn ==1) write(6,*) 'JUNG: ITT_CRM=',ITT_CRM,iyr,rjday

      CALL INI_TENDENCY (channel)
      
!     Individual CRM integration
      CALL AB_3D (N1, N2, ITT_CRM, channel)    ! remove ch_num

      if (DEBUG .AND. channel%num_chn ==1) write(6,*) 'JUNG: a4 AB_3D'
      
      IF (ITT_CRM.eq.itt_max) MAKE_AVG=.TRUE.

!     Determine eddy components
      CALL EDDY_PREPARE (channel)

!     Calculate the mean eddy effcts
      CALL CAL_EDDY_TRANS (itt_max,MAKE_AVG,channel)   ! remove ch_num
    
      CALL CAL_TENDENCY (itt_max,MAKE_AVG,channel)

      IF (physics) THEN
!      Calculate the mean sfc fluxes
       CALL CAL_SFC_EFT (itt_max,MAKE_AVG,channel)
      ENDIF

      IF (MAKE_AVG) CALL CAL_FEEDBACK(channel)
      if (DEBUG .AND. channel%num_chn ==1) write(6,*) 'JUNG: a4 CAL_FEEDBACK MAKE_AVG=',MAKE_AVG  

!-----------------------------------
      ENDDO   ! CRM time marching
!-----------------------------------

    END SUBROUTINE crm_prediction

!===================================================================================
    SUBROUTINE CRM_OUT_TEND (channel,vGCM_tend)
!   Arrange the CRM tendencies following the vGCM_tend data shape.
!   Interpolate theta at CRM layers to nLevel and convert the theta to Temp. at nLevel
!   INPUT: _E1 (total CRM effects)
!===================================================================================
    USE parmsld, only: ntracer,nk2,nk2_ext
    USE constld, only: physics

    type(channel_t),   intent(in)    :: channel       ! CRM data for one channel
    type(vGCM_tend_t), intent(inout) :: vGCM_tend(4)  ! vGCM tendency

    ! Local
    INTEGER, PARAMETER :: id_tend = 0                 ! handling tendency data

    INTEGER :: num_seg,mi1,mj1,nt,J,K

!*******************************
    DO num_seg = 1, 4
!*******************************
    mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
    mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

!#ifdef JUNG_TEST
!--------------------------------------------------
   IF (mi1.GT.mj1) THEN
!-------------------------------------------------- x-channel segment

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2_ext,id_tend,  &
                  channel%seg(num_seg)%TH_E1,                  &  ! theta tendency
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dT)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QV_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQV)

     if (physics) then
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QC_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQC)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QI_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQI)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QR_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQR)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QS_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQS)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QG_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQG)
     endif

     DO NT = 1,ntracer
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QT_E1(:,:,nt),          &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQT(:,:,nt))
     ENDDO

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%U_E1,                   &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dU)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%V_E1,                   &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dV)

!--------------------------------------------------
   ELSE
!-------------------------------------------------- y-channel segment

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2_ext,id_tend,  &
                  channel%seg(num_seg)%TH_E1,                  &  ! theta tendency
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dT)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QV_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQV)

     if (physics) then
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QC_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQC)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QI_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQI)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QR_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQR)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QS_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQS)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QG_E1,                  &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQG)
     endif

     DO NT = 1,ntracer
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%QT_E1(:,:,nt),          &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dQT(:,:,nt))
     ENDDO

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%U_E1,                   &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dU)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_tend,      &
                  channel%seg(num_seg)%V_E1,                   &
                  channel%seg(num_seg)%ZM_G,                   &
                  vGCM_tend(num_seg)%dV)

!--------------------------------------------------
   ENDIF
!--------------------------------------------------

!     Convert theta-tendency to temperature-tendency
      DO K = 1, nLevel
       DO J = 1,nVGCM_seg
        vGCM_tend(num_seg)%dT(J,K) = vGCM_tend(num_seg)%dT(J,K)*channel%seg(num_seg)%PI_G(K,J)
       ENDDO
      ENDDO

!#else
!------------------------------ JUNG_DEBUG (remove)
#ifdef JUNG_TEST
!------------------------------
     vGCM_tend(num_seg)%dT(:,:)  = 0.0_r8
     
     vGCM_tend(num_seg)%dQV(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dQC(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dQI(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dQR(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dQS(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dQG(:,:) = 0.0_r8
     
     vGCM_tend(num_seg)%dQT(:,:,:) = 0.0_r8
     
     vGCM_tend(num_seg)%dU(:,:) = 0.0_r8
     vGCM_tend(num_seg)%dV(:,:) = 0.0_r8
       
!------------------------------ JUNG_DEBUG (remove)
#endif
!------------------------------

!*******************************
   ENDDO  ! num_seg
!*******************************

    END SUBROUTINE crm_out_tend

!===================================================================================
    SUBROUTINE CRM_OUT (channel,vGCM_out,time)
!   Arrange the CRM diagnostic data following the vGCM_out data shape
!   Interpolate theta at CRM layers to nLevel and convert the theta to Temp. at nLevel
!   INPUT: _E0 (net-size and time average over the CRM cycle)
!   Currently, it has only 2D fields.
!===================================================================================
    USE parmsld, only: ntracer,nk2,nk2_ext,nLevel,nVGCM_seg
    USE constld, only: zt,zz,physics,val_missing

    type(channel_t),     intent(in) :: channel      ! CRM data for one channel
    type(vGCM_out_t), intent(inout) :: vGCM_out(4)  ! vGCM out data

    REAL(kind=r8), optional, intent(in) :: time

    ! Local
    INTEGER, PARAMETER :: id_out = 1                ! handling state (field) data
    
    INTEGER :: Kval(nLevel,nVGCM_seg)

    INTEGER :: num_seg,mi1,mj1,nt,J,K

!*******************************
    DO num_seg = 1, 4
!*******************************
    mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
    mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment
    
    CALL vIndex_z (nVGCM_seg,nLevel,nk2_ext,channel%seg(num_seg)%ZM_G,zt,zz,Kval)                      

!--------------------------------------------------
   IF (mi1.GT.mj1) THEN
!-------------------------------------------------- x-channel segment

!---------------------
! 2D Variables
!---------------------
     if (physics) then
       CALL GET_X (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%SPREC_E0,  &
                   vGCM_out(num_seg)%SPREC)

       CALL GET_X (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WTH_E0,    &
                   vGCM_out(num_seg)%WTH)

       CALL GET_X (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WQV_E0,    &
                   vGCM_out(num_seg)%WQV)

       CALL GET_X (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%UW_E0,     &
                   vGCM_out(num_seg)%UW)

       CALL GET_X (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WV_E0,     &
                   vGCM_out(num_seg)%WV)
     endif

!---------------------
! 3D Variables
!---------------------
     
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%TH3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%T)

     ! Use a cubic-spline interpolation & constrain the positive values.  
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QV3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QV,positive=.true.)

     ! Use a linear interpolation for QC,QI,QR,QS,QG,QT to avoid a generation of negative value.
     if (physics) then
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QC3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QC,kpoint=kval)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QI3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QI,kpoint=kval)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QR3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QR,kpoint=kval)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QS3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QS,kpoint=kval)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QG3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QG,kpoint=kval)
     endif

     DO NT = 1,ntracer
     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QT3D_sa(:,:,nt), &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QT(:,:,nt),kpoint=kval)
     ENDDO

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%U3DX_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%U)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%U3DY_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%V)

     CALL GET_XL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%W3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%W,cv_type=1)

!--------------------------------------------------
   ELSE
!-------------------------------------------------- y-channel segment

!---------------------
! 2D Variables
!---------------------
     if (physics) then
       CALL GET_Y (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%SPREC_E0,  &
                   vGCM_out(num_seg)%SPREC)

       CALL GET_Y (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WTH_E0,    &
                   vGCM_out(num_seg)%WTH)

       CALL GET_Y (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WQV_E0,    &
                   vGCM_out(num_seg)%WQV)

       CALL GET_Y (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%UW_E0,     &
                   vGCM_out(num_seg)%UW)

       CALL GET_Y (channel%seg(num_seg)%NFACE,channel%seg(num_seg)%WV_E0,     &
                   vGCM_out(num_seg)%WV)
     endif

!---------------------
! 3D Variables
!---------------------

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%TH3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%T)

     ! Use a cubic-spline interpolation & constrain the positive values.  
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QV3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QV,positive=.true.)

     ! Use a linear interpolation for QC,QI,QR,QS,QG,QT to avoid a generation of negative value. 
     if (physics) then
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QC3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QC,kpoint=kval)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QI3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QI,kpoint=kval)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QR3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QR,kpoint=kval)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QS3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QS,kpoint=kval)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QG3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QG,kpoint=kval)
     endif

     DO NT = 1,ntracer
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%QT3D_sa(:,:,nt), &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%QT(:,:,nt),kpoint=kval)
     ENDDO

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%U3DX_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%U)

     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%U3DY_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%V)
                  
     CALL GET_YL (channel%seg(num_seg)%NFACE,nk2,id_out,channel%seg(num_seg)%W3D_sa, &
                  channel%seg(num_seg)%ZM_G,vGCM_out(num_seg)%W,cv_type=1)                  
                  
!--------------------------------------------------
   ENDIF
!--------------------------------------------------
!     Convert theta to temperature
      DO K = 1, nLevel
       DO J = 1,nVGCM_seg
        if (vGCM_out(num_seg)%T(J,K) .ne. val_missing) then
        vGCM_out(num_seg)%T(J,K) = vGCM_out(num_seg)%T(J,K)*channel%seg(num_seg)%PI_G(K,J)
        endif
       ENDDO
      ENDDO

!*******************************
   ENDDO  ! num_seg
!*******************************

    END SUBROUTINE crm_out

    SUBROUTINE GET_X (NFACE,E1,EOUT)
!   arrange the data for the vGCM grid structure

    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs

    REAL(kind=r8),DIMENSION(nVGCM_seg),INTENT(IN)  :: E1
    REAL(kind=r8),DIMENSION(nVGCM_seg),INTENT(OUT) :: EOUT

    ! LOCAL
    INTEGER I,NPO

    IF (NFACE .NE. 6) THEN
      DO I = 1,nVGCM_seg
        EOUT(I) = E1(I)
      ENDDO
    ELSE
      DO I = 1,nVGCM_seg
        npo = nVGCM_seg - I + 1
        EOUT(I) = E1(npo)
      ENDDO
    ENDIF

    END SUBROUTINE get_x

    SUBROUTINE GET_XL (NFACE,kdim,id_type,E1,pibar_G,EOUT,cv_type,kpoint,positive)
!   arrange the data for the vGCM grid structure with vertical layer dimension
!   Horizontal : -nhalo_vGCM:nhalo_vGCM, 1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM
!   Vertical   : nLevel (nLevel+1)

    USE constld, only: pibar,pibarz
    USE constld, only: zt,zz        

    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data
    INTEGER, INTENT(IN) :: id_type   ! data type (0: tendency data, 1: state data)
    
    REAL(kind=r8),DIMENSION(kdim,nVGCM_seg),  INTENT(IN) :: E1
    REAL(kind=r8),DIMENSION(nLevel,nVGCM_seg),INTENT(IN) :: pibar_G

    REAL(kind=r8),DIMENSION(nVGCM_seg,nLevel),INTENT(OUT) :: EOUT
    
    INTEGER, optional, INTENT(IN) :: cv_type   ! convert type (1: vvm interface to vGCM layers)
    INTEGER, optional, INTENT(IN) :: kpoint(nLevel,nVGCM_seg)
    LOGICAL, optional, INTENT(IN) :: positive    

    ! LOCAL
    REAL(kind=r8),DIMENSION(nLevel,nVGCM_seg) :: E1_NEW
    INTEGER I,K,NPO

    if (present(cv_type)) then
      DO I = 1,nVGCM_seg
       CALL VVM_i_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I))
      ENDDO    
    else
      if (present(kpoint)) then
      
      DO I = 1,nVGCM_seg
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I),kval=kpoint(:,I))
      ENDDO      
      
      else
      
      if (present(positive)) then
      DO I = 1,nVGCM_seg
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I), positive=positive)
      ENDDO      
      else
      DO I = 1,nVGCM_seg
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I))
      ENDDO
      endif
      
      endif
    
    endif
    
    IF (NFACE .NE. 6) THEN
      DO I = 1,nVGCM_seg
       DO K = 1, nLevel
        EOUT(I,K) = E1_new(K,I)
       ENDDO
      ENDDO
    ELSE
      DO I = 1,nVGCM_seg
       npo = nVGCM_seg - I + 1
       DO K = 1, nLevel
        EOUT(I,K) = E1_new(K,npo)
       ENDDO
      ENDDO
    ENDIF

    END SUBROUTINE get_xl

    SUBROUTINE GET_Y (NFACE,E1,EOUT)
!   arrange the data for the vGCM grid structure

    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs

    REAL(kind=r8),DIMENSION(nVGCM_seg),INTENT(IN)  :: E1
    REAL(kind=r8),DIMENSION(nVGCM_seg),INTENT(OUT) :: EOUT

    ! LOCAL
    INTEGER I,NPO

    IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO I = 1,nVGCM_seg
        npo = nVGCM_seg - I + 1
        EOUT(I) = E1(npo)
      ENDDO
    ELSE
      DO I = 1,nVGCM_seg
        EOUT(I) = E1(I)
      ENDDO
    ENDIF

    END SUBROUTINE get_y

    SUBROUTINE GET_YL (NFACE,kdim,id_type,E1,pibar_G,EOUT,cv_type,kpoint,positive)
!   arrange the data for the vGCM grid structure with vertical layer dimension
!   Horizontal : -nhalo_vGCM:nhalo_vGCM, 1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM
!   Vertical   : nLevel (nLevel+1)

    USE constld, only: pibar,pibarz
    USE constld, only: zt,zz    
    
    INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
    INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data
    INTEGER, INTENT(IN) :: id_type   ! data type (0: tendency data, 1: state data)
    
    REAL(kind=r8),DIMENSION(kdim,nVGCM_seg),  INTENT(IN) :: E1
    REAL(kind=r8),DIMENSION(nLevel,nVGCM_seg),INTENT(IN) :: pibar_G

    REAL(kind=r8),DIMENSION(nVGCM_seg,nLevel),INTENT(OUT) :: EOUT

    INTEGER, optional, INTENT(IN) :: cv_type   ! convert type (1: vvm interface to vGCM layers)
    INTEGER, optional, INTENT(IN) :: kpoint(nLevel,nVGCM_seg)
    LOGICAL, optional, INTENT(IN) :: positive 
         
    ! LOCAL
    REAL(kind=r8),DIMENSION(nLevel,nVGCM_seg) :: E1_NEW
    INTEGER I,K,NPO
        
    if (present(cv_type)) then
      DO I = 1,nVGCM_seg
       CALL VVM_i_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I))
      ENDDO    
    else
      if (present(kpoint)) then
      
      DO I = 1,nVGCM_seg
      
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I),kval=kpoint(:,I))                      
      ENDDO      
      
      else
      
      if (present(positive)) then
      DO I = 1,nVGCM_seg
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I),positive=positive)
      ENDDO      
      else
      DO I = 1,nVGCM_seg
       CALL VVM_l_to_vGCM_l (kdim,nLevel,id_type,E1(:,I),zt(1:kdim),zz(1:kdim),  &
                             pibar_G(:,I),E1_NEW(:,I))
      ENDDO
      endif
      
      endif
      
    endif

    IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO I = 1,nVGCM_seg
       npo = nVGCM_seg - I + 1
       DO K = 1, nLevel
        EOUT(I,K) = E1_new(K,npo)
       ENDDO
      ENDDO
    ELSE
      DO I = 1,nVGCM_seg
       DO K = 1, nLevel
        EOUT(I,K) = E1_new(K,I)
       ENDDO
      ENDDO
    ENDIF

    END SUBROUTINE get_yl

!===========================================================================================
      SUBROUTINE VVM_layers_to_vGCM_layers (kdim,kmG,id_type, &
                                            val_C,exnerl_C,exner_C,exnerl_G,val_G)
!===========================================================================================
!     Vertical interpolation of VVM tendency to that on GCM grids
!     Val_G, exnerl, val: downward indexing (top to bottom)
!     Val_C: upward indexing (bottom to top)

      USE constld, only: val_missing
      USE utils,   only: spline,splint

      INTEGER, INTENT(IN) :: kdim      ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG       ! vertical dimension of VAL_G (nLevel)
      INTEGER, INTENT(IN) :: id_type   ! data type (0: tendency data, 1: state data)

      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: val_C      ! target variable of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: exnerl_C   ! exner ftn of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: exner_C    ! exner ftn of VVM (at interface)

      REAL(kind=r8), DIMENSION(kmG), INTENT(IN)  :: exnerl_G   ! exner ftn of GCM (at layer)
      REAL(kind=r8), DIMENSION(kmG), INTENT(OUT) :: val_G      ! target variable of GCM (at mid-layer)
       
      ! Local
      REAL(kind=r8), DIMENSION(kdim+1) :: exnerl,val,y2
      REAL(kind=r8) :: yp1,ypkm
      INTEGER k,kbeg

      ! Reverse data
      do k=2,kdim
       exnerl(k) = exnerl_C(kdim-k+2)
       val(k) = val_C(kdim-k+2)
      enddo

      exnerl(1) = exner_C(kdim)
      exnerl(kdim+1) = exner_C(1)

      val(1) = val(2)
      val(kdim+1) = val(kdim)
       
!---------------
      k=1
      do while (exnerl_G(k) .lt. exnerl(1))
      k = k + 1
      enddo
      kbeg = k
!---------------
      
      yp1  = 0.5_r8*((val(2)-val(1))/(exnerl(2)-exnerl(1))   &
                    +(val(3)-val(2))/(exnerl(3)-exnerl(2)))

      ypkm = 0.5_r8*((val(kdim+1)-val(kdim))/(exnerl(kdim+1)-exnerl(kdim))  &
                    +(val(kdim)-val(kdim-1))/(exnerl(kdim)-exnerl(kdim-1)))

      call spline(exnerl,val,kdim+1,yp1,ypkm,y2)

      do k=kbeg,kmG
       call splint(exnerl,val,y2,kdim+1,exnerl_G(k),val_G(k))
      enddo

      if (kbeg .gt. 1) then
!     For the GCM grids above the (active) CRM vertical domain.

      IF (id_type .EQ. 0) THEN
        ! tendency data
        do k=1,kbeg-1
         val_G(k) = 0.0_r8
        enddo
      ELSE
        ! field data
        do k=1,kbeg-1
         val_G(k) = val_missing
        enddo
      ENDIF

      endif

      END SUBROUTINE VVM_layers_to_vGCM_layers

!===========================================================================================
      Subroutine vIndex_z (JmG,kmG,kdim,hgtl_G,hgtl_C,hgt_C,Kval)
!===========================================================================================      
!     Find the vertical grid index for a linear interpolation.
 
!     hgtl_G       : downward indexing (top to bottom) 
!     hgtl_C, hgt_C: upward indexing (bottom to top)

      USE utils, only: indexr
      
      INTEGER, INTENT(IN) :: JmG    ! horizontal dimension of VAL_G (nVGCM_seg)  
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2)
      
      REAL(kind=r8), DIMENSION(kmG,JmG), INTENT(IN) :: hgtl_G   ! height of GCM (at mid-layer) 
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgtl_C   ! height of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgt_C    ! height of VVM (at mid-layer)
      
      INTEGER, DIMENSION(kmG,JmG), INTENT(OUT)  :: kval     

      ! LOCAL
      LOGICAL LF
      REAL(kind=r8), DIMENSION(kdim+1) :: hgtl
      INTEGER k,j
      
      do k=2,kdim
       hgtl(k) = hgtl_C(k)
      enddo
      
      hgtl(1)      = hgt_C(1)
      hgtl(kdim+1) = hgt_C(kdim) 
      
      do j=1,JmG
      do k=1,kmG 
       kval(k,j) = INDEXR (hgtl_G(k,j),kdim+1,hgtl,LF)  
      enddo     
      enddo 
            
      END SUBROUTINE vindex_z   

!===========================================================================================
      SUBROUTINE VVM_l_to_vGCM_l (kdim,kmG,id_type, &
                                  val_C,hgtl_C,hgt_C,hgtl_G,val_G,kval,positive)
!===========================================================================================
!     Vertical interpolation of VVM data at mid-layer to that on GCM mid-layers
!     Val_G, hgtl_G: downward indexing (top to bottom)
!     Val_C, hgt_C : upward indexing (bottom to top)

      USE constld, only: val_missing
      USE utils,   only: spline,splint,fintrp

      INTEGER, INTENT(IN) :: kdim      ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG       ! vertical dimension of VAL_G (nLevel)
      INTEGER, INTENT(IN) :: id_type   ! data type (0: tendency data, 1: state data)

      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: val_C    ! target variable of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: hgtl_C   ! height of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: hgt_C    ! height of VVM (at interface)

      REAL(kind=r8), DIMENSION(kmG), INTENT(IN)  :: hgtl_G   ! height of GCM (at layer)
      REAL(kind=r8), DIMENSION(kmG), INTENT(OUT) :: val_G    ! target variable of GCM (at mid-layer)
       
      INTEGER, optional, INTENT(IN) :: kval(kmG) 
      LOGICAL, OPTIONAL, INTENT(IN) :: positive
       
      ! Local
      REAL(kind=r8), parameter :: fac = 1.0_r8
      REAL(kind=r8), DIMENSION(kdim+1) :: hgtl,val,y2
      REAL(kind=r8) :: yp1,ypkm
      INTEGER k,kbeg,k1,k2

      do k=2,kdim
       hgtl(k) = hgtl_C(k)
       val(k) = val_C(k)
      enddo

      hgtl(1) = hgt_C(1)
      hgtl(kdim+1) = hgt_C(kdim)

!     val(1) = val(2)
      val(1) = val(2) + (hgtl(1)-hgtl(2))*fac*(val(2)-val(3))/(hgtl(2)-hgtl(3))
      val(kdim+1) = val(kdim)
      
!---------------
      k=1
      do while (hgtl_G(k) .gt. hgtl(kdim+1))
      k = k + 1
      enddo
      kbeg = k
!---------------

      IF (PRESENT(kval)) THEN
      
        ! linear interpolation
        do k=kbeg,kmG
         K1 = kval(k)
         K2 = K1 + 1
         val_G(k) = FINTRP(1,hgtl_G(k),hgtl(K1),val(K1),hgtl(K2),val(K2))
        enddo
        
      ELSE
      
        ! cubic-spline interpolation
        yp1  = 0.5_r8*((val(2)-val(1))/(hgtl(2)-hgtl(1))   &
                      +(val(3)-val(2))/(hgtl(3)-hgtl(2)))

        ypkm = 0.5_r8*((val(kdim+1)-val(kdim))/(hgtl(kdim+1)-hgtl(kdim))  &
                      +(val(kdim)-val(kdim-1))/(hgtl(kdim)-hgtl(kdim-1)))

        call spline(hgtl,val,kdim+1,yp1,ypkm,y2)

        do k=kbeg,kmG
         call splint(hgtl,val,y2,kdim+1,hgtl_G(k),val_G(k))
        enddo
        
      ENDIF
      
      if (present(positive)) then
        do k=kbeg,kmG
         val_G(k) = MAX(0.0_r8,val_G(k))
        enddo
      endif
      
      if (kbeg .gt. 1) then
!     For the GCM grids above the (active) CRM vertical domain.

      IF (id_type .EQ. 0) THEN
        ! tendency data
        do k=1,kbeg-1
         val_G(k) = 0.0_r8
        enddo
      ELSE
        ! field data
        do k=1,kbeg-1
         val_G(k) = val_missing
        enddo
      ENDIF

      endif

      END SUBROUTINE VVM_l_to_vGCM_l

!===========================================================================================
      SUBROUTINE VVM_i_to_vGCM_l (kdim,kmG,id_type,val_C,hgt_C,hgtl_G,val_G)
!===========================================================================================
!     Vertical interpolation of VVM data at interfaces to that on GCM mid-layers
!     Val_G, hgtl_G: downward indexing (top to bottom)
!     Val_C, hgt_C : upward indexing (bottom to top)

      USE constld, only: val_missing
      USE utils,   only: spline,splint

      INTEGER, INTENT(IN) :: kdim      ! vertical dimension of VAL_C (nk2 or nk2_ext)
      INTEGER, INTENT(IN) :: kmG       ! vertical dimension of VAL_G (nLevel)
      INTEGER, INTENT(IN) :: id_type   ! data type (0: tendency data, 1: state data)

      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: val_C    ! target variable of VVM (at interface)
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN) :: hgt_C    ! height of VVM (at interface)

      REAL(kind=r8), DIMENSION(kmG), INTENT(IN)  :: hgtl_G   ! height of GCM (at layer)
      REAL(kind=r8), DIMENSION(kmG), INTENT(OUT) :: val_G    ! target variable of GCM (at mid-layer)
       
      ! Local
      REAL(kind=r8), DIMENSION(kdim) :: hgt,val,y2
      REAL(kind=r8) :: yp1,ypkm
      INTEGER k,kbeg

      do k=1,kdim
       hgt(k) = hgt_C(k)
       val(k) = val_C(k)
      enddo

!---------------
      k=1
      do while (hgtl_G(k) .gt. hgt(kdim))
      k = k + 1
      enddo
      kbeg = k
!---------------

      yp1  = 0.5_r8*((val(2)-val(1))/(hgt(2)-hgt(1))   &
                    +(val(3)-val(2))/(hgt(3)-hgt(2)))

      ypkm = 0.5_r8*((val(kdim)-val(kdim-1))/(hgt(kdim)-hgt(kdim-1))  &
                    +(val(kdim-1)-val(kdim-2))/(hgt(kdim-1)-hgt(kdim-2)))

      call spline(hgt,val,kdim,yp1,ypkm,y2)

      do k=kbeg,kmG
       call splint(hgt,val,y2,kdim,hgtl_G(k),val_G(k))
      enddo

      if (kbeg .gt. 1) then
!     For the GCM grids above the (active) CRM vertical domain.

      IF (id_type .EQ. 0) THEN
        ! tendency data
        do k=1,kbeg-1
         val_G(k) = 0.0_r8
        enddo
      ELSE
        ! field data
        do k=1,kbeg-1
         val_G(k) = val_missing
        enddo
      ENDIF

      endif

      END SUBROUTINE VVM_i_to_vGCM_l
      
!===========================================================================================
      SUBROUTINE OUTCON (IFORT)
!===========================================================================================
      USE parmsld, only: nk2,nk2_ext
      USE constld

      INTEGER , INTENT(IN) :: ifort    ! unit number for output

      REAL (kind=r4), DIMENSION(nk2_ext) :: ZZ_s,ZT_s,FNZ_s,FNT_s
      REAL (kind=r4), DIMENSION(nk2_ext) :: RHO_s,PBAR_s,PIBAR_s,THBAR_s
      REAL (kind=r4), DIMENSION(nk2_ext) :: RHOZ_s,PBARZ_s,PIBARZ_s

      REAL (kind=r4), DIMENSION(nk2) :: FN1_s,FN2_s

      INTEGER K

      write(ifort,*)
      write(ifort,*)
      write(ifort,*) '***************************************'
      write(ifort,*) 'LOGICAL PARAMETERS'
      write(ifort,*) '***************************************'
      write(ifort,*)
      write(ifort,*) 'nosfx     = ', nosfx
      write(ifort,*) 'notopo    = ', notopo
      write(ifort,*) 'nomap     = ', nomap
      write(ifort,*) 'zconst    = ', zconst
      write(ifort,*) 'physics   = ', physics
      write(ifort,*) 'turbulence= ', turbulence
      write(ifort,*) 'radcode   = ', radcode
      write(ifort,*) 'microcode = ', microcode
      write(ifort,*) 'buoy      = ', buoy
      write(ifort,*) 'notherm   = ', notherm
      write(ifort,*) 'fconst    = ', fconst
      write(ifort,*) 'localp    = ', localp
      write(ifort,*)
      write(ifort,*)
      write(ifort,*) '***************************************'
      write(ifort,*) 'INTEGER PARAMETERS'
      write(ifort,*) '***************************************'
      write(ifort,*)
      write(ifort,*) 'itt_max   = ', itt_max
      write(ifort,*) 'itinit    = ', itinit
      write(ifort,*) 'itstop    = ', itstop
      write(ifort,*) 'nsflux    = ', nsflux
      write(ifort,*) 'nrad      = ', nrad
      write(ifort,*) 'niterw    = ', niterw
      write(ifort,*) 'niterxy   = ', niterxy
      write(ifort,*) 'niterxy_i = ', niterxy_init
      write(ifort,*) 'klow_gwd  = ', klow_gwd
      
      write(ifort,*) 'iver_wind = ', iver_wind          ! JUNG_temp
      
      write(ifort,*)
      write(ifort,*)
      write(ifort,*) '***************************************'
      write(ifort,*) 'REAL PARAMETERS'
      write(ifort,*) '***************************************'
      write(ifort,*)
      write(ifort,*) 'piv       = ', piv
      write(ifort,*) 'rearth    = ', rearth
      write(ifort,*) 'omega     = ', omega
      write(ifort,*) 'epsilo    = ', epsilo
      write(ifort,*) 'cappa     = ', cappa
      write(ifort,*) 'vk        = ', vk
      write(ifort,*) 'hlm       = ', hlm
      write(ifort,*) 'hlf       = ', hlf
      write(ifort,*) 'grav      = ', grav
      write(ifort,*) 'cp        = ', cp
      write(ifort,*) 'gas_const = ', gas_const
      write(ifort,*) 'wvg_const = ', wvg_const
      write(ifort,*) 'delta     = ', delta
      write(ifort,*) 'h2otrip   = ', h2otrip
      write(ifort,*) 't_melt    = ', t_melt
      write(ifort,*) 'cp_h2o    = ', cp_h2o
      write(ifort,*) 'cp_ice    = ', cp_ice
      write(ifort,*) 'cpwv      = ', cpwv
      write(ifort,*) 'sec_cday  = ', sec_cday
      write(ifort,*) 'sec_sday  = ', sec_sday
      write(ifort,*) 'pstd      = ', pstd
      write(ifort,*) 'rho_h2o   = ', rho_h2o
      write(ifort,*) 'pzero     = ', pzero
      write(ifort,*) 'hls       = ', hls
      write(ifort,*) 'rad2deg   = ', rad2deg
      write(ifort,*) 'dt        = ', dt
      write(ifort,*) 'dt_gcm    = ', dt_gcm
      write(ifort,*) 'dz        = ', dz
      write(ifort,*) 'dzsq      = ', dzsq
      write(ifort,*) 'dz_ext    = ', dz_ext
      write(ifort,*) 'dx        = ', dx
      write(ifort,*) 'dy        = ', dy
      write(ifort,*) 'dxsq      = ', dxsq
      write(ifort,*) 'dysq      = ', dysq
      write(ifort,*) 'dxdy      = ', dxdy
      write(ifort,*) 'dx_gcm    = ', dx_gcm
      write(ifort,*) 'dy_gcm    = ', dy_gcm
      write(ifort,*) 'zrsea     = ', zrsea
      write(ifort,*) 'zrland    = ', zrland
      write(ifort,*) 'sst       = ', sst
      write(ifort,*) 'psfc      = ', psfc
      write(ifort,*) 'pisfc     = ', pisfc
      write(ifort,*) 'aladv     = ', aladv
      write(ifort,*) 'rxtau_q   = ', rxtau_q
      write(ifort,*) 'rxtau_d   = ', rxtau_d
      write(ifort,*) 'wrxmu     = ', wrxmu
      write(ifort,*) 'dtrad     = ', dtrad
      write(ifort,*) 'crad      = ', crad
      write(ifort,*) 'gwdht     = ', gwdht
      write(ifort,*) 'dthmax    = ', dthmax
      write(ifort,*) 'z1pert    = ', z1pert
      write(ifort,*) 'z2pert    = ', z2pert
      write(ifort,*)
      write(ifort,*)
      write(ifort,*) '***************************************'
      write(ifort,*) 'BASIC PROFILES'
      write(ifort,*) '***************************************'
      write(ifort,*)
      write(ifort,*) '========================================================='
      write(ifort,*) ' K ,     ZZ(K)  ,    ZT(K)   ,   FNZ(K)  ,   FNT(K)   '
      write(ifort,*) '========================================================='

      do K = 1,nk2_ext
       ZZ_s(K)  = ZZ(K)
       ZT_s(K)  = ZT(K)
       FNZ_s(K) = FNZ(K)
       FNT_s(K) = FNT(K)
       write(ifort,'(i3,4F12.3)') K,ZZ_s(K),ZT_s(K),FNZ_s(K),FNT_s(K)
      enddo

      write(ifort,*)
      write(ifort,*) '========================================================='
      write(ifort,*) ' K ,    RHO(K)  ,  PBAR(K)  ,  PIBAR(K)  , THBAR(K) '
      write(ifort,*) '========================================================='

      do K = 1,nk2_ext
       RHO_s(K)   = RHO(K)
       THBAR_s(K) = THBAR(K)
       PBAR_s(K)  = PBAR(K)
       PIBAR_s(K) = PIBAR(K)
       write(ifort,'(i3,4F12.3)') K,RHO_s(K),PBAR_s(K),PIBAR_s(K),THBAR_s(K)
      enddo

      write(ifort,*)
      write(ifort,*) '========================================================='
      write(ifort,*) ' K ,   RHOZ(K)  ,  PBARZ(K) ,  PIBARZ(K)  '
      write(ifort,*) '========================================================='

      do K = 1,nk2_ext
       RHOZ_s(K) = RHOZ(K)
       PBARZ_s(K) = PBARZ(K)
       PIBARZ_s(K) = PIBARZ(K)
       write(ifort,'(i3,3F12.3)') K,RHOZ_s(K),PBARZ_s(K),PIBARZ_s(K)
      enddo

      write(ifort,*)
      write(ifort,*) '========================================================='
      write(ifort,*) ' K ,     FN1(K)   ,  FN2(K)  '
      write(ifort,*) '========================================================='

      do K = 1,nk2
       FN1_s(K) = FN1(K)
       FN2_s(K) = FN2(K)
       write(ifort,'(i3,2F12.3)') K,FN1_s(K),FN2_s(K)
      enddo

   END SUBROUTINE outcon

END MODULE q3d_interface_module

! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
 module mpas_init_atm_gocart
 use mpas_derived_types
 use mpas_kind_types
 use mpas_dmpar
 use mpas_pool_routines

 use init_atm_read_met
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_vinterp
 use mpas_atmphys_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: init_atm_gocart


!mpas_init_atm_gocart contains the subroutines needed for the interpolation of climatological
!monthly-averaged hygroscopic ("water friendly") and nonhygroscopic ("ice friendly") aerosols
!used in the Thompson parameterization of cloud microphysics with CCN and IN nucleation (use
!physics option mp_thompson_aerorosols=true).
!Laura D. Fowler (laura@ucar.edu) / 2016-05-27.


 contains


!=================================================================================================================
 subroutine init_atm_gocart(block,mesh,configs,diag,state)
!=================================================================================================================

!input arguments:
 type (mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(block_type),intent(inout),target:: block 
 type(mpas_pool_type),intent(inout)   :: mesh
 type(mpas_pool_type),intent(inout)   :: state

!local variables and pointers:
 character (len=StrKIND),pointer:: config_start_time
 character(len=StrKIND):: filename_gocart
 character(len=StrKIND):: initial_date,mess

 logical:: lexist

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine init_atm_gocart:'

!inquire if the GOCART input file exists:
 filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.dat"

 inquire(file=filename_gocart,exist=lexist)
 if(lexist) then

    call mpas_pool_get_config(configs,'config_start_time',config_start_time)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART data to the MPAS mesh:  
    call init_hinterp_gocart(block,mesh)

    !--- interpolation of the monthly-averaged GOCART data to the initial date, and vertical interpolation to
    !    the MPAS levels:
    initial_date = trim(config_start_time)
    call init_vinterp_gocart(initial_date,mesh,diag,state)
 else
    write(mess,*) 'QNWFA_QNIFA_Monthly_Sigma.dat was not found in local directory:'
    call physics_message(mess)
    write(mess,*) 'nwfa and nifa are set to zero and not interpolated from climatological data.'
    call physics_message(mess)
 endif

 write(0,*) '--- end subroutine init_atm_gocart.'
 write(0,*)

 end subroutine init_atm_gocart

!=================================================================================================================
 subroutine init_vinterp_gocart(initial_date,mesh,diag,state)
!=================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: initial_date
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: mesh
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 integer,pointer:: nCells,nGocartLevels,nVertLevels
 integer,pointer:: index_nifa,index_nwfa
 integer:: iCell,k,kk

 real(kind=RKIND),dimension(:,:),pointer  :: nifa,nwfa,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: nifa_clim,nwfa_clim,pwif_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: nifa_int,nwfa_int,pwif_int,sorted_arr

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine init_vinterp_gocart:'

 call mpas_pool_get_dimension(mesh,'nCells'       ,nCells       )
 call mpas_pool_get_dimension(mesh,'nGocartLevels',nGocartLevels)
 call mpas_pool_get_dimension(mesh,'nVertLevels'  ,nVertLevels  )

 call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
    
 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_array(mesh,'nifa_clim',nifa_clim)
 call mpas_pool_get_array(mesh,'nwfa_clim',nwfa_clim)
 call mpas_pool_get_array(mesh,'pwif_clim',pwif_clim)

 call mpas_pool_get_array(state,'scalars',scalars)
 nifa => scalars(index_nifa,:,:)
 nwfa => scalars(index_nwfa,:,:)

 if(.not.allocated(nifa_int)  ) allocate(nifa_int(nGocartLevels,nCells))
 if(.not.allocated(nwfa_int)  ) allocate(nwfa_int(nGocartLevels,nCells))
 if(.not.allocated(pwif_int)  ) allocate(pwif_int(nGocartLevels,nCells))
 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nGocartLevels))

!--- interpolation of the monthly-averaged GOCART data to the initial date, and vertical interpolation to the
!    MPAS levels:
 do k = 1, nGocartLevels
    call monthly_interp_to_date(nCells,initial_date,nifa_clim(:,k,:),nifa_int(k,:))
    call monthly_interp_to_date(nCells,initial_date,nwfa_clim(:,k,:),nwfa_int(k,:))
    call monthly_interp_to_date(nCells,initial_date,pwif_clim(:,k,:),pwif_int(k,:))
 enddo

 do iCell = 1, nCells
    !--- nifa:
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    do k = 1, nGocartLevels
       kk = nGocartLevels + 1 -k
       sorted_arr(1,kk) = pwif_int(k,iCell)
       sorted_arr(2,kk) = nifa_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       nifa(k,iCell) = vertical_interp(target_p,nGocartLevels-1, &
                                sorted_arr(:,1:nGocartLevels-1),order=1,extrap=0)
       if(target_p >= pwif_int(1,iCell) .and. k < nVertLevels) nifa(k,iCell) = nifa(k+1,iCell)
    enddo

    !--- nwfa:
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    do k = 1, nGocartLevels
       kk = nGocartLevels + 1 -k
       sorted_arr(1,kk) = pwif_int(k,iCell)
       sorted_arr(2,kk) = nwfa_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       nwfa(k,iCell) = vertical_interp(target_p,nGocartLevels-1, &
                                sorted_arr(:,1:nGocartLevels-1),order=1,extrap=0)
       if(target_p >= pwif_int(1,iCell) .and. k < nVertLevels) nwfa(k,iCell) = nwfa(k+1,iCell)
    enddo
 enddo

 if(allocated(nifa_int)  ) deallocate(nifa_int  )
 if(allocated(nwfa_int)  ) deallocate(nwfa_int  )
 if(allocated(pwif_int)  ) deallocate(pwif_int  )
 if(allocated(sorted_arr)) deallocate(sorted_arr)
 
 write(0,*) '--- enter subroutine init_vinterp_gocart:'

 end subroutine init_vinterp_gocart

!=================================================================================================================
 subroutine init_hinterp_gocart(block,mesh)
!=================================================================================================================

!inout arguments:
 type(block_type),intent(inout),target:: block 
 type (mpas_pool_type),intent(inout)  :: mesh

!local variables:
 type(dm_info),pointer:: dminfo
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 character(len=StrKIND):: filename_gocart
 logical:: have_landmask

 integer,pointer:: nCells
 integer:: i,j
 integer:: iCell,istatus,k,masked,nmonths,nInterpPoints
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:,:),pointer:: nifa_clim,nwfa_clim,pwif_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: destField3d

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine init_hinterp_gocart:'

 dminfo => block%domain%dminfo

 filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.dat"

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_array(mesh,'nifa_clim',nifa_clim)
 call mpas_pool_get_array(mesh,'nwfa_clim',nwfa_clim)
 call mpas_pool_get_array(mesh,'pwif_clim',pwif_clim)

!open intermediate file:
 istatus = 0
 call read_gocart_init(trim(filename_gocart),istatus)
 if(istatus /= 0) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Error opening gocart file ', trim(filename_gocart)
    write(0,*) '********************************************************************************'
    call mpas_dmpar_abort(dminfo)
 else
    write(0,*) 'Processing file ',trim(filename_gocart)
 end if

!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field % field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field%nx+3,field%ny))

       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
!      write(0,*) 'minval, maxval of LANDSEA = ', minval(maskslab), maxval(maskslab)
    end if
    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do
 call read_met_close()

 if(.not. have_landmask) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Landsea mask not available from the surface file '
    write(0,*) '********************************************************************************'
 end if


!read gocart data:
 istatus = 0
 call read_gocart_init(trim(filename_gocart),istatus)
 if(istatus /= 0) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Error opening gocart file ', trim(filename_gocart)
    write(0,*) '********************************************************************************'
    call mpas_dmpar_abort(dminfo)
 endif
 call read_next_met_field(field, istatus)

!horizontally interpolate GOCART data:
 do while(istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = 0.0

    mask_array => landmask

    if(index(field % field, 'QNIFA_JAN') /= 0 .or. &
       index(field % field, 'QNIFA_FEB') /= 0 .or. &
       index(field % field, 'QNIFA_MAR') /= 0 .or. &
       index(field % field, 'QNIFA_APR') /= 0 .or. &
       index(field % field, 'QNIFA_MAY') /= 0 .or. &
       index(field % field, 'QNIFA_JUN') /= 0 .or. &
       index(field % field, 'QNIFA_JUL') /= 0 .or. &
       index(field % field, 'QNIFA_AUG') /= 0 .or. &
       index(field % field, 'QNIFA_SEP') /= 0 .or. &
       index(field % field, 'QNIFA_OCT') /= 0 .or. &
       index(field % field, 'QNIFA_NOV') /= 0 .or. &
       index(field % field, 'QNIFA_DEC') /= 0 .or. &
       index(field % field, 'QNWFA_JAN') /= 0 .or. &
       index(field % field, 'QNWFA_FEB') /= 0 .or. &
       index(field % field, 'QNWFA_MAR') /= 0 .or. &
       index(field % field, 'QNWFA_APR') /= 0 .or. &
       index(field % field, 'QNWFA_MAY') /= 0 .or. &
       index(field % field, 'QNWFA_JUN') /= 0 .or. &
       index(field % field, 'QNWFA_JUL') /= 0 .or. &
       index(field % field, 'QNWFA_AUG') /= 0 .or. &
       index(field % field, 'QNWFA_SEP') /= 0 .or. &
       index(field % field, 'QNWFA_OCT') /= 0 .or. &
       index(field % field, 'QNWFA_NOV') /= 0 .or. &
       index(field % field, 'QNWFA_DEC') /= 0 .or. &
       index(field % field, 'P_WIF_JAN') /= 0 .or. &
       index(field % field, 'P_WIF_FEB') /= 0 .or. &
       index(field % field, 'P_WIF_MAR') /= 0 .or. &
       index(field % field, 'P_WIF_APR') /= 0 .or. &
       index(field % field, 'P_WIF_MAY') /= 0 .or. &
       index(field % field, 'P_WIF_JUN') /= 0 .or. &
       index(field % field, 'P_WIF_JUL') /= 0 .or. &
       index(field % field, 'P_WIF_AUG') /= 0 .or. &
       index(field % field, 'P_WIF_SEP') /= 0 .or. &
       index(field % field, 'P_WIF_OCT') /= 0 .or. &
       index(field % field, 'P_WIF_NOV') /= 0 .or. &
       index(field % field, 'P_WIF_DEC') /= 0) then

       !
       !set up projection:
       !
       call map_init(proj)
   
       if(field%iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON,proj, &
                       latinc = real(field%deltalat,RKIND), &
                       loninc = real(field%deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1   = real(field%startlat,RKIND), &
                       lon1   = real(field%startlon,RKIND))
       elseif(field%iproj == PROJ_GAUSS) then
          call map_set(PROJ_GAUSS,proj, &
                       nlat = nint(field%deltalat), &
                       loninc = 360.0_RKIND / real(field%nx,RKIND), &
                       lat1 = real(field%startlat,RKIND), &
                       lon1 = real(field%startlon,RKIND))
       endif
       
       !
       !horizontally interpolate field at level k:
       !
       if(index(field%field,'QNIFA_JAN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_JAN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'QNIFA_FEB') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_FEB at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'QNIFA_MAR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_MAR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'QNIFA_APR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_APR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'QNIFA_MAY') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_MAY at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'QNIFA_JUN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_JUN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'QNIFA_JUL') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_JUL at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'QNIFA_AUG') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_AUG at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'QNIFA_SEP') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_SEP at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'QNIFA_OCT') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_OCT at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'QNIFA_NOV') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_NOV at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'QNIFA_DEC') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNIFA_DEC at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 12
       elseif(index(field%field,'QNWFA_JAN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_JAN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'QNWFA_FEB') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_FEB at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'QNWFA_MAR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_MAR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'QNWFA_APR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_APR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'QNWFA_MAY') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_MAY at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'QNWFA_JUN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_JUN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'QNWFA_JUL') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_JUL at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'QNWFA_AUG') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_AUG at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'QNWFA_SEP') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_SEP at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'QNWFA_OCT') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_OCT at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'QNWFA_NOV') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_NOV at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'QNWFA_DEC') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating QNWFA_DEC at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 12
       elseif(index(field%field,'P_WIF_JAN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_JAN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'P_WIF_FEB') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_FEB at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'P_WIF_MAR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_MAR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'P_WIF_APR') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_APR at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'P_WIF_MAY') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_MAY at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'P_WIF_JUN') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_JUN at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'P_WIF_JUL') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_JUL at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'P_WIF_AUG') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_AUG at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'P_WIF_SEP') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_SEP at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'P_WIF_OCT') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_OCT at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'P_WIF_NOV') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_NOV at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'P_WIF_DEC') /= 0) then
          k = field%xlvl
          write(0,*) 'Interpolating P_WIF_DEC at ',k
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 12
       endif

       allocate(rslab(-2:field%nx+3,field%ny))
       rslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       rslab(0,1:field%ny)  = field%slab(field%nx  ,1:field%ny)
       rslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       rslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       rslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       rslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       rslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)

       do iCell = 1, nInterpPoints
          if(mask_array(iCell) /= masked) then
             lat = latPoints(iCell)*DEG_PER_RAD
             lon = lonPoints(iCell)*DEG_PER_RAD
             call latlon_to_ij(proj,lat,lon,x,y)
             if(x < 0.5) then
                lon = lon + 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             elseif(x > real(field%nx,kind=RKIND)+ 0.5) then
                lon = lon - 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             endif

             if(maskval /= -1.0) then
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                                  interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                                  interp_list,1)
             endif
          else
             destField3d(nmonths,k,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo

 call read_met_close()

 write(0,*) '--- end subroutine init_hinterp_gocart:'

 end subroutine init_hinterp_gocart
 
!=================================================================================================================
 subroutine read_gocart_init(fg_source,istatus)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: fg_source

!output arguments:
 integer,intent(out):: istatus

!local variables:
 logical:: is_used
 integer:: io_status

!-----------------------------------------------------------------------------------------------------------------

 istatus = 0

 do input_unit = 10, 100
    inquire(unit=input_unit,opened=is_used)
    if (.not. is_used) exit
 end do
 if(input_unit > 100) write(0,*) 'Error: In read_gocart_init(), couldn''t find an available Fortran unit.'
 open(unit=input_unit, file=trim(fg_source), status='old', form='unformatted', iostat=io_status)

 if (io_status > 0) istatus = 1

 end subroutine read_gocart_init

!=================================================================================================================
 end module mpas_init_atm_gocart
!=================================================================================================================

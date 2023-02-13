module conditional_diag_restart

  implicit none
  private

  public :: cnd_diag_init_restart
  public :: cnd_diag_write_restart
  public :: cnd_diag_read_restart

contains
  !==================================================================================================
  subroutine cnd_diag_init_restart( dimids, hdimcnt, pver, pver_id, pverp_id,            &! in
                                    file, cnd_metric_desc,  cnd_flag_desc,               &! inout
                                    cnd_qoi_val_desc, cnd_qoi_old_desc, cnd_qoi_inc_desc )! inout
  !------------------------------------------------------------------------------------------------
  ! Purpose: add variables to the restart file for conditional sampling and diagnostics.
  ! History: First version by Hui Wan, PNNL, 2021-04
  !------------------------------------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use conditional_diag, only: cnd_diag_info
  use pio,              only: file_desc_t, var_desc_t, pio_def_var, pio_double

  integer, intent(in) :: dimids(4)
  integer, intent(in) :: hdimcnt
  integer, intent(in) :: pver
  integer, intent(in) :: pver_id
  integer, intent(in) :: pverp_id

  type(file_desc_t), intent(inout) :: file

  type(var_desc_t),allocatable,intent(inout) :: cnd_metric_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_flag_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_val_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_inc_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_old_desc(:,:)

  ! Local variables

  integer :: ierr, ndims
  integer :: dimids_local(4)
  character(len=256) :: pname  !variable name in restart file

  integer :: ncnd, nchkpt, nqoi, icnd, ichkpt, iqoi, nver
  character(len=*),parameter :: subname = 'cnd_diag_init_restart'

  if (cnd_diag_info%ncnd <= 0 ) return

  !-------------------------------------------------------------------------------
  ! Copy dimension IDs to a local array.
  !-------------------------------------------------------------------------------
  ! Because the various metrics, flags and QoIs in the conditional diagnostics 
  ! data structre can have different sizes in the vertical dimension,
  ! the last element of the local variable dimids_local might get different values
  ! during this subroutine. The use of a local variable here ensures that 
  ! the array dimids in the calling subroutine remains intact.
  !-------------------------------------------------------------------------------
  dimids_local(:) = dimids(:)

  !-------------------------------------------------------------------------------
  ! Allocate memeory for the variable description arrays
  !-------------------------------------------------------------------------------
  ! (Question: would it be better to allocate and deallocate at the beginning
  ! and end of each run instead of each time step?)

  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  allocate( cnd_metric_desc(ncnd) )
  allocate( cnd_flag_desc(ncnd) )

  if (nqoi>0) then
     allocate( cnd_qoi_old_desc(       ncnd,nqoi) )
     allocate( cnd_qoi_val_desc(nchkpt,ncnd,nqoi) )
     allocate( cnd_qoi_inc_desc(nchkpt,ncnd,nqoi) )
  end if

  !-----------------------------------------------------------
  ! Add the metrics and corresponding flags to restart file
  !-----------------------------------------------------------
  do icnd = 1,ncnd

     nver = cnd_diag_info%metric_nver(icnd)

     ! Dimension information

     if ( nver == 1) then
        ndims = hdimcnt

     else if ( nver == pver ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pver_id

     else if ( nver == pver+1 ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pverp_id

     else
        call endrun(subname//': check cnd_diag_info%metric_nver')
     end if

     ! Add the metric variable

     write(pname,'(a,i2.2,a)') 'cnd',icnd,'_metric'
     ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_metric_desc(icnd))

     ! Add the flag variable

     write(pname,'(a,i2.2,a)') 'cnd',icnd,'_flag'
     ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims),   cnd_flag_desc(icnd))

  end do

  !--------------------------------
  ! Add the QoIs to restart file
  !--------------------------------
  do iqoi = 1,nqoi

     ! Dimension information

     nver = cnd_diag_info%qoi_nver_save(iqoi)

     if ( nver == 1) then
        ndims = hdimcnt

     else if ( nver == pver ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pver_id

     else if ( nver == pver+1 ) then
        ndims = hdimcnt+1
        dimids_local(ndims) = pverp_id

     else
        call endrun(subname//': check cnd_diag_info%qoi_nver_save')
     end if

     ! Add to restart file the variables containing QoIs at various checkpoints 

     if (cnd_diag_info%l_output_state) then
        do icnd = 1,ncnd
         do ichkpt = 1,nchkpt
            write(pname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_',ichkpt
            ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_qoi_val_desc(ichkpt,icnd,iqoi))
         end do
        end do
     end if

     ! Add to restart file the variables containing increments at various checkpoints 

     if (cnd_diag_info%l_output_incrm) then

        do icnd = 1,ncnd
         do ichkpt = 1,nchkpt
            write(pname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_inc',ichkpt
            ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_qoi_inc_desc(ichkpt,icnd,iqoi))
         end do
        end do

        ! Add to restart file the variable containing the "old" value of the field 

        do icnd = 1,ncnd
           write(pname,'(2(a,i2.2),a)') 'cnd',icnd, '_qoi',iqoi, '_old'
           ierr = pio_def_var(File, trim(pname), pio_double, dimids_local(1:ndims), cnd_qoi_old_desc(icnd,iqoi))
        end do

     end if

  end do !iqoi

  end subroutine cnd_diag_init_restart
  !=========================================================

  subroutine cnd_diag_write_restart( phys_diag, begchunk, endchunk,         &! in
                                     physgrid, file_hdimsizes, file_nhdims, &! in
                                     pcols, chunk_ncols, fillvalue,         &! in
                                     file, cnd_metric_desc,  cnd_flag_desc, &! inout
                                     cnd_qoi_val_desc, cnd_qoi_inc_desc,    &! inout
                                     cnd_qoi_old_desc                       )! inout
  !------------------------------------------------------------------------------------------------
  ! Purpose: add variables to the restart file for conditional sampling and diagnostics.
  ! History: First version by Hui Wan, PNNL, 2021-04
  !------------------------------------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use conditional_diag, only: cnd_diag_info, cnd_diag_t
  use pio,              only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_write_darray
  use cam_grid_support, only: cam_grid_get_decomp, cam_grid_write_dist_array
  use shr_kind_mod,     only: r8 => shr_kind_r8

  integer, intent(in) :: begchunk, endchunk
  type(cnd_diag_t),intent(in) :: phys_diag(begchunk:endchunk)

  integer, intent(in) :: physgrid
  integer, intent(in) :: file_nhdims                  ! number of horizontal dimensions in restart file 
  integer, intent(in) :: file_hdimsizes(file_nhdims)  ! horizontal dimension sizes in restart file

  integer, intent(in) :: pcols
  integer, intent(in) :: chunk_ncols(begchunk:endchunk)
  real(r8),intent(in) :: fillvalue

  type(file_desc_t), intent(inout) :: file

  type(var_desc_t),allocatable,intent(inout) :: cnd_metric_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_flag_desc(:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_val_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_inc_desc(:,:,:)
  type(var_desc_t),allocatable,intent(inout) :: cnd_qoi_old_desc(:,:)

  ! Local variables

  integer :: file_dims(3)   ! dimension sizes in restart file, local variable
  integer :: arry_dims(3)   ! dimension sizes of array holding values to be written out, local variable

  type(io_desc_t), pointer :: iodesc

  integer :: ierr
  integer :: lchnk,ncol
  integer :: ncnd, nchkpt, nqoi, icnd, ichkpt, iqoi, nver
  character(len=*),parameter :: subname = 'cnd_write_init_restart'

  real(r8):: tmpfield_2d_1(pcols, begchunk:endchunk)
  real(r8):: tmpfield_2d_2(pcols, begchunk:endchunk)

  real(r8), allocatable :: tmpfield_3d_1(:,:,:)
  real(r8), allocatable :: tmpfield_3d_2(:,:,:)


  if (cnd_diag_info%ncnd <= 0 ) return

  !---------------------------------------------------------------
  ! Gather dimension info and save in local variables
  !---------------------------------------------------------------
  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  file_dims(1:file_nhdims) = file_hdimsizes(1:file_nhdims)

  !-------------------
  ! metrics and flags
  !-------------------
  do icnd = 1,ncnd

     nver = cnd_diag_info%metric_nver(icnd)

     if (nver==1) then

        !--------------------------------------------
        ! get iodesc needed by pio_write_darray calls

        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1

        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_dims(1:file_nhdims), pio_double, iodesc) ! 4xin, out

        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out

        tmpfield_2d_1 = fillvalue
        tmpfield_2d_2 = fillvalue

        do lchnk = begchunk, endchunk
           ncol = chunk_ncols(lchnk) 
           tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)% metric(:ncol,1)
           tmpfield_2d_2(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%   flag(:ncol,1)
        end do

        call pio_write_darray(File, cnd_metric_desc(icnd), iodesc, tmpfield_2d_1, ierr)
        call pio_write_darray(File,   cnd_flag_desc(icnd), iodesc, tmpfield_2d_2, ierr)

     else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_write_dist_array calls

        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d_1(pcols,nver,begchunk:endchunk) )
        allocate( tmpfield_3d_2(pcols,nver,begchunk:endchunk) )

        tmpfield_3d_1 = fillvalue
        tmpfield_3d_2 = fillvalue
        
        !----------------------------------------------------------------
        ! pack metric and flag values into tmp arrays and write them out

        do lchnk = begchunk, endchunk
           ncol = chunk_ncols(lchnk)
           tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)% metric(:ncol,:)
           tmpfield_3d_2(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%   flag(:ncol,:)
        end do

        call cam_grid_write_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d_1, cnd_metric_desc(icnd))
        call cam_grid_write_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d_2,   cnd_flag_desc(icnd))

        deallocate(tmpfield_3d_1)
        deallocate(tmpfield_3d_2)

     end if 
  end do

  !----------------------------------------------------------
  ! Conditionally sampled QoIs: "old", snapshots, and "inc"
  !----------------------------------------------------------
  do iqoi = 1,nqoi

     nver = cnd_diag_info%qoi_nver_save(iqoi)

     if (nver==1) then

        tmpfield_2d_1 = fillvalue

        !--------------------------------------------
        ! get iodesc needed by pio_write_darray calls
        !--------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1
 
        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_dims(1:file_nhdims), pio_double, iodesc) ! 4xin, out

        !-------------------------------------------------------
        ! pack field values into tmp arrays and write them out
        !-------------------------------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% val(:ncol,1,ichkpt)
              end do

              call pio_write_darray(File, cnd_qoi_val_desc(ichkpt,icnd,iqoi), iodesc, tmpfield_2d_1, ierr)

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !-------------------------------------------------------------------
        ! pack increment and old values into tmp arrays and write them out
        !-------------------------------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           do icnd = 1,ncnd

              ! increments corresponding to various checkpoints 

              do ichkpt = 1,nchkpt

                 do lchnk = begchunk, endchunk
                    ncol = chunk_ncols(lchnk)
                    tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% inc(:ncol,1,ichkpt)
                 end do

                 call pio_write_darray(File, cnd_qoi_inc_desc(ichkpt,icnd,iqoi), iodesc, tmpfield_2d_1, ierr)

              end do !ichkpt

              ! the "old" values
 
              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_2d_1(:ncol,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)%old(:ncol,1)
              end do

              call pio_write_darray(File, cnd_qoi_old_desc(icnd,iqoi), iodesc, tmpfield_2d_1, ierr)

           end do !icnd

        end if !cnd_diag_info%l_output_incrm


    else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_write_dist_array calls
        !----------------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d_1(pcols,nver,begchunk:endchunk) )
        tmpfield_3d_1 = fillvalue
        
        !-------------------------------------------------------
        ! pack field values into tmp arrays and write them out
        !-------------------------------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% val(:ncol,:,ichkpt)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), &
                                              tmpfield_3d_1, cnd_qoi_val_desc(ichkpt,icnd,iqoi)            )

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !-------------------------------------------------------------------
        ! pack increment and old values into tmp arrays and write them out
        !-------------------------------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           ! increments associated with various checkpoints 

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% inc(:ncol,:,ichkpt)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), &
                                              tmpfield_3d_1, cnd_qoi_inc_desc(ichkpt,icnd,iqoi)            )

           end do
           end do

           ! the "old" values

           do icnd = 1,ncnd

              do lchnk = begchunk, endchunk
                 ncol = chunk_ncols(lchnk)
                 tmpfield_3d_1(:ncol,:,lchnk) = phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% old(:ncol,:)
              end do

              call cam_grid_write_dist_array( File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), &
                                              tmpfield_3d_1, cnd_qoi_old_desc(icnd,iqoi)             )

           end do

        end if !cnd_diag_info%l_output_incrm

        !----------
        ! Clean up
        !----------
        deallocate(tmpfield_3d_1)

    end if 
  end do

  !-------------------------------------------------------------------------
  ! Done writing variables for restart. Dealocate description info arrays.
  ! (Question: would it be better to allocate and deallocate at the beginning
  ! and end of each run instead of each time step?)
  !-------------------------------------------------------------------------
  deallocate( cnd_metric_desc )
  deallocate( cnd_flag_desc )

  if (nqoi>0) then
     deallocate( cnd_qoi_old_desc )
     deallocate( cnd_qoi_val_desc )
     deallocate( cnd_qoi_inc_desc )
  end if
 
  end subroutine cnd_diag_write_restart

  !======================================================================
  subroutine cnd_diag_read_restart( phys_diag, begchunk, endchunk,         &! in
                                    physgrid, file_hdimsizes, file_nhdims, &! in
                                    pcols, fillvalue, file                 )! in, in, inout 
  !------------------------------------------------------------------------------------------------
  ! Purpose: read variables for conditional sampling and budget analysis from restart file
  ! History: First version by Hui Wan, PNNL, 2021-04
  !------------------------------------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use conditional_diag, only: cnd_diag_info, cnd_diag_t
  use pio,              only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_read_darray, pio_inq_varid
  use cam_grid_support, only: cam_grid_get_decomp, cam_grid_read_dist_array
  use shr_kind_mod,     only: r8 => shr_kind_r8

  integer, intent(in) :: begchunk, endchunk
  type(cnd_diag_t),intent(inout) :: phys_diag(begchunk:endchunk)

  integer, intent(in) :: physgrid
  integer, intent(in) :: file_nhdims                  ! number of horizontal dimensions in restart file 
  integer, intent(in) :: file_hdimsizes(file_nhdims)  ! horizontal dimension sizes in restart file

  integer, intent(in) :: pcols
  real(r8),intent(in) :: fillvalue

  type(file_desc_t), intent(inout) :: file

  ! Local variables

  integer :: file_dims(3)   ! dimension sizes in restart file, local variable
  integer :: arry_dims(3)   ! dimension sizes of array holding values to be read in, local variable

  type(io_desc_t), pointer :: iodesc
  type(var_desc_t) :: vardesc

  integer :: ierr
  integer :: lchnk
  integer :: ncnd, nchkpt, nqoi, icnd, ichkpt, iqoi, nver

  real(r8) :: tmpfield_2d(pcols, begchunk:endchunk)
  real(r8), allocatable :: tmpfield_3d(:,:,:)

  character(len=256) :: varname  !variable name in restart file

  character(len=*),parameter :: subname = 'cnd_read_init_restart'

  !----------------------------------------

  if (cnd_diag_info%ncnd <= 0 ) return

  !---------------------------------------------------------------
  ! Gather dimension info and save in local variables
  !---------------------------------------------------------------
  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  file_dims(1:file_nhdims) = file_hdimsizes(1:file_nhdims)

  !-----------------
  ! metric and flag
  !-----------------
  do icnd = 1,ncnd

     nver = cnd_diag_info%metric_nver(icnd)

     if (nver==1) then

        tmpfield_2d = fillvalue
        !--------------------------------------------
        ! get iodesc needed by pio_read_darray calls
        !--------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1

        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_hdimsizes(1:file_nhdims), pio_double, iodesc)

        !----------------------------
        ! read and unpack the metric 
        !----------------------------
        write(varname,'(a,i2.2,a)') 'cnd',icnd,'_metric'
        ierr = pio_inq_varid(File, trim(varname), vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield_2d, ierr)

        do lchnk = begchunk,endchunk
           phys_diag(lchnk)%cnd(icnd)% metric(:,1) = tmpfield_2d(:,lchnk)
        end do

        !----------------------------
        ! read and unpack the flags
        !----------------------------
        write(varname,'(a,i2.2,a)') 'cnd',icnd,'_flag'
        ierr = pio_inq_varid(File, trim(varname), vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield_2d, ierr)

        do lchnk = begchunk,endchunk
           phys_diag(lchnk)%cnd(icnd)% flag(:,1) = tmpfield_2d(:,lchnk)
        end do

     else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_read_dist_array calls
        !----------------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d(pcols,nver,begchunk:endchunk) )
        tmpfield_3d = fillvalue
        
        !-------------------------------
        ! read and unpack metric values 
        !-------------------------------
        write(varname,'(a,i2.2,a)') 'cnd',icnd,'_metric'
        ierr = pio_inq_varid(File, trim(varname), vardesc)

        call cam_grid_read_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d, vardesc)

        do lchnk = begchunk, endchunk
           phys_diag(lchnk)%cnd(icnd)% metric(:,:) = tmpfield_3d(:,:,lchnk)
        end do

        !-------------------------------
        ! read and unpack flag values 
        !-------------------------------
        write(varname,'(a,i2.2,a)') 'cnd',icnd,'_flag'
        ierr = pio_inq_varid(File, trim(varname), vardesc)

        call cam_grid_read_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d, vardesc)

        do lchnk = begchunk, endchunk
           phys_diag(lchnk)%cnd(icnd)% flag(:,:) = tmpfield_3d(:,:,lchnk)
        end do

        !----------
        ! clean up
        !----------
        deallocate(tmpfield_3d)

     end if 
  end do

  !---------------------------------------------------------
  ! Conditionally sampled QoIs: "old", snapsots, and "inc"
  !---------------------------------------------------------
  do iqoi = 1,nqoi

     nver = cnd_diag_info%qoi_nver_save(iqoi)

     if (nver==1) then

        tmpfield_2d = fillvalue

        !--------------------------------------------
        ! get iodesc needed by pio_read_darray calls
        !--------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = endchunk - begchunk + 1
 
        call cam_grid_get_decomp(physgrid, arry_dims(1:2), file_dims(1:file_nhdims), pio_double, iodesc) ! 4xin, out

        !----------------------------------------------
        ! read and unpack QoIs at various checkpoints 
        !----------------------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              write(varname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_',ichkpt
              ierr = pio_inq_varid(File, trim(varname), vardesc)
              call pio_read_darray(File, vardesc, iodesc, tmpfield_2d, ierr)

              do lchnk = begchunk,endchunk
                 phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% val(:,1,ichkpt) = tmpfield_2d(:,lchnk)
              end do

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !-------------------------------------------
        ! read and unpack increments and old values 
        !-------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           do icnd = 1,ncnd

              ! increments corresponding to various checkpoints

              do ichkpt = 1,nchkpt

                 write(varname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_inc',ichkpt
                 ierr = pio_inq_varid(File, trim(varname), vardesc)
                 call pio_read_darray(File, vardesc, iodesc, tmpfield_2d, ierr)

                 do lchnk = begchunk,endchunk
                    phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% inc(:,1,ichkpt) = tmpfield_2d(:,lchnk)
                 end do

              end do !ichkpt

              ! the "old" values
 
              write(varname, '(2(a,i2.2),a)') 'cnd',icnd, '_qoi',iqoi, '_old'
              ierr = pio_inq_varid(File, trim(varname), vardesc)
              call pio_read_darray(File, vardesc, iodesc, tmpfield_2d, ierr)

              do lchnk = begchunk,endchunk
                 phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% old(:,1) = tmpfield_2d(:,lchnk)
              end do

           end do !icnd

        end if !cnd_diag_info%l_output_incrm

    else ! nver > 1

        !----------------------------------------------------
        ! prepare input for cam_grid_read_dist_array calls
        !----------------------------------------------------
        arry_dims(1) = pcols
        arry_dims(2) = nver
        arry_dims(3) = endchunk - begchunk + 1

        file_dims(file_nhdims+1) = nver

        allocate( tmpfield_3d(pcols,nver,begchunk:endchunk) )
        tmpfield_3d = fillvalue
        
        !------------------------------
        ! read and unpack field values 
        !------------------------------
        if (cnd_diag_info%l_output_state) then

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              write(varname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_',ichkpt
              ierr = pio_inq_varid(File, trim(varname), vardesc)

              call cam_grid_read_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d, vardesc)

              do lchnk = begchunk,endchunk
                 phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% val(:,:,ichkpt) = tmpfield_3d(:,:,lchnk)
              end do

           end do
           end do

        end if !cnd_diag_info%l_output_state

        !-------------------------------------------
        ! read and unpack increments and old values 
        !-------------------------------------------
        if (cnd_diag_info%l_output_incrm) then

           ! increments associated with various checkpoints 

           do icnd = 1,ncnd
           do ichkpt = 1,nchkpt

              write(varname,'(3(a,i2.2))') 'cnd',icnd, '_qoi',iqoi, '_inc',ichkpt
              ierr = pio_inq_varid(File, trim(varname), vardesc)

              call cam_grid_read_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d, vardesc)

              do lchnk = begchunk,endchunk
                 phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% inc(:,:,ichkpt) = tmpfield_3d(:,:,lchnk)
              end do

           end do
           end do

           ! the "old" values

           do icnd = 1,ncnd

              write(varname, '(2(a,i2.2),a)') 'cnd',icnd, '_qoi',iqoi, '_old'
              ierr = pio_inq_varid(File, trim(varname), vardesc)

              call cam_grid_read_dist_array(File, physgrid, arry_dims(1:3), file_dims(1:file_nhdims+1), tmpfield_3d, vardesc)

              do lchnk = begchunk,endchunk
                 phys_diag(lchnk)%cnd(icnd)%qoi(iqoi)% old(:,:) = tmpfield_3d(:,:,lchnk)
              end do

           end do

        end if !cnd_diag_info%l_output_incrm

        !----------
        ! Clean up
        !----------
        deallocate(tmpfield_3d)

    end if 
  end do
 
  end subroutine cnd_diag_read_restart

end module conditional_diag_restart

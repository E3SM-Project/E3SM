module restart_dynamics

! Write and read dynamics fields from the restart file.  For exact restart
! it is necessary to write all element data, including duplicate columns,
! to the file.  The namelist option, se_write_restart_unstruct, is
! available to write just the unique columns to the restart file using the
! same unstructured grid used by the history and initial files.  This
! results in the introduction of a roundoff size difference on restart, but
! writes the fields in the unstructured grid format which is easier to
! modify if the user desires to introduce perturbations or other
! adjustments into the run.  The restart file containing the unstructured
! grid format may also be used for an initial run.

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: iam

use constituents,     only: cnst_name
use dyn_grid,         only: timelevel, fvm, elem, edgebuf
use dyn_comp,         only: dyn_import_t, dyn_export_t, dyn_init, write_restart_unstruct
use hycoef,           only: init_restart_hycoef, write_restart_hycoef, &
                            hyai, hybi, ps0
use ref_pres,         only: ptop_ref

use pio,              only: pio_global, pio_unlimited, pio_offset_kind, pio_double, &
                            pio_seterrorhandling, pio_bcast_error, pio_noerr, &
                            file_desc_t, var_desc_t, io_desc_t, &
                            pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                            pio_def_dim, pio_def_var,  &
                            pio_enddef, &
                            pio_initdecomp, pio_freedecomp, pio_setframe, &
                            pio_put_att, pio_put_var, pio_write_darray, &
                            pio_get_att, pio_read_darray

use cam_pio_utils,    only: pio_subsystem, cam_pio_handle_error
use cam_grid_support, only: cam_grid_header_info_t, cam_grid_id, cam_grid_write_attr, &
                            cam_grid_write_var, cam_grid_get_decomp, cam_grid_dimensions, &
                            max_hcoordname_len, cam_grid_get_dim_names
use ncdio_atm,        only: infld

use infnan,           only: isnan
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

use parallel_mod,     only: par
use thread_mod,       only: horz_num_threads
use control_mod,      only: qsplit
use dimensions_mod,   only: np, npsq, ne, nlev, qsize, nelemd, nc, ntrac
use dof_mod,          only: UniquePoints
use element_mod,      only: element_t
use time_mod,         only: tstep, TimeLevel_Qdp

use edge_mod,         only: initEdgeBuffer, edgeVpack, edgeVunpack, FreeEdgeBuffer
use edgetype_mod,     only: EdgeBuffer_t
use bndry_mod,        only: bndry_exchange

use fvm_control_volume_mod, only: fvm_struct, n0_fvm

implicit none
private
save

public :: &
   init_restart_dynamics,  &
   write_restart_dynamics, &
   read_restart_dynamics

! these variables are module data so they can be shared between the
! file definition and write phases
type(var_desc_t)              :: psdry_desc, udesc, vdesc, tdesc
type(var_desc_t), allocatable :: qdesc_dp(:)
type(var_desc_t)              :: dp_fvm_desc
type(var_desc_t), pointer     :: c_fvm_desc(:)

integer, private :: nelem_tot = -1 ! Correct total number of elements

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine init_nelem_tot()
  use spmd_utils, only: mpicom, MPI_INTEGER, MPI_SUM

  integer :: ierr

  if (nelem_tot < 0) then
    call MPI_Allreduce(nelemd, nelem_tot, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
  end if
end subroutine init_nelem_tot

subroutine init_restart_dynamics(file, dyn_out)

   ! Define dimensions, variables, attributes for restart file.

   ! This is not really an "init" routine.  It is called before
   ! write_restart_dynamics every time an restart is written.

   ! arguments
   type(file_desc_t),  intent(inout) :: file
   type(dyn_export_t), intent(in)    :: dyn_out

   ! local variables
   integer :: i
   integer :: vdimids(2)
   integer :: nlev_dimid
   integer :: ncol_dimid
   integer :: ncol_fvm_dimid
   integer :: time_dimid

   integer :: ierr, err_handling

   integer :: grid_id
   type(cam_grid_header_info_t) :: info

   !----------------------------------------------------------------------------

   call init_nelem_tot()
   call init_restart_hycoef(file, vdimids)
   nlev_dimid = vdimids(1)

   call pio_seterrorhandling(File, pio_bcast_error, err_handling)

   ierr = PIO_Def_Dim(File, 'time', PIO_UNLIMITED, time_dimid)

   ! GLL restart fields

   ! number of columns written to restart depends on whether all columns in the
   ! element structures are written, or just the unique columns (unstructured grid)
   if (write_restart_unstruct) then
      grid_id = cam_grid_id('GLL')
      call cam_grid_write_attr(File, grid_id, info)
      ncol_dimid = info%get_hdimid(1)
   else
      ierr = PIO_Def_Dim(File,'nenpnp', nelem_tot*np*np, ncol_dimid)
      ierr = PIO_Put_Att(File, PIO_GLOBAL, 'ne', ne)
      ierr = PIO_Put_Att(File, PIO_GLOBAL, 'np', np)
   end if

   ierr = PIO_Def_Var(File, 'PSDRY', pio_double, (/ncol_dimid, time_dimid/), psdry_desc)
   ierr = PIO_Def_Var(File, 'U', pio_double, (/ncol_dimid, nlev_dimid, time_dimid/), Udesc)
   ierr = PIO_Def_Var(File, 'V', pio_double, (/ncol_dimid, nlev_dimid, time_dimid/), Vdesc)
   ierr = PIO_Def_Var(File, 'T', pio_double, (/ncol_dimid, nlev_dimid, time_dimid/), Tdesc)

   allocate(qdesc_dp(qsize))

   do i=1,qsize
      ierr = PIO_Def_Var(File,"dp"//trim(cnst_name(i)), pio_double, &
                         (/ncol_dimid, nlev_dimid, time_dimid/), Qdesc_dp(i))
   end do

   ! CSLAM restart fields

   if (ntrac > 0) then

      grid_id = cam_grid_id('FVM')
      call cam_grid_write_attr(File, grid_id, info)
      ncol_fvm_dimid = info%get_hdimid(1)

      ierr = PIO_Def_Var(File, 'dp_fvm', pio_double, &
         (/ncol_fvm_dimid, nlev_dimid, time_dimid/), dp_fvm_desc)

      allocate(c_fvm_desc(ntrac))
      do i = 1, ntrac
         ierr = PIO_Def_Var(File, trim(cnst_name(i))//"_fvm", pio_double, &
            (/ncol_fvm_dimid, nlev_dimid, time_dimid/), c_fvm_desc(i))
      end do

   end if

   call pio_seterrorhandling(File, err_handling)

end subroutine init_restart_dynamics

!=========================================================================================

subroutine write_restart_dynamics(File, dyn_out)

   type(file_desc_t), intent(inout) :: File
   type(dyn_export_t), intent(in)   :: dyn_out

   ! local variables
   integer(pio_offset_kind), parameter :: t_idx = 1

   type(element_t),  pointer :: elem(:)
   type(fvm_struct), pointer :: fvm(:)

   integer :: tl, tlqdp
   integer :: i, ie, ii, j, k, m
   integer :: ierr

   integer :: grid_id
   integer :: grid_dimlens(2)



   integer :: array_lens(3)
   integer :: file_lens(2)
   type(io_desc_t), pointer :: iodesc3d_fvm
   real(r8),    allocatable :: buf3d(:,:,:)



   character(len=*), parameter :: sub = 'write_restart_dynamics'
   !----------------------------------------------------------------------------

   call write_restart_hycoef(File)

   tl = timelevel%n0
   call TimeLevel_Qdp(timelevel, qsplit, tlQdp)

   if (iam .lt. par%nprocs) then
      elem => dyn_out%elem
      fvm => dyn_out%fvm
   else
      allocate (elem(0), fvm(0))
   endif

   ! write fields on GLL grid

   if (write_restart_unstruct) then
      call write_unstruct()
   else
      call write_elem()
   end if

   ! write CSLAM fields

   if (ntrac > 0) then

      grid_id = cam_grid_id('FVM')

      ! write coords for FVM grid
      call cam_grid_write_var(File, grid_id)

      call cam_grid_dimensions(grid_id, grid_dimlens)
      allocate(buf3d(nc*nc,nlev,nelemd))
      array_lens = (/nc*nc, nlev, nelemd/)
      file_lens  = (/grid_dimlens(1), nlev/)
      call cam_grid_get_decomp(grid_id, array_lens, file_lens, pio_double, iodesc3d_fvm)

      do ie = 1, nelemd
         do k = 1, nlev
            ii = 1
            do j = 1, nc
               do i = 1, nc
                  buf3d(ii,k,ie) = fvm(ie)%dp_fvm(i,j,k,n0_fvm)
                  ii = ii + 1
               end do
            end do
         end do
      end do
      call PIO_Setframe(file, dp_fvm_desc, t_idx)
      call PIO_Write_Darray(file, dp_fvm_desc, iodesc3d_fvm, buf3d, ierr)

      do m = 1, ntrac
         do ie = 1, nelemd
            do k = 1, nlev
               ii = 1
               do j = 1, nc
                  do i = 1, nc
                     buf3d(ii,k,ie) = fvm(ie)%c(i,j,k,m,n0_fvm)
                     ii = ii + 1
                  end do
               end do
            end do
         end do
         call PIO_Setframe(file, c_fvm_desc(m), t_idx)
         call PIO_Write_Darray(file, c_fvm_desc(m), iodesc3d_fvm, buf3d, ierr)
      end do

      deallocate(c_fvm_desc)
      deallocate(buf3d)
      ! should this call be made on a pointer?
      !call pio_freedecomp(File, iodesc3d_fvm)

   end if

   if (iam >= par%nprocs) then
      deallocate(elem, fvm)
   endif

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine write_elem()

   ! local variables
   integer          :: i, ie, j, k
   integer          :: ierr
   integer, pointer :: ldof(:)

   type(io_desc_t)  :: iodesc2d, iodesc3d

   real(kind=r8), pointer :: var3d(:,:,:,:), var2d(:,:,:)
   !----------------------------------------------------------------------------

   ldof => get_restart_decomp(elem, 1)
   call PIO_InitDecomp(pio_subsystem, pio_double, (/nelem_tot*np*np/), ldof, iodesc2d)
   deallocate(ldof)

   ldof => get_restart_decomp(elem, nlev)
   call PIO_InitDecomp(pio_subsystem, pio_double, (/nelem_tot*np*np,nlev/), ldof, iodesc3d)
   deallocate(ldof)

   allocate(var2d(np,np,nelemd))
   allocate(var3d(np,np,nelemd,nlev))

   !$omp parallel do num_threads(horz_num_threads) private(ie, j, i)
   do ie = 1, nelemd
      do j = 1, np
         do i = 1, np
            var2d(i,j,ie) = elem(ie)%state%psdry(i,j)
         end do
      end do
   end do
   call PIO_Setframe(File, psdry_desc, t_idx)
   call PIO_Write_Darray(File, psdry_desc, iodesc2d, var2d, ierr)

   !$omp parallel do num_threads(horz_num_threads) private(ie, k, j, i)
   do ie = 1, nelemd
      do k = 1, nlev
         do j = 1, np
            do i = 1, np
               var3d(i,j,ie,k) = elem(ie)%state%V(i,j,1,k,tl)
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Udesc, t_idx)
   call PIO_Write_Darray(File, Udesc, iodesc3d, var3d, ierr)

   !$omp parallel do num_threads(horz_num_threads) private(ie, k, j, i)
   do ie = 1, nelemd
      do k = 1, nlev
         do j = 1, np
            do i = 1, np
               var3d(i,j,ie,k) = elem(ie)%state%V(i,j,2,k,tl)
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Vdesc, t_idx)
   call PIO_Write_Darray(File, Vdesc, iodesc3d, var3d, ierr)

   !$omp parallel do num_threads(horz_num_threads) private(ie, k, j, i)
   do ie = 1, nelemd
      do k = 1, nlev
         do j = 1, np
            do i = 1, np
               var3d(i,j,ie,k) = elem(ie)%state%T(i,j,k,tl)
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Tdesc, t_idx)
   call PIO_Write_Darray(File, Tdesc, iodesc3d, var3d, ierr)

   do m = 1, qsize

      !$omp parallel do num_threads(horz_num_threads) private(ie, k, j, i)
      do ie = 1, nelemd
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  var3d(i,j,ie,k) = elem(ie)%state%Qdp(i,j,k,m,tlQdp)
               end do
            end do
         end do
      end do
      call PIO_Setframe(File, Qdesc_dp(m), t_idx)
      call PIO_Write_Darray(File, Qdesc_dp(m), iodesc3d, var3d, ierr)

   end do

   deallocate(var2d)
   deallocate(var3d)
   deallocate(qdesc_dp)

   call pio_freedecomp(File, iodesc2d)
   call pio_freedecomp(File, iodesc3d)

end subroutine write_elem

!-------------------------------------------------------------------------------

subroutine write_unstruct()

   ! local variables
   integer          :: i, ie, ii, j, k
   integer          :: ierr

   integer :: array_lens_3d(3), array_lens_2d(2)
   integer :: file_lens_2d(2), file_lens_1d(1)

   type(io_desc_t), pointer :: iodesc
   real(r8),    allocatable :: var2d(:,:), var3d(:,:,:)
   !----------------------------------------------------------------------------

   grid_id = cam_grid_id('GLL')

   ! write coordinate variables for unstructured GLL grid
   call cam_grid_write_var(File, grid_id)

   ! create map for distributed write
   call cam_grid_dimensions(grid_id, grid_dimlens)

   ! create map for distributed write of 2D fields
   array_lens_2d = (/npsq, nelemd/)
   file_lens_1d  = (/grid_dimlens(1)/)
   call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc)

   allocate(var2d(npsq,nelemd))

   do ie = 1, nelemd
      ii = 1
      do j = 1, np
         do i = 1, np
            var2d(ii,ie) = elem(ie)%state%psdry(i,j)
            ii = ii + 1
         end do
      end do
   end do
   call PIO_Setframe(File, psdry_desc, t_idx)
   call PIO_Write_Darray(File, psdry_desc, iodesc, var2d, ierr)

   nullify(iodesc)
   deallocate(var2d)

   ! create map for distributed write of 3D fields
   array_lens_3d = (/npsq, nlev, nelemd/)
   file_lens_2d  = (/grid_dimlens(1), nlev/)
   call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc)

   allocate(var3d(npsq,nlev,nelemd))

   do ie = 1, nelemd
      do k = 1, nlev
         ii = 1
         do j = 1, np
            do i = 1, np
               var3d(ii,k,ie) = elem(ie)%state%V(i,j,1,k,tl)
               ii = ii + 1
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Udesc, t_idx)
   call PIO_Write_Darray(File, Udesc, iodesc, var3d, ierr)

   do ie = 1, nelemd
      do k = 1, nlev
         ii = 1
         do j = 1, np
            do i = 1, np
               var3d(ii,k,ie) = elem(ie)%state%V(i,j,2,k,tl)
               ii = ii + 1
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Vdesc, t_idx)
   call PIO_Write_Darray(File, Vdesc, iodesc, var3d, ierr)

   do ie = 1, nelemd
      do k = 1, nlev
         ii = 1
         do j = 1, np
            do i = 1, np
               var3d(ii,k,ie) = elem(ie)%state%T(i,j,k,tl)
               ii = ii + 1
            end do
         end do
      end do
   end do
   call PIO_Setframe(File, Tdesc, t_idx)
   call PIO_Write_Darray(File, Tdesc, iodesc, var3d, ierr)

   do m = 1, qsize

      !$omp parallel do num_threads(horz_num_threads) private(ie, k, j, i)
      do ie = 1, nelemd
         do k = 1, nlev
            ii = 1
            do j = 1, np
               do i = 1, np
                  var3d(ii,k,ie) = elem(ie)%state%Qdp(i,j,k,m,tlQdp)
                  ii = ii + 1
               end do
            end do
         end do
      end do
      call PIO_Setframe(File, Qdesc_dp(m), t_idx)
      call PIO_Write_Darray(File, Qdesc_dp(m), iodesc, var3d, ierr)

   end do

   deallocate(var3d)
   deallocate(qdesc_dp)

end subroutine write_unstruct

!-------------------------------------------------------------------------------

end subroutine write_restart_dynamics

!=========================================================================================

subroutine read_restart_dynamics(File, dyn_in, dyn_out)

   ! arguments
   type(File_desc_t), intent(inout) :: File
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out

   ! local variables
   integer(pio_offset_kind), parameter :: t_idx = 1

   integer :: tl, tlQdp
   integer :: i, ie, ii, k, m, j
   integer :: ierr, err_handling
   integer :: fne, fnp, fnlev, fnc
   integer :: hdim_len, ncols_fvm

   integer :: nlev_dimid
   integer :: ncol_dimid
   integer :: ncol_fvm_dimid

   type(var_desc_t) :: udesc
   type(var_desc_t) :: vdesc
   type(var_desc_t) :: tdesc
   type(var_desc_t) :: psdry_desc
   type(var_desc_t), allocatable :: qdesc_dp(:)

   integer :: grid_id
   integer :: grid_dimlens(2)
   character(len=max_hcoordname_len) :: dimname1, dimname2

   real(r8),    allocatable :: var3d_fvm(:,:,:)

   logical :: readvar

   character(len=*), parameter :: sub = 'read_restart_dynamics'
   !----------------------------------------------------------------------------

   ! Note1: the hybrid coefficients are read from the same location as for an
   !        initial run (e.g., dyn_grid_init).

   ! Note2: the dyn_in and dyn_out objects are not associated with the elem and fvm
   !        objects until dyn_init is called.  Until the restart is better integrated
   !        into dyn_init we just access elem and fvm directly from the dyn_grid
   !        module.

   tl = timelevel%n0
   call TimeLevel_Qdp(timelevel, qsplit, tlQdp)
   call init_nelem_tot()

   call pio_seterrorhandling(File, pio_bcast_error, err_handling)

   ! some checks that the restart contains the same grid as the running model.

   ierr = PIO_Get_Att(File, PIO_GLOBAL, 'ne', fne)
   ierr = PIO_Get_Att(File, PIO_GLOBAL, 'np', fnp)
   if (ne /= fne .or. np /= fnp) then
      write(iulog,*) 'Restart file np or ne does not match model. np (file, model):', &
                     fnp, np, ' ne (file, model) ', fne, ne
      call endrun(sub//': Restart file np or ne does not match model.')
   end if

   ierr = PIO_Inq_DimID(File, 'lev', nlev_dimid)
   ierr = PIO_Inq_dimlen(File, nlev_dimid, fnlev)
   if (nlev /= fnlev) then
      write(iulog,*) 'Restart file nlev does not match model. nlev (file, namelist):', &
                     fnlev, nlev
      call endrun(sub//': Restart file nlev does not match model.')
   end if

   ! variable descriptors of required dynamics fields
   ierr = PIO_Inq_varid(File, 'U',     udesc)
   call cam_pio_handle_error(ierr, sub//': cannot find U')
   ierr = PIO_Inq_varid(File, 'V',     Vdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find V')
   ierr = PIO_Inq_varid(File, 'T',     tdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find T')
   ierr = PIO_Inq_varid(File, 'PSDRY', psdry_desc)
   call cam_pio_handle_error(ierr, sub//': cannot find PSDRY')
   allocate(qdesc_dp(qsize))
   do m = 1, qsize
      ierr = PIO_Inq_varid(File, "dp"//trim(cnst_name(m)), Qdesc_dp(m))
      call cam_pio_handle_error(ierr, sub//': cannot find dp'//trim(cnst_name(m)))
   end do

   ! check whether the restart fields on the GLL grid contain unique columns
   ! or the element structure (nenpnp = nelem_tot*np*np columns)

   ierr = PIO_Inq_DimID(File, 'nenpnp', ncol_dimid)
   if (ierr == pio_noerr) then

      call read_elem()

   else

      call read_unstruct()

   end if

   deallocate(qdesc_dp)

   ! recompute dp3d from psdry
   do ie = 1, nelemd
      do k = 1, nlev
         elem(ie)%state%dp3d(:,:,k,tl) = ((hyai(k+1) - hyai(k))*ps0) + &
                              ((hybi(k+1) - hybi(k))*elem(ie)%state%psdry(:,:))
      end do
   end do

   ! Seems like this initialization should be done somewhere else.
   do ie = 1, nelemd
      elem(ie)%derived%fM = 0._r8
      elem(ie)%derived%fT = 0._r8
      elem(ie)%derived%fQ = 0._r8
   end do

   ! read cslam fields

   if (ntrac > 0) then

      ! Checks that file and model dimensions agree.

      ierr = PIO_Get_Att(File, PIO_GLOBAL, 'nc', fnc)
      if (nc /= fnc) then
         write(iulog,*) 'Restart file nc does not match model. nc (file, model):',fnc,nc,&
             ' ne (file, model) ', fne, ne
         call endrun(sub//': Restart file nc does not match model.')
      end if

      ierr = PIO_Inq_DimID(File, 'ncol_fvm', ncol_fvm_dimid)
      call cam_pio_handle_error(ierr, sub//': cannot find ncol_fvm')
      ierr = PIO_Inq_dimlen(File, ncol_fvm_dimid, ncols_fvm)

      grid_id = cam_grid_id('FVM')
      call cam_grid_dimensions(grid_id, grid_dimlens)

      if (ncols_fvm /= grid_dimlens(1)) then
         write(iulog,*) 'Restart file ncol_fvm does not match model. ncols_fvm (file, model):',&
                        ncols_fvm, grid_dimlens(1)
         call endrun(sub//': Restart file ncols_fvm does not match model.')
      end if

      allocate(var3d_fvm(nc*nc,nlev,nelemd))
      var3d_fvm = 0._r8

      ! dp_fvm
      call infld('dp_fvm', file, 'ncol_fvm', 'lev', 1, nc**2, 1, nlev, 1, nelemd, &
                 var3d_fvm, readvar, gridname='FVM', timelevel=int(t_idx))
      do ie = 1, nelemd
         do k = 1, nlev
            ii = 1
            do j = 1, nc
               do i = 1, nc
                  fvm(ie)%dp_fvm(i,j,k,n0_fvm) = var3d_fvm(ii,k,ie)
                  ii = ii + 1
               end do
            end do
         end do
      end do

      ! tracers
      do m = 1, ntrac
         call infld(trim(cnst_name(m))//'_fvm', file, 'ncol_fvm', 'lev', &
                    1, nc**2, 1, nlev, 1, nelemd, var3d_fvm, readvar,    &
                    gridname='FVM', timelevel=int(t_idx))
         do ie = 1, nelemd
            do k = 1, nlev
               ii = 1
               do j = 1, nc
                  do i = 1, nc
                     fvm(ie)%c(i,j,k,m,n0_fvm) = var3d_fvm(ii,k,ie)
                     ii = ii + 1
                  end do
               end do
            end do
         end do
      end do

      ! compute dry surface pressure (a derived quantity)
      do ie = 1, nelemd
         do j = 1, nc
            do i = 1, nc
               fvm(ie)%psc(i,j) = sum(fvm(ie)%dp_fvm(i,j,:,n0_fvm)) +  ptop_ref
            end do
         end do
      end do

   end if

   call pio_seterrorhandling(File, err_handling)

   call dyn_init(dyn_in, dyn_out)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine read_elem()

   ! local variables
   integer :: ierr
   integer :: ncol
   integer :: i, ie, ii, j, k, m

   integer, pointer :: ldof(:)

   type(io_desc_t) :: iodesc2d, iodesc3d
   real(r8), allocatable :: var3d(:), var2d(:)

   character(len=*), parameter :: sub='read_elem'
   !----------------------------------------------------------------------------

   ierr = PIO_Inq_dimlen(File, ncol_dimid, ncol)
   call cam_pio_handle_error(ierr, sub//': reading nenpnp')

   ldof => get_restart_decomp(elem, 1)
   call PIO_InitDecomp(pio_subsystem, pio_double, (/ncol/), ldof, iodesc2d)
   deallocate(ldof)

   ldof => get_restart_decomp(elem, nlev)
   call PIO_InitDecomp(pio_subsystem, pio_double, (/ncol,nlev/), ldof, iodesc3d)
   deallocate(ldof)

   allocate(var3d(ncol*nlev), var2d(ncol))
   var2d = 0._r8
   var3d = 0._r8

   call pio_setframe(File, psdry_desc, t_idx)
   call pio_read_darray(File, psdry_desc, iodesc2d, var2d, ierr)
   call cam_pio_handle_error(ierr, sub//': reading PSDRY')
   ii = 0
   do ie = 1, nelemd
      do j = 1, np
         do i = 1, np
            ii = ii + 1
            elem(ie)%state%psdry(i,j) = var2d(ii)
         end do
      end do
   end do

   call pio_setframe(File, udesc, t_idx)
   call pio_read_darray(File, udesc, iodesc3d, var3d, ierr)
   call cam_pio_handle_error(ierr, sub//': reading U')
   ii = 0
   do k = 1, nlev
      do ie = 1, nelemd
         do j = 1, np
            do i = 1, np
               ii = ii + 1
               elem(ie)%state%v(i,j,1,k,tl) = var3d(ii)
            end do
         end do
      end do
   end do

   call pio_setframe(File, vdesc, t_idx)
   call pio_read_darray(File, vdesc, iodesc3d, var3d, ierr)
   call cam_pio_handle_error(ierr, sub//': reading V')
   ii = 0
   do k = 1, nlev
      do ie = 1, nelemd
         do j = 1, np
            do i = 1, np
               ii = ii + 1
               elem(ie)%state%v(i,j,2,k,tl) = var3d(ii)
            end do
         end do
      end do
   end do

   call pio_setframe(File, tdesc, t_idx)
   call pio_read_darray(File, tdesc, iodesc3d, var3d, ierr)
   call cam_pio_handle_error(ierr, sub//': reading T')
   ii = 0
   do k = 1, nlev
      do ie = 1, nelemd
         do j = 1, np
            do i = 1, np
               ii = ii + 1
               elem(ie)%state%T(i,j,k,tl) = var3d(ii)
            end do
         end do
      end do
   end do

   do m = 1, qsize
      call pio_setframe(File, qdesc_dp(m), t_idx)
      call pio_read_darray(File, qdesc_dp(m), iodesc3d, var3d, ierr)
      call cam_pio_handle_error(ierr, sub//': reading dp'//trim(cnst_name(m)))
      ii = 0
      do k = 1, nlev
         do ie = 1, nelemd
            do j = 1, np
               do i = 1, np
                  ii = ii + 1
                  elem(ie)%state%Qdp(i,j,k,m,tlQdp) = var3d(ii)
               end do
            end do
         end do
      end do
   end do

   deallocate(var3d, var2d)

end subroutine read_elem

!-------------------------------------------------------------------------------

subroutine read_unstruct()

   ! local variables
   integer :: grid_id

   integer :: i, ie, ii, j, kptr, m

   real(r8), allocatable :: dbuf2(:,:)         ! (npsq,nelemd)
   real(r8), allocatable :: dbuf3(:,:,:)       ! (npsq,nlev,nelemd)

   type(EdgeBuffer_t) :: edge

   character(len=*), parameter :: sub='read_unstruct'
   !----------------------------------------------------------------------------

   ! The name of the unstructured grid dimension is in the grid object
   ! since the coordinate date was written to the restart file using
   ! that object.
   grid_id = cam_grid_id('GLL')
   call cam_grid_get_dim_names(grid_id, dimname1, dimname2)

   allocate(dbuf2(npsq,nelemd))

   call read_2d('PSDRY', dbuf2)
   do ie = 1, nelemd
      ii = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%psdry(i,j) = dbuf2(ii,ie)
            ii = ii + 1
         end do
      end do
   end do

   deallocate(dbuf2)

   allocate(dbuf3(npsq,nlev,nelemd))

   call read_3d('U', dbuf3)
   do ie = 1, nelemd
      ii = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%v(i,j,1,:,tl) = dbuf3(ii,:,ie)
            ii = ii + 1
         end do
      end do
   end do

   call read_3d('V', dbuf3)
   do ie = 1, nelemd
      ii = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%v(i,j,2,:,tl) = dbuf3(ii,:,ie)
            ii = ii + 1
         end do
      end do
   end do

   call read_3d('T', dbuf3)
   do ie = 1, nelemd
      ii = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%T(i,j,:,tl) = dbuf3(ii,:,ie)
            ii = ii + 1
         end do
      end do
   end do

   do m = 1, qsize

      call read_3d('dp'//trim(cnst_name(m)), dbuf3)
      do ie = 1, nelemd
         ii = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%Qdp(i,j,:,m,tlQdp) = dbuf3(ii,:,ie)
               ii = ii + 1
            end do
         end do
      end do

   end do

   deallocate(dbuf3)

   ! boundary exchange
   if (iam < par%nprocs) then
      call initEdgeBuffer(par, edge, elem, (3+qsize)*nlev + 1 )
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVpack(edge, elem(ie)%state%psdry(:,:), 1, kptr, ie)
      kptr = kptr + 1
      call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,tl), 2*nlev, kptr, ie)
      kptr = kptr + (2 * nlev)
      call edgeVpack(edge, elem(ie)%state%T(:,:,:,tl), nlev, kptr, ie)
      kptr = kptr + nlev
      call edgeVpack(edge, elem(ie)%state%Qdp(:,:,:,:,tlQdp), nlev*qsize, kptr, ie)
   end do
   if (iam < par%nprocs) then
      call bndry_exchange(par, edge, location='read_restart_dynamics::read_ustruct')
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edge, elem(ie)%state%psdry(:,:), 1, kptr, ie)
      kptr = kptr + 1
      call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,tl), 2*nlev, kptr, ie)
      kptr = kptr + (2 * nlev)
      call edgeVunpack(edge, elem(ie)%state%T(:,:,:,tl), nlev, kptr, ie)
      kptr = kptr + nlev
      call edgeVunpack(edge, elem(ie)%state%Qdp(:,:,:,:,tlQdp), nlev*qsize, kptr, ie)
   end do

   if (iam < par%nprocs) then
      call FreeEdgeBuffer(edge)
   end if

end subroutine read_unstruct

!-------------------------------------------------------------------------------

subroutine read_2d(fieldname, buffer)

   character(len=*),  intent(in)    :: fieldname
   real(r8),          intent(inout) :: buffer(:, :)

   logical :: found
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), file, dimname1, 1, npsq, 1, nelemd, buffer, &
              found, gridname='GLL', timelevel=int(t_idx))
   if (.not. found) then
      call endrun('read_restart_dynamics: read_unstruct: read_2d: Could not find ' // &
                   trim(fieldname))
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   where (isnan(buffer)) buffer = 0.0_r8

end subroutine read_2d

!-------------------------------------------------------------------------------

subroutine read_3d(fieldname, buffer)

   character(len=*),  intent(in)    :: fieldname
   real(r8),          intent(inout) :: buffer(:,:,:)

   logical :: found
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), file, dimname1, 'lev', 1, npsq, 1, nlev, &
              1, nelemd, buffer, found, gridname='GLL', timelevel=int(t_idx))
   if (.not. found) then
      call endrun('read_restart_dynamics: read_unstruct: read_3d: Could not find ' // &
                   trim(fieldname))
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   where (isnan(buffer)) buffer = 0.0_r8

end subroutine read_3d

!-------------------------------------------------------------------------------
end subroutine read_restart_dynamics

!=========================================================================================
! Private
!=========================================================================================

function get_restart_decomp(elem, lev) result(ldof)

   ! Get the integer mapping of a variable in the dynamics decomp in memory.
   ! The canonical ordering is as on the file. A 0 value indicates that the
   ! variable is not on the file (eg halo or boundary values)

   type(element_t), intent(in) :: elem(:)
   integer,         intent(in) :: lev
   integer,         pointer    :: ldof(:)

   integer :: i, j, k, ie
   !----------------------------------------------------------------------------

   allocate(ldof(nelemd*np*np*lev))

   j = 1
   do k = 1, lev
      do ie = 1, nelemd
         do i = 1, np*np
            ldof(j) = (elem(ie)%GlobalID-1)*np*np + (k-1)*nelem_tot*np*np + i
            j = j + 1
         end do
      end do
   end do

end function get_restart_decomp

!=========================================================================================

function get_restart_decomp_fvm(elem, lev) result(ldof)

   type(element_t), intent(in) :: elem(:)
   integer,         intent(in) :: lev
   integer,         pointer    :: ldof(:)

   integer :: i, j, k, ie
   !----------------------------------------------------------------------------

   allocate(ldof(nelemd*nc*nc*lev))

   j = 1
   do k = 1, lev
      do ie = 1, nelemd
         do i = 1, nc*nc
            ldof(j) = (elem(ie)%GlobalID-1)*nc*nc + (k-1)*nelem_tot*nc*nc + i
            j = j + 1
         end do
      end do
   end do

end function get_restart_decomp_fvm

!=========================================================================================

end module restart_dynamics

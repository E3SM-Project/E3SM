module prep_aoflux_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: CXX => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_xao, num_inst_frc, num_inst_ocn
  use seq_comm_mct,     only: CPLID, logunit
  use seq_comm_mct,     only : mbofxid ! iMOAB id for mpas ocean migrated mesh to coupler pes, just for xao flux calculations
#ifdef MOABDEBUG
  use seq_comm_mct,     only : mbox2id ! used only for debugging ocn and mct
#endif
  use seq_comm_mct,     only : mbaxid ! iMOAB app id for atm on cpl pes
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_comm_mct, only : num_moab_exports
  use seq_infodata_mod, only: seq_infodata_getdata, seq_infodata_type
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: atm, ocn

  ! use m_List               ,only: mct_list_nitem         => nitem
  ! use mct_mod ! for mct_list_nitem 

  use iso_c_binding


  implicit none
  private ! except
  save

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_aoflux_init

  public :: prep_aoflux_calc_xao_ox
  public :: prep_aoflux_calc_xao_ax

  public :: prep_aoflux_get_xao_ox
  public :: prep_aoflux_get_xao_ax

  ! these are to expose the artificial arrays created for setting moab tag
  ! these are the transpose of the AVs for fluxes;
  public :: prep_aoflux_get_xao_omct
  public :: prep_aoflux_get_xao_amct

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! attribute vectors
  type(mct_aVect), pointer :: xao_ox(:)   ! Atm-ocn fluxes, ocn grid, cpl pes
  type(mct_aVect), pointer :: xao_ax(:)   ! Atm-ocn fluxes, atm grid, cpl pes

  ! allocate xao_omct, but use lsize_o, size of the local mct ocn gsmap (and AVs)
  real(r8) ,  private, pointer :: xao_omct(:,:) ! atm-ocn fluxes, ocn grid, mct local sizes 
  real(r8) ,  private, pointer :: xao_amct(:,:) ! atm-ocn fluxes, atm grid, mct local sizes 

  ! seq_comm_getData variables
  logical :: iamroot_CPLID                ! .true. => CPLID masterproc
  integer :: mpicom_CPLID                 ! MPI cpl communicator

  ! seq_infodata_getData variables
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_aoflux_init (infodata)

    !---------------------------------------------------------------
    ! Description
    ! Initialize atm/ocn flux component and compute ocean albedos
    ! module variables
    !
    use iMOAB, only : iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage, iMOAB_GetMeshInfo
#ifdef MOABDEBUG
    use iMOAB, only : iMOAB_WriteMesh ! for writing debug file
#endif
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    !
    ! Local Variables
    integer                     :: exi,ierr
    integer                     :: lsize_o
    integer                     :: lsize_a
    integer                     :: tagtype, numco, tagindex, ent_type
    character(CXX)              :: tagname
    integer                     :: size_list ! for number of tags 
    real(r8),    allocatable    :: tagValues(:) ! used for setting some default tags
    integer                     :: arrSize ! for the size of tagValues
    type(mct_list)              :: temp_list  ! used to count the number of strings / fields 
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! used to find out the size of local mesh, in moab
#ifdef MOABDEBUG
    character*100 outfile, wopts ! for writing debug file
#endif
    character(CS)      :: aoflux_grid ! grid for atm ocn flux calc
    type(mct_avect) , pointer   :: a2x_ax
    type(mct_avect) , pointer   :: o2x_ox
    character(*)    , parameter :: subname = '(prep_aoflux_init)'

    !---------------------------------------------------------------

    call seq_infodata_getdata(infodata,  &
         aoflux_grid=aoflux_grid)

    call seq_comm_getdata(CPLID, &
         mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

    a2x_ax => component_get_c2x_cx(atm(1))
    if (associated(a2x_ax)) then
       lsize_a = mct_aVect_lsize(a2x_ax)
    else
       lsize_a = 0
    end if

    o2x_ox => component_get_c2x_cx(ocn(1))
    if (associated(o2x_ox)) then
       lsize_o = mct_aVect_lsize(o2x_ox)
    else
       lsize_o = 0
    end if

    allocate(xao_ax(num_inst_xao))
    do exi = 1,num_inst_xao
       call mct_aVect_init(xao_ax(exi), rList=seq_flds_xao_fields, lsize=lsize_a)
       call mct_aVect_zero(xao_ax(exi))
    end do
    allocate(xao_ox(num_inst_xao))
    do exi = 1,num_inst_xao
       call mct_aVect_init(xao_ox(exi), rList=seq_flds_xao_fields, lsize=lsize_o)
       call mct_aVect_zero(xao_ox(exi))
    enddo

! define flux tags on the moab ocean mesh, second copy of ocean mesh on coupler
    if (mbofxid .ge. 0 ) then ! //
      !add the normalization tag
       tagname = trim(seq_flds_xao_fields)//":norm8wt"//C_NULL_CHAR
       tagtype = 1 ! dense, double
       numco = 1
       ierr = iMOAB_DefineTagStorage(mbofxid, tagname, tagtype, numco, tagindex )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' error in defining tags on ocn phys mesh on cpl '
          call shr_sys_abort(subname//' ERROR in defining tags on ocn phys mesh on cpl')
       endif

       ! make it zero
       ! first form a list and get size.
       call mct_list_init(temp_list ,seq_flds_xao_fields)
       size_list=mct_list_nitem (temp_list) + 1 ! 1 more for the normalization tag
       call mct_list_clean(temp_list)
       ! find out the number of local elements in moab mesh
       ierr  = iMOAB_GetMeshInfo ( mbofxid, nvert, nvise, nbl, nsurf, nvisBC ); ! could be different of lsize_o
      ! local size of vertices is different from lsize_o
      ! nvsise(1) is the number of primary elements locally 
       arrSize = nvise(1) * size_list ! there are size_list tags that need to be zeroed out
       allocate(tagValues(arrSize) )
       ent_type = 1 ! cell type
       tagValues = 0._r8
       ierr = iMOAB_SetDoubleTagStorage ( mbofxid, tagname, arrSize , ent_type, tagValues(1))
       deallocate(tagValues)
       if (ierr .ne. 0) then
         write(logunit,*) subname,' error in zeroing out xao_fields  '
         call shr_sys_abort(subname//' ERROR in zeroing out xao_fields in init ')
       endif

       allocate(xao_omct(lsize_o, size_list)) ! the transpose of xao_ox(size_list, lsize_o) 
       xao_omct = 0._r8
#ifdef MOABDEBUG
       ! create for debugging the tags on mbox2id (mct grid on coupler)
       ierr = iMOAB_DefineTagStorage(mbox2id, tagname, tagtype, numco, tagindex )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' error in defining tags on ocn mct mesh on cpl '
          call shr_sys_abort(subname//' ERROR in defining tags on ocn mct mesh on cpl')
       endif
       ent_type = 0 ! cell type, this is point cloud mct
       arrSize = lsize_o * size_list
       ierr = iMOAB_SetDoubleTagStorage ( mbox2id, tagname, arrSize , ent_type, xao_omct )
       if (ierr .ne. 0) then
         write(logunit,*) subname,' error in zeroing out xao_fields on mct instance ocn '
         call shr_sys_abort(subname//' ERROR in zeroing out xao_fields on mct instance ocn ')
       endif
       !deallocate(xao_omct)
        ! debug out file
      outfile = 'o_flux.h5m'//C_NULL_CHAR
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(mbofxid, outfile, wopts)
 
      if (ierr .ne. 0) then
         write(logunit,*) subname,' error in writing o_flux mesh '
         call shr_sys_abort(subname//' ERROR in writing o_flux mesh ')
      endif
       ! debug out file
      outfile = 'ox_mct.h5m'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(mbox2id, outfile, wopts)
 
      if (ierr .ne. 0) then
         write(logunit,*) subname,' error in writing ox_mct mesh with 0 values '
         call shr_sys_abort(subname//' ERROR in writing ox_mct mesh ')
      endif
#endif
    endif

! define atm-ocn flux tags on the moab atm mesh
    if (mbaxid .ge. 0 ) then ! //
       tagname = trim(seq_flds_xao_fields)//C_NULL_CHAR
       tagtype = 1 ! dense, double
       numco = 1
       ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco, tagindex )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' error in defining tags on atm phys mesh on cpl '
          call shr_sys_abort(subname//' ERROR in defining tags on atm phys mesh on cpl')
       endif
       ! make it zero
       ! first form a list 
       call mct_list_init(temp_list ,seq_flds_xao_fields)
       size_list=mct_list_nitem (temp_list)
       call mct_list_clean(temp_list)
       ! find out the number of local elements in moab mesh
       ierr  = iMOAB_GetMeshInfo ( mbaxid, nvert, nvise, nbl, nsurf, nvisBC ); ! could be different of lsize_o
      ! local size of vertices is different from lsize_o
       arrSize = nvise(1) * size_list ! there are size_list tags that need to be zeroed out
       allocate(tagValues(arrSize) )
       ent_type = 1 ! cell type now, not a point cloud anymore
       tagValues = 0._r8
       ierr = iMOAB_SetDoubleTagStorage ( mbaxid, tagname, arrSize , ent_type, tagValues)
       deallocate(tagValues)
       if (ierr .ne. 0) then
         write(logunit,*) subname,' error in zeroing out xao_fields  '
         call shr_sys_abort(subname//' ERROR in zeroing out xao_fields in init ')
       endif
       allocate(xao_amct(lsize_a, size_list)) ! the transpose of xao_ax(size_list, lsize_a) 
    endif

  end subroutine prep_aoflux_init

  !================================================================================================

  subroutine prep_aoflux_calc_xao_ax(fractions_ox, flds, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create xao_ox
    !
    ! Uses
    use prep_atm_mod, only: prep_atm_get_mapper_So2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fo2a
    use prep_atm_mod, only: prep_atm_get_mapper_Sof2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fof2a
#ifdef MOABDEBUG
    use iMOAB, only :  iMOAB_WriteMesh
#endif
    !
    ! Arguments
    type(mct_aVect) , intent(in)    :: fractions_ox(:)
    character(len=*), intent(in)    :: flds
    character(len=*), intent(in)    :: timer
    !
    ! Local Variables
    type(seq_map)   , pointer :: mapper_So2a
    type(seq_map)   , pointer :: mapper_Fo2a
    type(seq_map)   , pointer :: mapper_Sof2a
    type(seq_map)   , pointer :: mapper_Fof2a
    integer :: exi, efi
    character(*), parameter :: subname = '(prep_aoflux_calc_xao_ax)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
#ifdef MOABDEBUG
    character*50             :: outfile, wopts, lnum
    integer :: ierr
#endif
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (trim(flds) == 'albedos') then
       do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1

          mapper_Sof2a => prep_atm_get_mapper_Sof2a()
          call seq_map_map(mapper_Sof2a, xao_ox(exi), xao_ax(exi), &
               fldlist=seq_flds_xao_albedo, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       enddo
    end if

    if (trim(flds) == 'states_and_fluxes') then
       do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1

          mapper_Sof2a => prep_atm_get_mapper_Sof2a()
          call seq_map_map(mapper_Sof2a, xao_ox(exi), xao_ax(exi), &
               fldlist=seq_flds_xao_states, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')

          mapper_Fof2a => prep_atm_get_mapper_Fof2a()
          call seq_map_map(mapper_Fof2a, xao_ox(exi), xao_ax(exi),&
               fldlist=seq_flds_xao_fluxes, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       enddo
    end if
#ifdef MOABDEBUG
! albedos is called second so wait until then to write
    if (trim(flds) == 'albedos') then
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
         write(lnum,"(I0.2)")num_moab_exports
      if(mbaxid > 0 ) then
            ! projections on atm
         outfile = 'FlxAlb2Atm'//trim(lnum)//'.h5m'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing ocean to atm projection'
            call shr_sys_abort(subname//' ERROR in writing ocean to atm projection')
         endif
      endif
      if(mbofxid > 0) then
         outfile = 'FlxAlb2Ocn'//trim(lnum)//'.h5m'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(mbofxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing ocean to atm projection'
            call shr_sys_abort(subname//' ERROR in writing ocean to atm projection')
         endif
      endif
    end if
#endif
    call t_drvstopf  (trim(timer))

  end subroutine prep_aoflux_calc_xao_ax

  !================================================================================================

  subroutine prep_aoflux_calc_xao_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create xao_ox
    !
    ! Uses
    use prep_ocn_mod, only: prep_ocn_get_mapper_Fa2o
    !
    ! Arguments
    character(len=*), intent(in)    :: timer
    !
    ! Local Variables
    type(seq_map), pointer :: mapper_Fa2o
    integer :: exi
    character(*), parameter :: subname = '(prep_aoflux_calc_xao_ox)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    ! this mapping has to be done with area overlap mapping for all fields
    ! due to the masking of the xao_ax data and the fact that a2oS is bilinear

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do exi = 1,num_inst_xao
       !       if (iamroot_CPLID .and. exi == 1) then
       !          write(logunit,F00) 'Calling map_atm2ocn_mct for mapping xao_ax to xao_ox'
       !       end if

       mapper_Fa2o => prep_ocn_get_mapper_Fa2o()
       call seq_map_map(mapper_Fa2o, xao_ax(exi), xao_ox(exi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_aoflux_calc_xao_ox

  !================================================================================================

  function prep_aoflux_get_xao_ox()
    type(mct_aVect), pointer :: prep_aoflux_get_xao_ox(:)
    prep_aoflux_get_xao_ox => xao_ox(:)
  end function prep_aoflux_get_xao_ox

  function prep_aoflux_get_xao_ax()
    type(mct_aVect), pointer :: prep_aoflux_get_xao_ax(:)
    prep_aoflux_get_xao_ax => xao_ax(:)
  end function prep_aoflux_get_xao_ax

  function prep_aoflux_get_xao_omct()
    real(r8), pointer :: prep_aoflux_get_xao_omct(:,:)
    prep_aoflux_get_xao_omct => xao_omct
  end function prep_aoflux_get_xao_omct

  function prep_aoflux_get_xao_amct()
    real(r8), pointer :: prep_aoflux_get_xao_amct(:,:)
    prep_aoflux_get_xao_amct => xao_amct
  end function prep_aoflux_get_xao_amct

end module prep_aoflux_mod

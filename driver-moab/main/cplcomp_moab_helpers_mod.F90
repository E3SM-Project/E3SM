module cplcomp_moab_helpers_mod

  use shr_kind_mod, only: CXX => SHR_KIND_CXX
  use shr_sys_mod, only: shr_sys_abort
  use iso_c_binding

  implicit none
  private

  public :: moab_register_app
  public :: moab_define_global_id_tag
  public :: moab_define_double_tag
  public :: moab_define_tags        ! iMOAB_DefineTagStorage on a coupler app, for every component instance
  public :: moab_set_tags_all_inst  ! iMOAB_SetDoubleTagStorage of the same values for every component instance
  public :: moab_load_mesh
  public :: moab_send_mesh
  public :: moab_receive_mesh
  public :: moab_free_sender_buffers

contains

  subroutine moab_register_app(appname, mpicom, id_join, appid, subctx)
    use iMOAB, only: iMOAB_RegisterApplication

    character(len=*), intent(in) :: appname
    integer, intent(in) :: mpicom
    integer, intent(in) :: id_join
    integer, intent(out) :: appid
    character(len=*), intent(in) :: subctx

    integer :: ierr

    ierr = iMOAB_RegisterApplication(trim(appname)//C_NULL_CHAR, mpicom, id_join, appid)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR cannot register app '//trim(appname))
    end if
  end subroutine moab_register_app

  subroutine moab_define_global_id_tag(appid, subctx)
    use iMOAB, only: iMOAB_DefineTagStorage

    integer, intent(in) :: appid
    character(len=*), intent(in) :: subctx

    integer :: ierr
    integer :: tagindex
    integer :: tagtype
    integer :: numco
    character(CXX) :: tagname

    tagtype = 0
    numco = 1
    tagname = 'GLOBAL_ID'//C_NULL_CHAR

    ierr = iMOAB_DefineTagStorage(appid, tagname, tagtype, numco, tagindex)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR adding GLOBAL_ID tag')
    end if
  end subroutine moab_define_global_id_tag

  subroutine moab_define_double_tag(appid, tagname_in, subctx)

    integer, intent(in) :: appid
    character(len=*), intent(in) :: tagname_in
    character(len=*), intent(in) :: subctx

    integer :: ierr
    integer :: tagindex
    integer :: tagtype
    integer :: numco
    character(CXX) :: tagname

    tagtype = 1
    numco = 1
    tagname = trim(tagname_in)//C_NULL_CHAR

    ierr = moab_define_tags(appid, tagname, tagtype, numco, tagindex)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR defining tag '//trim(tagname_in))
    end if
  end subroutine moab_define_double_tag

  !---------------------------------------------------------------------------
  ! Same arguments and return code as iMOAB_DefineTagStorage. Use it for tags on
  ! coupler applications: besides the given (instance 1) tags it defines the
  ! tags of instances 2..num_inst_max (see seq_comm_mbinst_tags), so that
  ! Get/SetDoubleTagStorage can address the data of every instance on the one
  ! coupler mesh. With a single instance it is identical to iMOAB_DefineTagStorage.
  ! tagindex returns the index of the first instance 1 tag.
  !---------------------------------------------------------------------------
  integer function moab_define_tags(appid, tagname, tagtype, numco, tagindex) result(ierr)
    use iMOAB, only: iMOAB_DefineTagStorage
    use seq_comm_mct, only: num_inst_max, seq_comm_mbinst_tags

    integer, intent(in) :: appid
    character(len=*), intent(in) :: tagname   ! colon-separated names, C_NULL_CHAR terminated
    integer, intent(in) :: tagtype, numco
    integer, intent(out) :: tagindex

    integer :: k, tagindex_k
    character(CXX) :: tagname_k

    ierr = iMOAB_DefineTagStorage(appid, tagname, tagtype, numco, tagindex)
    if (ierr /= 0) return
    do k = 2, num_inst_max
       tagname_k = seq_comm_mbinst_tags(tagname, k)
       ierr = iMOAB_DefineTagStorage(appid, tagname_k, tagtype, numco, tagindex_k)
       if (ierr /= 0) return
    enddo
  end function moab_define_tags

  !---------------------------------------------------------------------------
  ! Same arguments and return code as iMOAB_SetDoubleTagStorage. Sets the same
  ! values on the given (instance 1) tags and on the tags of instances
  ! 2..num_inst_max; used for initial values of coupler fields.
  !---------------------------------------------------------------------------
  integer function moab_set_tags_all_inst(appid, tagname, arrsize, ent_type, values) result(ierr)
    use iMOAB, only: iMOAB_SetDoubleTagStorage
    use seq_comm_mct, only: num_inst_max, seq_comm_mbinst_tags
    use shr_kind_mod, only: R8 => SHR_KIND_R8

    integer, intent(in) :: appid
    character(len=*), intent(in) :: tagname   ! colon-separated names, C_NULL_CHAR terminated
    integer, intent(in) :: arrsize, ent_type
    real(R8), intent(in) :: values(*)

    integer :: k
    character(CXX) :: tagname_k

    do k = 1, num_inst_max
       tagname_k = seq_comm_mbinst_tags(tagname, k)
       ierr = iMOAB_SetDoubleTagStorage(appid, tagname_k, arrsize, ent_type, values)
       if (ierr /= 0) return
    enddo
  end function moab_set_tags_all_inst

  subroutine moab_load_mesh(appid, infile, ropts, nghlay, subctx)
    use iMOAB, only: iMOAB_LoadMesh

    integer, intent(in) :: appid
    character(len=*), intent(in) :: infile
    character(len=*), intent(in) :: ropts
    integer, intent(in) :: nghlay
    character(len=*), intent(in) :: subctx

    integer :: ierr

    ierr = iMOAB_LoadMesh(appid, trim(infile)//C_NULL_CHAR, trim(ropts)//C_NULL_CHAR, nghlay)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR loading mesh from '//trim(infile))
    end if
  end subroutine moab_load_mesh

  subroutine moab_send_mesh(appid, mpicom_join, mpigrp_cplid, id_join, partmethod, subctx)
    use iMOAB, only: iMOAB_SendMesh

    integer, intent(in) :: appid
    integer, intent(in) :: mpicom_join
    integer, intent(in) :: mpigrp_cplid
    integer, intent(in) :: id_join
    integer, intent(in) :: partmethod
    character(len=*), intent(in) :: subctx

    integer :: ierr

    ierr = iMOAB_SendMesh(appid, mpicom_join, mpigrp_cplid, id_join, partmethod)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR sending mesh')
    end if
  end subroutine moab_send_mesh

  subroutine moab_receive_mesh(appid, mpicom_join, mpigrp_old, id_old, subctx)
    use iMOAB, only: iMOAB_ReceiveMesh

    integer, intent(in) :: appid
    integer, intent(in) :: mpicom_join
    integer, intent(in) :: mpigrp_old
    integer, intent(in) :: id_old
    character(len=*), intent(in) :: subctx

    integer :: ierr

    ierr = iMOAB_ReceiveMesh(appid, mpicom_join, mpigrp_old, id_old)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR receiving mesh')
    end if
  end subroutine moab_receive_mesh

  subroutine moab_free_sender_buffers(appid, context_id, subctx)
    use iMOAB, only: iMOAB_FreeSenderBuffers

    integer, intent(in) :: appid
    integer, intent(in) :: context_id
    character(len=*), intent(in) :: subctx

    integer :: ierr

    ierr = iMOAB_FreeSenderBuffers(appid, context_id)
    if (ierr /= 0) then
      call shr_sys_abort(trim(subctx)//' ERROR freeing sender buffers')
    end if
  end subroutine moab_free_sender_buffers

end module cplcomp_moab_helpers_mod

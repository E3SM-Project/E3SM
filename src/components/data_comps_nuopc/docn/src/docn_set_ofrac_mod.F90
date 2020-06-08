module docn_set_ofrac_mod

  use ESMF
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_pio_mod      , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_methods_mod , only : dshr_fldbun_getfldptr,  chkerr
  use pio

  implicit none
  private ! except

  public :: docn_set_ofrac

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine docn_set_ofrac(exportState, model_mesh, model_meshfile, model_maskfile, &
       model_nxg, model_nyg, compid, ocn_fraction, rc)

    ! -------------------------------------
    ! Determine ocean fraction
    ! -------------------------------------

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(ESMF_Mesh)        , intent(in)    :: model_mesh
    character(len=*)       , intent(in)    :: model_meshfile 
    character(len=*)       , intent(in)    :: model_maskfile 
    integer                , intent(in)    :: model_nxg
    integer                , intent(in)    :: model_nyg
    integer                , intent(in)    :: compid 
    real(r8)               , intent(out)   :: ocn_fraction(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer             :: numOwnedElements ! number of elements owned by this PET
    type(ESMF_DistGrid) :: distGrid         ! mesh distGrid
    type(ESMF_Array)    :: elemMaskArray
    integer, pointer    :: imask(:)
    type(file_desc_t)   :: pioid
    type(var_desc_t)    :: varid
    type(io_desc_t)     :: pio_iodesc
    integer             :: io_type          ! pio info
    integer             :: io_format        ! pio info
    integer             :: rcode
    integer, pointer    :: model_gindex(:)  ! model global index spzce
    type(iosystem_desc_t), pointer :: pio_subsystem => null() ! pio info
    integer :: n
    character(len=*), parameter :: subname='(docn_set_ofrac): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    
    call ESMF_MeshGet(model_mesh, numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(model_maskfile) /= trim(model_meshfile)) then

       ! Read in the ocean fraction from the input namelist ocean mask file and assume 'frac' name on domain file
       allocate(model_gindex(numOwnedElements))
       call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=model_gindex, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       pio_subsystem => shr_pio_getiosys(compid)
       io_type       =  shr_pio_getiotype(compid)
       io_format     =  shr_pio_getioformat(compid)
       rcode = pio_openfile(pio_subsystem, pioid, io_type, trim(model_maskfile), pio_nowrite)
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid, 'frac', varid)
       if ( rcode /= PIO_NOERR ) then
          call shr_sys_abort(' ERROR: variable frac not found in file '//trim(model_maskfile))
       end if
       call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
       call pio_initdecomp(pio_subsystem, pio_double, (/model_nxg, model_nyg/), model_gindex, pio_iodesc)
       call pio_read_darray(pioid, varid, pio_iodesc, ocn_fraction, rcode)
       call pio_closefile(pioid)
       call pio_freedecomp(pio_subsystem, pio_iodesc)
       deallocate(model_gindex)

    else

       ! Obtain the ocean fraction from the mask values in the ocean mesh file
       allocate(imask(numOwnedElements))
       elemMaskArray = ESMF_ArrayCreate(distGrid, imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! the following call sets the varues of imask
       call ESMF_MeshGet(model_mesh, elemMaskArray=elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! now set the fraction as just the real mask
       ocn_fraction(:) = real(imask(:), kind=r8)
       call ESMF_ArrayDestroy(elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(imask)

    end if

  end subroutine docn_set_ofrac

end module docn_set_ofrac_mod

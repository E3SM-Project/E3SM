module GetGlobalValuesMod

  !-----------------------------------------------------------------------
  ! Obtain and Write Global Index information
  !-----------------------------------------------------------------------
  implicit none
  private

  ! PUBLIC MEMBER FUNCTIONS:

  public :: GetGlobalIndex
  public :: GetGlobalWrite
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  integer function GetGlobalIndex(decomp_index, clmlevel)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given clmlevel
    !
    ! Uses:
    use clmtype    , only: nameg, namel, namec, namep
    use decompMod  , only: bounds_type, get_clmlevel_gsmap, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varctl , only: iulog
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use mct_mod
    !
    ! Arguments 
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    type(bounds_type)             :: bounds_proc   ! processor bounds
    type(mct_gsMap),pointer       :: gsmap         ! global seg map
    integer, pointer,dimension(:) :: gsmap_ordered ! gsmap ordered points
    integer                       :: beg_index     ! beginning proc index for clmlevel
    !----------------------------------------------------------------

    call get_proc_bounds(bounds_proc)

    if (trim(clmlevel) == nameg) then
       beg_index = bounds_proc%begg
    else if (trim(clmlevel) == namel) then
       beg_index = bounds_proc%begl
    else if (trim(clmlevel) == namec) then
       beg_index = bounds_proc%begc
    else if (trim(clmlevel) == namep) then
       beg_index = bounds_proc%begp
    else
       call shr_sys_abort('clmlevel of '//trim(clmlevel)//' not supported' // &
            errmsg(__FILE__, __LINE__))
    end if

    call get_clmlevel_gsmap(clmlevel=trim(clmlevel), gsmap=gsmap)
    call mct_gsmap_op(gsmap, iam, gsmap_ordered)
    GetGlobalIndex = gsmap_ordered(decomp_index - beg_index + 1)

  end function GetGlobalIndex

  !-----------------------------------------------------------------------
  subroutine GetGlobalWrite(decomp_index, clmlevel)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write global index information for input local indices
    !
    use clmtype           
    use shr_sys_mod , only: shr_sys_flush
    use shr_sys_mod , only: shr_sys_abort
    use shr_log_mod , only: errMsg => shr_log_errMsg
    use clm_varctl  , only: iulog
    !
    ! Arguments:
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    integer :: igrc, ilun, icol, ipft 
    !-----------------------------------------------------------------------

    if (trim(clmlevel) == nameg) then

       igrc = decomp_index
       write(iulog,*)'local  gridcell index = ',igrc
       write(iulog,*)'global gridcell index = ',GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,*)'gridcell longitude    = ',grc%londeg(igrc)
       write(iulog,*)'gridcell latitude     = ',grc%latdeg(igrc)

    else if (trim(clmlevel) == namel) then

       ilun = decomp_index
       igrc = lun%gridcell(ilun)
       write(iulog,*)'local  landunit index = ',ilun
       write(iulog,*)'global landunit index = ',GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,*)'global gridcell index = ',GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,*)'gridcell longitude    = ',grc%londeg(igrc)
       write(iulog,*)'gridcell latitude     = ',grc%latdeg(igrc)
       write(iulog,*)'landunit type         = ',lun%itype(decomp_index)

    else if (trim(clmlevel) == namec) then

       icol = decomp_index
       ilun = col%landunit(icol)
       igrc = col%gridcell(icol)
       write(iulog,*)'local  column   index = ',icol
       write(iulog,*)'global column   index = ',GetGlobalIndex(decomp_index=icol, clmlevel=namec)
       write(iulog,*)'global landunit index = ',GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,*)'global gridcell index = ',GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,*)'gridcell longitude    = ',grc%londeg(igrc)
       write(iulog,*)'gridcell latitude     = ',grc%latdeg(igrc)
       write(iulog,*)'column   type         = ',col%itype(icol)
       write(iulog,*)'landunit type         = ',lun%itype(ilun)
   
    else if (trim(clmlevel) == namep) then

       ipft = decomp_index
       icol = pft%column(ipft)
       ilun = pft%landunit(ipft)
       igrc = pft%gridcell(ipft)
       write(iulog,*)'local  pft      index = ',ipft
       write(iulog,*)'global pft      index = ',GetGlobalIndex(decomp_index=ipft, clmlevel=namep)
       write(iulog,*)'global column   index = ',GetGlobalIndex(decomp_index=icol, clmlevel=namec)
       write(iulog,*)'global landunit index = ',GetGlobalIndex(decomp_index=ilun, clmlevel=namel)
       write(iulog,*)'global gridcell index = ',GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,*)'gridcell longitude    = ',grc%londeg(igrc)
       write(iulog,*)'gridcell latitude     = ',grc%latdeg(igrc)
       write(iulog,*)'pft      type         = ',pft%itype(ipft)
       write(iulog,*)'column   type         = ',col%itype(icol)
       write(iulog,*)'landunit type         = ',lun%itype(ilun)

    else		       
       call shr_sys_abort('clmlevel '//trim(clmlevel)//'not supported '//errmsg(__FILE__, __LINE__))

    end if

    call shr_sys_flush(iulog)

  end subroutine GetGlobalWrite

end module GetGlobalValuesMod

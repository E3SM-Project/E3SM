!===============================================================================
! SVN $Id: shr_map_mod.F90 35318 2012-03-08 23:40:50Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_130528/shr/shr_map_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_map_mod -- generic map data type and associated methods
!
! !DESCRIPTION:
!    Generic map data type and associated methods
!    \newline
!    This module supports mapping of fields from one grid to another.
!    A general datatype, shr\_map\_mapType, stores the mapping information
!    set in shr\_map\_mapSet.  shr\_map\_mapData then allows this mapping
!    to be applied to an input array to generate the output array.
!    \newline
!    The mapType has several flags that give the user various options
!    for setting the mapping
!      type: [remap,fill]
!         remap - mapping of data between different grids, primarily
!            for the active grid area
!         fill - mapping of data on the same grid, primarily to fill missing
!            areas, copy data, or set the array to a spval.
!      algo: [copy,bilinear,nn,nnoni,nnonj,spval]
!         copy - copy data from one array to another using indexing
!         bilinear - bilinear remapping using 4 corner points
!         nn - nearest neighbor, set value to nn value
!         nnoni - nearest neighbor using i, search for nearest neighbor in the
!            i direction first, then j
!         nnonj - nearest neighbor using j, search for nearest neighbor in the
!            j direction first, then i
!         spval - set values to the spval
!      mask: [srcmask,dstmask,nomask,bothmask]
!         srcmask - use only src points with mask = true in mapping
!         dstmask - map only to dst points where mask = true
!         nomask - ignore both src and dst mask in mapping
!         bothmask - use both src and dst mask in mapping (srcmask and dstmask)
!      vect: [scalar,vector]
!         scalar - fields are scalar type (default)
!         vector - fields are vector type, operates only on 2 fields to 2 fields
!    NOTE: Not all combinatations are unique and not all combinations are valid
!    \newline
!    The above settings are put into the maptype using shr\_map\_put.  Public
!    parameters are available to users to set the switches.  The first three
!    switches must be set then the mapSet method can be called.  After the
!    mapSet method is called, the mapData method can be used.
!    \newline
!    A Note on Subroutine Arguments:
!    Lat, lon, and mask arguments in these routines are 2d (nx,ny) 
!    Array arguments are 2d (nf,nxy), number of fields by grid point
!    \newline
!    General Usage:
!       type(shr\_map\_mapType) :: mymap
!       call shr\_map\_put(mymap,'type','remap')
!       call shr\_map\_put(mymap,shr\_map\_fs\_algo,shr\_map\_fs\_bilinear)
!       call shr\_map\_put(mymap,shr\_map\_fs\_mask,'bothmask')
!       call shr\_map\_put(mymap,shr\_map\_fs\_vect,'scalar')
!       call shr\_map\_mapSet(mymap,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,rc=rCode)
!       call shr\_map\_mapData(Asrc,Adst,mymap)
!       \newline
!       call shr\_map\_mapSet(mymap,Xs,Ys,Ms,Xd,Yd,Md,name='fillnnoni',type='fill',algo='nnoni',mask='dstmask',rc=rc)
!       call shr\_map\_mapData(Asrc,Adst,mymap)
!       \newline
!       call shr\_map\_mapData(Ain,Aout,Xs,Ys,Ms,Xd,Yd,Md,type='remap',algo='nn',mask='dstmask',rc)
!
! !REMARKS:
!     nn needs a faster algorithm
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

module shr_map_mod

! !USES:

  use shr_const_mod
  use shr_kind_mod
  use shr_sys_mod
  use shr_timer_mod
  use shr_log_mod, only: s_loglev  => shr_log_Level
  use shr_log_mod, only: s_logunit => shr_log_Unit

  implicit none
  private

! !PUBLIC TYPES:

  public :: shr_map_maptype           ! shr_map datatype

  type shr_map_mapType                       ! like mct sparsematrix datatype
    private
    character(SHR_KIND_CS)       :: name 
    character(SHR_KIND_CS)       :: type 
    character(SHR_KIND_CS)       :: algo
    character(SHR_KIND_CS)       :: mask
    character(SHR_KIND_CS)       :: vect
    integer(SHR_KIND_IN)         :: nsrc     ! grid size or src
    integer(SHR_KIND_IN)         :: ndst     ! grid size of dst
    integer(SHR_KIND_IN)         :: nwts     ! number of total weights
    real(SHR_KIND_R8)   ,pointer :: xsrc(:)  ! longitude, for vector, rad
    real(SHR_KIND_R8)   ,pointer :: ysrc(:)  ! latitude , for vector, rad
    real(SHR_KIND_R8)   ,pointer :: xdst(:)  ! longitude, for vector, rad
    real(SHR_KIND_R8)   ,pointer :: ydst(:)  ! latitude , for vector, rad
    real(SHR_KIND_R8)   ,pointer :: wgts(:)  ! weights
    integer(SHR_KIND_IN),pointer :: isrc(:)  ! input grid index
    integer(SHR_KIND_IN),pointer :: idst(:)  ! output grid index
    character(SHR_KIND_CS)       :: fill     ! string to check if filled
    character(SHR_KIND_CS)       :: init     ! initialization of dst array
  end type shr_map_mapType

! PUBLIC MEMBER FUNCTIONS:

  public :: shr_map_checkInit         ! check whether map type is set
  public :: shr_map_checkFilled       ! check whether map wts are set
  public :: shr_map_put               ! put stuff into the datatype
  public :: shr_map_get               ! get stuff out of the datatype
  public :: shr_map_mapSet            ! compute weights in map
  public :: shr_map_mapData           ! map data
  public :: shr_map_listValidOpts     ! list valid options
  public :: shr_map_print             ! print map datatype info
  public :: shr_map_clean             ! clean map datatype
  public :: shr_map_setAbort          ! set abort flag for shr_map
  public :: shr_map_setDebug          ! set debug level for shr_map
  public :: shr_map_setDopole         ! set dopole flag

! PUBLIC DATA MEMBERS:

  !--- Field Strings (fldStr) ---
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_name  = 'name'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_type  = 'type'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_algo  = 'algo'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_mask  = 'mask'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_vect  = 'vect'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_nwts  = 'nwts'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_nsrc  = 'nsrc'
  character(SHR_KIND_CS),public,parameter :: shr_map_fs_ndst  = 'ndst'

  !--- "type" options ---
  character(len=*),public,parameter :: shr_map_fs_fill     = 'fill    '
  character(len=*),public,parameter :: shr_map_fs_cfill    = 'cfill   '
  character(len=*),public,parameter :: shr_map_fs_remap    = 'remap   '

  !--- "algorithm" options ---
  character(len=*),public,parameter :: shr_map_fs_copy     = 'copy    '
  character(len=*),public,parameter :: shr_map_fs_bilinear = 'bilinear'
  character(len=*),public,parameter :: shr_map_fs_nn       = 'nn      '
  character(len=*),public,parameter :: shr_map_fs_nnoni    = 'nnoni   '
  character(len=*),public,parameter :: shr_map_fs_nnonj    = 'nnonj   '
  character(len=*),public,parameter :: shr_map_fs_spval    = 'spval   '

  !--- "mask" options ---
  character(len=*),public,parameter :: shr_map_fs_srcmask  = 'srcmask '
  character(len=*),public,parameter :: shr_map_fs_dstmask  = 'dstmask '
  character(len=*),public,parameter :: shr_map_fs_nomask   = 'nomask  '
  character(len=*),public,parameter :: shr_map_fs_bothmask = 'bothmask'

  !--- "vect" options ---
  character(len=*),public,parameter :: shr_map_fs_scalar   = 'scalar  '
  character(len=*),public,parameter :: shr_map_fs_vector   = 'vector  '

  !--- other public parameters ---
  character(SHR_KIND_CS),public,parameter :: shr_map_setTru = 'TRUE map'
  character(SHR_KIND_CS),public,parameter :: shr_map_setFal = 'FALSE m '
  integer(SHR_KIND_IN)  ,public,parameter :: shr_map_ispval = -99
  real(SHR_KIND_R8)     ,public,parameter :: shr_map_spval  = shr_const_spval

!EOP

  !--- Must update these if anything above changes ---
  integer(SHR_KIND_IN),public,parameter :: shr_map_fs_ntype = 3
  character(len=*),public,parameter :: &
  shr_map_fs_types(shr_map_fs_ntype) = (/shr_map_fs_fill, &
                                         shr_map_fs_cfill, &
                                         shr_map_fs_remap /)

  integer(SHR_KIND_IN),public,parameter :: shr_map_fs_nalgo = 6
  character(len=*),public,parameter :: &
  shr_map_fs_algos(shr_map_fs_nalgo) = (/shr_map_fs_copy, &
                                         shr_map_fs_bilinear, &
                                         shr_map_fs_nn, &
                                         shr_map_fs_nnoni, &
                                         shr_map_fs_nnonj, &
                                         shr_map_fs_spval  /)

  integer(SHR_KIND_IN),public,parameter :: shr_map_fs_nmask = 4
  character(len=*),public,parameter :: &
  shr_map_fs_masks(shr_map_fs_nmask) = (/shr_map_fs_srcmask, &
                                         shr_map_fs_dstmask, &
                                         shr_map_fs_nomask , &
                                         shr_map_fs_bothmask /)

  integer(SHR_KIND_IN),public,parameter :: shr_map_fs_nvect = 2
  character(len=*),public,parameter :: &
  shr_map_fs_vects(shr_map_fs_nvect) = (/shr_map_fs_scalar, &
                                         shr_map_fs_vector /)

  interface shr_map_put ; module procedure &
    shr_map_putCS, &
    shr_map_putR8, &
    shr_map_putIN
  end interface

  interface shr_map_get ; module procedure &
    shr_map_getCS, &
    shr_map_getR8, &
    shr_map_getIN, &
    shr_map_getAR
  end interface

  interface shr_map_mapSet ; module procedure &
    shr_map_mapSet_global, &
    shr_map_mapSet_dest
  end interface

  interface shr_map_mapData ; module procedure &
    shr_map_mapDatam, &
    shr_map_mapDatanm
  end interface

  logical,save                   :: doabort = .true.
  logical,save                   :: dopole  = .true.   ! for bilinear
  integer(SHR_KIND_IN),save      :: debug = 0
  character(SHR_KIND_CS),parameter :: fillstring = 'mapisfilled'
  character(SHR_KIND_CS),parameter :: inispval   = 'spval'
  character(SHR_KIND_CS),parameter :: initcopy   = 'copy'
  real(SHR_KIND_R8)   ,parameter :: c0  = 0._SHR_KIND_R8
  real(SHR_KIND_R8)   ,parameter :: c1  = 1._SHR_KIND_R8
  real(SHR_KIND_R8)   ,parameter :: c2  = 2._SHR_KIND_R8
  real(SHR_KIND_R8)   ,parameter :: eps = 1.0e-12_SHR_KIND_R8
  real(SHR_KIND_R8)   ,parameter :: pi  = shr_const_pi

!===============================================================================
contains
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_checkInit -- returns init state of map
!
! !DESCRIPTION:
!     Returns init state of map.  shr\_map\_checkInit is true
!     if the type, algo, and mask are set to valid values.
!     \newline
!     test = shr\_map\_checkInit(map)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_map_checkInit(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType),intent(in) :: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkInit') "

!-------------------------------------------------------------------------------

  if (shr_map_checkFldStrOpt(shr_map_fs_type,map%type) .and. &
      shr_map_checkFldStrOpt(shr_map_fs_algo,map%algo) .and. &
      shr_map_checkFldStrOpt(shr_map_fs_mask,map%mask)) then
    shr_map_checkInit = .true.
  else
    shr_map_checkInit = .false.
  endif

end function shr_map_checkInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_checkFilled -- returns fill state of map
!
! !DESCRIPTION:
!     Returns fill state of map.  shr\_map\_checkFilled is true
!     if the number of weights are greater than zero in map
!     and if the wgts, isrc, and idst arrays have been allocated to
!     that size.
!     \newline
!     test = shr\_map\_checkFilled(map)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_map_checkFilled(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType),intent(in) :: map

!EOP

  !--- local ---
  integer(SHR_KIND_IN) :: nwts

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkFilled') "

!-------------------------------------------------------------------------------

  shr_map_checkFilled = .false.

  nwts = map%nwts
  if (map%fill == fillstring .and. nwts >= 0) then
     if (size(map%wgts) == nwts .and. size(map%isrc) == nwts &
                                .and. size(map%idst) == nwts ) then
        shr_map_checkFilled = .true.
     endif
  endif

end function shr_map_checkFilled

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_checkFldStr -- checks fldstr for validity
!
! !DESCRIPTION:
!     Returns true if fldstr is valid (ie. 'type','algo','mask')
!     \newline
!     test = shr\_map\_checkFldStr('type')
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_map_checkFldStr(fldStr)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*)   :: fldStr

!XXEOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkFldStr') "

!-------------------------------------------------------------------------------

  shr_map_checkFldStr = .false.

  if     (trim(fldStr) == trim(shr_map_fs_type).or. &
          trim(fldStr) == trim(shr_map_fs_name).or. &
          trim(fldStr) == trim(shr_map_fs_algo).or. &
          trim(fldStr) == trim(shr_map_fs_mask).or. &
          trim(fldStr) == trim(shr_map_fs_vect).or. &
          trim(fldStr) == trim(shr_map_fs_nsrc).or. &
          trim(fldStr) == trim(shr_map_fs_ndst).or. &
          trim(fldStr) == trim(shr_map_fs_nwts)) then
     shr_map_checkFldStr = .true.
  endif

end function shr_map_checkFldStr

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_checkFldStrOpt -- checks cval for validity with fldstr
!
! !DESCRIPTION:
!     Returns true if cval is valid for fldstr (ie. 'type,remap','algo,bilinear',
!     'mask,srcmask')
!     \newline
!     test = shr\_map\_checkFldStrOpt(shr_map_fs_algo,'bilinear')
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_map_checkFldStrOpt(fldStr,cval)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),intent(in) :: fldStr
  character(*),intent(in) :: cval

!XXEOP

  !--- local ---
  integer(SHR_KIND_IN)   :: n

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkFldStrOpt') "

!-------------------------------------------------------------------------------

  shr_map_checkFldStrOpt = .false.

  if (.not.shr_map_checkFldStr(fldStr)) return

  if (trim(fldStr) == trim(shr_map_fs_name)) then
       shr_map_checkFldStrOpt = .true.
  elseif (trim(fldStr) == trim(shr_map_fs_type)) then
    do n = 1,shr_map_fs_ntype
      if (trim(cval) == trim(shr_map_fs_types(n))) shr_map_checkFldStrOpt = .true.
    enddo
  elseif (trim(fldStr) == trim(shr_map_fs_algo)) then
    do n = 1,shr_map_fs_nalgo
      if (trim(cval) == trim(shr_map_fs_algos(n))) shr_map_checkFldStrOpt = .true.
    enddo
  elseif (trim(fldStr) == trim(shr_map_fs_mask)) then
    do n = 1,shr_map_fs_nmask
      if (trim(cval) == trim(shr_map_fs_masks(n))) shr_map_checkFldStrOpt = .true.
    enddo
  elseif (trim(fldStr) == trim(shr_map_fs_vect)) then
    do n = 1,shr_map_fs_nvect
      if (trim(cval) == trim(shr_map_fs_vects(n))) shr_map_checkFldStrOpt = .true.
    enddo
  endif

end function shr_map_checkFldStrOpt

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_getCS -- get string from map
!
! !DESCRIPTION:
!     one of the shr\_map\_get methods for chars
!     returns value cval for input fldstr in map
!     \newline
!     call shr\_map\_get(mymap,shr\_map\_fs\_type,cval)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_getCS(map,fldStr,cval)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  character(*)          ,intent(out):: cval

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_getCS') "

!-------------------------------------------------------------------------------

  cval = shr_map_setFal
  if (.not.shr_map_checkFldStr(fldStr)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
    return
  endif

  if     (trim(fldStr) == trim(shr_map_fs_name)) then
    cval = map%name
  elseif (trim(fldStr) == trim(shr_map_fs_type)) then
    cval = map%type
  elseif (trim(fldStr) == trim(shr_map_fs_algo)) then
    cval = map%algo
  elseif (trim(fldStr) == trim(shr_map_fs_mask)) then
    cval = map%mask
  elseif (trim(fldStr) == trim(shr_map_fs_vect)) then
    cval = map%vect
  else
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
  endif

end subroutine shr_map_getCS

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_getIN -- get integer from map
!
! !DESCRIPTION:
!     one of the shr\_map\_get methods for integers
!     returns value ival for input fldstr in map
!     \newline
!     call shr\_map\_get(mymap,shr\_map\_fs\_nwts,ival)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_getIN(map,fldStr,ival)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  integer(SHR_KIND_IN)  ,intent(out):: ival

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_getIN') "

!-------------------------------------------------------------------------------

  ival = shr_map_ispval
  if (.not.shr_map_checkFldStr(fldStr)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
    return
  endif

  if     (trim(fldStr) == trim(shr_map_fs_nwts)) then
    ival = map%nwts
  elseif (trim(fldStr) == trim(shr_map_fs_nsrc)) then
    ival = map%nsrc
  elseif (trim(fldStr) == trim(shr_map_fs_ndst)) then
    ival = map%ndst
  else
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
  endif

end subroutine shr_map_getIN

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_getR8 -- get real from map
!
! !DESCRIPTION:
!     one of the shr\_map\_get methods for reals
!     returns value rval for input fldstr in map
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_getR8(map,fldStr,rval)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  real(SHR_KIND_R8)     ,intent(out):: rval

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_getR8') "

!-------------------------------------------------------------------------------

  rval = shr_map_spval
  if (.not.shr_map_checkFldStr(fldStr)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
    return
  endif

  call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))

end subroutine shr_map_getR8

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_getAR -- get arrays from map
!
! !DESCRIPTION:
!     one of the shr\_map\_get methods for arrays
!     returns value ival for input fldstr in map
!     \newline
!     call shr\_map\_get(mymap,idst,isrc,wgts)
!
! !REVISION HISTORY:
!     2009-Jul-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_getAR(map,isrc,idst,wgts)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(in) :: map
  integer(SHR_KIND_IN),pointer,optional :: isrc(:)
  integer(SHR_KIND_IN),pointer,optional :: idst(:)
  real   (SHR_KIND_R8),pointer,optional :: wgts(:)

!EOP

  !--- local ---
  integer(SHR_KIND_IN) :: nwts

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_getAR') "

!-------------------------------------------------------------------------------

  nwts = map%nwts

  if (present(isrc)) then
     if (size(isrc) < nwts) then
        call shr_sys_abort(subName//' ERROR is isrc size')
     endif
     isrc(1:nwts) = map%isrc(1:nwts)
  endif

  if (present(idst)) then
     if (size(idst) < nwts) then
        call shr_sys_abort(subName//' ERROR is idst size')
     endif
     idst(1:nwts) = map%idst(1:nwts)
  endif

  if (present(wgts)) then
     if (size(wgts) < nwts) then
        call shr_sys_abort(subName//' ERROR is wgts size')
     endif
     wgts(1:nwts) = map%wgts(1:nwts)
  endif

end subroutine shr_map_getAR

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_putCS -- put char to map
!
! !DESCRIPTION:
!     one of the shr\_map\_put methods for chars
!     puts value cval for input fldstr in map
!     verify is optional argument that check validity and will
!     call abort if cval is not valid option for fldstr.
!     \newline
!     call shr\_map\_put(mymap,shr\_map\_fs\_algo,shr\_map\_fs\_bilinear)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_putCS(map,fldStr,cval,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  character(*)          ,intent(in) :: cval
  logical,optional      ,intent(in) :: verify     ! check if string is valid

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_putCS') "

!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify = verify
  if (lverify .and. .not.shr_map_checkFldStrOpt(fldStr,cval)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr)//' '//trim(cval))
    return
  endif

  if (trim(fldStr) == trim(shr_map_fs_name)) then
    map%name = cval
  elseif (trim(fldStr) == trim(shr_map_fs_type)) then
    map%type = cval
  elseif (trim(fldStr) == trim(shr_map_fs_algo)) then
    map%algo = cval
  elseif (trim(fldStr) == trim(shr_map_fs_mask)) then
    map%mask = cval
  elseif (trim(fldStr) == trim(shr_map_fs_vect)) then
    map%vect = cval
  else
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
  endif

end subroutine shr_map_putCS

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_putIN -- put integer to map
!
! !DESCRIPTION:
!     one of the shr\_map\_put methods for integers
!     puts value ival for input fldstr in map
!     verify is optional argument that check validity and will
!     call abort if ival is not valid option for fldstr.
!     \newline
!     call shr\_map\_put(mymap,shr\_map\_fs\_nwts,-1)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_putIN(map,fldStr,ival,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  integer(SHR_KIND_IN)  ,intent(in) :: ival
  logical,optional      ,intent(in) :: verify     ! check if string is valid

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_putIN') "
  character(*),parameter :: F01     = "('(shr_map_putIN) ',a,i8) "

!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify = verify
  if (lverify .and. .not.shr_map_checkFldStr(fldStr)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
    return
  endif

  if     (trim(fldStr) == trim(shr_map_fs_nwts)) then
    map%nwts = ival
  elseif (trim(fldStr) == trim(shr_map_fs_nsrc)) then
    map%nsrc = ival
  elseif (trim(fldStr) == trim(shr_map_fs_ndst)) then
    map%ndst = ival
  else
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
  endif

end subroutine shr_map_putIN

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_putR8 -- put real to map
!
! !DESCRIPTION:
!     one of the shr\_map\_put methods for reals
!     puts value rval for input fldstr in map
!     verify is optional argument that check validity and will
!     call abort if rval is not valid option for fldstr.
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_putR8(map,fldStr,rval,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  real(SHR_KIND_R8)     ,intent(in) :: rval
  logical,optional      ,intent(in) :: verify     ! check if string is valid

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_putR8') "

!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify = verify
  if (lverify .and. .not.shr_map_checkFldStr(fldStr)) then
    call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))
    return
  endif

  call shr_map_abort(subName//' ERROR illegal fldStr '//trim(fldStr))

end subroutine shr_map_putR8

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_print -- write map to stdout
!
! !DESCRIPTION:
!     Write map info to stdout
!     \newline
!     call shr\_map\_print(mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_print(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(in) :: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_print') "
  character(*),parameter :: F00     = "('(shr_map_print) ',a) "
  character(*),parameter :: F01     = "('(shr_map_print) ',a,2l2) "
  character(*),parameter :: F02     = "('(shr_map_print) ',a,i8) "
  character(*),parameter :: F03     = "('(shr_map_print) ',a,3i8) "
  character(*),parameter :: F04     = "('(shr_map_print) ',a,2i8) "
  character(*),parameter :: F05     = "('(shr_map_print) ',a,2e20.13) "

  if (s_loglev > 0) then
     write(s_logunit,*) ' '
     write(s_logunit,F01) '   name : '//trim(map%name),shr_map_checkInit(map),shr_map_checkFilled(map)
     write(s_logunit,F00) '   type : '//trim(map%type)
     write(s_logunit,F00) '   algo : '//trim(map%algo)
     write(s_logunit,F00) '   mask : '//trim(map%mask)
     write(s_logunit,F00) '   vect : '//trim(map%vect)
     write(s_logunit,F04) '   gsiz : ',map%nsrc,map%ndst
     write(s_logunit,F05) '   xsrc : ',minval(map%xsrc),maxval(map%xsrc)
     write(s_logunit,F05) '   ysrc : ',minval(map%ysrc),maxval(map%ysrc)
     write(s_logunit,F05) '   xdst : ',minval(map%xdst),maxval(map%xdst)
     write(s_logunit,F05) '   ydst : ',minval(map%ydst),maxval(map%ydst)
     write(s_logunit,F02) '   nwts : ',map%nwts
     write(s_logunit,F03) '   wsiz : ',size(map%wgts),size(map%isrc),size(map%idst)
     write(s_logunit,F00) '   init : '//trim(map%init)

     call shr_sys_flush(s_logunit)
  endif

end subroutine shr_map_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_listValidOpts -- list the valid switches for map
!
! !DESCRIPTION:
!     Lists the valid switches for map, informational only
!     \newline
!     call shr\_map\_listValidOpts()
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_listValidOpts()

  implicit none

! !INPUT/OUTPUT PARAMETERS:

!EOP

  !--- local ---
  integer(SHR_KIND_IN) :: n

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_listValidOpts') "
  character(*),parameter :: F00     = "('(shr_map_listValidOpts) ',a) "

!-------------------------------------------------------------------------------

  if (s_loglev > 0) then
     write(s_logunit,F00) ':'
     write(s_logunit,F00) '  '//trim(shr_map_fs_name)//' : any character string'
     do n = 1,shr_map_fs_ntype
        write(s_logunit,F00) '  '//trim(shr_map_fs_type)//' : '//trim(shr_map_fs_types(n))
     enddo
     do n = 1,shr_map_fs_nalgo
        write(s_logunit,F00) '  '//trim(shr_map_fs_algo)//' : '//trim(shr_map_fs_algos(n))
     enddo
     do n = 1,shr_map_fs_nmask
        write(s_logunit,F00) '  '//trim(shr_map_fs_mask)//' : '//trim(shr_map_fs_masks(n))
     enddo
     do n = 1,shr_map_fs_nvect
        write(s_logunit,F00) '  '//trim(shr_map_fs_vect)//' : '//trim(shr_map_fs_vects(n))
     enddo
     call shr_sys_flush(s_logunit)
  endif

end subroutine shr_map_listValidOpts

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_clean -- cleans map
!
! !DESCRIPTION:
!     Cleans map by resetting switches, deallocating arrays
!     \newline
!     call shr\_map\_clean(mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_clean(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map

!EOP

  !--- local ---
  integer :: rc

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_clean') "
  character(*),parameter :: F00     = "('(shr_map_clean) ',a) "

!-------------------------------------------------------------------------------

  map%fill = ' '
  map%init = ' '
  call shr_map_put(map,shr_map_fs_name,shr_map_setFal,verify=.false.)
  call shr_map_put(map,shr_map_fs_type,shr_map_setFal,verify=.false.)
  call shr_map_put(map,shr_map_fs_algo,shr_map_setFal,verify=.false.)
  call shr_map_put(map,shr_map_fs_mask,shr_map_setFal,verify=.false.)
  call shr_map_put(map,shr_map_fs_mask,shr_map_setFal,verify=.false.)
  call shr_map_put(map,shr_map_fs_nwts,shr_map_ispval)
  call shr_map_put(map,shr_map_fs_nsrc,shr_map_ispval)
  call shr_map_put(map,shr_map_fs_ndst,shr_map_ispval)
  deallocate(map%xsrc,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map wgts'
  deallocate(map%ysrc,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map wgts'
  deallocate(map%xdst,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map wgts'
  deallocate(map%ydst,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map wgts'
  deallocate(map%wgts,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map wgts'
  deallocate(map%isrc,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map isrc'
  deallocate(map%idst,stat=rc)
  if (rc > 0.and.debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'Warning: unable to deallocate map idst'

end subroutine shr_map_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_mapSet_global -- Compute mapping weights
!
! !DESCRIPTION:
!     Compute mapping weights based on setting in map.  Fill the
!     weights in the map.  Currently supported maps and action:
!        fill :copy  = copy array by index, mask switch used
!        fill :spval = copy array, fill with spval, mask switch not used
!        fill :nn*   = copy array, fill with nnval, mask switch not used
!        remap:copy  = copy array by index, mask switch used
!        remap:spval = sets array to spval, mask switch used
!        remap:bil*  = bilinear interpolation, mask switch used
!        remap:nn*   = sets array to nnval, mask switch used
!     \newline
!     Requirements for input grids:
!       Xsrc,Ysrc must be regular lat/lon grid, monotonically increasing,
!          can be degrees or radians
!       Xdst,Ydst are arbitrary list of lats/lons, must be same units as src
!       Msrc,Mdst have nonzero for active grid point, zero for non-active
!       src and dst must be the grid for type = fill
!     Grids are check for validity
!     \newline
!     call shr\_map\_mapSet(mymap,Xs,Ys,Ms,Xd,Yd,Md)
!     \newline
!     call shr\_map\_mapSet(mymap,Xs,Ys,Ms,Xd,Yd,Md,algo='bilinear')
!
! !REMARKS
!     If bothmask or srcmask is used with remap and some algorithms, active
!     dst grid points can have invalid values.  A report is produced after
!     weights are calculated and this information will be detailed.
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_mapSet_global(map,Xsrc,Ysrc,Msrc,Xdst_in,Ydst,Mdst,name,type,algo,mask,vect,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map      ! map
  real(SHR_KIND_R8)     ,intent(in) :: Xsrc(:,:)  ! lon of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Ysrc(:,:)  ! lat of src grid
  integer(SHR_KIND_IN)  ,intent(in) :: Msrc(:,:)  ! mask of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Xdst_in(:,:)  ! lon of dst grid
  real(SHR_KIND_R8)     ,intent(in) :: Ydst(:,:)  ! lat of dst grid
  integer(SHR_KIND_IN)  ,intent(in) :: Mdst(:,:)  ! mask of dst grid
  character(*) ,optional,intent(in) :: name       ! name
  character(*) ,optional,intent(in) :: type       ! type
  character(*) ,optional,intent(in) :: algo       ! algo
  character(*) ,optional,intent(in) :: mask       ! mask
  character(*) ,optional,intent(in) :: vect       ! vect
  integer(SHR_KIND_IN),optional,intent(out) :: rc ! error code

!EOP

  !--- local ---
  integer(SHR_KIND_IN)   :: nis,njs,nid,njd
  integer(SHR_KIND_IN)   :: nwts,n,n1,n2,ncnt,i,j,inn,jnn
  integer(SHR_KIND_IN)   :: irc,lrc
  real(SHR_KIND_R8) :: rmin,rmax                     ! min/max value
  real(SHR_KIND_R8) :: cang                          ! circle angle, deg or rad
  real(SHR_KIND_R8),allocatable    :: Xdst(:,:)      ! lon of dst grid, wrapped as needed

  integer(SHR_KIND_IN)             :: pmax           ! max num of wgts in pti...
  integer(SHR_KIND_IN)             :: ptot,ptot2     ! max num of wgts in lis...
  integer(SHR_KIND_IN)             :: pnum           ! num of wgts set in pti...
  integer(SHR_KIND_IN),allocatable :: pti(:)         ! i index for wgts
  integer(SHR_KIND_IN),allocatable :: ptj(:)         ! j index for wgts
  real(SHR_KIND_R8)   ,allocatable :: ptw(:)         ! weights for pti,ptj

  integer(SHR_KIND_IN),allocatable :: lis(:)         ! tmp src/dst index
  integer(SHR_KIND_IN),allocatable :: lid(:)         ! tmp src/dst index
  real(SHR_KIND_R8)   ,allocatable :: lwt(:)         ! tmp wgt array
  real(SHR_KIND_R8)   ,allocatable :: sum(:)         ! tmp sum array
  integer(SHR_KIND_IN),allocatable :: ltmp(:)        ! tmp src/dst index, for resize
  real(SHR_KIND_R8)   ,allocatable :: lwtmp(:)       ! tmp wgt array, for resize

  character(len=8) :: units    ! radians or degrees

  logical :: masksrc       ! local var to turn on masking using src mask
  logical :: maskdst       ! local var to turn on masking using dst mask
  logical :: maskdstbysrc  ! local var to turn on masking using src mask for
                           ! dst array, especially for fill
  logical :: renorm        ! local var to turn on renormalization

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_mapSet_global') "
  character(*),parameter :: F00     = "('(shr_map_mapSet_global) ',a) "
  character(*),parameter :: F01     = "('(shr_map_mapSet_global) ',a,l2) "
  character(*),parameter :: F02     = "('(shr_map_mapSet_global) ',a,2i8) "
  character(*),parameter :: F03     = "('(shr_map_mapSet_global) ',a,2e20.13) "

!-------------------------------------------------------------------------------

  lrc = 0
  if (present(rc)) rc = lrc

  if (present(name)) call shr_map_put(map,shr_map_fs_name,name)
  if (present(type)) call shr_map_put(map,shr_map_fs_type,type,verify=.true.)
  if (present(algo)) call shr_map_put(map,shr_map_fs_algo,algo,verify=.true.)
  if (present(mask)) call shr_map_put(map,shr_map_fs_mask,mask,verify=.true.)
  if (present(vect)) call shr_map_put(map,shr_map_fs_vect,vect,verify=.true.)
  map%init = inispval

  if (.NOT.shr_map_checkInit(map)) then
    call shr_map_abort(subName//' ERROR map not initialized')
  endif

  !--- is lat/lon degrees or radians? ---
  cang = 360._SHR_KIND_R8
  units = 'degrees'
  if (shr_map_checkRad(Ysrc)) then
    cang=c2*pi
    units = 'radians'
  endif

  nis = size(Xsrc,1)
  njs = size(Xsrc,2)
  nid = size(Xdst_in,1)
  njd = size(Xdst_in,2)

  !--- shift Xdst by 2pi to range of Xsrc as needed ---
  allocate(Xdst(nid,njd))
  rmin = minval(Xsrc)
  rmax = maxval(Xsrc)
  do j=1,njd
  do i=1,nid
    Xdst(i,j) = Xdst_in(i,j)
    do while ((Xdst(i,j) < rmin .and. Xdst(i,j)+cang <= rmax).or.  &
              (Xdst(i,j) > rmax .and. Xdst(i,j)-cang >= rmin))
      if (Xdst(i,j) < rmin) then
        Xdst(i,j) = Xdst(i,j) + cang
      elseif (Xdst(i,j) > rmax) then
        Xdst(i,j) = Xdst(i,j) - cang
      else
        call shr_sys_abort(subName//' ERROR in Xdst wrap')
      endif
    enddo
  enddo
  enddo
  
  call shr_map_checkGrids_global(Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,map,lrc)

  map%nwts = 0
  map%nsrc = nis*njs
  map%ndst = nid*njd

! deallocate(map%xsrc,stat=irc) ! this used to be a safe way to delloc when necessary,
! deallocate(map%ysrc,stat=irc) ! but do nothing when pointers were undefined or
! deallocate(map%xdst,stat=irc) ! un-associated, in Oct 2005, undefined ptrs started
! deallocate(map%ydst,stat=irc) ! causing seg-faults on bluesky (B. Kauffman)
  allocate(map%xsrc(nis*njs))
  allocate(map%ysrc(nis*njs))
  allocate(map%xdst(nid*njd))
  allocate(map%ydst(nid*njd))
  do j=1,njs
  do i=1,nis
    call shr_map_2dto1d(n1,nis,njs,i,j)
    map%xsrc(n1) = Xsrc(i,j)*c2*pi/cang
    map%ysrc(n1) = Ysrc(i,j)*c2*pi/cang
  enddo
  enddo
  do j=1,njd
  do i=1,nid
    call shr_map_2dto1d(n1,nid,njd,i,j)
    map%xdst(n1) = Xdst(i,j)*c2*pi/cang
    map%ydst(n1) = Ydst(i,j)*c2*pi/cang
  enddo
  enddo

  masksrc = .false.
  maskdstbysrc = .false.
  maskdst = .false.
  renorm  = .true.

  if (trim(map%type) /= trim(shr_map_fs_fill) .and. &
      trim(map%type) /= trim(shr_map_fs_cfill)) then
    if (trim(map%mask) == trim(shr_map_fs_bothmask) .or. &
        trim(map%mask) == trim(shr_map_fs_srcmask)) masksrc = .true.
    if (trim(map%mask) == trim(shr_map_fs_bothmask) .or. &
        trim(map%mask) == trim(shr_map_fs_dstmask)) maskdst = .true.
  endif
  if (trim(map%algo) == trim(shr_map_fs_spval)) then
    masksrc = .false.
    renorm = .false.
  endif

  if (debug > 1) then
    if (s_loglev > 0) write(s_logunit,*) ' '
    call shr_map_print(map)
  endif

  if (lrc /= 0) then
    if (present(rc)) rc = lrc
    return
  endif

  if (trim(map%algo) == trim(shr_map_fs_bilinear)) then
    if (dopole) then
      pmax = nis+2   ! possible for high lat points
      ptot = 4*nid*njd  ! start with bilinear estimate
    else
      pmax = 4   ! bilinear with 4 wts/map
      ptot = 4*nid*njd
    endif
  else
    pmax = 1     ! nn with 1 wts/map
    ptot = 1*nid*njd
  endif
  allocate(lis(ptot))
  allocate(lid(ptot))
  allocate(lwt(ptot))
  allocate(pti(pmax))
  allocate(ptj(pmax))
  allocate(ptw(pmax))

  !--- full array copy is default ---
  nwts = nid*njd
  do n=1,nwts
    lid(n) = n
    lis(n) = mod(n-1,nis*njs)+1
    lwt(n) = c1
  enddo

  !--- index copy anytime algo = copy ---
  if (trim(map%algo) == trim(shr_map_fs_copy)) then
    map%init = initcopy
    ! just use copy default

  !--- for fill ---
  elseif (trim(map%type) == trim(shr_map_fs_fill) .or. &
          trim(map%type) == trim(shr_map_fs_cfill)) then
    map%init = initcopy
    if (trim(map%algo) == trim(shr_map_fs_spval)) then
      maskdstbysrc = .true.
    elseif (trim(map%algo) == trim(shr_map_fs_nn)) then
      do n=1,nwts
        call shr_map_1dto2d(lis(n),nis,njs,i,j)
        if (Msrc(i,j) == 0) then
          call shr_map_findnn(Xsrc(i,j),Ysrc(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnoni)) then
      do n=1,nwts
        call shr_map_1dto2d(lis(n),nis,njs,i,j)
        if (Msrc(i,j) == 0) then
          call shr_map_findnnon('i',Xsrc(i,j),Ysrc(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnonj)) then
      do n=1,nwts
        call shr_map_1dto2d(lis(n),nis,njs,i,j)
        if (Msrc(i,j) == 0) then
          call shr_map_findnnon('j',Xsrc(i,j),Ysrc(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    else
      call shr_map_abort(subName//' ERROR: unsupported map option combo')
    endif

  !--- for remap ---
  elseif (trim(map%type) == trim(shr_map_fs_remap)) then
    map%init = inispval
    if (trim(map%algo) == trim(shr_map_fs_spval)) then
      nwts = 0
    elseif (trim(map%algo) == trim(shr_map_fs_nn)) then
      do n=1,nwts
        call shr_map_1dto2d(lid(n),nid,njd,i,j)
        call shr_map_findnn(Xdst(i,j),Ydst(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnoni)) then
      do n=1,nwts
        call shr_map_1dto2d(lid(n),nid,njd,i,j)
        call shr_map_findnnon('i',Xdst(i,j),Ydst(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnonj)) then
      do n=1,nwts
        call shr_map_1dto2d(lid(n),nid,njd,i,j)
        call shr_map_findnnon('j',Xdst(i,j),Ydst(i,j),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_bilinear)) then
      nwts = 0
      do n=1,nid*njd
        call shr_map_1dto2d(n,nid,njd,i,j)
        call shr_map_getWts(Xdst(i,j),Ydst(i,j),Xsrc,Ysrc,pti,ptj,ptw,pnum,units)
        if (nwts + pnum > size(lwt)) then
           !--- resize lis, lid, lwt.  ptot is old size, ptot2 is new size
           ptot = size(lwt)
           ptot2 = ptot + max(ptot/2,pnum*10)
           if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) 'resize wts ',ptot,ptot2
           allocate(ltmp(ptot))
           ltmp(1:nwts) = lis(1:nwts)
           deallocate(lis)
           allocate(lis(ptot2))
           lis(1:nwts) = ltmp(1:nwts)
           ltmp(1:nwts) = lid(1:nwts)
           deallocate(lid)
           allocate(lid(ptot2))
           lid(1:nwts) = ltmp(1:nwts)
           deallocate(ltmp)
           allocate(lwtmp(ptot))
           lwtmp(1:nwts) = lwt(1:nwts)
           deallocate(lwt)
           allocate(lwt(ptot2))
           lwt(1:nwts) = lwtmp(1:nwts)
           deallocate(lwtmp)
        endif
        do n1 = 1,pnum
          nwts = nwts + 1
          lid(nwts) = n
          call shr_map_2dto1d(lis(nwts),nis,njs,pti(n1),ptj(n1))
          lwt(nwts) = ptw(n1)
        enddo
      enddo
    else
      call shr_map_abort(subName//' ERROR: unsupported map option combo')
      if (present(rc)) rc = 1
      return
    endif
  else
    call shr_map_abort(subName//' ERROR: unsupported map option combo')
    if (present(rc)) rc = 1
    return
  endif

!--- compress weights and copy to map ---
  !--- remove 1:1 copies if initcopy
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'init: ',map%init
  if (map%init == initcopy .and. &
      trim(map%type) /= trim(shr_map_fs_cfill)) then
    ncnt = 0
    do n=1,nwts
      if (lid(n) == lis(n) .and. abs(lwt(n)-c1) < eps) then
        ! skipit
      else
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove dst grid points ---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'maskdst: ',maskdst
  if (maskdst) then
    ncnt = 0
    do n=1,nwts
      call shr_map_1dto2d(lid(n),nid,njd,i,j)
      if (Mdst(i,j) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove dst grid points based on src mask---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'maskdstbysrc: ',maskdstbysrc
  if (maskdstbysrc) then
    ncnt = 0
    do n=1,nwts
      call shr_map_1dto2d(lid(n),nid,njd,i,j)
      if (Msrc(i,j) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points by src, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove src grid points ---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'masksrc: ',masksrc
  if (masksrc) then
    ncnt = 0
    do n=1,nwts
      call shr_map_1dto2d(lis(n),nis,njs,i,j)
      if (Msrc(i,j) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm src grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif
  
  !--- renormalize wgts to 1.0 ---
  allocate(sum(nid*njd))
  !--- sum weights for dst grid points ---
  sum(:) = c0
  do n=1,nwts
    sum(lid(n)) = sum(lid(n)) + lwt(n)
  enddo
  !--- print min/max sum ---
  rmin = maxval(sum)
  rmax = minval(sum)
  do n=1,nid*njd
    if (sum(n) > eps) then
      rmin = min(rmin,sum(n))
      rmax = max(rmax,sum(n))
    endif
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F03) 'sum wts min/max ',rmin,rmax
  !--- renormalize so sum on destination is always 1.0 for active dst points
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'renorm: ',renorm
  if (renorm) then
    do n=1,nwts
      if (sum(lid(n)) > eps) then
        lwt(n) = lwt(n) / sum(lid(n))
      endif
    enddo
    !--- sum weights for dst grid points ---
    sum(:) = c0
    do n=1,nwts
      sum(lid(n)) = sum(lid(n)) + lwt(n)
    enddo
    !--- print min/max sum ---
    rmin = maxval(sum)
    rmax = minval(sum)
    do n=1,nid*njd
      if (sum(n) > eps) then
        rmin = min(rmin,sum(n))
        rmax = max(rmax,sum(n))
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F03) 'sum wts min/max ',rmin,rmax
  endif

  map%nwts = nwts
! deallocate(map%idst,stat=irc)
! deallocate(map%isrc,stat=irc)
! deallocate(map%wgts,stat=irc)
  allocate(map%idst(nwts))
  allocate(map%isrc(nwts))
  allocate(map%wgts(nwts))
  do n=1,nwts
    map%idst(n) = lid(n)
    map%isrc(n) = lis(n)
    map%wgts(n) = lwt(n)
  enddo

  deallocate(Xdst)

  deallocate(lis)
  deallocate(lid)
  deallocate(lwt)
  deallocate(sum)

  deallocate(pti)
  deallocate(ptj)
  deallocate(ptw)

  map%fill = fillstring
  call shr_map_checkWgts_global(Msrc,Mdst,map)

  if (present(rc)) rc = lrc

end subroutine shr_map_mapSet_global  

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_mapSet_dest -- Compute mapping weights
!
! !DESCRIPTION:
!     Compute mapping weights based on setting in map.  Fill the
!     weights in the map.  Currently supported maps and action:
!        fill :copy  = copy array by index, mask switch used
!        fill :spval = copy array, fill with spval, mask switch not used
!        fill :nn*   = copy array, fill with nnval, mask switch not used
!        remap:copy  = copy array by index, mask switch used
!        remap:spval = sets array to spval, mask switch used
!        remap:bil*  = bilinear interpolation, mask switch used
!        remap:nn*   = sets array to nnval, mask switch used
!     \newline
!     Requirements for input grids:
!       Xsrc,Ysrc must be regular lat/lon grid, monotonically increasing
!          or decreasing, can be degrees or radians
!       Xdst,Ydst are arbitrary list of lats/lons, must be same units as src
!       Msrc,Mdst have nonzero for active grid point, zero for non-active
!       src and dst must be the grid for type = fill
!     Grids are check for validity
!     \newline
!     call shr\_map\_mapSet(mymap,Xs,Ys,Ms,Xd,Yd,Md)
!     \newline
!     call shr\_map\_mapSet(mymap,Xs,Ys,Ms,Xd,Yd,Md,algo='bilinear')
!
! !REMARKS
!     If bothmask or srcmask is used with remap and some algorithms, active
!     dst grid points can have invalid values.  A report is produced after
!     weights are calculated and this information will be detailed.
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_mapSet_dest(map,Xsrc,Ysrc,Msrc,Xdst_in,Ydst,Mdst,ndst,Idst,name,type,algo,mask,vect,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_map_mapType) ,intent(inout):: map      ! map
  real(SHR_KIND_R8)     ,intent(in) :: Xsrc(:,:)  ! lon of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Ysrc(:,:)  ! lat of src grid
  integer(SHR_KIND_IN)  ,intent(in) :: Msrc(:,:)  ! mask of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Xdst_in(:) ! lon of dst grid
  real(SHR_KIND_R8)     ,intent(in) :: Ydst(:)    ! lat of dst grid
  integer(SHR_KIND_IN)  ,intent(in) :: Mdst(:)    ! mask of dst grid
  integer(SHR_KIND_IN)  ,intent(in) :: ndst       ! global size of dst
  integer(SHR_KIND_IN)  ,intent(in) :: Idst(:)    ! global index of dst grid
  character(*) ,optional,intent(in) :: name       ! name
  character(*) ,optional,intent(in) :: type       ! type
  character(*) ,optional,intent(in) :: algo       ! algo
  character(*) ,optional,intent(in) :: mask       ! mask
  character(*) ,optional,intent(in) :: vect       ! vect
  integer(SHR_KIND_IN),optional,intent(out) :: rc ! error code

!EOP

  !--- local ---
  integer(SHR_KIND_IN)   :: nis,njs,nid,njd
  integer(SHR_KIND_IN)   :: nwts,n,n1,n2,ncnt,i,j,inn,jnn
  integer(SHR_KIND_IN)   :: irc,lrc
  real(SHR_KIND_R8) :: rmin,rmax                     ! min/max value
  real(SHR_KIND_R8) :: cang                          ! circle angle, deg or rad
  real(SHR_KIND_R8),allocatable    :: Xdst(:)        ! lon of dst grid, wrapped as needed

  integer(SHR_KIND_IN)             :: pmax           ! max num of wgts in pti...
  integer(SHR_KIND_IN)             :: ptot,ptot2     ! max num of wgts in lis...
  integer(SHR_KIND_IN)             :: pnum           ! num of wgts set in pti...
  integer(SHR_KIND_IN),allocatable :: pti(:)         ! i index for wgts
  integer(SHR_KIND_IN),allocatable :: ptj(:)         ! j index for wgts
  real(SHR_KIND_R8)   ,allocatable :: ptw(:)         ! weights for pti,ptj

  integer(SHR_KIND_IN),allocatable :: lis(:)         ! tmp src/dst index
  integer(SHR_KIND_IN),allocatable :: lid(:)         ! tmp src/dst index
  real(SHR_KIND_R8)   ,allocatable :: lwt(:)         ! tmp wgt array
  real(SHR_KIND_R8)   ,allocatable :: sum(:)         ! tmp sum array
  integer(SHR_KIND_IN),allocatable :: ltmp(:)        ! tmp src/dst index, for resize
  real(SHR_KIND_R8)   ,allocatable :: lwtmp(:)       ! tmp wgt array, for resize

  character(len=8) :: units    ! radians or degrees

  logical :: masksrc       ! local var to turn on masking using src mask
  logical :: maskdst       ! local var to turn on masking using dst mask
  logical :: maskdstbysrc  ! local var to turn on masking using src mask for
                           ! dst array, especially for fill
  logical :: renorm        ! local var to turn on renormalization

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_mapSet_dest') "
  character(*),parameter :: F00     = "('(shr_map_mapSet_dest) ',a) "
  character(*),parameter :: F01     = "('(shr_map_mapSet_dest) ',a,l2) "
  character(*),parameter :: F02     = "('(shr_map_mapSet_dest) ',a,2i8) "
  character(*),parameter :: F03     = "('(shr_map_mapSet_dest) ',a,2e20.13) "

!-------------------------------------------------------------------------------

  write(s_logunit,F00) 'ERROR this routine is not validated'
  call shr_sys_abort(subName//' ERROR subroutine not validated')

  lrc = 0
  if (present(rc)) rc = lrc

  if (present(name)) call shr_map_put(map,shr_map_fs_name,name)
  if (present(type)) call shr_map_put(map,shr_map_fs_type,type,verify=.true.)
  if (present(algo)) call shr_map_put(map,shr_map_fs_algo,algo,verify=.true.)
  if (present(mask)) call shr_map_put(map,shr_map_fs_mask,mask,verify=.true.)
  if (present(vect)) call shr_map_put(map,shr_map_fs_vect,vect,verify=.true.)
  map%init = inispval

  if (.NOT.shr_map_checkInit(map)) then
    call shr_map_abort(subName//' ERROR map not initialized')
  endif

  !--- is lat/lon degrees or radians? ---
  cang = 360._SHR_KIND_R8
  units = 'degrees'
  if (shr_map_checkRad(Ysrc)) then
    cang=c2*pi
    units = 'radians'
  endif

  nis = size(Xsrc,1)
  njs = size(Xsrc,2)
  nid = size(Xdst_in,1)
  njd = 1

  !--- shift Xdst by 2pi to range of Xsrc as needed ---
  allocate(Xdst(nid))
  rmin = minval(Xsrc)
  rmax = maxval(Xsrc)
  do i=1,nid
    Xdst(i) = Xdst_in(i)
    do while ((Xdst(i) < rmin .and. Xdst(i)+cang <= rmax).or.  &
              (Xdst(i) > rmax .and. Xdst(i)-cang >= rmin))
      if (Xdst(i) < rmin) then
        Xdst(i) = Xdst(i) + cang
      elseif (Xdst(i) > rmax) then
        Xdst(i) = Xdst(i) - cang
      else
        call shr_sys_abort(subName//' ERROR in Xdst wrap')
      endif
    enddo
  enddo
  
  call shr_map_checkGrids_dest(Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,map,lrc)

  map%nwts = 0
  map%nsrc = nis*njs
  map%ndst = ndst

! deallocate(map%xsrc,stat=irc) ! this used to be a safe way to delloc when necessary,
! deallocate(map%ysrc,stat=irc) ! but do nothing when pointers were undefined or
! deallocate(map%xdst,stat=irc) ! un-associated, in Oct 2005, undefined ptrs started
! deallocate(map%ydst,stat=irc) ! causing seg-faults on bluesky (B. Kauffman)
  allocate(map%xsrc(nis*njs))
  allocate(map%ysrc(nis*njs))
  allocate(map%xdst(nid*njd))
  allocate(map%ydst(nid*njd))
  do j=1,njs
  do i=1,nis
    call shr_map_2dto1d(n1,nis,njs,i,j)
    map%xsrc(n1) = Xsrc(i,j)*c2*pi/cang
    map%ysrc(n1) = Ysrc(i,j)*c2*pi/cang
  enddo
  enddo
  do i=1,nid
    map%xdst(i) = Xdst(i)*c2*pi/cang
    map%ydst(i) = Ydst(i)*c2*pi/cang
  enddo

  masksrc = .false.
  maskdstbysrc = .false.
  maskdst = .false.
  renorm  = .true.

  if (trim(map%type) /= trim(shr_map_fs_fill) .and. &
      trim(map%type) /= trim(shr_map_fs_cfill)) then
    if (trim(map%mask) == trim(shr_map_fs_bothmask) .or. &
        trim(map%mask) == trim(shr_map_fs_srcmask)) masksrc = .true.
    if (trim(map%mask) == trim(shr_map_fs_bothmask) .or. &
        trim(map%mask) == trim(shr_map_fs_dstmask)) maskdst = .true.
  endif
  if (trim(map%algo) == trim(shr_map_fs_spval)) then
    masksrc = .false.
    renorm = .false.
  endif

  if (debug > 1) then
    if (s_loglev > 0) write(s_logunit,*) ' '
    call shr_map_print(map)
  endif

  if (lrc /= 0) then
    if (present(rc)) rc = lrc
    return
  endif

  if (trim(map%algo) == trim(shr_map_fs_bilinear)) then
    if (dopole) then
      pmax = nis+2   ! possible for high lat points
      ptot = 4*nid*njd  ! start with bilinear estimate
    else
      pmax = 4   ! bilinear with 4 wts/map
      ptot = 4*nid*njd
    endif
  else
    pmax = 1     ! nn with 1 wts/map
    ptot = 1*nid*njd
  endif
  allocate(lis(ptot))
  allocate(lid(ptot))
  allocate(lwt(ptot))
  allocate(pti(pmax))
  allocate(ptj(pmax))
  allocate(ptw(pmax))

  !--- full array copy is default ---
  nwts = nid*njd
  do n=1,nwts
    lid(n) = Idst(n)
    lis(n) = Idst(n)
    lwt(n) = c1
  enddo

  !--- index copy anytime algo = copy ---
  if (trim(map%algo) == trim(shr_map_fs_copy)) then
    map%init = initcopy
    ! just use copy default

  !--- for fill ---
  elseif (trim(map%type) == trim(shr_map_fs_fill) .or. &
          trim(map%type) == trim(shr_map_fs_cfill)) then
    map%init = initcopy
    if (trim(map%algo) == trim(shr_map_fs_spval)) then
      maskdstbysrc = .true.
    elseif (trim(map%algo) == trim(shr_map_fs_nn)) then
      do n=1,nwts
        if (Mdst(n) == 0) then
          call shr_map_findnn(Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnoni)) then
      do n=1,nwts
        if (Mdst(n) == 0) then
          call shr_map_findnnon('i',Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnonj)) then
      do n=1,nwts
        if (Mdst(n) == 0) then
          call shr_map_findnnon('j',Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
          call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
        endif
      enddo
    else
      call shr_map_abort(subName//' ERROR: unsupported map option combo')
    endif

  !--- for remap ---
  elseif (trim(map%type) == trim(shr_map_fs_remap)) then
    map%init = inispval
    if (trim(map%algo) == trim(shr_map_fs_spval)) then
      nwts = 0
    elseif (trim(map%algo) == trim(shr_map_fs_nn)) then
      do n=1,nwts
        call shr_map_findnn(Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnoni)) then
      do n=1,nwts
        call shr_map_findnnon('i',Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_nnonj)) then
      do n=1,nwts
        call shr_map_findnnon('j',Xdst(n),Ydst(n),Xsrc,Ysrc,Msrc,inn,jnn)
        call shr_map_2dto1d(lis(n),nis,njs,inn,jnn)
      enddo
    elseif (trim(map%algo) == trim(shr_map_fs_bilinear)) then
      nwts = 0
      do n=1,nid*njd
        call shr_map_getWts(Xdst(n),Ydst(n),Xsrc,Ysrc,pti,ptj,ptw,pnum,units)
        if (nwts + pnum > size(lwt)) then
           !--- resize lis, lid, lwt.  ptot is old size, ptot2 is new size
           ptot = size(lwt)
           ptot2 = ptot + max(ptot/2,pnum*10)
           if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) 'resize wts ',ptot,ptot2
           allocate(ltmp(ptot))
           ltmp(1:nwts) = lis(1:nwts)
           deallocate(lis)
           allocate(lis(ptot2))
           lis(1:nwts) = ltmp(1:nwts)
           ltmp(1:nwts) = lid(1:nwts)
           deallocate(lid)
           allocate(lid(ptot2))
           lid(1:nwts) = ltmp(1:nwts)
           deallocate(ltmp)
           allocate(lwtmp(ptot))
           lwtmp(1:nwts) = lwt(1:nwts)
           deallocate(lwt)
           allocate(lwt(ptot2))
           lwt(1:nwts) = lwtmp(1:nwts)
           deallocate(lwtmp)
        endif
        do n1 = 1,pnum
          nwts = nwts + 1
          lid(nwts) = Idst(n)
          call shr_map_2dto1d(lis(nwts),nis,njs,pti(n1),ptj(n1))
          lwt(nwts) = ptw(n1)
        enddo
      enddo
    else
      call shr_map_abort(subName//' ERROR: unsupported map option combo')
      if (present(rc)) rc = 1
      return
    endif
  else
    call shr_map_abort(subName//' ERROR: unsupported map option combo')
    if (present(rc)) rc = 1
    return
  endif

!--- compress weights and copy to map ---
  !--- remove 1:1 copies if initcopy
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'init: ',map%init
  if (map%init == initcopy .and. &
      trim(map%type) /= trim(shr_map_fs_cfill)) then
    ncnt = 0
    do n=1,nwts
      if (lid(n) == lis(n) .and. abs(lwt(n)-c1) < eps) then
        ! skipit
      else
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove dst grid points ---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'maskdst: ',maskdst
  if (maskdst) then
    ncnt = 0
    do n=1,nwts
      if (Mdst(n) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove dst grid points based on src mask---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'maskdstbysrc: ',maskdstbysrc
  if (maskdstbysrc) then
    ncnt = 0
    do n=1,nwts
      call shr_map_1dto2d(lid(n),nis,njs,i,j)
      if (Msrc(i,j) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm dst grid points by src, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif

  !--- remove src grid points ---
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'masksrc: ',masksrc
  if (masksrc) then
    ncnt = 0
    do n=1,nwts
      call shr_map_1dto2d(lis(n),nis,njs,i,j)
      if (Msrc(i,j) /= 0) then
        ncnt = ncnt+1
        lid(ncnt) = lid(n)
        lis(ncnt) = lis(n)
        lwt(ncnt) = lwt(n)
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F02) ' rm src grid points, nwts old/new = ',nwts,ncnt
    nwts = ncnt
  endif
  
  !--- renormalize wgts to 1.0 ---
  allocate(sum(ndst))
  !--- sum weights for dst grid points ---
  sum(:) = c0
  do n=1,nwts
    sum(lid(n)) = sum(lid(n)) + lwt(n)
  enddo
  !--- print min/max sum ---
  rmin = maxval(sum)
  rmax = minval(sum)
  do n=1,ndst
    if (sum(n) > eps) then
      rmin = min(rmin,sum(n))
      rmax = max(rmax,sum(n))
    endif
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F03) 'sum wts min/max ',rmin,rmax
  !--- renormalize so sum on destination is always 1.0 for active dst points
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) 'renorm: ',renorm
  if (renorm) then
    do n=1,nwts
      if (sum(lid(n)) > eps) then
        lwt(n) = lwt(n) / sum(lid(n))
      endif
    enddo
    !--- sum weights for dst grid points ---
    sum(:) = c0
    do n=1,nwts
      sum(lid(n)) = sum(lid(n)) + lwt(n)
    enddo
    !--- print min/max sum ---
    rmin = maxval(sum)
    rmax = minval(sum)
    do n=1,nid*njd
      if (sum(n) > eps) then
        rmin = min(rmin,sum(n))
        rmax = max(rmax,sum(n))
      endif
    enddo
    if (debug > 1 .and. s_loglev > 0) write(s_logunit,F03) 'sum wts min/max ',rmin,rmax
  endif

  map%nwts = nwts
! deallocate(map%idst,stat=irc)
! deallocate(map%isrc,stat=irc)
! deallocate(map%wgts,stat=irc)
  allocate(map%idst(nwts))
  allocate(map%isrc(nwts))
  allocate(map%wgts(nwts))
  do n=1,nwts
    map%idst(n) = lid(n)
    map%isrc(n) = lis(n)
    map%wgts(n) = lwt(n)
  enddo

  deallocate(Xdst)

  deallocate(lis)
  deallocate(lid)
  deallocate(lwt)
  deallocate(sum)

  deallocate(pti)
  deallocate(ptj)
  deallocate(ptw)

  map%fill = fillstring
!!  call shr_map_checkWgts_dest(Msrc,Mdst,map)

  if (present(rc)) rc = lrc

end subroutine shr_map_mapSet_dest  

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_mapDatam -- maps arrays using input map
!
! !DESCRIPTION:
!     Maps arrays using preset map
!     \newline
!     call shr\_map\_mapData(Ain,Aout,mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_mapDatam(arrsrc,arrdst,map)
  !--- map arrsrc to arrdst, each array is dimension (fields,grid index) ---

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  real(SHR_KIND_R8)     ,intent(in) :: arrsrc(:,:)  ! src array(fields,grid)
  real(SHR_KIND_R8)     ,intent(out):: arrdst(:,:)  ! dst array(fields,grid)
  type(shr_map_mapType) ,intent(in) :: map      ! map

!EOP

  !--- local ---
  integer(SHR_KIND_IN) :: n,n2            ! counters
  integer(SHR_KIND_IN) :: indi,indo       ! array indices, in/out
  real(SHR_KIND_R8)    :: wgt             ! value of weight
  integer(SHR_KIND_IN) :: nfi,nfo         ! number of fields in array, in/out
  integer(SHR_KIND_IN) :: nsi,nso         ! size of grid in array, in/out
  real(SHR_KIND_R8)    :: theta           ! angle difference
  integer(SHR_KIND_IN),save      :: t0=-1,t1,t2,t3,t4,t5 ! timers
  integer(SHR_KIND_IN),parameter :: timing=0        ! turn timers off/on (0/1)
  logical,pointer      :: initnew(:)      ! mask for initialization

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_mapDatam') "
  character(*),parameter :: F00     = "('(shr_map_mapDatam) ',a) "
  character(*),parameter :: F01     = "('(shr_map_mapDatam) ',a,2i8) "

!-------------------------------------------------------------------------------

  if (timing>0 .and. t0 == -1) then
     call shr_timer_get(t0,subName//"everything")
     call shr_timer_get(t1,subName//"initial checks")
     call shr_timer_get(t2,subName//"dst to spval")
     call shr_timer_get(t4,subName//"map vector")
     call shr_timer_get(t5,subName//"map scalar")
  end if

  if (timing>0) call shr_timer_start(t0)
  if (timing>0) call shr_timer_start(t1)

  !--- get number of fields ---
  nfi = size(arrsrc,1)
  nfo = size(arrdst,1)

  !--- check number of fields ---
  if (nfi /= nfo) then
    write(s_logunit,F01) ' field numbers dont match ',nfi,nfo
    call shr_map_abort(subName//' ERROR number of fields')
  endif

  !--- check two fields for vector ---
  if (trim(map%vect) == trim(shr_map_fs_vector).and.(nfi /= 2)) then
    write(s_logunit,F01) ' vector mapping, must map only two fields',nfi,nfo
    call shr_map_abort(subName//' ERROR vector mapping fields not two')
  endif

  !--- check that map is set ---
  if (.not.shr_map_checkFilled(map)) then
    write(s_logunit,F00) ' map is not filled'
    call shr_map_abort(subName//' ERROR map is not filled')
  endif

  !--- get size of grid ---
  nsi = size(arrsrc,2)
  nso = size(arrdst,2)

  !--- check size of grid ---
  if (nsi /= map%nsrc) then
    write(s_logunit,F01) ' src grid size doesnt match ',nsi,map%nsrc
    call shr_map_abort(subName//' ERROR src grid size')
  endif
  if (nso /= map%ndst) then
    write(s_logunit,F01) ' dst grid size doesnt match ',nso,map%ndst
    call shr_map_abort(subName//' ERROR dst grid size')
  endif

  if (timing>0) call shr_timer_stop(t1)
  if (timing>0) call shr_timer_start(t2)

  allocate(initnew(1:nso))
  initnew = .true.
  !--- set arrdst to spval, all points, default ---
  if (map%init == inispval) then
    arrdst = shr_map_spval
  elseif (map%init == initcopy) then
    if (nsi /= nso) then
      write(s_logunit,F01) ' initcopy has nsi ne nso ',nsi,nso
      call shr_map_abort(subName//' ERROR initcopy size')
    else
      do n = 1,nsi
      do n2 = 1,nfo
        arrdst(n2,n) = arrsrc(n2,n)
      enddo
      enddo
    endif
  else
    write(s_logunit,F00) ' map%init illegal '//trim(map%init)
    call shr_map_abort(subName//' ERROR map init')
  endif

  if (timing>0) call shr_timer_stop(t2)

  !--- generate output array ---
  if (trim(map%vect) == trim(shr_map_fs_vector)) then
    if (timing>0) call shr_timer_start(t4)
    do n=1,map%nwts
      indi = map%isrc(n)
      indo = map%idst(n)
      wgt  = map%wgts(n)
      theta = map%xdst(indo) - map%xsrc(indi)
      if (initnew(indo)) then
        initnew(indo) = .false.
        arrdst(1,indo) = wgt*( arrsrc(1,indi)*cos(theta)  &
                              +arrsrc(2,indi)*sin(theta))
        arrdst(2,indo) = wgt*(-arrsrc(1,indi)*sin(theta)  &
                              +arrsrc(2,indi)*cos(theta))
      else
        arrdst(1,indo) = arrdst(1,indo) + wgt*( arrsrc(1,indi)*cos(theta)  &
                                               +arrsrc(2,indi)*sin(theta))
        arrdst(2,indo) = arrdst(2,indo) + wgt*(-arrsrc(1,indi)*sin(theta)  &
                                               +arrsrc(2,indi)*cos(theta))
      endif
    enddo
    if (timing>0) call shr_timer_stop(t4)
  else
    if (timing>0) call shr_timer_start(t5)
    do n=1,map%nwts
      indi = map%isrc(n)
      indo = map%idst(n)
      wgt  = map%wgts(n)
      if (initnew(indo)) then
        initnew(indo) = .false.
        do n2 = 1,nfo
          arrdst(n2,indo) = arrsrc(n2,indi)*wgt
        enddo
      else
        do n2 = 1,nfo
          arrdst(n2,indo) = arrdst(n2,indo) + arrsrc(n2,indi)*wgt
        enddo
      endif
    enddo
    if (timing>0) call shr_timer_stop(t5)
  endif

  deallocate(initnew)

  if (timing>0) call shr_timer_stop(t0)

end subroutine shr_map_mapDatam

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_mapDatanm -- maps arrays without map
!
! !DESCRIPTION:
!     Maps arrays, don't save the map
!     \newline
!     call shr\_map\_mapData(Ain,Aout,Xs,Ys,Ms,Xd,Yd,Md,name,type,algo,vect,rc)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_mapDatanm(arrsrc,arrdst,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,name,type,algo,mask,vect,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  !--- map arrsrc to arrdst, each array is dimension (fields,grid index) ---
  real(SHR_KIND_R8)     ,intent(in) :: arrsrc(:,:)  ! src array(fields,grid)
  real(SHR_KIND_R8)     ,intent(out):: arrdst(:,:)  ! dst array(fields,grid)
  real(SHR_KIND_R8)     ,intent(in) :: Xsrc(:,:)  ! lon of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Ysrc(:,:)  ! lat of src grid
  integer(SHR_KIND_IN)  ,intent(in) :: Msrc(:,:)  ! mask of src grid
  real(SHR_KIND_R8)     ,intent(in) :: Xdst(:,:)  ! lon of dst grid
  real(SHR_KIND_R8)     ,intent(in) :: Ydst(:,:)  ! lat of dst grid
  integer(SHR_KIND_IN)  ,intent(in) :: Mdst(:,:)  ! mask of dst grid
  character(*)          ,intent(in) :: name       ! name
  character(*)          ,intent(in) :: type       ! type
  character(*)          ,intent(in) :: algo       ! algo
  character(*)          ,intent(in) :: mask       ! mask
  character(*) ,optional,intent(in) :: vect       ! vect
  integer(SHR_KIND_IN),optional,intent(out) :: rc ! error code

!EOP

  !--- local ---
  type(shr_map_mapType) :: map
  integer(SHR_KIND_IN)  :: lrc

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_mapDatanm') "
  character(*),parameter :: F00     = "('(shr_map_mapDatanm) ',a) "

!-------------------------------------------------------------------------------

  lrc = 0

  call shr_map_put(map,shr_map_fs_name,name,verify=.false.)
  call shr_map_put(map,shr_map_fs_type,type,verify=.true.)
  call shr_map_put(map,shr_map_fs_algo,algo,verify=.true.)
  call shr_map_put(map,shr_map_fs_mask,mask,verify=.true.)
  if (present(vect)) then
    call shr_map_put(map,shr_map_fs_vect,vect,verify=.true.)
  else
    call shr_map_put(map,shr_map_fs_vect,'scalar',verify=.true.)
  endif
  call shr_map_mapSet(map,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,rc=lrc)
  call shr_map_mapData(arrsrc,arrdst,map)

  call shr_map_clean(map)

  if (present(rc)) rc = lrc

end subroutine shr_map_mapDatanm

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_setAbort -- local interface to abort routine
!
! !DESCRIPTION:
!     Set doabort flag for shr_map methods, true = call shr\_sys\_abort,
!     false = write error message and continue
!     \newline
!     call shr\_map\_abort(subName//' ERROR: illegal option')
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_setAbort') "
  character(*),parameter :: F00     = "('(shr_map_setAbort) ',a) "

!-------------------------------------------------------------------------------

  doabort = flag

end subroutine shr_map_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_setDebug -- set local debug level
!
! !DESCRIPTION:
!     Set debug level for shr_map methods, 0 = production
!     \newline
!     call shr\_map\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_setDebug') "
  character(*),parameter :: F00     = "('(shr_map_setDebug) ',a) "

!-------------------------------------------------------------------------------

  debug = iflag

end subroutine shr_map_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_setDopole -- set dopole flag
!
! !DESCRIPTION:
!     set dopole flag
!     \newline
!     call shr\_map\_setDopole(flag)
!
! !REVISION HISTORY:
!     2009-Jun-22 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_setDopole(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical, intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_setDopole') "
  character(*),parameter :: F00     = "('(shr_map_setDopole) ',a) "

  dopole = flag

end subroutine shr_map_setDopole

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_abort -- local interface to abort routine
!
! !DESCRIPTION:
!     Local interface to abort routine.  Depending on local flag, abort,
!     either calls shr\_sys\_abort or writes abort message and continues.
!     \newline
!     call shr\_map\_abort(subName//' ERROR: illegal option')
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(in) :: string

!XXEOP

  !--- local ---
  character(shr_kind_CL) :: lstring

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_abort') "
  character(*),parameter :: F00     = "('(shr_map_abort) ',a) "

!-------------------------------------------------------------------------------

  lstring = ''
  if (present(string)) lstring = string

  if (doabort) then
    call shr_sys_abort(lstring)
  else
    write(s_logunit,F00) trim(lstring)
  endif

end subroutine shr_map_abort

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_checkGrids_global -- local routine to check mapSet grids
!
! !DESCRIPTION:
!     Local method to check grid arguments in shr\_map\_mapSet
!     \newline
!     call shr\_map\_checkGrids_global(Xs,Ys,Ms,Xd,Yd,Md,mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_checkGrids_global(Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,map,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:
  real(SHR_KIND_R8)    ,intent(in) :: Xsrc(:,:)   ! src lat
  real(SHR_KIND_R8)    ,intent(in) :: Ysrc(:,:)   ! src lon
  integer(SHR_KIND_IN) ,intent(in) :: Msrc(:,:)   ! src mask
  real(SHR_KIND_R8)    ,intent(in) :: Xdst(:,:)   ! dst lat
  real(SHR_KIND_R8)    ,intent(in) :: Ydst(:,:)   ! dst lon
  integer(SHR_KIND_IN) ,intent(in) :: Mdst(:,:)   ! dst mask
  type(shr_map_mapType),intent(in) :: map     ! map
  integer(SHR_KIND_IN),optional,intent(out) :: rc ! error code

!XXEOP

  !--- local ---
  integer(SHR_KIND_IN) :: i,j,nis,njs,nid,njd,ncnt
  logical              :: error,flag

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkGrids_global') "
  character(*),parameter :: F00     = "('(shr_map_checkGrids_global) ',a) "
  character(*),parameter :: F01     = "('(shr_map_checkGrids_global) ',a,2i8) "
  character(*),parameter :: F02     = "('(shr_map_checkGrids_global) ',a,4i8) "
  character(*),parameter :: F03     = "('(shr_map_checkGrids_global) ',a,2g20.13) "
  character(*),parameter :: F04     = "('(shr_map_checkGrids_global) ',a,i8,a,i8) "
  character(*),parameter :: F05     = "('(shr_map_checkGrids_global) ',a,i8,2g20.13) "
  character(*),parameter :: F06     = "('(shr_map_checkGrids_global) ',a,2i8,2g20.13) "

!-------------------------------------------------------------------------------

  error = .false.
  if (present(rc)) rc = 0

  !--- get size of X arrays
  nis = size(Xsrc,1)
  njs = size(Xsrc,2)
  nid = size(Xdst,1)
  njd = size(Xdst,2)

  !--- check array size consistency for src and dst
  if (size(Ysrc,1) /= nis) then
    write(s_logunit,F01) 'ERROR Xsrc,Ysrc i-dim mismatch',nis,size(Ysrc,1)
    error = .true.
  endif
  if (size(Ysrc,2) /= njs) then
    write(s_logunit,F01) 'ERROR Xsrc,Ysrc j-dim mismatch',njs,size(Ysrc,2)
    error = .true.
  endif
  if (size(Msrc,1) /= nis) then
    write(s_logunit,F01) 'ERROR Xsrc,Msrc i-dim mismatch',nis,size(Msrc,1)
    error = .true.
  endif
  if (size(Msrc,2) /= njs) then
    write(s_logunit,F01) 'ERROR Xsrc,Msrc j-dim mismatch',njs,size(Msrc,2)
    error = .true.
  endif
  if (size(Ydst,1) /= nid) then
    write(s_logunit,F01) 'ERROR Xdst,Ydst i-dim mismatch',nid,size(Ydst,1)
    error = .true.
  endif
  if (size(Ydst,2) /= njd) then
    write(s_logunit,F01) 'ERROR Xdst,Ydst j-dim mismatch',njd,size(Ydst,2)
    error = .true.
  endif
  if (size(Mdst,1) /= nid) then
    write(s_logunit,F01) 'ERROR Xdst,Mdst i-dim mismatch',nid,size(Mdst,1)
    error = .true.
  endif
  if (size(Mdst,2) /= njd) then
    write(s_logunit,F01) 'ERROR Xdst,Mdst j-dim mismatch',njd,size(Mdst,2)
    error = .true.
  endif

  !--- fill type must have same grid size on src and dst ---
  if (trim(map%type) == trim(shr_map_fs_fill) .or. &
      trim(map%type) == trim(shr_map_fs_cfill)) then
    if (nis*njs /= nid*njd) then
      write(s_logunit,F02) 'ERROR: fill type, src/dst sizes ',nis*njs,nid*njd
      error = .true.
    endif
  endif

  !--- write min/max or X, Y and M count ---
  if (debug > 1 .and. s_loglev > 0) then
    write(s_logunit,F03) ' Xsrc min/max ',minval(Xsrc),maxval(Xsrc)
    write(s_logunit,F03) ' Ysrc min/max ',minval(Ysrc),maxval(Ysrc)
    write(s_logunit,F03) ' Xdst min/max ',minval(Xdst),maxval(Xdst)
    write(s_logunit,F03) ' Ydst min/max ',minval(Ydst),maxval(Ydst)
  endif

  ncnt = 0
  do j=1,njs
  do i=1,nis
    if (Msrc(i,j) == 0) ncnt = ncnt + 1
  enddo
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F04) ' Msrc mask T ',nis*njs-ncnt,' of ',nis*njs

  ncnt = 0
  do j=1,njd
  do i=1,nid
    if (Mdst(i,j) == 0) ncnt = ncnt + 1
  enddo
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F04) ' Mdst mask T ',nid*njd-ncnt,' of ',nid*njd

  if (trim(map%algo) == trim(shr_map_fs_bilinear)) then

    !--- check that Xsrc is monotonically increasing for bilinear ---
    flag = .false.
    i = 1
    do while (i < nis .and. .not.flag)
      if (((Xsrc(nis,1) > Xsrc(1,1)) .and. (Xsrc(i+1,1) <=  Xsrc(i,1))) .or. &
          ((Xsrc(nis,1) < Xsrc(1,1)) .and. (Xsrc(i+1,1) >=  Xsrc(i,1)))) then
        write(s_logunit,F05) 'ERROR Xsrc not monotonic ',i,Xsrc(i+1,1),Xsrc(i,1)
        flag = .true.
        error = .true.
      endif
      i = i+1
    enddo

    !--- check that Ysrc is monotonically increasing for bilinear ---
    flag = .false.
    j = 1
    do while (j < njs .and. .not.flag)
      if (((Ysrc(njs,1) > Ysrc(1,1)) .and. (Ysrc(1,j+1) <=  Ysrc(1,j))) .or. &
          ((Ysrc(njs,1) < Ysrc(1,1)) .and. (Ysrc(1,j+1) >=  Ysrc(1,j)))) then
        write(s_logunit,F05) 'ERROR Ysrc not monotonic ',i,Ysrc(1,j+1),Ysrc(1,j)
        flag = .true.
        error = .true.
      endif
      j = j+1
    enddo

    !--- check that Xsrc and Ysrc are regular lat/lon grids for bilinear
    flag = .false.
    i = 1
    do while (i < nis .and. .not.flag)
      j = 2
      do while (j < njs .and. .not.flag)
        if (abs(Xsrc(i,j)-Xsrc(i,1)) > eps) then
          write(s_logunit,F06) ' ERROR Xsrc not regular lat,lon ',i,j, &
              Xsrc(i,j),Xsrc(1,j)
          flag = .true.
          error = .true.
        endif
        j = j+1
      enddo
      i = i+1
    enddo

    flag = .false.
    j = 1
    do while (j < njs .and. .not.flag)
      i = 2
      do while (i < nis .and. .not.flag)
        if (abs(Ysrc(i,j)-Ysrc(1,j)) > eps) then
          write(s_logunit,F06) ' ERROR Ysrc not regular lat,lon ',i,j, &
            Ysrc(i,j),Ysrc(1,j)
          flag = .true.
          error = .true.
        endif
        i = i+1
      enddo
      j = j+1
    enddo
  endif

  if (error) then
    call shr_map_abort(subName//' ERROR ')
    if (present(rc)) rc = 1
  endif

end subroutine shr_map_checkGrids_global

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_checkGrids_dest -- local routine to check mapSet grids
!
! !DESCRIPTION:
!     Local method to check grid arguments in shr\_map\_mapSet
!     \newline
!     call shr\_map\_checkGrids_dest(Xs,Ys,Ms,Xd,Yd,Md,mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_checkGrids_dest(Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,map,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:
  real(SHR_KIND_R8)    ,intent(in) :: Xsrc(:,:)   ! src lat
  real(SHR_KIND_R8)    ,intent(in) :: Ysrc(:,:)   ! src lon
  integer(SHR_KIND_IN) ,intent(in) :: Msrc(:,:)   ! src mask
  real(SHR_KIND_R8)    ,intent(in) :: Xdst(:)     ! dst lat
  real(SHR_KIND_R8)    ,intent(in) :: Ydst(:)     ! dst lon
  integer(SHR_KIND_IN) ,intent(in) :: Mdst(:)     ! dst mask
  type(shr_map_mapType),intent(in) :: map     ! map
  integer(SHR_KIND_IN),optional,intent(out) :: rc ! error code

!XXEOP

  !--- local ---
  integer(SHR_KIND_IN) :: i,j,nis,njs,nid,njd,ncnt
  logical              :: error,flag

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkGrids_dest') "
  character(*),parameter :: F00     = "('(shr_map_checkGrids_dest) ',a) "
  character(*),parameter :: F01     = "('(shr_map_checkGrids_dest) ',a,2i8) "
  character(*),parameter :: F02     = "('(shr_map_checkGrids_dest) ',a,4i8) "
  character(*),parameter :: F03     = "('(shr_map_checkGrids_dest) ',a,2g20.13) "
  character(*),parameter :: F04     = "('(shr_map_checkGrids_dest) ',a,i8,a,i8) "
  character(*),parameter :: F05     = "('(shr_map_checkGrids_dest) ',a,i8,2g20.13) "
  character(*),parameter :: F06     = "('(shr_map_checkGrids_dest) ',a,2i8,2g20.13) "

!-------------------------------------------------------------------------------

  error = .false.
  if (present(rc)) rc = 0

  !--- get size of X arrays
  nis = size(Xsrc,1)
  njs = size(Xsrc,2)
  nid = size(Xdst,1)
  njd = 1

  !--- check array size consistency for src and dst
  if (size(Ysrc,1) /= nis) then
    write(s_logunit,F01) 'ERROR Xsrc,Ysrc i-dim mismatch',nis,size(Ysrc,1)
    error = .true.
  endif
  if (size(Ysrc,2) /= njs) then
    write(s_logunit,F01) 'ERROR Xsrc,Ysrc j-dim mismatch',njs,size(Ysrc,2)
    error = .true.
  endif
  if (size(Msrc,1) /= nis) then
    write(s_logunit,F01) 'ERROR Xsrc,Msrc i-dim mismatch',nis,size(Msrc,1)
    error = .true.
  endif
  if (size(Msrc,2) /= njs) then
    write(s_logunit,F01) 'ERROR Xsrc,Msrc j-dim mismatch',njs,size(Msrc,2)
    error = .true.
  endif
  if (size(Ydst,1) /= nid) then
    write(s_logunit,F01) 'ERROR Xdst,Ydst i-dim mismatch',nid,size(Ydst,1)
    error = .true.
  endif
  if (size(Mdst,1) /= nid) then
    write(s_logunit,F01) 'ERROR Xdst,Mdst i-dim mismatch',nid,size(Mdst,1)
    error = .true.
  endif

!---  tcraig, can't check this with dest mapset ---
!  !--- fill type must have same grid size on src and dst ---
!  if (trim(map%type) == trim(shr_map_fs_fill) .or. &
!      trim(map%type) == trim(shr_map_fs_cfill)) then
!    if (nis*njs /= nid*njd) then
!      write(s_logunit,F02) 'ERROR: fill type, src/dst sizes ',nis*njs,nid*njd
!      error = .true.
!    endif
!  endif

  !--- write min/max or X, Y and M count ---
  if (debug > 1 .and. s_loglev > 0) then
    write(s_logunit,F03) ' Xsrc min/max ',minval(Xsrc),maxval(Xsrc)
    write(s_logunit,F03) ' Ysrc min/max ',minval(Ysrc),maxval(Ysrc)
    write(s_logunit,F03) ' Xdst min/max ',minval(Xdst),maxval(Xdst)
    write(s_logunit,F03) ' Ydst min/max ',minval(Ydst),maxval(Ydst)
  endif

  ncnt = 0
  do j=1,njs
  do i=1,nis
    if (Msrc(i,j) == 0) ncnt = ncnt + 1
  enddo
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F04) ' Msrc mask T ',nis*njs-ncnt,' of ',nis*njs

  ncnt = 0
  do i=1,nid
    if (Mdst(i) == 0) ncnt = ncnt + 1
  enddo
  if (debug > 1 .and. s_loglev > 0) write(s_logunit,F04) ' Mdst mask T ',nid*njd-ncnt,' of ',nid*njd

  if (trim(map%algo) == trim(shr_map_fs_bilinear)) then

    !--- check that Xsrc is monotonically increasing for bilinear ---
    flag = .false.
    i = 1
    do while (i < nis .and. .not.flag)
      if (Xsrc(i+1,1) <=  Xsrc(i,1)) then
        write(s_logunit,F05) 'ERROR Xsrc not increasing ',i,Xsrc(i+1,1),Xsrc(i,1)
        flag = .true.
        error = .true.
      endif
      i = i+1
    enddo

    !--- check that Ysrc is monotonically increasing for bilinear ---
    flag = .false.
    j = 1
    do while (j < njs .and. .not.flag)
      if (Ysrc(1,j+1) <=  Ysrc(1,j)) then
        write(s_logunit,F05) 'ERROR Ysrc not increasing ',i,Ysrc(1,j+1),Ysrc(1,j)
        flag = .true.
        error = .true.
      endif
      j = j+1
    enddo

    !--- check that Xsrc and Ysrc are regular lat/lon grids for bilinear
    flag = .false.
    i = 1
    do while (i < nis .and. .not.flag)
      j = 2
      do while (j < njs .and. .not.flag)
        if (abs(Xsrc(i,j)-Xsrc(i,1)) > eps) then
          write(s_logunit,F06) ' ERROR Xsrc not regular lat,lon ',i,j, &
              Xsrc(i,j),Xsrc(1,j)
          flag = .true.
          error = .true.
        endif
        j = j+1
      enddo
      i = i+1
    enddo

    flag = .false.
    j = 1
    do while (j < njs .and. .not.flag)
      i = 2
      do while (i < nis .and. .not.flag)
        if (abs(Ysrc(i,j)-Ysrc(1,j)) > eps) then
          write(s_logunit,F06) ' ERROR Ysrc not regular lat,lon ',i,j, &
            Ysrc(i,j),Ysrc(1,j)
          flag = .true.
          error = .true.
        endif
        i = i+1
      enddo
      j = j+1
    enddo
  endif

  if (error) then
    call shr_map_abort(subName//' ERROR ')
    if (present(rc)) rc = 1
  endif

end subroutine shr_map_checkGrids_dest

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_checkWgts_global -- checks weights
!
! !DESCRIPTION:
!     Checks weights in map for validity
!     \newline
!     call shr\_map\_checkWgts_global(Ms,Md,mymap)
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_checkWgts_global(Msrc,Mdst,map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:
  integer(SHR_KIND_IN) ,intent(in) :: Msrc(:,:)  ! src mask
  integer(SHR_KIND_IN) ,intent(in) :: Mdst(:,:)  ! dst mask
  type(shr_map_mapType),intent(in) :: map        ! map

!XXEOP

  !--- local ---
  integer(SHR_KIND_IN) :: i,j,nis,njs,nid,njd,n
  integer(SHR_KIND_IN) :: ic1,ic2,ic3,ic4,ic5     ! counters
  logical              :: error
  real(SHR_KIND_R8),allocatable :: Csrc(:,:)
  real(SHR_KIND_R8),allocatable :: Cdst(:,:)

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_checkWgts_global') "
  character(*),parameter :: F00     = "('(shr_map_checkWgts_global) ',a) "
  character(*),parameter :: F01     = "('(shr_map_checkWgts_global) ',a,i8) "
  character(*),parameter :: F02     = "('(shr_map_checkWgts_global) ',a,3i8) "
  character(*),parameter :: F03     = "('(shr_map_checkWgts_global) ',a,i8,a) "

!-------------------------------------------------------------------------------

  error = .false.

  if (debug > 0) call shr_map_print(map)

  if (map%nwts < 1) then
    if (s_loglev > 0) write(s_logunit,F00) 'WARNING map size is zero'
  endif

  if (size(map%wgts) /= map%nwts .or. &
      size(map%isrc) /= map%nwts .or. &
      size(map%idst) /= map%nwts) then
    call shr_map_abort(subName//'ERROR sizes inconsistent')
  endif

  !--- get size of X arrays
  nis = size(Msrc,1)
  njs = size(Msrc,2)
  nid = size(Mdst,1)
  njd = size(Mdst,2)

  allocate(Csrc(nis,njs))
  allocate(Cdst(nid,njd))

  Csrc = c0
  Cdst = c0

  do n = 1,map%nwts
    call shr_map_1dto2d(map%isrc(n),nis,njs,i,j)
    Csrc(i,j) = c1
    call shr_map_1dto2d(map%idst(n),nid,njd,i,j)
    Cdst(i,j) = Cdst(i,j) + map%wgts(n)
  enddo

  ic1 = 0
  ic2 = 0
  ic3 = 0
  ic4 = 0
  ic5 = 0
  do j=1,njs
  do i=1,nis
    if (Msrc(i,j) /= 0) then   ! live src pt
      if (abs(Csrc(i,j)-c1) < eps) then
        ic1 = ic1 + 1          ! in use
      else
        ic2 = ic2 + 1          ! not used
      endif
    else                       ! dead src pt
      if (abs(Csrc(i,j)-c1) < eps) then
        ic3 = ic3 + 1          ! in use
      else
        ic5 = ic5 + 1          ! not used
      endif
    endif
  enddo
  enddo
! if (ic3 > 0) error = .true.
  if (debug > 0 .and. s_loglev > 0) then
    write(s_logunit,F01) ' total number of SRC points           : ',nis*njs
    write(s_logunit,F01) ' wgts from SRC TRUE  points; used     : ',ic1
    write(s_logunit,F01) ' wgts from SRC TRUE  points; not used : ',ic2
    write(s_logunit,F01) ' wgts from SRC FALSE points; used     : ',ic3
    write(s_logunit,F01) ' wgts from SRC FALSE points; not used : ',ic5
  endif

  ic1 = 0
  ic2 = 0
  ic3 = 0
  ic4 = 0
  ic5 = 0
  do j=1,njd
  do i=1,nid
    if (Mdst(i,j) /= 0) then   ! wgts should sum to one
      if (abs(Cdst(i,j)-c1) < eps) then
        ic1 = ic1 + 1          ! wgts sum to one
      else
        ic2 = ic2 + 1          ! invalid wgts
      endif
    else                       ! wgts should sum to one or zero
      if (abs(Cdst(i,j)-c1) < eps) then
        ic3 = ic3 + 1          ! wgts sum to one
      elseif (abs(Cdst(i,j)) < eps) then
        ic4 = ic4 + 1          ! wgts sum to zero
      else
        ic5 = ic5 + 1          ! invalid wgts
      endif
    endif
  enddo
  enddo
! if (ic2 > 0) error = .true.
! if (ic5 > 0) error = .true.
  if (debug > 0 .and. s_loglev > 0) then
    write(s_logunit,F01) ' total number of DST points           : ',nid*njd
    write(s_logunit,F01) ' sum wgts for DST TRUE  points; one   : ',ic1
    if (ic2 > 0) then
      write(s_logunit,F03) ' sum wgts for DST TRUE  points; not   : ',ic2,' **-WARNING-**'
    else
      write(s_logunit,F01) ' sum wgts for DST TRUE  points; not   : ',ic2
    endif
    write(s_logunit,F01) ' sum wgts for DST FALSE points; one   : ',ic3
    write(s_logunit,F01) ' sum wgts for DST FALSE points; zero  : ',ic4
    write(s_logunit,F01) ' sum wgts for DST FALSE points; not   : ',ic5
  endif

  deallocate(Csrc)
  deallocate(Cdst)

  if (error) call shr_map_abort(subName//' ERROR invalid weights')

end subroutine shr_map_checkWgts_global

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_map_getWts -- local code that sets weights for a point
!
! !DESCRIPTION:
!     Local code that sets weights for a point.  Executes searches
!     and computes weights.  For bilinear remap for example.
!
! !REMARKS:
!     Assumes Xsrc,Ysrc are regular lat/lon grids, monotonicallly increasing
!        on constant latitude and longitude lines.  
!     Assumes Xdst,Ydst,Xsrc,Ysrc are all either radians or degrees
!
! !REVISION HISTORY:
!     2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_getWts(Xdst,Ydst,Xsrc,Ysrc,pti,ptj,ptw,pnum,units)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

!XXEOP

  real(SHR_KIND_R8)   ,intent(in) :: Xdst,Ydst
  real(SHR_KIND_R8)   ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
  integer(SHR_KIND_IN),intent(out):: pti(:),ptj(:)
  real(SHR_KIND_R8)   ,intent(out):: ptw(:)
  integer(SHR_KIND_IN),intent(out):: pnum
  character(len=*),optional,intent(in) :: units

  !--- local ---
  integer(SHR_KIND_IN) :: isize,jsize   ! array sizes
  integer(SHR_KIND_IN) :: i,j           ! indices
  integer(SHR_KIND_IN) :: n             ! do loop counter
  integer(SHR_KIND_IN) :: il,ir         ! index of i left/mid/right
  integer(SHR_KIND_IN) :: jl,ju         ! index of j lower/mid/upper
  integer(SHR_KIND_IN) :: pmax          ! size of pti,ptj,ptw
  real(SHR_KIND_R8)    :: xsl,xsr       ! value of Xsrc, left/right
  real(SHR_KIND_R8)    :: ysl,ysu       ! value of Ysrc, left/right
  real(SHR_KIND_R8)    :: xd,yd         ! value of Xdst,Ydst
  real(SHR_KIND_R8)    :: dx,dy,dx1,dy1 ! some d_lengths for weights calc
  real(SHR_KIND_R8)    :: csize         ! circle angle/radians
  real(SHR_KIND_R8)    :: rmin,rmax     ! min/max
  real(SHR_KIND_R8)    :: cpole         ! the r8 lat value of the pole
  integer(SHR_KIND_IN) :: pole          ! 0=no, 1=north, 2=south

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_getWts') "
  character(*),parameter :: F00     = "('(shr_map_getWts) ',a) "
  character(*),parameter :: F02     = "('(shr_map_getWts) ',a,4g20.13) "
  character(*),parameter :: F03     = "('(shr_map_getWts) ',a,2g20.13) "
  character(*),parameter :: F04     = "('(shr_map_getWts) ',a,4i8) "
  character(*),parameter :: F05     = "('(shr_map_getWts) ',a,3g20.13) "

!-------------------------------------------------------------------------------

  pmax = size(pti,1)

  !--- is lat/lon degrees or radians?  needed for X wraparound ---
  if (present(units)) then
    if (trim(units) == 'degrees') then
       csize = 360._SHR_KIND_R8
    elseif (trim(units) == 'radians') then
       csize = c2*pi
    else
       call shr_sys_abort(subName//' ERROR in optional units = '//trim(units))
    endif
  else
     csize = 360._SHR_KIND_R8
     if (shr_map_checkRad(Ysrc)) csize = c2*pi
  endif

  isize = size(Xsrc,1)
  jsize = size(Xsrc,2)
  pti = 0
  ptj = 0
  ptw = c0

  cpole = csize/(c2*c2)

  xd = Xdst
  yd = Ydst

  if (yd >  cpole + 1.0e-3 .or. &
      yd < -cpole - 1.0e-3) then
     write(s_logunit,*) trim(subname),' ERROR: yd outside bounds ',yd
     call shr_map_abort(subName//' ERROR yd outside 90 degree bounds')
  endif
  if (yd >  cpole) yd =  cpole
  if (yd < -cpole) yd = -cpole

  call shr_map_find4corners(Xdst,yd,Xsrc,Ysrc,il,ir,jl,ju)

  !--- bilinear ---
  pnum = 4
  pole = 0
  xsl = Xsrc(il,1)
  xsr = Xsrc(ir,1)
  ysl = Ysrc(1,jl)
  ysu = Ysrc(1,ju)

  if (Xdst < Xsrc(1,1) .or. Xdst > Xsrc(isize,1)) then
    xsl = mod(Xsrc(il,1),csize)
    xsr = mod(Xsrc(ir,1),csize)
    xd  = mod(Xdst       ,csize)
    if (xsl > xd) xsl = xsl - csize
    if (xsr < xd) xsr = xsr + csize
  endif

  if (yd > Ysrc(1,jsize)) then
    if (dopole) then
      pnum = isize+2
      pole = 1
    endif
    ysu =  cpole
  elseif (yd < Ysrc(1,1)) then
    if (dopole) then
      pnum = isize+2
      pole = 2
    endif
    ysl = -cpole
  endif

  !--- compute dx1,dy1; distance from src(1) to dst
  dx  = (xsr-xsl)
  dy  = (ysu-ysl)
  dx1 = ( xd-xsl)
  dy1 = ( yd-Ysl)

  if (dx1 > dx .and. dx1-dx < 1.0e-7 ) dx1 = dx
  if (dy1 > dy .and. dy1-dy < 1.0e-7 ) dy1 = dy

  if (dx <= c0 .or. dy <= c0 .or. dx1 > dx .or. dy1 > dy) then
    write(s_logunit,*) ' '
    write(s_logunit,F02) 'ERROR in dx,dy: ',dx1,dx,dy1,dy
    write(s_logunit,F03) '   dst: ',Xdst,Ydst
    write(s_logunit,F04) '   ind: ',il,ir,jl,ju
    write(s_logunit,F02) '   dis: ',dx1,dx,dy1,dy
    write(s_logunit,F05) '   x3 : ',xsl,xd,xsr
    write(s_logunit,F05) '   y3 : ',ysl,yd,ysu
    write(s_logunit,*) ' '
    call shr_map_abort(subName//' ERROR in dx,dy calc')
    stop
    return
  endif

  dx1 = dx1 / dx
  dy1 = dy1 / dy

  if (pnum > pmax) then
    call shr_sys_abort(subName//' ERROR pti not big enough')
  endif

  if (pole == 0) then       ! bilinear

    pti(1) = il
    pti(2) = ir
    pti(3) = il
    pti(4) = ir

    ptj(1) = jl
    ptj(2) = jl
    ptj(3) = ju
    ptj(4) = ju

    ptw(1) = (c1-dx1)*(c1-dy1)
    ptw(2) = (   dx1)*(c1-dy1)
    ptw(3) = (c1-dx1)*(   dy1)
    ptw(4) = (   dx1)*(   dy1)

  elseif (pole == 1) then   ! north pole

    pti(1) = il
    pti(2) = ir

    ptj(1) = jl
    ptj(2) = jl

    ptw(1) = (c1-dx1)*(c1-dy1)
    ptw(2) = (   dx1)*(c1-dy1)

    do n=1,isize
      pti(2+n) = n
      ptj(2+n) = ju
      ptw(2+n) = (dy1)/real(isize,SHR_KIND_R8)
    enddo

  elseif (pole == 2) then   ! south pole

    pti(1) = il
    pti(2) = ir

    ptj(1) = ju
    ptj(2) = ju

    ptw(1) = (c1-dx1)*(   dy1)
    ptw(2) = (   dx1)*(   dy1)

    do n=1,isize
      pti(2+n) = n
      ptj(2+n) = jl
      ptw(2+n) = (c1-dy1)/real(isize,SHR_KIND_R8)
    enddo

  else

    write(s_logunit,F00) ' ERROR illegal pnum situation '
    call shr_map_abort(subName//' ERROR illegal pnum situation')

  endif

end subroutine shr_map_getWts

!===============================================================================

subroutine shr_map_find4corners(Xdst,Ydst,Xsrc,Ysrc,il,ir,jl,ju)

! finds 4 corner points surrounding dst in src
! returns left, right, lower, and upper i and j index

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  real(SHR_KIND_R8)   ,intent(in) :: Xdst,Ydst
  real(SHR_KIND_R8)   ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
  integer(SHR_KIND_IN),intent(out):: il,ir,jl,ju

  !--- local ---
  integer(SHR_KIND_IN) :: isize,jsize
  integer(SHR_KIND_IN) :: im,jm

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_find4corners') "
  character(*),parameter :: F00     = "('(shr_map_find4corners) ',a,2i8) "

!-------------------------------------------------------------------------------

  isize = size(Xsrc,1)
  jsize = size(Xsrc,2)

  if (Xsrc(isize,1) > Xsrc(1,1)) then
    ! increasing Xsrc
    if (Xdst < Xsrc(1,1) .or. Xdst > Xsrc(isize,1)) then
      il = isize
      ir = 1
    else
      !--- find i index where Xsrc(i) <=  Xdst < Xsrc(i+1) ---
      il = 1
      ir = isize
      do while (ir-il > 1)
        im = (ir+il)/2
        if (Xdst >= Xsrc(im,1)) then
          il = im
        else
          ir = im
        endif
      enddo
    endif
  else
    ! decreasing Xsrc
    if (Xdst > Xsrc(1,1) .or. Xdst < Xsrc(isize,1)) then
      il = 1
      ir = isize
    else
      !--- find i index where Xsrc(i) >  Xdst >= Xsrc(i+1) ---
      il = isize
      ir = 1
      do while (il-ir > 1)
        im = (ir+il)/2
        if (Xdst >= Xsrc(im,1)) then
          il = im
        else
          ir = im
        endif
      enddo
    endif
  endif

  if (Ysrc(1,jsize) > Ysrc(1,1)) then
    ! increasing Ysrc
    if (Ydst > Ysrc(1,jsize)) then
      jl = jsize
      ju = jsize
    elseif (Ydst < Ysrc(1,1)) then
      jl = 1
      ju = 1
    else
      !--- find j index where Ysrc(j) <=  Ydst < Ysrc(j+1) ---
      jl = 1
      ju = jsize
      do while (ju-jl > 1)
        jm = (ju+jl)/2
        if (Ydst >= Ysrc(1,jm)) then
          jl = jm
        else
          ju = jm
        endif
      enddo
    endif
  else
    ! decreasing Ysrc
    if (Ydst < Ysrc(1,jsize)) then
      jl = jsize
      ju = jsize
    elseif (Ydst > Ysrc(1,1)) then
      jl = 1
      ju = 1
    else
      !--- find j index where Ysrc(j) <=  Ydst < Ysrc(j+1) ---
      jl = jsize
      ju = 1
      do while (jl-ju > 1)
        jm = (ju+jl)/2
        if (Ydst >= Ysrc(1,jm)) then
          jl = jm
        else
          ju = jm
        endif
      enddo
    endif
  endif

end subroutine shr_map_find4corners

!===============================================================================

subroutine shr_map_findnn(Xdst,Ydst,Xsrc,Ysrc,Msrc,inn,jnn)

! finds point in src nearest to dst, returns inn,jnn src index
! searches using Msrc active points only

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  real(SHR_KIND_R8)   ,intent(in) :: Xdst,Ydst
  real(SHR_KIND_R8)   ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
  integer(SHR_KIND_IN),intent(in) :: Msrc(:,:)
  integer(SHR_KIND_IN),intent(out):: inn,jnn

  !--- local ---
  integer(SHR_KIND_IN) :: isize,jsize
  integer(SHR_KIND_IN) :: i,j
  real(SHR_KIND_R8)    :: dnn,dist

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_findnn') "
  character(*),parameter :: F00     = "('(shr_map_findnn) ',a,2i8) "

!-------------------------------------------------------------------------------

  isize = size(Xsrc,1)
  jsize = size(Xsrc,2)

  inn = -1
  jnn = -1
  dnn = -1._SHR_KIND_R8
  do j=1,jsize
  do i=1,isize
    if (Msrc(i,j) /= 0) then
      dist = shr_map_finddist(Xdst,Ydst,Xsrc(i,j),Ysrc(i,j))
      if (dist < dnn .or. inn < 0) then
        dnn = dist
        inn = i
        jnn = j
      endif
    endif
  enddo
  enddo

end subroutine shr_map_findnn

!===============================================================================

subroutine shr_map_findnnon(dir,Xdst,Ydst,Xsrc,Ysrc,Msrc,inn,jnn)

! finds point in src nearest to dst searching i dir first
! returns inn,jnn src index
! searches using Msrc active points only

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*)        ,intent(in) :: dir
  real(SHR_KIND_R8)   ,intent(in) :: Xdst,Ydst
  real(SHR_KIND_R8)   ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
  integer(SHR_KIND_IN),intent(in) :: Msrc(:,:)
  integer(SHR_KIND_IN),intent(out):: inn,jnn

  !--- local ---
  integer(SHR_KIND_IN) :: isize,jsize
  integer(SHR_KIND_IN) :: il,ir,jl,ju
  integer(SHR_KIND_IN) :: n,i,j
  integer(SHR_KIND_IN) :: is,js
  integer(SHR_KIND_IN) :: i2,j2
  real(SHR_KIND_R8)    :: dnn,dist,ds

  !--- formats ---
  character(*),parameter :: subName = "('shr_map_findnnon') "
  character(*),parameter :: F00     = "('(shr_map_findnnon) ',a,2i8) "

!-------------------------------------------------------------------------------

  isize = size(Xsrc,1)
  jsize = size(Xsrc,2)

  !--- find 4 corner points
  call shr_map_find4corners(Xdst,Ydst,Xsrc,Ysrc,il,ir,jl,ju)

  !--- find closest of 4 corner points to dst, set that to is,js
  is = il
  js = jl
  ds = shr_map_finddist(Xdst,Ydst,Xsrc(il,jl),Ysrc(il,jl))
  dist = shr_map_finddist(Xdst,Ydst,Xsrc(ir,jl),Ysrc(ir,jl))
  if (dist < ds) then
    is = ir
    js = jl
    ds = dist
  endif
  dist = shr_map_finddist(Xdst,Ydst,Xsrc(il,ju),Ysrc(il,ju))
  if (dist < ds) then
    is = il
    js = ju
    ds = dist
  endif
  dist = shr_map_finddist(Xdst,Ydst,Xsrc(ir,ju),Ysrc(ir,ju))
  if (dist < ds) then
    is = ir
    js = ju
    ds = dist
  endif

  inn = -1
  jnn = -1
  dnn = -1._SHR_KIND_R8
  i2 = 0
  j2 = 0

  if (trim(dir) == 'i') then
    !--- search biased over i ---
    do while (inn < 0 .and. j2 < jsize)
    do n=1,2
      if (n == 1) j = min(js + j2,jsize)
      if (n == 2) j = max(js - j2,1)
      do i=1,isize
        if (Msrc(i,j) /= 0) then
          dist = shr_map_finddist(Xdst,Ydst,Xsrc(i,j),Ysrc(i,j))
          if (dist < dnn .or. inn < 0) then
            dnn = dist
            inn = i
            jnn = j
          endif
        endif
      enddo
    enddo
    j2 = j2 + 1
    enddo
  elseif (trim(dir) == 'j') then
    !--- search biased over j ---
    do while (inn < 0 .and. i2 < isize)
    do n=1,2
      if (n == 1) i = min(is + i2,isize)
      if (n == 2) i = max(is - i2,1)
      do j=1,jsize
        if (Msrc(i,j) /= 0) then
          dist = shr_map_finddist(Xdst,Ydst,Xsrc(i,j),Ysrc(i,j))
          if (dist < dnn .or. inn < 0) then
            dnn = dist
            inn = i
            jnn = j
          endif
        endif
      enddo
    enddo
    i2 = i2 + 1
    enddo
  else
    call shr_map_abort(subName//' ERROR illegal dir '//trim(dir))
  endif

end subroutine shr_map_findnnon

!===============================================================================

real(SHR_KIND_R8) function shr_map_finddist(Xdst,Ydst,Xsrc,Ysrc)

! x,y distance computation

  implicit none
  real(SHR_KIND_R8),intent(in) :: Xdst,Ydst,Xsrc,Ysrc
  character(*),parameter :: subName = "('shr_map_finddist') "

!-------------------------------------------------------------------------------

  shr_map_finddist = sqrt((Ydst-Ysrc)**2 + (Xdst-Xsrc)**2)

end function shr_map_finddist

!===============================================================================

logical function shr_map_checkRad(Grid)

! check if grid is rad or degree

  implicit none
  real(SHR_KIND_R8),intent(in) :: Grid(:,:)
  character(*),parameter :: subName = "('shr_map_checkRad') "
  real(SHR_KIND_R8) :: rmin,rmax

!-------------------------------------------------------------------------------

  shr_map_checkRad = .false.
  rmin = minval(Grid)
  rmax = maxval(Grid)
  if ((rmax - rmin) < 1.01_SHR_KIND_R8*c2*pi) shr_map_checkRad = .true.

end function shr_map_checkRad

!===============================================================================

subroutine shr_map_1dto2d(gid,ni,nj,i,j)

! convert from a 1d index system to a 2d index system
! gid is 1d index; ni,nj are 2d grid size; i,j are local 2d index

  implicit none
  integer(SHR_KIND_IN),intent(in) :: gid,ni,nj
  integer(SHR_KIND_IN),intent(out):: i,j
  character(*),parameter :: subName = "('shr_map_1dto2d') "
  character(*),parameter :: F01     = "('(shr_map_1dto2d) ',a,3i8)" 

!-------------------------------------------------------------------------------

  if (gid < 1 .or. gid > ni*nj) then
    write(s_logunit,F01) 'ERROR: illegal gid ',gid,ni,nj
    call shr_map_abort(subName//' ERROR')
  endif
  j = (gid-1)/ni+1
  i = mod(gid-1,ni)+1

end subroutine shr_map_1dto2d

!===============================================================================

subroutine shr_map_2dto1d(gid,ni,nj,i,j)

! convert from a 2d index system to a 1d index system
! gid is 1d index; ni,nj are 2d grid size; i,j are local 2d index

  implicit none
  integer(SHR_KIND_IN),intent(in) :: ni,nj,i,j
  integer(SHR_KIND_IN),intent(out):: gid
  character(*),parameter :: subName = "('shr_map_2dto1d') "
  character(*),parameter :: F01     = "('(shr_map_2dto1d) ',a,4i8)" 

!-------------------------------------------------------------------------------

  if (i < 1 .or. i > ni .or. j < 1 .or. j > nj) then
    write(s_logunit,F01) 'ERROR: illegal i,j ',i,ni,j,nj
    call shr_map_abort(subName//' ERROR')
  endif
  gid = (j-1)*ni + i

end subroutine shr_map_2dto1d

!===============================================================================
!===============================================================================
end module shr_map_mod
  

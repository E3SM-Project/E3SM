!===============================================================================
! SVN $Id: shr_alarm_mod.F90 239 2006-02-08 19:02:33Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk/shr/shr_alarm_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_ncio_mod.F90 --- Module to handle reading and writing of NetCDF data
!
! !DESCRIPTION:
!
!     A module to handle the reading and writing of scalar data to NetCDF files for the
!  purpose of storing scalar netCDF data in either: integer, character, or logical format.
!
! Typical usage:
!
! ! Initialize a descripVar
! call shr_ncio_descripInit( DescripVar(n), Name="T", LongName=                   &
!                            "Globally averaged Temperature of surface",       &
!                            Units="K",   RealR8Data=.true., RealR8Fill=1.e+36 )
! ! Read in the file
! FileType="glob avg T file"
! call shr_ncio_open( NCFileName, MasterTask, FileType, ncId, exists )
! call shr_ncio_descripRead( ncId, nVars, mpicom, MasterTask, DescripVar )
! ! Do any other NetCDF reading on ncID....
!...
! call shr_ncio_close( ncId, MasterTask, NCFilename, type=FileType )
! ! Return the name of the descripVar
! name = shr_ncio_descripName( DescripVar(n) )
! ! Get the data inside the variable
! glob_t = shr_ncio_descripGetRealR8( DescripVar(n) )
! ! Now set the scalar data of the variable
! call shr_ncio_descripSet( DescriptVar(n), Name="T", RealR8Data=glob_avg_t )
! ! Write the scalar data to a file
! FileType="glob avg T file"
! call shr_ncio_open( NCFileName, MasterTask, FileType, ncId, exists )
! call shr_ncio_descripWrite( ncId, nVars, mpicom, MasterTask, exists, DescripVar )
! ! Do any other NetCDF writing on ncID....
!...
! call shr_ncio_close( ncId, MasterTask, NCFilename, FileType )
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Dec-20 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_ncio_mod

! !USES:

   use shr_kind_mod,      only: SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_IN, &
                                SHR_KIND_R8, SHR_KIND_I8, &
                                SHR_KIND_R8
   use shr_sys_mod,       only: shr_sys_flush, shr_sys_getenv
   use shr_const_mod,     only: shr_const_spval
   use netcdf

   implicit none

   private  ! default private

! !PUBLIC TYPES:

   public :: shr_ncio_descripType     ! NetCDF data description type

! !PUBLIC MEMBER FUNCTIONS

   public :: shr_ncio_descripSetDefault ! Set list of descrip variables to default settings
   public :: shr_ncio_descripInit       ! Setup variable and variable Names and types
   public :: shr_ncio_descripPutData    ! Set a variable with data
   public :: shr_ncio_descripName       ! Get the Name from the ncio variable
   public :: shr_ncio_descripGetString  ! Get the string value from the ncio variable
   public :: shr_ncio_descripGetInteger ! Get the integer value from the ncio variable
   public :: shr_ncio_descripGetLogical ! Get the logical value from the ncio variable
   public :: shr_ncio_descripGetRealR8  ! Get the real value from the ncio variable
   public :: shr_ncio_descripRead       ! Read in file
   public :: shr_ncio_descripWrite      ! Write out to file
   public :: shr_ncio_open              ! Open NetCDF file for reading/writing/appending
   public :: shr_ncio_close             ! Close a NetCDF file
   public :: shr_ncio_setDebug          ! Set debug level of printing
   public :: shr_ncio_setAbort          ! Set flag if should abort on error or not

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

! Private member functions:

   private :: shr_ncio_logical2Int ! Convert logical into an integer
   private :: shr_ncio_int2Logical ! Convert integer back into a logical
   private :: shr_ncio_abort       ! Local abort method

   !--- Private data to use for data that is not specified yet ---
   integer,           parameter :: shr_ncio_integerFill = -99999
   character(len=*),  parameter :: shr_ncio_stringFill  = '----NOTSET----'
   real(SHR_KIND_R8), parameter :: shr_ncio_realR8Fill  = shr_const_spval
   integer,           parameter :: maxDims              = 1

   !--- Public NetCDF description type ---
   type shr_ncio_descripType
       private    ! Opaque type
       character(SHR_KIND_CS) :: Name = shr_ncio_stringFill        ! Variable Name
       integer                :: XType = shr_ncio_integerFill      ! Variable type
       integer                :: id = shr_ncio_integerFill         ! variable id number
       integer                :: nDims                             ! number of dimensions
       integer                :: DimSizes(maxDims)                 ! Sizes of dimensions
       character(SHR_KIND_CS) :: DimNames(maxDims)                 ! Names of dimensions
       character(SHR_KIND_CS) :: Units = shr_ncio_stringFill       ! Units of variable
       character(SHR_KIND_CL) :: StringData = shr_ncio_stringFill  ! String scalar to write
       character(SHR_KIND_CL) :: ListDescrips = shr_ncio_stringFill! List of descriptions
       integer                :: IntegerData = shr_ncio_integerFill! Integer scalar to write
       real(SHR_KIND_R8)      :: RealR8Data  = shr_ncio_realR8Fill ! Real scalar to write
       integer                :: IntegerFill = shr_ncio_integerFill! Integer Fill-value
       real(SHR_KIND_R8)      :: RealR8Fill  = shr_ncio_realR8Fill ! Real Fill-value
       logical                :: LogicalData                       ! Logical scalar to write
       integer, pointer       :: ListValues(:)                     ! List values
       character(SHR_KIND_CL) :: LongName = shr_ncio_stringFill    ! Long Name of variable
   end type shr_ncio_descripType

   !--- Private data to signify which type of data will be written ---
   integer, parameter :: shr_ncio_integerDataValue = 1   ! Data is integer
   integer, parameter :: shr_ncio_stringDataValue  = 2   ! Data is string
   integer, parameter :: shr_ncio_logicalDataValue = 3   ! Data is logical
   integer, parameter :: shr_ncio_realR8DataValue  = 4   ! Data is real r8
   logical, save :: doAbort = .true.                     ! If abort on error or not
   integer, save :: debugLevel = 1                       ! Debug level
   character(len=*), parameter :: shrCharacterDimName = "shr_character"

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripSetDefault -- Set ncio_descrip variables to defaults
!
! !DESCRIPTION:
!
! Set an array of shr_ncio_descrip variables to the default settings. This
! allows you to re-initialize them to different values if you so choose, It's
! also important on compilers that do NOT recognize structure initialization
! (as above) such as on PGI.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_ncio_descripSetDefault( NVars, DescripVars )

! !USES:

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,                    intent(IN)    :: nVars              ! Number of variables
  type(shr_ncio_descripType), intent(INOUT) :: DescripVars(nVars) ! Output description variables

!EOP

   !----- local -----
   character(len=*), parameter :: subName = '(shr_ncio_descripSetDefault) '
   integer                     :: i      ! index

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   do i = 1, NVars
       DescripVars(i)%Name         = shr_ncio_stringFill    ! Variable Name
       DescripVars(i)%XType        = shr_ncio_integerFill   ! Variable type
       DescripVars(i)%id           = shr_ncio_integerFill   ! variable id number
       DescripVars(i)%Units        = shr_ncio_stringFill    ! Units of variable
       DescripVars(i)%StringData   = shr_ncio_stringFill    ! String scalar to write
       DescripVars(i)%ListDescrips = shr_ncio_stringFill    ! List of descriptions
       DescripVars(i)%IntegerData  = shr_ncio_integerFill   ! Integer scalar to write
       DescripVars(i)%RealR8Data   = shr_ncio_realR8Fill    ! Real scalar to write
       DescripVars(i)%IntegerFill  = shr_ncio_integerFill   ! Integer Fill-value
       DescripVars(i)%RealR8Fill   = shr_ncio_realR8Fill    ! Real Fill-value
       DescripVars(i)%LongName     = shr_ncio_stringFill    ! Long Name of variable
   end do

END SUBROUTINE shr_ncio_descripSetDefault

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripInit -- Set ncio_descrip variable Names
!
! !DESCRIPTION:
!
! Initial setup of a ncio descrip variable. Set logical data flags either to: 
!   StringData, ! IntegerData, LogicalData, or RealR8Data. Can optionally set 
! the LongName, Units, and fillvalue (only for RealR8 or Integer data). Note, 
! if you select IntegerData for the data-type, you can only set the IntergerFill 
! value, NOT the RealR8Fill value and vica-versa. If you pick IntegerData
! as the data-type you can also select ListDescrips and ListIntValues with 
! the list of valid values and a description for the meaning of each of those 
! values. In this case, attributes giving the values and the matching 
! description of that value will be added to the output NetCDF file.
! For example, if you have integer data to describe the vegetation type, with 
! a list of 12 different vegetation types, you give a colen delimited list to 
! describe each type, as well as ListIntValues to give the corresponding value 
! that goes with each description in the list.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_ncio_descripInit( DescripVar, Name, LongName, Units,    &
                                 StringData, IntegerData, LogicalData, &
                                 RealR8Data, IntegerFill,              &
                                 RealR8Fill, ListDescrips,             &
                                 ListIntValues, nDims, dimSizes )

! !USES:

  use shr_string_mod,    only: shr_string_listIsValid, shr_string_listGetNum

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(INOUT) :: DescripVar! Output description variable
  character(len=*), intent(IN) :: Name                   ! Short name of data
  character(len=*), intent(IN), optional :: LongName     ! Longname of data
  character(len=*), intent(IN), optional :: Units        ! Physical units of data
  logical, intent(IN), optional, target :: StringData    ! Flag if data is string-type
  logical, intent(IN), optional, target :: IntegerData   ! Flag if data is integer-type
  logical, intent(IN), optional, target :: LogicalData   ! Flag if data is logical-type
  logical, intent(IN), optional, target :: RealR8Data    ! Flag if data is string-type
  real(SHR_KIND_R8), intent(IN), optional :: RealR8Fill  ! Fill value for real-R8 data
  integer, intent(IN), optional, target :: IntegerFill   ! Fill value for integer data
  ! --- For integer-list, description of each value as colen delimited list ---
  character(len=*), intent(IN), optional :: ListDescrips 
  ! --- For integer-list, integer values that go with each description ----
  integer, intent(IN), optional :: ListIntValues(:)      
  integer, intent(IN), optional :: nDims                 ! Number of dimensions of variable
  integer, intent(IN), optional :: dimSizes(:)           ! Sizes of each dimension

!EOP

   !----- local -----
   character(len=*), parameter :: subName = '(shr_ncio_descripInit) '
   logical :: NotSet        ! Type not set yet or not
   integer :: nList         ! Number of values in list

!-------------------------------------------------------------------------------
! Notes:
!
! Currently only implemented for scalar data, nDims and dimSizes are not fully 
! tested.
!-------------------------------------------------------------------------------

  ! --- Check if setting more than one data-type ------
  if ( (PRESENT(StringData)  .and. PRESENT(IntegerData)) .or. &
       (PRESENT(StringData)  .and. PRESENT(LogicalData)) .or. &
       (PRESENT(StringData)  .and. PRESENT(RealR8Data )) .or. &
       (PRESENT(IntegerData) .and. PRESENT(LogicalData)) .or. &
       (PRESENT(IntegerData) .and. PRESENT(RealR8Data )) .or. &
       (PRESENT(LogicalData) .and. PRESENT(RealR8Data )) )then
     call shr_ncio_abort( subName//': can not set more than one variable type' )
  end if
  if ( PRESENT(nDims) .or. PRESENT(dimSizes) )then
     call shr_ncio_abort( subName//': currently nDims and dimSizes are NOT implemented' )
  end if
  if (  PRESENT(nDims) .and. .not. PRESENT(dimSizes)  &
  .or.  (.not. PRESENT(nDims) .and. PRESENT(dimSizes) ) )then
     call shr_ncio_abort( subName//': if nDims set so must also dimSizes' )
  else if ( PRESENT(nDims) .and. PRESENT(dimSizes) )then
     if ( nDims /= size(dimSizes) )then
        call shr_ncio_abort( subName//': if nDims NOT set to same dimensionality as dimSizes' )
     end if
     if ( nDims > maxDims )then
        call shr_ncio_abort( subName//': nDims > size of max dimensions' )
     end if
     DescripVar%nDims       = nDims
     DescripVar%dimSizes(:) = dimSizes(:nDims)
  end if

  ! --- Name, check if name already defined ------
  if ( DescripVar%Name == shr_ncio_stringFill )then
     DescripVar%Name = Name
  ! --- If name already defined, but you are asking for a different name -- abort
  else if ( trim(Name) /= trim(DescripVar%Name) )then
     call shr_ncio_abort( subName//': bad Name: '//trim(Name)// &
                         ' sent to variable already defined as:'// &
                         trim(DescripVar%Name) )
  ! --- If name is the same, as previously defined name -- silently continue ------
  end if
  ! --- Long-name and units attributes ------
  if ( PRESENT(LongName) ) DescripVar%LongName = LongName
  if ( PRESENT(Units)     ) DescripVar%Units   = Units
  NotSet = .true.
  ! --- For String Data type ------
  if ( PRESENT(StringData) )then
     NotSet = .false.
     DescripVar%XType      = shr_ncio_stringDataValue
     if ( .not. PRESENT(nDims) )then
        DescripVar%nDims       = 1
        DescripVar%dimSizes(1) = len(DescripVar%StringData)
     end if
     DescripVar%dimNames(1) = shrCharacterDimName
     if ( .not. PRESENT(Units) ) DescripVar%Units = "string"
  end if
  ! --- For Integer Data type ------
  if ( PRESENT(IntegerData) )then
     if ( .not. NotSet ) call shr_ncio_abort( subName//': trying to define '// &
                                             'to more than one data-type' )
     NotSet = .false.
     DescripVar%XType       = shr_ncio_integerDataValue
     if ( .not. PRESENT(nDims) )then
        DescripVar%nDims       = 0
     end if
     if ( PRESENT(IntegerFill) )then
        DescripVar%IntegerFill = IntegerFill
        DescripVar%IntegerData = IntegerFill
     end if
     ! --- For an integer list ------
     if ( PRESENT(ListDescrips) )then
        if ( .not. PRESENT(ListIntValues) )then
           call shr_ncio_abort( subName//': setting ListDescrips without '// &
                               'setting ListIntValues' )
        end if
        if ( .not. shr_string_listIsValid( ListDescrips ) )then
           call shr_ncio_abort( subName//': ListDescrips is not a valid '// &
                               'list of descriptions' )
        end if
        nList = shr_string_listGetNum( ListDescrips )
        if ( size(ListIntValues) /= nList )then
           call shr_ncio_abort( subName//': number of list descriptions '// &
                      'inconsistent with number of list integer values' )
        end if
        allocate( DescripVar%ListValues(nList) )
        DescripVar%ListDescrips  = ListDescrips
        DescripVar%ListValues(:) = ListIntValues(:)
     else
        if ( PRESENT(ListIntValues) )then
           call shr_ncio_abort( subName//': setting ListIntValues '// &
                               'without setting ListDescrips' )
        end if
     end if
  else
     if ( PRESENT(IntegerFill) )then
        call shr_ncio_abort( subName//': setting integer FillValue '// &
                            'without setting IntegerData' )
     end if
     if ( PRESENT(ListDescrips) .or. PRESENT(ListIntValues) )then
        call shr_ncio_abort( subName//': setting ListDescrips or '// &
                     'ListIntValues without setting IntegerData' )
     end if
  end if
  ! --- Logical data ------
  if ( PRESENT(LogicalData) )then
     if ( .not. NotSet ) call shr_ncio_abort( subName//': trying '// &
                            'to define to more than one data-type' )
     NotSet = .false.
     DescripVar%XType       = shr_ncio_logicalDataValue
     if ( .not. PRESENT(nDims) )then
        DescripVar%nDims       = 0
     end if
     if ( .not. PRESENT(Units) ) DescripVar%Units = &
                                    "logical flag (0=false)"
  end if
  ! --- Real-R8 data ------
  if ( PRESENT(RealR8Data) )then
     if ( .not. NotSet ) call shr_ncio_abort( &
         subName//': trying to define to more than one data-type' )
     NotSet = .false.
     DescripVar%XType      = shr_ncio_realR8DataValue
     if ( .not. PRESENT(nDims) )then
        DescripVar%nDims       = 0
     end if
     if ( PRESENT(RealR8Fill) )then
        DescripVar%RealR8Fill = RealR8Fill
        DescripVar%RealR8Data = RealR8Fill
     end if
  else
     if ( PRESENT(RealR8Fill) )then
        call shr_ncio_abort( subName//': setting realr8 FillValue '// &
                            'without setting RealR8Data' )
     end if
  end if
  if ( NotSet ) call shr_ncio_abort( subName//': called without giving '// &
                                    'a value to set' )

END SUBROUTINE shr_ncio_descripInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripPutData -- Set a ncio_descrip variable type
!
! !DESCRIPTION:
!
! Set the data of a ncio descrip variable. Data is set to either the string, 
! integer, logical or Real-R8 data value input. Can NOT specify more than one 
! data-type. Check that name agrees with name on variable type, or else abort.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_ncio_descripPutData( DescripVar, Name, StringData, &
                                    IntegerData, LogicalData, RealR8Data )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(INOUT) :: DescripVar      ! descrip variable to set
  character(len=*),           intent(IN)    :: Name            ! Name of this variable
  character(len=*), optional, intent(IN), target :: StringData ! String-data value
  integer, target,  optional, intent(IN)    :: IntegerData     ! Integer-data value
  logical, target,  optional, intent(IN)    :: LogicalData     ! Logical-data value
  real(SHR_KIND_R8), optional,intent(IN), target :: RealR8Data ! Real-R8-data value

!EOP

  !----- local -----
   character(len=*), parameter :: subName = '(shr_ncio_descripPutData) '
   logical :: NotSet      ! Flag if data is not set yet

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Check if too many data values are set (can only be one!) ----- 
  if ( (PRESENT(StringData)  .and. PRESENT(IntegerData)) .or. &
       (PRESENT(StringData)  .and. PRESENT(LogicalData)) .or. &
       (PRESENT(StringData)  .and. PRESENT(RealR8Data )) .or. &
       (PRESENT(IntegerData) .and. PRESENT(LogicalData)) .or. &
       (PRESENT(IntegerData) .and. PRESENT(RealR8Data )) .or. &
       (PRESENT(LogicalData) .and. PRESENT(RealR8Data )) )then
     call shr_ncio_abort( subName//': can not set more than one variable type' )
  end if
  ! --- Make sure name is set and agrees with the input ----- 
  if ( DescripVar%Name == shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': descrip variable has not been defined '// &
                         'yet with shr_ncio_descripInit' )
  else if ( trim(Name) /= trim(DescripVar%Name) )then
     call shr_ncio_abort( subName//': bad Name: '//trim(Name)// &
                         ' sent to variable already defined as:'// &
                         trim(DescripVar%Name) )
  end if
  NotSet = .true.
  !-----------------------------------------------------------------------------
  ! ----------- Test for each type of data ------------
  ! Make sure data sent in is consistent with expected type
  ! And track if data was set, so can abort if no data was set
  !-----------------------------------------------------------------------------
  ! --- If string data -------
  if ( PRESENT(StringData) )then
     NotSet = .false.
     if ( DescripVar%XType /= shr_ncio_stringDataValue )then
        call shr_ncio_abort( subName//': Setting descrip variable to string '// &
                            'which is wrong variable type' )
     end if
     if ( len_trim(StringData) > len(DescripVar%StringData) )then
        call shr_ncio_abort( subName//': Length of input string data longer '// &
                            'than storage size of DescripVar type' )
     end if
     DescripVar%StringData(1:len(DescripVar%StringData)) = ' '
     DescripVar%StringData(:len_trim(StringData)) = trim(StringData)
  end if
  ! --- If integer data -------
  if ( PRESENT(IntegerData) )then
     NotSet = .false.
     if ( DescripVar%XType /= shr_ncio_integerDataValue )then
        call shr_ncio_abort( subName//': Setting descrip variable to integer '// &
                            'which is wrong variable type' )
     end if
     DescripVar%IntegerData = IntegerData
  end if
  ! --- If logical data -------
  if ( PRESENT(LogicalData) )then
     NotSet = .false.
     if ( DescripVar%XType /= shr_ncio_logicalDataValue )then
        call shr_ncio_abort( subName//': Setting descrip variable to logical '// &
                            'which is wrong variable type' )
     end if
     DescripVar%LogicalData = LogicalData
  end if
  ! --- If real R8 data -------
  if ( PRESENT(RealR8Data) )then
     NotSet = .false.
     if ( DescripVar%XType /= shr_ncio_realR8DataValue )then
        call shr_ncio_abort( subName//': Setting descrip variable to realr8 '// &
                            'which is wrong variable type' )
     end if
     DescripVar%RealR8Data = RealR8Data
  end if
  ! --- If no data was set abort with an error -------
  if ( NotSet ) call shr_ncio_abort( subName//': called without giving a '// &
                                    'value to set' )

END SUBROUTINE shr_ncio_descripPutData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripName -- Return the Name of this ncio variable
!
! !DESCRIPTION:
!
! Returns Name of the variable
!
! !INTERFACE: ------------------------------------------------------------------

FUNCTION shr_ncio_descripName( DescripVar )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(IN) :: DescripVar
  character(len=SHR_KIND_CS) :: shr_ncio_descripName

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  shr_ncio_descripName = trim(DescripVar%Name)

END FUNCTION shr_ncio_descripName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripGetString -- Return the string value of this ncio variable
!
! !DESCRIPTION:
!
! Returns the string value from this variable
!
! !INTERFACE: ------------------------------------------------------------------

FUNCTION shr_ncio_descripGetString( DescripVar )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(IN) :: DescripVar    ! NetCDF description variable
  character(len=SHR_KIND_CL) :: shr_ncio_descripGetString ! Returned string value of data

!EOP

  !----- local -----
  character(len=*), parameter :: subName = '(shr_ncio_descripGetString) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- If this structure hasn't been initialized yet -- abort with an error -------
  if ( DescripVar%Name == shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': input shr_ncio description variable '// &
                         'has not been initializated yet: ' )
  end if
  ! --- If this isn't a string type of variable -- abort with an error -------
  if ( DescripVar%XType /= shr_ncio_stringDataValue )then
     call shr_ncio_abort( subName//': trying to get a string from a '// &
                         'different variable type' )
  end if
  shr_ncio_descripGetString(:) = ' '
  shr_ncio_descripGetString(1:len_trim(DescripVar%StringData)) = &
                      trim(DescripVar%StringData)
  ! --- If data not set yet abort with an error -------
  if ( shr_ncio_descripGetString== shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': Returned string has not been set yet' )
  end if

END FUNCTION shr_ncio_descripGetString

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripGetInteger -- Return the integer value of this ncio variable
!
! !DESCRIPTION:
!
! Returns the integer value from this variable
!
! !INTERFACE: ------------------------------------------------------------------

integer FUNCTION shr_ncio_descripGetInteger( DescripVar )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(IN) :: DescripVar ! Input NetCDF description variable

!EOP

   !----- local -----
  character(len=*), parameter :: subName = '(shr_ncio_descripGetInteger) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- If this structure hasn't been initialized yet -- abort with an error -------
  if ( DescripVar%Name == shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': input shr_ncio description '// &
                         'variable has not been initializated yet: ' )
  end if
  ! --- If this isn't a integer  type of variable -- abort with an error -------
  if ( DescripVar%XType /= shr_ncio_integerDataValue )then
     call shr_ncio_abort( subName//': trying to get an integer from a '// &
                         'different variable type' )
  end if
  ! --- If data not set yet abort with an error -------
  shr_ncio_descripGetInteger = DescripVar%IntegerData
  if ( shr_ncio_descripGetInteger == shr_ncio_integerFill  )then
     call shr_ncio_abort( subName//': Returned integer has not been set yet' )
  end if

END FUNCTION shr_ncio_descripGetInteger

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripGetLogical -- Return the logical  value of this ncio variable
!
! !DESCRIPTION:
!
! Returns the logical value from this variable
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_ncio_descripGetLogical( DescripVar )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(IN) :: DescripVar ! Input NetCDF description variable

!EOP

  !----- local -----
  character(len=*), parameter :: subName = '(shr_ncio_descripGetLogical) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- If this structure hasn't been initialized yet -- abort with an error ------
  if ( DescripVar%Name == shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': input shr_ncio description '// &
                         'variable has not been initializated yet: ' )
  end if
  ! --- If this isn't a logical type of variable -- abort with an error ------
  if ( DescripVar%XType /= shr_ncio_logicalDataValue )then
     call shr_ncio_abort( subName//': trying to get an logical from a '// &
                         'different variable type' )
  end if
  shr_ncio_descripGetLogical = DescripVar%LogicalData

END FUNCTION shr_ncio_descripGetLogical

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripGetRealR8 -- Return the realr8 value of this ncio variable
!
! !DESCRIPTION:
!
! Returns the realr8 value from this variable
!
! !INTERFACE: ------------------------------------------------------------------

real(SHR_KIND_R8) FUNCTION shr_ncio_descripGetRealR8( DescripVar )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_ncio_descripType), intent(IN) :: DescripVar ! Input NetCDF description variable

!EOP

  !----- local -----
  character(len=*), parameter :: subName = '(shr_ncio_descripGetRealR8) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- If this structure hasn't been initialized yet -- abort with an error ------
  if ( DescripVar%Name == shr_ncio_stringFill )then
     call shr_ncio_abort( subName//': input shr_ncio description variable '// &
                         'has not been initializated yet: ' )
  end if
  ! --- If this isn't a real r8 type of variable -- abort with an error ------
  if ( DescripVar%XType /= shr_ncio_realR8DataValue )then
     call shr_ncio_abort( subName//': trying to get an real from a '// &
                         'different variable type' )
  end if
  ! --- If data not set yet abort with an error ------
  shr_ncio_descripGetRealR8 = DescripVar%RealR8Data
  if ( shr_ncio_descripGetRealR8 == shr_ncio_realR8Fill  )then
     call shr_ncio_abort( subName//': Returned realr8  has not been set yet' )
  end if

END FUNCTION shr_ncio_descripGetRealR8

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripRead  -- Read in netCDF restart file
!
! !DESCRIPTION:
!
! Read in restart file information from netCDF input restart file
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_ncio_descripRead( ncId, nVars, prefix, mpicom, MasterTask, var )

! !USES:

  use shr_string_mod, only: shr_string_lastIndex

  use shr_mpi_mod,    only: shr_mpi_bcast

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,          intent(IN)           :: ncId       ! NetCDF file unit id
  integer,          intent(IN)           :: nVars      ! Number of variables
  character(len=*), intent(IN), optional :: prefix     ! Prefix in front of variable Names
  integer,          intent(IN), optional :: MPICom     ! MPI communicator
  logical,          intent(IN), optional :: MasterTask ! If Master PE task or not
  type(shr_ncio_descripType), intent(INOUT) :: var(:)  ! NetCDF description structure
  !----- local -----
  character(len=*), parameter :: subName = "(shr_ncio_descripRead) "
  logical :: MasterTask2                  ! If Master PE task or not
  integer :: rCode                        ! Return code
  integer :: nDims                        ! Number of dimensions on file
  integer :: VarDimIds(maxDims)           ! Array of dimension id's for this variable
  integer :: dimIds(maxDims)              ! Dimension id number
  integer :: type                         ! Type of variable
  integer :: i, dim                       ! Index
  integer :: nChars                       ! Number of characters in string
  integer :: n                            ! Character index
  character(SHR_KIND_CL) :: StringData    ! String value read in
  integer                :: IntegerData   ! Integer value read in
  real(SHR_KIND_R8)      :: RealR8Data    ! Real value to read in
  character(len=SHR_KIND_CS) :: prefixUse ! Prefix in front of variable Names

!-------------------------------------------------------------------------------
! Notes:
!   Currently this interface can NOT be used to read in non-scalar data.
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if

  if ( present(prefix) )then
     prefixUse = prefix
  else
     prefixUse = ""
  end if
  if ( MasterTask2 )then
     !-------------------------------------------------------------------------
     ! Loop through variables
     !-------------------------------------------------------------------------
     do i = 1, nVars
        if ( debugLevel > 1 ) write(6,*) 'Read variable: ', trim(var(i)%Name)
        ! --- If name not set, hasn't been initialized -- abort with an error -----
        if ( trim(var(i)%Name) == shr_ncio_stringFill )then
           write(6,'(a,i3,a,i3)') 'variable number = ', i, ' of ', nVars
           call shr_ncio_abort( subName//': variable name not defined -- '// &
                               'DescripVarSet not called' )
        end if
        ! --- Get variable ID -----
        rcode = nf90_inq_varid(ncId,trim(prefixUse)//trim(var(i)%Name),var(i)%id )
        call shr_ncio_abort( subName//': variable '// trim(var(i)%Name)//' not found', rcode )
        ! --- Get type of variable -----
        rcode = nf90_inquire_variable(ncId, var(i)%id, nDims=nDims, &
                                      XType=type )
        call shr_ncio_abort( subName// ': error on inquiry of '//var(i)%Name, rcode )
        if (nDims /= var(i)%nDims )then
           write(6,'(a,a,a,i4,a,i4)') 'Number of dimensions for variable', &
           trim(var(i)%name), ' :', nDims, ' expected:', var(i)%nDims
           call shr_ncio_abort( subName//': '//var(i)%Name// &
                               ' dimension size different than expected' )
        end if
        ! --- Read in the dimension names for this variable ----------------------
        do dim = 1, nDims
           rcode = nf90_inq_dimId(ncId, var(i)%dimNames(dim), dimIds(dim) )
           call shr_ncio_abort( subName// ': error gettting dimension', rcode )
        end do
        ! --- Check that dimension ids are correct -------------------------------
        if ( nDims > 0 )then
           rcode = nf90_inquire_variable(ncId, var(i)%id, dimIds=VarDimIds )
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                               ': error on inquiry of dimIds'//var(i)%Name )
           if ( any(dimIds /= VarDimIds) ) call shr_ncio_abort( subName// &
                                ': dimIds not correct' )
        end if
        ! --- If variable string data type -----
        if (var(i)%XType == shr_ncio_stringDataValue )then
           if ( type /= NF90_CHAR ) call shr_ncio_abort( subName//': '// &
                                        var(i)%Name//' not proper type' )
           rcode = nf90_get_att(ncId, var(i)%id, "nChars", nChars)
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                             ': error getting nChars from : '//var(i)%Name )
           rcode = nf90_get_var( ncId, varid=var(i)%id, values=StringData, &
                                 count=(/nChars/) )
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                      ': error on get of '//var(i)%Name )
           var(i)%StringData(:) = ' '
           var(i)%StringData(:nChars) = StringData(:nChars)
        ! --- If integer data type --------
        else if (var(i)%XType == shr_ncio_integerDataValue )then
           if ( type /= NF90_INT ) call shr_ncio_abort( subName//': '// &
                                      var(i)%Name//' not proper type' )
           rcode = nf90_get_var( ncId, var(i)%id, IntegerData)
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                          ': error on get of '//var(i)%Name )
           var(i)%IntegerData = IntegerData
        ! --- If real r8 data type --------
        else if (var(i)%XType == shr_ncio_realR8DataValue )then
           if ( type /= NF90_DOUBLE) call shr_ncio_abort( subName//': '// &
                                        var(i)%Name//' not proper type' )
           rcode = nf90_get_var( ncId, var(i)%id, RealR8Data )
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                      ': error on get of '//var(i)%Name )
           var(i)%RealR8Data = RealR8Data
        ! --- If logical data type --------
        else if (var(i)%XType == shr_ncio_logicalDataValue )then
           if ( type /= NF90_INT ) call shr_ncio_abort( subName//': '// &
                                       var(i)%Name//' not proper type' )
           rcode = nf90_get_var( ncId, var(i)%id, IntegerData)
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                     ': error on get of '//var(i)%Name )
           var(i)%LogicalData = shr_ncio_logical2Int( IntegerData )
        ! --- Otherwise data type is unknown abort with an error --------
        else
           call shr_ncio_abort( subName//': only integer, logical, '// &
                 'real-r8 or character can not be read for variable ' &
                  //var(i)%Name )
        end if
     
     end do
  end if

  ! --- If MPI broadcast to other MPI tasks -----------------------------------
  if ( present(mpicom) )then
     if ( debugLevel > 1 ) write(6,*) 'Broadcast variables to all tasks'
     do i = 1, nVars
        if (      var(i)%XType == shr_ncio_stringDataValue )then
           call shr_mpi_bcast( var(i)%StringData, mpicom )
        else if ( var(i)%XType == shr_ncio_integerDataValue )then
           call shr_mpi_bcast( var(i)%IntegerData, mpicom )
        else if ( var(i)%XType == shr_ncio_realR8DataValue )then
           call shr_mpi_bcast( var(i)%RealR8Data, mpicom )
        else if ( var(i)%XType == shr_ncio_logicalDataValue )then
           call shr_mpi_bcast( var(i)%LogicalData, mpicom )
        else
           call shr_ncio_abort( subName//': invalid ncio type: ' )
        end if
     end do
  end if


END SUBROUTINE shr_ncio_descripRead

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_descripWrite -- Write out netCDF restart file information
!
! !DESCRIPTION:
!
! Write out restart file information to netCDF restart file
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_ncio_descripWrite( ncId, nVars, prefix, mpicom, MasterTask, exists, &
                           var )

! !USES:

  use shr_string_mod,    only: shr_string_lastindex, shr_string_listGetName, &
                               shr_string_listGetNum

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,          intent(IN)           :: ncId      ! NetCDF file id
  integer,          intent(IN)           :: nVars     ! Number of variables
  character(len=*), intent(IN), optional :: prefix    ! Prefix in front of variable Names
  integer,          intent(IN), optional :: MPICom    ! MPI communicator
  logical,          intent(IN), optional :: MasterTask! If Master PE task or not
  logical,          intent(IN)           :: exists    ! If datafile exists already or not
  type(shr_ncio_descripType), intent(INOUT) :: var(:) ! NetCDF description structure

!EOP 

  !----- local -----
  character(len=*), parameter :: subName = "(shr_ncio_descripWrite) "
  logical :: MasterTask2                  ! If Master PE task or not
  integer :: rcode                        ! NetCDF return code
  integer :: dimIds(maxDims)              ! Dimension id numbers
  integer :: i                            ! Variable index
  integer :: dim                          ! Dimension index
  logical :: NotSet                       ! If variables set yet or not
  logical :: DimSet                       ! If dimension set yet or not
  integer :: n                            ! character index
  integer :: type                         ! data type to write to file
  integer :: list                         ! List index number
  integer :: nList                        ! number of items in list
  character(len=SHR_KIND_CL) :: name      ! List description name
  character(len=SHR_KIND_CS) :: prefixUse ! Prefix in front of variable Names
  character(len=*), parameter :: F00= &
  "(a,' input: nciD=',i3,' nvars=',i3, ' prefix=',a,' master=',l1,' exists=',l1)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if

  if ( present(prefix) )then
     prefixUse = prefix
  else
     prefixUse = ""
  end if
  !-----------------------------------------------------------------------------
  if ( debugLevel > 2 ) write(6,F00) subName, ncid, nVars, trim(prefixUse), &
                                     MasterTask2, exists
  if ( MasterTask2 )then
     ! --- If file already exists, go to re-define mode -----
     if ( exists )then
         rcode = nf90_redef(ncId)
        call shr_ncio_abort( subName// ': error on redefine output NetCDF file', rcode )
     end if
     !--------------------------------------------------------------------------
     ! Define variables
     !--------------------------------------------------------------------------
     do i = 1, nVars
        NotSet = .false.
        if ( debugLevel > 1 ) write(6,*) 'Write variable: ', trim(var(i)%Name)
        ! --- If file exists, check if given variable is already defined -------
        if ( exists )then
           rcode = nf90_inq_varid(ncId, Name=trim(prefixUse)//trim(var(i)%Name), &
                                  varid=var(i)%id )
           if (rcode == nf90_enotvar )then
              NotSet = .true.
           else if (rcode /= nf90_noerr)then
              call shr_ncio_abort( subName//': error getting variable id' )
           end if
        end if
        ! --- If new file or variable NOT defined, define it now ---------------
        if ( .not. exists .or. NotSet )then
           !--------------------------------------------------------------------
           ! Define dimensions for this variable if not already defined
           !--------------------------------------------------------------------
           do dim = 1, var(i)%nDims
              DimSet = .false.
              rcode = nf90_inq_dimId(ncId, Name=var(i)%dimNames(dim), &
                                     dimId=dimIds(dim))
              if ( rcode == nf90_ebaddim )then
                 DimSet = .true.
              else if (rcode /= nf90_noerr)then
                 call shr_ncio_abort( subName//': error getting correct '// &
                      'dimension id for :'//trim(var(i)%dimNames(dim)) )
              end if
              if ( DimSet )then
                 rcode = nf90_def_dim(ncId, Name=var(i)%dimNames(dim), &
                                      len=SHR_KIND_CL, dimId=dimIds(dim) )
                  call shr_ncio_abort( subName//': error writing dimension', &
                                       rcode )
              end if
           end do
           if ( var(i)%XType == shr_ncio_stringDataValue )then
              type = NF90_CHAR
           else if ( var(i)%XType == shr_ncio_integerDataValue .or. &
                     var(i)%XType == shr_ncio_logicalDataValue )then
              type = NF90_INT
           else if ( var(i)%XType == shr_ncio_realR8DataValue )then
              type = NF90_DOUBLE
           else
              call shr_ncio_abort( subName// &
                            ': error on variable definition: '// &
                            trim(var(i)%Name) )
           end if
           
           if ( var(i)%NDims > 0 )then
              rcode = nf90_def_var( ncId, trim(prefixUse)//trim(var(i)%Name), &
                                    type, dimIds=dimIds, varid=var(i)%id )
           else
              rcode = nf90_def_var( ncId, trim(prefixUse)//trim(var(i)%Name), &
                                    type, var(i)%id )
           end if
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                            ': error on variable definition: '//   &
                            trim(var(i)%Name) )
        end if
        !-----------------------------------------------------------------------
        ! Attributes on variables
        !-----------------------------------------------------------------------
        rcode = nf90_put_att( ncId, var(i)%id, "long_name", var(i)%LongName )
        if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                              ': error on putting LongName '// &
                                              'attribute on variable: '//      &
                                              trim(var(i)%Name) )
        rcode = nf90_put_att( ncId, var(i)%id, "units", var(i)%Units )
        if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                    ': error on putting Units attribute '// &
                                    'on variable: '//trim(var(i)%Name) )
        if ( var(i)%XType == shr_ncio_stringDataValue )then
          rcode = nf90_put_att( ncId, var(i)%id, "nChars", &
                                len_trim(var(i)%StringData) )
          if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                        ': error on putting nChars '// &
                                        'attribute on variable: '//    &
                                        trim(var(i)%Name) )
        end if
        ! --- Fill value for integer data -------
        if ( var(i)%IntegerFill /= shr_ncio_integerFill )then
          rcode = nf90_put_att( ncId, var(i)%id, "_FillValue", &
                                var(i)%IntegerFill )
          if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                        ': error on putting _FillValue '//  &
                                        'integer attribute on variable: '// &
                                         var(i)%Name )
          rcode = nf90_put_att( ncId, var(i)%id, "missing_value", &
                                var(i)%IntegerFill )
          if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                        ': error on putting missing_value'// &
                                        'integer attribute on variable: '//  &
                                        var(i)%Name )
        end if
        ! --- Fill value for real-r8 data -------
        if ( var(i)%RealR8Fill /= shr_ncio_realR8Fill )then
          rcode = nf90_put_att( ncId, var(i)%id, "_FillValue", &
                                var(i)%RealR8Fill )
          if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                        ': error on putting _FillValue '// &
                                        'RealR8 attribute on variable: '// &
                                        var(i)%Name )
          rcode = nf90_put_att( ncId, var(i)%id, "missing_value", &
                                var(i)%RealR8Fill )
          if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                        ': error on putting missing_value'// &
                                        'RealR8 attribute on variable: '//   &
                                        var(i)%Name )
        end if
        ! --- List description attributes -------------------------------------
        if ( var(i)%ListDescrips /= shr_ncio_stringFill )then
           rcode = nf90_put_att( ncId, var(i)%id, "type", "Integer list" )
           if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                  ': error on putting '// &
                                  ' type string attribute on variable: '// &
                                  var(i)%Name )
           nList = shr_string_listGetNum( var(i)%ListDescrips )
           do list = 1, nList
              call shr_string_listGetName( var(i)%ListDescrips, list, &
                                           name, rcode )
              rcode = nf90_put_att( ncId, var(i)%id, trim(name), &
                                    var(i)%ListValues(list) )
              if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                              ': error on putting ' &
                              //name//' integer attribute on variable: '// &
                              var(i)%Name )
           end do
        end if
     end do
     ! --- Take it out of definition mode ------
     rcode = nf90_enddef( ncId )
     if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                           ': error ending definition mode' )
     !--------------------------------------------------------------------------
     ! Write out variables
     !--------------------------------------------------------------------------
     do i = 1, nVars
        if (      var(i)%XType == shr_ncio_stringDataValue )then
           rcode = nf90_put_var( ncId, var(i)%id, trim(var(i)%StringData), &
                                 count=(/len_trim(var(i)%StringData)/) )
        else if ( var(i)%XType == shr_ncio_integerDataValue )then
           rcode = nf90_put_var( ncId, var(i)%id, var(i)%IntegerData )
        else if ( var(i)%XType == shr_ncio_realR8DataValue )then
           rcode = nf90_put_var( ncId, var(i)%id, var(i)%RealR8Data )
        else if ( var(i)%XType == shr_ncio_logicalDataValue )then
           rcode = nf90_put_var( ncId, var(i)%id, shr_ncio_int2Logical( &
                                 var(i)%LogicalData ) )
        else
           call shr_ncio_abort( subName//': invalid ncio type: ' )
        end if
        if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                                 ': error writing variable: '  &
                                 //var(i)%Name )
     end do
  end if
END SUBROUTINE shr_ncio_descripWrite

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_open -- Open a NetCDF file for reading, writing or appending
!
! !DESCRIPTION:
!
! Open NetCDF file for reading, writing or appending. If reading, or appending 
! abort if file does not exist or on error. If appending
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncio_open( NCFileName, MasterTask, FileType, ncId, exists, &
                          writing, appending )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=*), intent(IN) :: NCFileName  ! Full path to input fileName to read
  logical, intent(IN) :: MasterTask           ! If Master PE task or not
  character(len=*), intent(IN) :: FileType    ! Description of file-type
  integer, intent(OUT) :: ncId                ! NetCDF file unit
  logical, intent(OUT) :: exists              ! If file exists or not
  logical, intent(IN), optional :: writing    ! If should open for writing
  logical, intent(IN), optional :: appending  ! If should open for appending

!EOP

  !----- local -----
  character(len=*), parameter :: subName = "(shr_ncio_open) "
  integer :: rCode                       ! Return code
  logical :: writing2                    ! If should open for writing
  logical :: appending2                  ! If should open for appending

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Only open file on MasterTask -----
  if ( MasterTask ) then
     if ( present(writing) )then
        writing2 = writing
     else
        writing2 = .false.
     end if
     if ( present(appending) )then
        appending2 = appending
     else
        appending2 = .false.
     end if
     if ( appending2 .and. writing2 )then
        call shr_ncio_abort( subName//': can NOT set both appending and writing'// &
                            ' option!' )
     end if
   
     if ( debugLevel > 0 ) write(6,*) 'Open NetCDF file FileType: '//trim(FileType)// &
                                      ' Filename: ', trim(NCFileName)
     inquire( file = trim(NCFileName), exist = exists )
     ! If file exists and reading but not appending
     if (      exists .and. (.not. writing2) .and. (.not. appending2) )then
        if ( debugLevel > 1 ) write(6,*) 'File exists open for reading not appending: '
        rCode = nf90_open( NCFileName, nf90_nowrite, ncId )
     ! If file exists and appending
     else if ( exists .and. appending2 )then
        if ( debugLevel > 1 ) write(6,*) 'File exists open for appending: '
        rCode = nf90_open( NCFileName, nf90_noclobber, ncId )
     ! If file exists and writing
     else if ( exists .and. writing2 )then
        if ( debugLevel > 1 ) write(6,*) 'File exists open for writing: '
        rCode = nf90_open( NCFileName, nf90_write, ncId )
     ! If file does NOT exist and writing
     else if ( (.not. exists) .and. writing2 )then
        if ( debugLevel > 1 ) write(6,*) 'File does NOT exist open for writing: '
        rCode = nf90_create( NCFileName, nf90_noclobber, ncId )
     ! If file does NOT exist and want to read or append -- flag an error
     else if ( .not. exists .and. ((.not. writing2) .or. appending2) )then
        call shr_ncio_abort( subName//': input file does not exist -- can '// &
                            'NOT open for reading!' )
     end if
     ! Check error code from above options...
     call shr_ncio_abort( subName//': error opening : '//NCFileName, rcode )
  end if

end subroutine shr_ncio_open

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_close -- Close a NetCDF file
!
! !DESCRIPTION:
!
! Close a NetCDF file
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncio_close( ncId, MasterTask, NCFilename, type )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer, intent(IN) :: ncId                ! NetCDF file unit
  logical, intent(IN) :: MasterTask          ! If Master PE task or not
  character(len=*), optional :: NCFileName   ! Filename to close
  character(len=*), optional :: type         ! Type of file to close

!EOP

  !----- local -----
  character(len=*), parameter :: subName = "(shr_ncio_close) "
  character(len=SHR_KIND_CL) :: FileName ! Filename to close
  character(len=SHR_KIND_CL) :: FileType ! Description of file to close
  integer :: rCode                       ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Only close file on MasterTask -----
  if ( MasterTask )then
     if ( PRESENT(NCFileName) )then
       FileName = NCFileName
     else
       FileName = " "
     end if
     if ( PRESENT(type) )then
       FileType = type
     else
       FileType = " "
     end if
     rCode = nf90_close( ncId )
     if (rcode /= nf90_noerr) call shr_ncio_abort( subName// &
                              ': error closing '//trim(FileType)//' file: '// &
                              trim(FileName) )
  end if
end subroutine shr_ncio_close

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_setDebug -- set debug level
!
! !DESCRIPTION:
!
! Set debug printing level...
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncio_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer, intent(in) :: iflag    ! Flag value to set printing level to

!EOP

  !--- local ---
  character(*),parameter :: subName =   "(shr_ncio_setDebug)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

  debugLevel = iflag

end subroutine shr_ncio_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_setAbort -- set abort level
!
! !DESCRIPTION:
!
! Set abort level...
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncio_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical, intent(in) :: flag    ! Flag value to set abort level to

!EOP

  !--- local ---
  character(*),parameter :: subName =   "(shr_ncio_setAbort)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

  doAbort = flag

end subroutine shr_ncio_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_logical2Int -- Return logical value based on integer input
!
! !DESCRIPTION:
!
! Convert integer to logical. Private to this module.
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_ncio_logical2Int( int_input )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer, intent(IN) :: int_input      ! Integer input value to convert to logical

!EOP

  !----- local -----
   character(len=*), parameter :: subName = '(shr_ncio_logical2Int) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (      int_input == 0 )then
      shr_ncio_logical2Int = .false.
   else if ( int_input == 1 )then
      shr_ncio_logical2Int = .true.
   else
      call shr_ncio_abort( subName//': bad input to shr_ncio_logical2Int: ' )
   end if

END FUNCTION shr_ncio_logical2Int

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_int2Logical -- Return integer value based on logical input
!
! !DESCRIPTION:
!
! Convert logical to integer. Private to this module
!
! !INTERFACE: ------------------------------------------------------------------

integer FUNCTION shr_ncio_int2Logical( log_input )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical, intent(IN) :: log_input   ! input logical value to convert to integer

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( .not. log_input )then
      shr_ncio_int2Logical = 0
   else
      shr_ncio_int2Logical = 1
   end if

END FUNCTION shr_ncio_int2Logical

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncio_abort -- abort on error
!
! !DESCRIPTION:
!
! Private routine to abort on error
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncio_abort( string, rcode)

! !USES:

   use shr_sys_mod,       only: shr_sys_abort

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(IN) :: string  ! abort message
  integer,     optional,intent(IN) :: rcode   ! NetCDF error code

!EOP

  !--- local ---
  character(SHR_KIND_CL) :: lstring
  character(*),parameter :: subName =   "(shr_ncio_abort)"
  character(*),parameter :: F00     = "('(shr_ncio_abort) ',a)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

  if (present(rcode))then
     if ( rcode == nf90_noerr )then
        return
     else
        write(6,'(a,a,i3)') subname, ' : NetCDF error code = ', rcode
     end if
  end if
  lstring = ''
  if (present(string)) lstring = string

  if (doAbort) then
    call shr_sys_abort(lstring)
  else
    write(6,F00) ' no abort:'//trim(lstring)
  endif

end subroutine shr_ncio_abort

!===============================================================================
!===============================================================================

END MODULE shr_ncio_mod

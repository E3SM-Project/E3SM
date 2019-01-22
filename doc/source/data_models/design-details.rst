.. _design-details:

================
 Design Details
================

----------------------
Data Model Performance
----------------------

There are two primary costs associated with strdata: reading data and spatially mapping data.
Time interpolation is relatively cheap in the current implementation. 
As much as possible, redundant operations are minimized.
Fill and mapping weights are generated at initialization and saved. 
The upper and lower bound mapped input data is saved between time steps to reduce mapping costs in cases where data is time interpolated more often than new data is read.
If the input data timestep is relatively small (for example, hourly data as opposed to daily or monthly data) the cost of reading input data can be quite large. 
Also, there can be significant variation in cost of the data model over the coarse of the run, for instance, when new inputdata must be read and interpolated, although it's relatively predictable.
The present implementation doesn't support changing the order of operations, for instance, time interpolating the data before spatial mapping. 
Because the present computations are always linear, changing the order of operations will not fundamentally change the results.
The present order of operations generally minimizes the mapping cost for typical data model use cases.

----------------------
Data Model Limitations
----------------------

There are several limitations in both options and usage within the data models at the present time.
Spatial interpolation can only be performed from a two-dimensional latitude-longitude input grid. 
The target grid can be arbitrary but the source grid must be able to be described by simple one-dimensional lists of longitudes and latitudes, although they don't have to be equally spaced.

----------------------
IO Through Data Models
----------------------

At the present time, data models can only read netcdf data, and IO is handled through either standard netcdf interfaces or through the PIO library using either netcdf or pnetcdf.
If standard netcdf is used, global fields are read and then scattered one field at a time. 
If PIO is used, then data will be read either serially or in parallel in chunks that are approximately the global field size divided by the number of IO tasks.
If pnetcdf is used through PIO, then the pnetcdf library must be included during the build of the model. 
The pnetcdf path and option is hardwired into the ``Macros.make`` file for the specific machine.
To turn on ``pnetcdf`` in the build, make sure the ``Macros.make`` variables ``PNETCDF_PATH``, ``INC_PNETCDF``, and ``LIB_PNETCDF`` are set and that the PIO ``CONFIG_ARGS`` sets the ``PNETCDF_PATH`` argument. 

Beyond just the option of selecting IO with PIO, several namelist variables are available to help optimize PIO IO performance.
Those are **TODO** - list these. 
The total mpi tasks that can be used for IO is limited to the total number of tasks used by the data model.
Often though, using fewer IO tasks results in improved performance. 
In general, [io_root + (num_iotasks-1)*io_stride + 1] has to be less than the total number of data model tasks.
In practice, PIO seems to perform optimally somewhere between the extremes of 1 task and all tasks, and is highly machine and problem dependent.

-------------
Restart Files
-------------
Restart files are generated automatically by the data models based on a flag sent from the driver.
The restart files must meet the CIME naming convention and an ``rpointer`` file is generated at the same time. 
An ``rpointer`` file is a *restart pointer* file which contains the name of the most recently created restart file. 
Normally, if restart files are read, the restart filenames are specified in the ``rpointer`` file. 
Optionally though, there are namelist variables such as ``restfilm`` to specify the restart filenames via namelist. If those namelist variables are set, the ``rpointer`` file will be ignored. 

In most cases, no restart file is required for the data models to restart exactly.
This is because there is no memory between timesteps in many of the data model science modes. 
If a restart file is required, it will be written automatically and then must be used to continue the previous run.

There are separate stream restart files that only exist for performance reasons. 
A stream restart file contains information about the time axis of the input streams. 
This information helps reduce the startup costs associated with reading the input dataset time axis information. 
If a stream restart file is missing, the code will restart without it but may need to reread data from the input data files that would have been stored in the stream restart file. 
This will take extra time but will not impact the results.

.. _data-structures:

---------------
Data Structures
---------------

The data models all use three fundamental routines.

- $CIMEROOT/src/utils/shr_dmodel_mod.F90

- $CIMEROOT/src/utils/shr_stream_mod.F90

- $CIMEROOT/src/utils/shr_strdata.F90

These routines contain three data structures that are leveraged by all the data model code.

The most basic type, ``shr_stream_fileType`` is contained in ``shr_stream_mod.F90`` and specifies basic information related to a given stream file.

.. code-block:: Fortran

   type shr_stream_fileType
      character(SHR_KIND_CL) :: name = shr_stream_file_null	! the file name
      logical                :: haveData = .false.		! has t-coord data been read in?
      integer  (SHR_KIND_IN) :: nt = 0				! size of time dimension
      integer  (SHR_KIND_IN),allocatable :: date(:)		! t-coord date: yyyymmdd
      integer  (SHR_KIND_IN),allocatable :: secs(:)		! t-coord secs: elapsed on date
   end type shr_stream_fileType

The following type, ``shr_stream_streamType`` contains information
that encapsulates the information related to all files specific to a
target stream. These are the list of files found in the ``domainInfo``
and ``fieldInfo`` blocks of the target stream description file (see the overview of the :ref:`stream_description_file`).

.. code-block:: Fortran

   type shr_stream_streamType
      !private                                    ! no public access to internal components
      !--- input data file names and data ---
      logical                   :: init           ! has stream been initialized?
      integer  (SHR_KIND_IN),pointer :: initarr(:) => null()! surrogate for init flag
      integer  (SHR_KIND_IN)    :: nFiles         ! number of data files
      character(SHR_KIND_CS)    :: dataSource     ! meta data identifying data source
      character(SHR_KIND_CL)    :: filePath       ! remote location of data files
      type(shr_stream_fileType), allocatable :: file(:) ! data specific to each file

      !--- specifies how model dates align with data dates ---
      integer(SHR_KIND_IN)      :: yearFirst      ! first year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearLast       ! last  year to use in t-axis (yyyymmdd)
      integer(SHR_KIND_IN)      :: yearAlign      ! align yearFirst with this model year
      integer(SHR_KIND_IN)      :: offset         ! offset in seconds of stream data
      character(SHR_KIND_CS)    :: taxMode        ! cycling option for time axis

      !--- useful for quicker searching ---
      integer(SHR_KIND_IN) :: k_lvd,n_lvd         ! file/sample of least valid date
      logical              :: found_lvd           ! T <=> k_lvd,n_lvd have been set
      integer(SHR_KIND_IN) :: k_gvd,n_gvd         ! file/sample of greatest valid date
      logical              :: found_gvd           ! T <=> k_gvd,n_gvd have been set

      !---- for keeping files open
      logical                 :: fileopen         ! is current file open
      character(SHR_KIND_CL)  :: currfile         ! current filename
      type(file_desc_t)       :: currpioid        ! current pio file desc

      !--- stream data not used by stream module itself ---
      character(SHR_KIND_CXX):: fldListFile       ! field list: file's  field names
      character(SHR_KIND_CXX):: fldListModel      ! field list: model's field names
      character(SHR_KIND_CL) :: domFilePath       ! domain file: file path of domain file
      character(SHR_KIND_CL) :: domFileName       ! domain file: name
      character(SHR_KIND_CS) :: domTvarName       ! domain file: time-dim var name
      character(SHR_KIND_CS) :: domXvarName       ! domain file: x-dim var name
      character(SHR_KIND_CS) :: domYvarName       ! domain file: y-dim var name
      character(SHR_KIND_CS) :: domZvarName       ! domain file: z-dim var name
      character(SHR_KIND_CS) :: domAreaName       ! domain file: area  var name
      character(SHR_KIND_CS) :: domMaskName       ! domain file: mask  var name

      character(SHR_KIND_CS) :: tInterpAlgo       ! Algorithm to use for time interpolation
      character(SHR_KIND_CL) :: calendar          ! stream calendar
   end type shr_stream_streamType

Finally, the ``shr_strdata_type`` is the heart of the CIME data
model implemenentation and contains information for all the streams
that are active for the target data model. The first part of the
``shr_strdata_type`` is filled in by the namelist values read in from the
namelist group (see the :ref:`stream data namelist section <shr-strdata-nml>`).

.. code-block:: Fortran

   type shr_strdata_type
     ! --- set by input namelist ---
    character(CL)  :: dataMode          ! flags physics options wrt input data
    character(CL)  :: domainFile        ! file   containing domain info
    character(CL)  :: streams (nStrMax) ! stream description file names
    character(CL)  :: taxMode (nStrMax) ! time axis cycling mode
    real(R8)       :: dtlimit (nStrMax) ! dt max/min limit
    character(CL)  :: vectors (nVecMax) ! define vectors to vector map
    character(CL)  :: fillalgo(nStrMax) ! fill algorithm
    character(CL)  :: fillmask(nStrMax) ! fill mask
    character(CL)  :: fillread(nStrMax) ! fill mapping file to read
    character(CL)  :: fillwrit(nStrMax) ! fill mapping file to write
    character(CL)  :: mapalgo (nStrMax) ! scalar map algorithm
    character(CL)  :: mapmask (nStrMax) ! scalar map mask
    character(CL)  :: mapread (nStrMax) ! regrid mapping file to read
    character(CL)  :: mapwrit (nStrMax) ! regrid mapping file to write
    character(CL)  :: tintalgo(nStrMax) ! time interpolation algorithm
    integer(IN)    :: io_type           ! io type, currently pnetcdf or netcdf

    !--- data required by cosz t-interp method, ---
    real(R8)     :: eccen   ! orbital eccentricity
    real(R8)     :: mvelpp  ! moving vernal equinox long
    real(R8)     :: lambm0  ! mean long of perihelion at vernal equinox (radians)
    real(R8)     :: obliqr  ! obliquity in degrees
    integer(IN)  :: modeldt ! data model dt in seconds (set to the coupling frequency)

    ! --- data model grid, public ---
    integer(IN)     :: nxg          ! data model grid lon size
    integer(IN)     :: nyg          ! data model grid lat size
    integer(IN)     :: nzg          ! data model grid vertical size
    integer(IN)     :: lsize        ! data model grid local size
    type(mct_gsmap) :: gsmap        ! data model grid global seg map
    type(mct_ggrid) :: grid         ! data model grid ggrid
    type(mct_avect) :: avs(nStrMax) ! data model stream attribute vectors

    ! --- stream specific arrays, stream grid ---
    type(shr_stream_streamType)    :: stream(nStrMax)
    type(iosystem_desc_t), pointer :: pio_subsystem => null()
    type(io_desc_t)    :: pio_iodesc(nStrMax)
    integer(IN)        :: nstreams          ! actual number of streams
    integer(IN)        :: strnxg(nStrMax)   ! stream grid lon sizes
    integer(IN)        :: strnyg(nStrMax)   ! stream grid lat sizes
    integer(IN)        :: strnzg(nStrMax)   ! tream grid global sizes
    logical            :: dofill(nStrMax)   ! true if stream grid is different from data model grid
    logical            :: domaps(nStrMax)   ! true if stream grid is different from data model grid
    integer(IN)        :: lsizeR(nStrMax)   ! stream local size of gsmapR on processor
    type(mct_gsmap)    :: gsmapR(nStrMax)   ! stream global seg map
    type(mct_rearr)    :: rearrR(nStrMax)   ! rearranger 
    type(mct_ggrid)    :: gridR(nStrMax)    ! local stream grid on processor
    type(mct_avect)    :: avRLB(nStrMax)    ! Read attrvect
    type(mct_avect)    :: avRUB(nStrMax)    ! Read attrvect
    type(mct_avect)    :: avFUB(nStrMax)    ! Final attrvect
    type(mct_avect)    :: avFLB(nStrMax)    ! Final attrvect
    type(mct_avect)    :: avCoszen(nStrMax) ! data assocaited with coszen time interp
    type(mct_sMatP)    :: sMatPf(nStrMax)   ! sparse matrix map for fill on stream grid
    type(mct_sMatP)    :: sMatPs(nStrMax)   ! sparse matrix map for mapping from stream to data model grid 
    integer(IN)        :: ymdLB(nStrMax)    ! lower bound time for stream
    integer(IN)        :: todLB(nStrMax)    ! lower bound time for stream
    integer(IN)        :: ymdUB(nStrMax)    ! upper bound time for stream
    integer(IN)        :: todUB(nStrMax)    ! upper bound time for stream
    real(R8)           :: dtmin(nStrMax)    
    real(R8)           :: dtmax(nStrMax)

    ! --- internal ---
    integer(IN)        :: ymd  ,tod
    character(CL)      :: calendar          ! model calendar for ymd,tod
    integer(IN)        :: nvectors          ! number of vectors
    integer(IN)        :: ustrm (nVecMax)
    integer(IN)        :: vstrm (nVecMax)
    character(CL)      :: allocstring
  end type shr_strdata_type


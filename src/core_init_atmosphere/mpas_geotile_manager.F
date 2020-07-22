module mpas_geotile_manager

    use iso_c_binding, only : c_float, c_char

    use mpas_constants, only : pii
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_log, only : mpas_log_write
    use mpas_derived_types, only : MPAS_LOG_ERR, MPAS_POOL_SILENT
    use mpas_pool_routines, only : mpas_pool_type, mpas_pool_destroy_pool, mpas_pool_create_pool
    use mpas_pool_routines, only : mpas_pool_add_config, mpas_pool_get_config
    use mpas_pool_routines, only : mpas_pool_get_error_level, mpas_pool_set_error_level
    use mpas_stack

    implicit none

    public :: mpas_geotile_mgr_type
    public :: mpas_geotile_type

    private

    type mpas_geotile_mgr_type
        type (mpas_pool_type), pointer :: pool
        type (tile_hash), dimension(:,:), pointer :: hash
        type (mpas_stack_type), pointer :: stack

        character (len=StrKIND) :: directory ! Path to the dataset directory
        character (len=StrKIND) :: index     ! Path the index file of the dataset directory

        integer :: nTileX ! Number of tiles in the X direction
        integer :: nTileY ! Number of tiles in the Y direction
        integer :: pixel_nx ! Total number of pixels in the x direction
        integer :: pixel_ny ! Total number of pixels in the y direction
    contains
        ! Public Procedures
        procedure, public :: init => mpas_geotile_mgr_init
        procedure, public :: finalize => mpas_geotile_mgr_finalize
        procedure, public :: get_tile => mpas_geotile_mgr_get_tile
        procedure, public :: pixel_to_latlon => mpas_geotile_mgr_pixel_to_latlon
        procedure, public :: latlon_to_pixel => mpas_geotile_mgr_latlon_to_pixel
        procedure, public :: tile_to_latlon => mpas_geotile_mgr_tile_to_latlon
        procedure, public :: push_neighbors => mpas_geotile_mgr_push_neighbors

        ! Stack Procedures
        procedure, public :: push_tile => mpas_geotile_mgr_push_tile
        procedure, public :: pop_tile => mpas_geotile_mgr_pop_tile
        procedure, public :: is_stack_empty => mpas_geotile_mgr_stack_is_empty

        ! Private Procedures
        procedure, private :: search_tile => mpas_geotile_mgr_search_tile
        procedure, private :: add_tile => mpas_geotile_mgr_add_tile
        procedure, private :: gen_filename => mpas_geotile_mgr_gen_tile_name
        procedure, private :: hash_to_ll => mpas_geotile_mgr_hash_to_latlon
    end type mpas_geotile_mgr_type


    type, extends(mpas_stack_payload_type) :: mpas_geotile_type
        real (c_float), dimension(:,:,:), pointer :: tile

        character (len=StrKIND) :: fname ! Path to the file that contains the data for this tile
        integer :: hash_x ! The x offset of this tile in the hash table
        integer :: hash_y ! The y offset of this tile in the hash table

        integer :: x, y ! The tiles range, in pixels
        logical :: is_processed = .false.
    end type mpas_geotile_type


    type tile_hash
        type(mpas_geotile_type), pointer :: ptr => null()
    end type tile_hash


    contains


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_init => init
    !
    !> \brief   Initialize a mpas_geotile_mgr class
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Initialize a geotile manager class by parsing the index file located
    !> within path and allocated needed data structures for static interpolation.
    !> Init should be called before calling any other mpas_geotile_mgr_type
    !> procedures. If path is not a directory or no index file is found in path,
    !> 1 will be returned. Upon success 0 will be returned.
    !>
    !> This function will also allocate the following variables in the pool attribute
    !> of this geotile manager instance if they are not found within the index file:
    !>     * tile_bdr = 0
    !>     * signed = 0 ! No
    !>     * scalefactor = 1.0_RKIND
    !>     * endian = "big"
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_init(mgr, path) result(ierr)

        use mpas_parse_geoindex, only : mpas_parse_index

        implicit none

        ! Input variables
        class (mpas_geotile_mgr_type) :: mgr
        character (len=*), intent(in) :: path

        ! Local variables
        character (len=StrKIND), pointer :: endian
        integer, pointer :: tile_nx ! Number of pixels in the x-direction for a single tile
        integer, pointer :: tile_ny ! Number of pixels in the y-direction for a single tile
        integer, pointer :: signed
        integer, pointer :: tile_bdr
        integer :: err_level
        real (kind=RKIND), pointer :: dx ! Grid spacing in the x-direction
        real (kind=RKIND), pointer :: dy ! Grid spacing in the y-direction
        real (kind=RKIND), pointer :: scalefactor
        logical :: res

        ! Return variable
        integer :: ierr

        ierr = 0

        mgr % directory = path

        ! Check to see if the index file exists in the directory
        inquire(file=trim(mgr % directory)//"index", exist=res)
        if (.not. res) then
            call mpas_log_write("Could not find an 'index' file in geotile directory: "//trim(path), messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif
        mgr % index = trim(mgr % directory)//"index"

        ! Create the pool for this geotile and call mpas_parse_index
        call mpas_pool_create_pool(mgr % pool)
        ierr = mpas_parse_index(mgr % index, mgr % pool)
        if (ierr /= 0) then
            call mpas_log_write("Error parsing geotile index file: "//trim(mgr % index), messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif

        signed => null()
        endian => null()
        scalefactor => null()
        tile_bdr => null()

        err_level = mpas_pool_get_error_level()
        call mpas_pool_set_error_level(MPAS_POOL_SILENT)

        call mpas_pool_get_config(mgr % pool, 'signed', signed)
        call mpas_pool_get_config(mgr % pool, 'endian', endian)
        call mpas_pool_get_config(mgr % pool, 'scale_factor', scalefactor)
        call mpas_pool_get_config(mgr % pool, 'tile_bdr', tile_bdr)


        !
        ! tile_bdr, signed, endian, and scale_factor all have default values, so if they
        ! are not present in the index file then set them as the default values, as
        ! reported in section 3-53 of the WRF-ARW User's Guide
        !
        if (.not. associated(endian)) then
            allocate(endian)
            endian = "big"
            call mpas_pool_add_config(mgr % pool, 'endian', endian)
        endif

        if (.not. associated(scalefactor)) then
            allocate(scalefactor)
            scalefactor = 1.0_RKIND
            call mpas_pool_add_config(mgr % pool, 'scale_factor', scalefactor)
        endif

        if (.not. associated(signed)) then
            allocate(signed)
            signed = 0 ! 0 = 'no', 1 = 'yes'
            call mpas_pool_add_config(mgr % pool, 'signed', signed)
        endif

        if (.not. associated(tile_bdr)) then
            allocate(tile_bdr)
            tile_bdr = 0
            call mpas_pool_add_config(mgr % pool, 'tile_bdr', tile_bdr)
        endif

        ! Reset the pool's error level
        call mpas_pool_set_error_level(err_level)

        call mpas_pool_get_config(mgr % pool, 'tile_x', tile_nx)
        call mpas_pool_get_config(mgr % pool, 'tile_y', tile_ny)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        ! Calculate the total number of pixels in x dir
        ! NOTE: This calculation assumes that a dataset is a global dataset and may
        ! not work correctly for non-global datasets
        mgr % pixel_nx = nint(360.0_RKIND / dx)
        mgr % pixel_ny = nint(180.0_RKIND / dy)

        ! Calculate the number of tiles in the x, y directions
        ! NOTE: This calculation assumes that a dataset is a global dataset and may
        ! not work correctly for non-global datasets
        mgr % nTileX = mgr % pixel_nx / tile_nx
        mgr % nTileY = mgr % pixel_ny / tile_ny

        ! Allocate hash table
        allocate(mgr % hash(0: mgr % nTileX, 0: mgr % nTileY))

        ! Mark the stack as empty
        mgr % stack => null()

    end function mpas_geotile_mgr_init


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_finalize => finalize
    !
    !> \brief   Free all memory used by the mpas_geotile_mgr_type
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Deallocated all memory used by this geotile_mgr_type and destroy the
    !> associated pool. After calling this function, none of the methods
    !> should be used, unless the class is reinitialized by recalling
    !> mpas_geotile_mgr_init.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_finalize(mgr) result(ierr)

        implicit none

        ! Input variable
        class (mpas_geotile_mgr_type) :: mgr

        ! Return variable
        integer :: ierr

        ! Local variable
        integer :: i
        integer :: j

        ierr = 0

        ! Loop through the hash table and deallocate any loaded tiles
        ! Then deallocate the hash table
        do i = 0, mgr % nTileX
            do j = 0, mgr % nTileY
                if (associated(mgr % hash(i, j) % ptr)) then
                    deallocate(mgr % hash(i, j) % ptr % tile)
                    deallocate(mgr % hash(i, j) % ptr)
                endif
            enddo
        enddo
        deallocate(mgr % hash)

        if (associated(mgr % hash)) then
            call mpas_log_write("Problem deallocating the geotile hash table", messageType=MPAS_LOG_ERR)
            ierr = -1
            return
        endif

        call mpas_pool_destroy_pool(mgr % pool)
        if (associated(mgr % pool)) then
            call mpas_log_write("Problem deallocating the geotile pool", messageType=MPAS_LOG_ERR)
            ierr = -1
            return
        endif

        call mpas_stack_free(mgr % stack)
        if (associated(mgr % stack)) then
            call mpas_log_write("Problem deallocating the stack", messageType=MPAS_LOG_ERR)
            ierr = -1
            return
        endif

    end function mpas_geotile_mgr_finalize


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_get_tile => get_tile
    !
    !> \brief   Return an array containing the values of a datatile
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Retrieve the datatile that contains the coordinate lat, lon of the dataset
    !> that this mpas_geotile_manager instance was initalized with. Both lat,
    !> lon should be in radians and lon should be in the range of -1/2 * pi to
    !> 1/2 * pi. lat values that are greater than 2.0 * pi or less than -2.0 * pi
    !> will be normalized to be between -pi and pi. Upon success 0 will be returned
    !> and tile will point to the mpas_geotile_type that holds the datatile which
    !> contains the coordinate lat, lon.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_get_tile(mgr, lat, lon, tile) result(ierr)
        implicit none

        ! Input variables
        class (mpas_geotile_mgr_type) :: mgr
        real (kind=RKIND), value :: lat
        real (kind=RKIND), value :: lon
        type (mpas_geotile_type), pointer :: tile

        ! Return variable
        integer :: ierr

        ierr = 0
        tile => null()

        ! Normalize longitude to be between -pi and pi
        call normalize_lon(lon)

        if (.not. mgr % search_tile(lat, lon, tile)) then
            ierr = mgr % add_tile(lat, lon, tile)
        endif

    end function mpas_geotile_mgr_get_tile


    !***********************************************************************
    !
    !  private function mpas_geotile_mgr_search_tile => search_tile
    !
    !> \brief   Search to see if a tile has already been loaded
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Private function that searches to see if the datatile that contains
    !> the coordinate lat, lon has already been loaded. If the datatile has been
    !> loaded, .true. will be returned and tile will point to the mpas_geotile_type
    !> that contains the datatile. If the datatile has not been loaded, .false.
    !> will be returned and tile will be unassociated.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_search_tile(mgr, lat, lon, tile) result(loaded)

        implicit none

        ! Input variables
        class (mpas_geotile_mgr_type) :: mgr
        real (kind=RKIND), value :: lat
        real (kind=RKIND), value :: lon
        type (mpas_geotile_type), pointer :: tile

        ! Return variable
        logical :: loaded

        ! Local variables
        integer, pointer :: tile_nx
        integer, pointer :: tile_ny
        character (len=StrKIND) :: fname
        integer :: x, y
        integer :: start_x
        integer :: start_y
        integer :: ierr

        loaded = .false.
        tile => null()

        call mpas_pool_get_config(mgr % pool, 'tile_x', tile_nx)
        call mpas_pool_get_config(mgr % pool, 'tile_y', tile_ny)

        !
        ! Using gen_filename, get the tiles start x and y pixel values of the tile
        !
        ierr = mgr % gen_filename(lat, lon, fname, start_x, start_y)
        if (ierr /= 0) then
            call mpas_log_write("Error generating filename", messageType=MPAS_LOG_ERR)
            return
        endif

        !
        ! Access the tile in the hash table (-1 here as the hash table is from
        ! 0:tile_nx, and 0:tile_ny).
        !
        x = (start_x - 1) / tile_nx
        if (x < 0 .or. x > size(mgr % hash, 1)) then
            return
        endif

        y = (start_y - 1) / tile_ny
        if (y < 0 .or. y > size(mgr % hash, 2)) then
            return
        endif

        tile => mgr % hash(x,y) % ptr
        if (associated(tile)) then
            loaded = .true.
        endif

    end function mpas_geotile_mgr_search_tile


    !***********************************************************************
    !
    !  private function mpas_geotile_mgr_add_tile => add_tile
    !
    !> \brief   Read in a datatile file and store a reference to it
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Read the datatile that contains the coordinate lat, lon. Open success,
    !> 0 will be returned and tile will point to the mpas_geotile_type which
    !> contains the coordiate lat, lon. Upon success a reference to that
    !> mpas_geotile_type will be placed into the hash table, which can later
    !> be searched via search_tile. On error, 1 will be returned and tile
    !> will be unassociated.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_add_tile(mgr, lat, lon, tile) result(ierr)

        use iso_c_binding, only : c_loc, c_ptr, c_int, c_float
        use mpas_c_interfacing, only : mpas_f_to_c_string

        implicit none

        interface
           subroutine read_geogrid(fname, rarray, nx, ny, nz, isigned, endian, &
                                   scalefactor, wordsize, status) bind(C)
              use iso_c_binding, only : c_char, c_int, c_float, c_ptr
              character (c_char), dimension(*), intent(in) :: fname
              type (c_ptr), value :: rarray
              integer (c_int), intent(in), value :: nx
              integer (c_int), intent(in), value :: ny
              integer (c_int), intent(in), value :: nz
              integer (c_int), intent(in), value :: isigned
              integer (c_int), intent(in), value :: endian
              real (c_float), intent(in), value :: scalefactor
              integer (c_int), intent(in), value :: wordsize
              integer (c_int), intent(inout) :: status
           end subroutine read_geogrid
        end interface

        ! Arguments
        class (mpas_geotile_mgr_type) :: mgr
        real (kind=RKIND), intent(in) :: lat
        real (kind=RKIND), intent(in) :: lon
        type (mpas_geotile_type), intent(inout), pointer :: tile
        integer :: ierr

        ! Local variables
        integer, pointer :: tile_nx, tile_ny, tile_nz
        integer, pointer :: tile_bdr
        integer, pointer :: wordsize
        integer :: start_x, start_y
        integer, pointer :: signed
        character (len=StrKIND), pointer :: endian
        real (kind=RKIND), pointer :: scalefactor
        logical :: res

        character (len=StrKIND) :: fname
        character (kind=c_char), dimension(StrKIND+1) :: c_fname
        integer (c_int) :: c_tile_nx, c_tile_ny, c_tile_nz
        integer (c_int) :: c_endian
        integer (c_int) :: c_wordsize
        integer (c_int) :: c_signed
        real (c_float) :: c_scalefactor
        integer (c_int) :: status
        type (c_ptr) :: c_tile_ptr

        ierr = 0

        tile => null()
        endian => null()
        scalefactor => null()
        signed => null()

        call mpas_pool_get_config(mgr % pool, 'tile_x', tile_nx)
        call mpas_pool_get_config(mgr % pool, 'tile_y', tile_ny)
        call mpas_pool_get_config(mgr % pool, 'tile_z', tile_nz)
        call mpas_pool_get_config(mgr % pool, 'tile_bdr', tile_bdr)
        call mpas_pool_get_config(mgr % pool, 'wordsize', wordsize)
        call mpas_pool_get_config(mgr % pool, 'signed', signed)
        call mpas_pool_get_config(mgr % pool, 'endian', endian)
        call mpas_pool_get_config(mgr % pool, 'scale_factor', scalefactor)

        c_tile_nx = tile_nx + 2 * tile_bdr ! The number of pixels in the x direction, including halo cells
        c_tile_ny = tile_ny + 2 * tile_bdr ! The number of pixels in the y direction, including halo cells
        c_tile_nz = tile_nz ! The number of pixels in the z direction
        c_wordsize = wordsize
        c_scalefactor = scalefactor
        c_signed = signed

        if (endian == "big") then
            c_endian = 0
        else if (endian == "little") then
            c_endian = 1
        endif

        !
        ! Determine the file that contains lat, lon.
        !
        ierr = mgr % gen_filename(lat, lon, fname, start_x, start_y)
        if (ierr /= 0) then
            call mpas_log_write("Error creating filename for coordinate: ($r, $r)", realArgs=(/lat, lon/), messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif

        !
        ! See if this file actually exists
        !
        fname = trim(mgr % directory)//trim(fname)
        inquire(file=trim(fname), exist=res)
        if (.not. res) then
            call mpas_log_write("This geotile file does not exist: "//trim(fname), messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif

        call mpas_f_to_c_string(fname, c_fname)

        !
        ! Allocate and read the tile
        !
        allocate(tile)
        allocate(tile % tile(tile_nx + (tile_bdr * 2), tile_ny + (tile_bdr * 2), tile_nz))
        c_tile_ptr = c_loc(tile % tile)
        call read_geogrid(c_fname, c_tile_ptr, c_tile_nx, c_tile_ny, c_tile_nz, c_signed, c_endian, &
                            c_scalefactor, c_wordsize, status)
        if (status /= 0) then
            call mpas_log_write("Error reading this geogrid file: "//trim(fname), messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif

        tile % fname = fname
        tile % hash_x = (start_x - 1) / tile_nx
        tile % hash_y = (start_y - 1) / tile_ny
        tile % x = start_x
        tile % y = start_y

        !
        ! Add the tile to the hash table
        !
        mgr % hash(tile % hash_x, tile % hash_y)  % ptr => tile

    end function mpas_geotile_mgr_add_tile


    !***********************************************************************
    !
    !  private function mpas_geotile_mgr_gen_tile_name => gen_filename
    !
    !> \brief   Generate the filename of the tile at lat, lon
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Generate the name of the file that contains the coordinate lat, lon
    !> (in radians) and optionally return start_x and start_y (the location
    !> of the first global pixel coordinate of tile). Warning: This function
    !> can produce filenames that may not exist (For lon less than -.5 * pi and
    !> greater than .5 * pi and lat less than -pi and greater than pi).
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_gen_tile_name(mgr, lat, lon, fname, start_x, start_y) result(ierr)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        real (kind=RKIND), value :: lat
        real (kind=RKIND), value :: lon
        character (len=StrKIND), intent(out) :: fname
        integer, intent(out), optional :: start_x
        integer, intent(out), optional :: start_y

        character (len=StrKIND), parameter :: fname_format = "(i5.5, '-', i5.5, '.', i5.5, '-', i5.5)"

        real (kind=RKIND), pointer :: dx
        real (kind=RKIND), pointer :: dy
        integer, pointer :: tile_nx
        integer, pointer :: tile_ny
        integer, dimension(2) :: x
        integer, dimension(2) :: y

        integer :: ierr
        ierr = 0

        call mpas_pool_get_config(mgr % pool, 'tile_x', tile_nx)
        call mpas_pool_get_config(mgr % pool, 'tile_y', tile_ny)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        ! Find the global pixel location that contains lat, lon
        call mgr % latlon_to_pixel(lat, lon, x(1), y(1))

        ! Calculate the range of this tile, which will be its filename
        x(1) = (x(1) - modulo(x(1), tile_nx)) + 1
        x(2) = x(1) + tile_nx - 1

        y(1) = (y(1) - modulo(y(1), tile_ny)) + 1
        y(2) = y(1) + tile_ny - 1

        write(fname, fname_format) x(1), x(2), y(1), y(2)

        if (present(start_x)) then
            start_x = x(1)
        endif
        if (present(start_y)) then
            start_y = y(1)
        endif

    end function mpas_geotile_mgr_gen_tile_name


    !***********************************************************************
    !
    !  public subroutine mpas_geotile_mgr_tile_to_latlon => tile_to_latlon
    !
    !> \brief   Translate a tile indices to latitude and longitude,
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Given a tile, translate the relative tile coordinates, i, j, of that
    !> tile to a latitude and longitude coordinate. Upon success, lat, lon
    !> will be in the range of -1/2 * pi to 1/2 * pi and 0 to 2.0 * pi, respectively.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_geotile_mgr_tile_to_latlon(mgr, tile, j, i, lat, lon)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        type (mpas_geotile_type), pointer, intent(in) :: tile
        integer, value :: j
        integer, value :: i
        real (kind=RKIND), intent(out) :: lat
        real (kind=RKIND), intent(out) :: lon

        integer, pointer :: tile_bdr
        real (kind=RKIND), pointer :: known_lon
        real (kind=RKIND), pointer :: known_lat
        real (kind=RKIND), pointer :: dx
        real (kind=RKIND), pointer :: dy
        integer :: ierr

        ierr = 0

        call mpas_pool_get_config(mgr % pool, 'tile_bdr', tile_bdr)
        call mpas_pool_get_config(mgr % pool, 'known_lat', known_lat)
        call mpas_pool_get_config(mgr % pool, 'known_lon', known_lon)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        lat = real((j - (tile_bdr + 1) + tile % y - 1), kind=RKIND) * dy + known_lat
        lon = real((i - (tile_bdr + 1) + tile % x - 1), kind=RKIND) * dx + known_lon

        call deg2Rad(lat)
        call deg2Rad(lon)

    end subroutine mpas_geotile_mgr_tile_to_latlon


    !***********************************************************************
    !
    !  public subroutine mpas_geotile_mgr_latlon_to_pixel => latlon_to_pixel
    !
    !> \brief   Translate a latitude, longitude coordinate to pixel coordinates
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Translate a latitude, longitude coordinate into global pixel coordinates. lat
    !> should be in the range of -.5 * pi to .5 * pi and lon should be between -pi
    !> and pi.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_geotile_mgr_latlon_to_pixel(mgr, lat, lon, pixel_x, pixel_y)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        real (kind=RKIND), value :: lat
        real (kind=RKIND), value :: lon
        integer, intent(out) :: pixel_x
        integer, intent(out) :: pixel_y

        integer, pointer :: tile_bdr
        real (kind=RKIND), pointer :: known_lon
        real (kind=RKIND), pointer :: known_lat
        real, pointer :: dx
        real, pointer :: dy
        integer :: ierr

        ierr = 0

        call mpas_pool_get_config(mgr % pool, 'tile_bdr', tile_bdr)
        call mpas_pool_get_config(mgr % pool, 'known_lon', known_lon)
        call mpas_pool_get_config(mgr % pool, 'known_lat', known_lat)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        call rad2Deg(lat)
        call rad2Deg(lon)

        pixel_x = nint((lon - known_lon) / dx)
        pixel_y = nint((lat - known_lat) / dy)

        if (pixel_x < 0) then
            pixel_x = pixel_x + mgr % pixel_nx
        endif

        if (pixel_y < 0) then
            pixel_y = 0
        else if (pixel_y >= mgr % pixel_ny) then
            pixel_y = mgr % pixel_ny - 1
        endif

    end subroutine mpas_geotile_mgr_latlon_to_pixel


    !***********************************************************************
    !
    !  public subroutine mpas_geotile_mgr_pixel_to_latlon => pixel_to_latlon
    !
    !> \brief   Translate a pixel coordinate to latitude longitude
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Translate a global pixel coordinate to a latitude, longitude coordinate
    !
    !-----------------------------------------------------------------------
    subroutine mpas_geotile_mgr_pixel_to_latlon(mgr, x, y, lat, lon)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        integer, value :: x
        integer, value :: y
        real (kind=RKIND), intent(out) :: lat
        real (kind=RKIND), intent(out) :: lon

        integer, pointer :: tile_bdr
        real (kind=RKIND), pointer :: known_lat
        real (kind=RKIND), pointer :: known_lon
        real (kind=RKIND), pointer :: known_x
        real (kind=RKIND), pointer :: known_y
        real, pointer :: dx
        real, pointer :: dy
        integer :: ierr

        ierr = 0

        call mpas_pool_get_config(mgr % pool, 'tile_bdr', tile_bdr)
        call mpas_pool_get_config(mgr % pool, 'known_lat', known_lat)
        call mpas_pool_get_config(mgr % pool, 'known_lon', known_lon)
        call mpas_pool_get_config(mgr % pool, 'known_x', known_x)
        call mpas_pool_get_config(mgr % pool, 'known_y', known_y)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        lat = y * dy + known_lat
        lon = x * dx + known_lon

        call deg2Rad(lat)
        call deg2Rad(lon)

    end subroutine mpas_geotile_mgr_pixel_to_latlon


    !***********************************************************************
    !
    !  private function mpas_geotile_mgr_hash_to_latlon => hash_to_ll
    !
    !> \brief   Find the lat, lon center from a hash entry
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Translate the index within the hash table into the latitude and longitude
    !> coordinate (in radians) of the center of the datatile at that index.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_geotile_mgr_hash_to_latlon(mgr, xHash, yHash, lat, lon)

        implicit none

        class(mpas_geotile_mgr_type) :: mgr
        integer, intent(in), value :: xHash
        integer, intent(in), value :: yHash
        real, intent(out) :: lat
        real, intent(out) :: lon

        integer, pointer :: tile_nx
        integer, pointer :: tile_ny
        real (kind=RKIND), pointer :: known_lat
        real (kind=RKIND), pointer :: known_lon
        real (kind=RKIND), pointer :: dx
        real (kind=RKIND), pointer :: dy

        integer :: x
        integer :: y

        call mpas_pool_get_config(mgr % pool, 'tile_x', tile_nx)
        call mpas_pool_get_config(mgr % pool, 'tile_y', tile_ny)
        call mpas_pool_get_config(mgr % pool, 'known_lat', known_lat)
        call mpas_pool_get_config(mgr % pool, 'known_lon', known_lon)
        call mpas_pool_get_config(mgr % pool, 'dx', dx)
        call mpas_pool_get_config(mgr % pool, 'dy', dy)

        x = (xHash * tile_nx) + (tile_nx / 2)
        y = (yHash * tile_ny) + (tile_ny / 2)

        lon = (real((x), kind=RKIND) * dx ) + known_lon
        lat = (real((y), kind=RKIND) * dy ) + known_lat

        call deg2Rad(lat)
        call deg2Rad(lon)

    end subroutine mpas_geotile_mgr_hash_to_latlon


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_push_neighbors => push_neighbors
    !
    !> \brief   Determine the tile nighbors and push them onto the stack
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Determine the neighbors of a tile and push them onto the stack. If the
    !> tile neighbors have not been loaded, via add_tile, then they will be.
    !> Upon success, the neighbors of a tile will be pushed onto the stack and
    !> 0 will be returned.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_push_neighbors(mgr, tile) result(ierr)

        implicit none

        class(mpas_geotile_mgr_type) :: mgr
        type (mpas_geotile_type), pointer, intent(in) :: tile

        integer :: ierr
        type (mpas_geotile_type), pointer :: neighbor
        real (kind=RKIND) :: lat
        real (kind=RKIND) :: lon
        integer :: xHash
        integer :: yHash

        ierr = 0

        ! Up
        ! Calculate the tile's hash coordinates
        neighbor => null()
        if (tile % hash_y + 1 > mgr % nTileY - 1) then
            xHash = modulo(tile % hash_x + (mgr % nTileX / 2), mgr % nTileX)
            yHash = tile % hash_y
        else
            xHash = tile % hash_x
            yHash = tile % hash_y + 1
        endif
        call mgr % hash_to_ll(xHash, yHash, lat, lon)

        ierr = mgr % get_tile(lat, lon, neighbor)
        if (ierr /= 0) then
            call mpas_log_write("There was a problem getting the up tile", messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif
        ierr = mgr % push_tile(neighbor)

        ! Down
        neighbor => null()
        if (tile % hash_y - 1 < 0) then
            xHash = modulo(tile % hash_x + (mgr % nTileX / 2), mgr % nTileX)
            yHash = tile % hash_y
        else
            xHash = tile % hash_x
            yHash = tile % hash_y - 1
        endif
        call mgr % hash_to_ll(xHash, yHash, lat, lon)

        ierr = mgr % get_tile(lat, lon, neighbor)
        if (ierr /= 0) then
            call mpas_log_write("There was a problem getting the down tile", messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif
        ierr = mgr % push_tile(neighbor)

        ! Right
        neighbor => null()
        if (tile % hash_x + 1 > mgr % nTileX - 1) then
            yHash = tile % hash_y
            xHash = 0
        else
            xHash = tile % hash_x + 1
            yHash = tile % hash_y
        endif
        call mgr % hash_to_ll(xHash, yHash, lat, lon)

        ierr = mgr % get_tile(lat, lon, neighbor)
        if (ierr /= 0) then
            call mpas_log_write("There was a problem getting the right tile", messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif
        ierr = mgr % push_tile(neighbor)

        ! Left
        neighbor => null()
        if (tile % hash_x -1 < 0) then
            xHash = mgr % nTileX - 1
            yHash = tile % hash_y
        else
            xHash = tile % hash_x - 1
            yHash = tile % hash_y
        endif
        call mgr % hash_to_ll(xHash, yHash, lat, lon)

        ierr = mgr % get_tile(lat, lon, neighbor)
        if (ierr /= 0) then
            call mpas_log_write("There was a problem getting the left tile", messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        endif
        ierr = mgr % push_tile(neighbor)

    end function mpas_geotile_mgr_push_neighbors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Stack wrappers and helper functions !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_push_tile => push_tile
    !
    !> \brief   Push a mpas_geotile_type onto the stack
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Wrapper subroutine for mpas_stack_push from mpas_stack.F. Aftering calling
    !> this subroutine, the tile pushed will be on the top of the stack associated
    !> with mpas_geotile_mgr instance (TODO: Is instance the correct term??) and pop_tile
    !> can be used to retrive the tile that was last pushed onto the stack.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_push_tile(mgr, tile) result(ierr)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        type (mpas_geotile_type), pointer :: tile
        integer :: ierr

        ierr = 0

        mgr % stack => mpas_stack_push(mgr % stack, tile)

    end function mpas_geotile_mgr_push_tile


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_pop_tile => pop_tile
    !
    !> \brief   Pop the top mpas_geotile_type off of the stack
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Retrieve and remove the last tile that was pushed onto the stack that
    !> is associated with this mpas_geotile_mgr instance (TODO: Is instance the correct term?). If the stack is empty,
    !> then tile will be unassociated.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_pop_tile(mgr) result(tile)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        class (mpas_stack_payload_type), pointer :: top
        type (mpas_geotile_type), pointer :: tile

        tile => null()

        if (mpas_stack_is_empty(mgr % stack)) then
            return
        endif

        top => mpas_stack_pop(mgr % stack)

        select type(top)
            type is(mpas_geotile_type)
                tile => top
                return
            class default
                ! Should not get here
        end select

    end function mpas_geotile_mgr_pop_tile


    !***********************************************************************
    !
    !  public function mpas_geotile_mgr_stack_is_empty => is_stack_empty
    !
    !> \brief   Return .true. if stack is empty and .false. otherwise
    !> \author  Miles A. Curry
    !> \date    02/20/2020
    !> \details
    !> Return .true. if the stack associated with this mpas_geotile_mgr instance (TODO: Is instance the correct term?)
    !> is empty, and .false. if it is not empty.
    !
    !-----------------------------------------------------------------------
    function mpas_geotile_mgr_stack_is_empty(mgr) result(empty)

        implicit none

        class (mpas_geotile_mgr_type) :: mgr
        logical :: empty

        empty = mpas_stack_is_empty(mgr % stack)

    end function mpas_geotile_mgr_stack_is_empty



    ! Convert radians to degrees
    subroutine rad2Deg(rad)

        implicit none
        real (kind=RKIND), intent(inout) :: rad

        rad = rad * (180.0_RKIND / pii)

    end subroutine rad2Deg


    ! Convert degrees to radians
    subroutine deg2Rad(deg)

        implicit none
        real (kind=RKIND), intent(inout) :: deg

        deg = deg * (pii / 180.0_RKIND)

    end subroutine deg2Rad


    ! Normalize logitude (in radians) to be between -pi and pi.
    subroutine normalize_lon(lon)

        implicit none
        real (kind=RKIND), intent(inout) :: lon

        if (lon > pii) then
            lon = lon - (2.0 * pii)
        else if (lon < -pii) then
            lon = lon + (2.0 * pii)
        endif

    end subroutine normalize_lon

end module mpas_geotile_manager

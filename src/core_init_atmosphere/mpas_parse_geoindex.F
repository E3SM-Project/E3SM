module mpas_parse_geoindex

   use mpas_log, only : mpas_log_write
   use mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_WARN
   use mpas_pool_routines

   implicit none

   private

   public :: mpas_parse_index

   contains

   !***********************************************************************
   !
   !  subroutine mpas_parse_index
   !
   !> \brief   Parse a geogrid's index file and put the results into an MPAS pool
   !> \author  Miles A. Curry
   !> \date    02/20/2020
   !> \details
   !> Parse an index file of a static data set into an MPAS pool, allocating
   !> each keyword=value pair into the pool with the pool member key being
   !> keyword, and the value being value.
   !>
   !> This function can parse index files with one keyword=value pair
   !> per line; a "#" at the start of a line, which will cause the line to be
   !> ignored; or an empty line containing only a newline/return character, which
   !> will also be ignored. Spaces or tabs before, between or after the
   !> keyword=value tokens are > ignored.
   !>
   !> If a line contains anything but the above valid syntaxes, a syntax
   !> error will raised and -1 will be returned.
   !>
   !> Case is ignored.
   !>
   !> The definitions of a keyword, which can found in section 3-53
   !> of the WRF-AWR User's Guide, will determine the corresponding type
   !> of that keyword. A keyword that has been assigned the wrong type
   !> will raise a type error and -1 will be returned.
   !>
   !> Keywords that are not handled explicitly by this function will produce
   !> a warning that the keyword is unrecognized.
   !
   !-----------------------------------------------------------------------
   function mpas_parse_index(path, geo_pool) result(ierr)

      use mpas_io_units

      implicit none
      ! Input Variables
      character (len=*), intent(in) :: path
      type (mpas_pool_type), intent(inout) :: geo_pool
      integer :: ierr

      ! Local Variables
      character (len=StrKIND) :: line, lhs, rhs
      character (len=StrKIND) :: read_err_msg, open_msg
      integer :: geo_unit
      integer :: open_stat, read_stat, line_read_stat
      integer :: i, k
      logical :: res

      character (len=StrKIND), pointer :: char_t
      integer :: iceiling, ifloor
      integer, pointer :: int_t
      real(kind=RKIND), pointer :: real_t

      ierr = 0

      inquire(file=trim(path), exist=res)
      if ( .not. res) then
         call mpas_log_write("Could not find or open the file at: "//trim(path), messageType=MPAS_LOG_ERR)
         ierr = -1
         return
      endif

      call mpas_new_unit(geo_unit)
      open_stat = 0
      open(geo_unit, FILE=trim(path), action='READ', iostat=open_stat, iomsg=open_msg)
      if (open_stat /= 0) then
         call mpas_release_unit(geo_unit)
         call mpas_log_write("Could not open 'index' file at:'"//trim(path)//"'", messageType=MPAS_LOG_ERR)
         call mpas_log_write(trim(open_msg), messageType=MPAS_LOG_ERR)
         ierr = -1
         return
      endif

      line_read_stat = 0
      read_stat = 0
      k = 1 ! Keep track of line numbers for error reporting
      read(geo_unit,'(a)', iostat=line_read_stat) line
      do while ( line_read_stat == 0 )
         line = lowercase(line)

         !
         ! If a blank or comment line is encountered, read the next line
         !
         if (line(1:1) == '#' .or. len_trim(line) == 0) then
            k = k + 1
            read(geo_unit,'(a)', iostat=line_read_stat) line
            cycle
         endif

         do i = 1, len(trim(line)), 1
            if (line(i:i) == '=') then
               lhs = adjustl(trim(line(1:i-1)))
               rhs = adjustl(trim(line(i+1:len(trim(line)))))
               exit
            endif
            ! If i is at the end of the string, and we haven't broken out of this loop,
            ! then we do not have a '=' present in this line, thus we have a syntax error
            if (i == len(trim(line))) then
               close(geo_unit)
               call mpas_release_unit(geo_unit)
               call mpas_log_write("Syntax error on line $i of index file: '"//trim(path)//"'", &
                                                                               intArgs=(/k/), messageType=MPAS_LOG_ERR)
               call mpas_log_write("Line $i: '"//trim(line)//"'", intArgs=(/k/), messageType=MPAS_LOG_ERR)
               ierr = -1
               return
            endif
         enddo

         !
         ! Strings
         !
         if (        trim(lhs) == 'type' &
                .or. trim(lhs) == 'projection' &
                .or. trim(lhs) == 'units' &
                .or. trim(lhs) == 'description' &
                .or. trim(lhs) == 'row_order' &
                .or. trim(lhs) == 'endian' &
                .or. trim(lhs) == 'mminlu'   ) then

            allocate(char_t)
            char_t = rhs
            call mpas_pool_add_config(geo_pool, trim(lhs), char_t)

         !
         ! Reals
         !
         else if (     trim(lhs) == 'dx' &
                  .or. trim(lhs) == 'dy' &
                  .or. trim(lhs) == 'known_x' &
                  .or. trim(lhs) == 'known_y' &
                  .or. trim(lhs) == 'known_lat' &
                  .or. trim(lhs) == 'known_lon' &
                  .or. trim(lhs) == 'scale_factor' &
                  .or. trim(lhs) == 'stdlon' &
                  .or. trim(lhs) == 'truelat1' &
                  .or. trim(lhs) == 'truelat2' &
                  .or. trim(lhs) == 'missing_value' ) then

            allocate(real_t)
            read(rhs, *, iostat=read_stat, iomsg=read_err_msg) real_t
            call mpas_pool_add_config(geo_pool, trim(lhs), real_t)

         !
         ! Integers
         !
         else if (     trim(lhs) == 'tile_x' &
                  .or. trim(lhs) == 'tile_y' &
                  .or. trim(lhs) == 'tile_z' &
                  .or. trim(lhs) == 'tile_z_start' &
                  .or. trim(lhs) == 'tile_z_end' &
                  .or. trim(lhs) == 'tile_bdr' &
                  .or. trim(lhs) == 'wordsize' &
                  .or. trim(lhs) == 'category_max' &
                  .or. trim(lhs) == 'category_min' &
                  .or. trim(lhs) == 'iswater' &
                  .or. trim(lhs) == 'islake' &
                  .or. trim(lhs) == 'isice' &
                  .or. trim(lhs) == 'isurban' &
                  .or. trim(lhs) == 'isoilwater' &
                  .or. trim(lhs) == 'filename_digits' ) then

            ! Because each compiler handles reporting type errors when transferring
            ! data in a read statement a little bit differently, we will have to type check
            ! integer values ourselves.
            allocate(real_t)
            read(rhs, *, iostat=read_stat, iomsg=read_err_msg) real_t
            iceiling = ceiling(real_t)
            ifloor = floor(real_t)
            if (iceiling /= ifloor) then
               close(geo_unit)
               call mpas_release_unit(geo_unit)
               call mpas_log_write("Type error while reading '"//trim(path)//"'.", messageType=MPAS_LOG_ERR)
               call mpas_log_write("Could not convert '"//trim(rhs)//"' to an integer on line $i: '"//trim(line)//"'", &
                                                                             intArgs=(/k/), messageType=MPAS_LOG_ERR)
               ierr = -1
               return
            endif

            allocate(int_t)
            int_t = int(real_t)
            deallocate(real_t)
            call mpas_pool_add_config(geo_pool, trim(lhs), int_t)

         !
         ! Booleans - Yes will be assigned 1, and no will be assigned to 0
         !
         else if (lhs == 'signed') then
            if (trim(rhs) == 'yes') then
               allocate(int_t)
               int_t = 1
               call mpas_pool_add_config(geo_pool, trim(lhs), int_t)
            else if (trim(rhs) == 'no') then
               allocate(int_t)
               int_t = 0
               call mpas_pool_add_config(geo_pool, trim(lhs), int_t)
            else
               read_stat = -1
               read_err_msg = "Logical was not correct type"
            endif
         else
            call mpas_log_write("Unrecognized keyword: '"//trim(lhs)//"' on line $i of '"//trim(path)//"'", intArgs=(/k/), &
                               messageType=MPAS_LOG_WARN)
         endif
         ! Since read gives us an error string in iomsg on a type error, we
         ! can handle all errors for any type in one place
         if ( read_stat /= 0) then
            close(geo_unit)
            call mpas_release_unit(geo_unit)
            call mpas_log_write("Type error on line $i of: '"//trim(path)//"'.", intArgs=(/k/), messageType=MPAS_LOG_ERR)
            call mpas_log_write(trim(read_err_msg)//": '"//trim(line)//"'", messageType=MPAS_LOG_ERR)
            ierr = -1
            return
         endif

         k = k + 1
         read(geo_unit,'(a)', iostat=line_read_stat) line
      enddo

      close(geo_unit)
      call mpas_release_unit(geo_unit)

   end function mpas_parse_index


   ! Returns a copy of 'str' in which all upper-case letters have been converted
   ! to lower-case letters.
   function lowercase(str) result(lowerStr)

       character(len=*), intent(in) :: str
       character(len=len(str)) :: lowerStr

       integer :: i
       integer, parameter :: offset = (iachar('a') - iachar('A'))


       do i=1,len(str)
           if (iachar(str(i:i)) >= iachar('A') .and. iachar(str(i:i)) <= iachar('Z')) then
              lowerStr(i:i) = achar(iachar(str(i:i)) + offset)
           else
              lowerStr(i:i) = str(i:i)
           end if
       end do

   end function lowercase


end module mpas_parse_geoindex

module mpas_stream_list

#define COMMA ,
#define LIST_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define LIST_WARN_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define LIST_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

    use mpas_kind_types, only : StrKIND
    use mpas_io_units, only : stderrUnit
    use mpas_derived_types


    contains


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_create
    !
    !> \brief Initialize a new MPAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Instantiates and initializes a stream_list type, to store all active streams.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_list_create(list, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), pointer :: list
        integer, intent(out), optional :: ierr


        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_create()')

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        allocate(list)

        list % nItems = 0
        nullify(list % head)

    end subroutine MPAS_stream_list_create !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_destroy
    !
    !> \brief Free all memory associated with an MPAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Destroys a stream list type, freeing all memory that was created as
    !>  part of the list.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_list_destroy(list, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), pointer:: list
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: node


        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_destroy()')

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        if (.not. associated(list)) return

        if (list % nItems == 0) then
            deallocate(list)
            return
        end if

        node => list % head
        do while (associated(node))
            list % head => list % head % next
            deallocate(node)
            node => list % head
        end do

        deallocate(list)

    end subroutine MPAS_stream_list_destroy !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_insert
    !
    !> \brief Add a stream to a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Adds a stream to the list of streams, first makes sure the stream doesn't exist in the stream list.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_list_insert(list, stream, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: list
        type (MPAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: node


        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_insert()')

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        nullify(stream % next)

        if (.not. associated(list % head)) then
            list % head => stream
        else
            node => list % head
            do while (associated(node % next))
                if (node % name == stream % name) then
                    if (present(ierr)) ierr = MPAS_STREAM_LIST_DUPLICATE
                    LIST_ERROR_WRITE('Found duplicate item '//trim(stream % name)//' in list.')
                    return
                end if
                node => node % next
            end do
            node % next => stream
        end if

        list % nItems = list % nItems + 1

    end subroutine MPAS_stream_list_insert !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_remove
    !
    !> \brief Remove a stream from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Removes a stream from the list of streams. Exits if the stream doesn't exist.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_list_remove(list, streamName, stream, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: list
        character (len=*), intent(in) :: streamName
        type (MPAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: node, nodePrev


        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_remove()')

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        ! Return if no streams exist in stream list
        if ( list % nItems == 0 ) then
            if (present(ierr)) ierr = MPAS_STREAM_LIST_NOT_FOUND
            LIST_ERROR_WRITE('Item '//trim(streamName)//' not found in list.')
            nullify(stream)
            return
        end if

        ! Check the head of the stream list
        node => list % head
        if (associated(node)) then
            if (node % name == streamName) then
                list % head => node % next
                stream => node
                list % nItems = list % nItems - 1
                return
            end if
        end if

        ! Loop through all streams until we find the one with StreamName
        nodePrev => node
        node => node % next
        do while (associated(node))
            if (node % name == streamName) then
                nodePrev % next => node % next
                stream => node
                list % nItems = list % nItems - 1
                return
            end if

            nodePrev => node
            node => node % next
        end do

        ! If the routine hasn't returned yet, the stream was not found. Return an error.
        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOT_FOUND
        LIST_ERROR_WRITE('Item '//trim(streamName)//' not found in list.')
        nullify(stream)

    end subroutine MPAS_stream_list_remove !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_query
    !
    !> \brief Get streams from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  This function allows one to search through a stream list for a set of
    !>  matching streams. In the first call to this function, the stream
    !>  argument should be a null pointer. When called a single time, this
    !>  function sets the stream argument to the first stream in the list that
    !>  matches the streamPattern argument. Upon subsequent calls, searching
    !>  begins at the first stream *after* the stream argument, to allow one to
    !>  find all streams that match streamPattern.  Each time a valid stream is
    !>  found, this function returns .true. and if no more streams are found, it
    !>  will return .false.
    !
    !-----------------------------------------------------------------------
    logical function MPAS_stream_list_query(list, streamPattern, stream, ierr) result(found) !{{{

        use iso_c_binding, only: c_char, c_int
        use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

        implicit none

        type (MPAS_stream_list_type), intent(in) :: list
        character (len=*), intent(in) :: streamPattern
        type (MPAS_stream_list_type), pointer :: stream
        integer, intent(out), optional :: ierr

        character(kind=c_char), dimension(StrKIND+1) :: c_pattern
        character(kind=c_char), dimension(StrKIND+1) :: c_test_string
        integer(kind=c_int) :: c_match

        type (MPAS_stream_list_type), pointer :: node

        interface
           subroutine check_regex_match(pattern, string, imatch) bind(c)
              use iso_c_binding, only: c_char, c_int
              character(kind=c_char), dimension(*), intent(in) :: pattern, string
              integer(kind=c_int), intent(out) :: imatch
           end subroutine check_regex_match
        end interface

        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_query()')

        call mpas_f_to_c_string(streamPattern, c_pattern)

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR
        found = .false.

        ! Return if no streams exist in stream list
        if ( list % nItems == 0 ) then
            LIST_DEBUG_WRITE(' -- No items matching '//trim(streamPattern)//' found in list.')
            nullify(stream)
            return
        end if

        if ( associated(stream) ) then
           if ( associated(stream % next) ) then
              node => stream % next
           else
              nullify(node)
           end if
        else
           node => list % head
        end if

        do while (associated(node))
            call mpas_f_to_c_string(node % name, c_test_string)

            call check_regex_match(c_pattern, c_test_string, c_match)

            if ( c_match == -1 ) then
                write(stderrUnit, *) 'ERROR: Regular expression matching failed.'
                write(stderrUnit, *) '       Pattern was: ', trim(streamPattern)
                write(stderrUnit, *) '       Test string was: ', trim(node % name)
            end if

            if ( c_match == 1 ) then
                found = .true.
                stream => node
                return
            end if
            node => node % next
        end do

        LIST_DEBUG_WRITE(' -- No items matching '//trim(streamPattern)//' found in list.')
        nullify(stream)

    end function MPAS_stream_list_query !}}}


    !-----------------------------------------------------------------------
    !  routine printlist
    !
    !> \brief Prints the contents of a list.
    !> \author Michael Duda
    !> \date   25 August 2014
    !> \details
    !>  Traverses a list, printing the 'name' component from each node.
    !
    !-----------------------------------------------------------------------
    subroutine printlist(list)!{{{

        implicit none

        type (MPAS_stream_list_type), intent(in) :: list

        type (MPAS_stream_list_type), pointer :: node
        integer :: i

        i = 1
        write(stderrUnit,*) '----------------------------'
        write(stderrUnit,*) 'List contains:'
        node => list % head
        do while (associated(node))
            write(stderrUnit,'(a,i3,a)') '    ', i, ') '//trim(node % name)
            i = i + 1
            node => node % next
        end do
        write(stderrUnit,*) '----------------------------'

    end subroutine printlist!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_list_length
    !
    !> \brief Returns the length of a stream list
    !> \author Michael Duda
    !> \date   25 August 2014
    !> \details
    !>  Returns the number of items stored in a stream list.
    !
    !-----------------------------------------------------------------------
    integer function MPAS_stream_list_length(list, ierr) result(nItems) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(in) :: list
        integer, intent(out), optional :: ierr


        LIST_DEBUG_WRITE(' -- Called MPAS_stream_list_length()')

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        nItems = list % nItems

    end function MPAS_stream_list_length !}}}


end module mpas_stream_list

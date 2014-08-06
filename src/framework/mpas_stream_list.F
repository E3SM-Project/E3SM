module mpas_stream_list

    use mpas_kind_types
    use mpas_grid_types
    use mpas_timekeeping
    use mpas_io_units
    use mpas_io_streams

    implicit none
    public

    integer, parameter :: MPAS_STREAM_LIST_NOERR = 0 ,&
                          MPAS_STREAM_LIST_FATAL = 1

    type mpas_stream_list_node_type
        character (len=StrKIND) :: name
        type (mpas_stream_type), pointer :: stream => null()
        type (mpas_stream_list_node_type), pointer :: next => null()
    end type mpas_stream_list_node_type

    type mpas_stream_list_type
        integer :: nStreams
        type (mpas_stream_list_node_type), pointer :: head => null()
    end type mpas_stream_list_type


    contains


    !-----------------------------------------------------------------------
    !  routine mpas_stream_list_init
    !
    !> \brief Initialize a new MPAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Instantiates and initializes a stream_list type, to store all active streams.
    !
    !-----------------------------------------------------------------------
    subroutine mpaS_stream_list_init(list, ierr)!{{{

        type (mpas_stream_list_type), pointer :: list
        integer, intent(out), optional :: ierr

        allocate(list)
        list % nStreams = 0
        nullify(list % head)

    end subroutine mpas_stream_list_init!}}}


    !-----------------------------------------------------------------------
    !  routine mpas_stream_list_finalize
    !
    !> \brief Free all memory associated with an MPAS stream list.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Destroys a stream list type, freeing all memory that was created as
    !>  part of the list.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_stream_list_finalize(list, ierr)!{{{

        type (mpas_stream_list_type), pointer:: list
        integer, intent(out), optional :: ierr

        type (mpas_stream_list_node_type), pointer :: node

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        if (list % nStreams > 0) then
            do while(associated(list % head))
                node => list % head

                if ( associated(node % next) ) then
                    list % head => node % next
                else
                    nullify(list % head)
                end if

                call mpas_closeStream(node % stream)
                deallocate(node)
            end do
        end if


    end subroutine mpas_stream_list_finalize!}}}


    !-----------------------------------------------------------------------
    !  routine mpas_stream_list_add_stream
    !
    !> \brief Add a stream to a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Adds a stream to the list of streams, first makes sure the stream doesn't exist in the stream list.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_stream_list_add_stream(list, stream, streamName, ierr)!{{{
        type (mpas_stream_list_type), pointer :: list
        type (mpas_stream_type), pointer :: stream
        character (len=*), intent(in) :: streamName
        integer, intent(out), optional :: ierr

        type (mpas_stream_list_node_type), pointer :: node, nodePrev

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        ! Check that no stream exists with the name already.
        node => list % head
        nodePrev => node
        do while(associated(node))
            if (trim(node % name) == trim(streamName)) then
                if (present(ierr)) ierr = MPAS_STREAM_LIST_FATAL
                write(stderrUnit, *) 'ERROR: mpas_stream_list_add_stream -- Stream already exists.'
                return
            end if

            nodePrev => node
            node => node % next
        end do

        allocate(node)
        node % name = streamName
        node % stream => stream
        nullify(node % next)

        nodePrev % next => node

        list % nStreams = list % nStreams + 1

    end subroutine mpas_stream_list_add_stream!}}}

    !-----------------------------------------------------------------------
    !  routine mpas_stream_list_remove_stream
    !
    !> \brief Remove a stream from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Removes a stream from the list of streams. Exits if the stream doesn't exist.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_stream_list_remove_stream(list, streamName, ierr)!{{{
        type (mpas_stream_list_type), pointer :: list
        character (len=*), intent(in) :: streamName
        integer, intent(out), optional :: ierr

        type (mpas_stream_list_node_type), pointer :: node, nodePrev

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        ! Check the head of the stream list
        if (associated(list % head)) then
            if (trim(list % head % name) == trim(streamName)) then
                node => list % head
                list % head => node % next
                call mpas_closeStream(node % stream)
                deallocate(node)
                list % nStreams = list % nStreams - 1
                return
            end if
        end if

        ! If the head was not the correct stream, start looking at the second stream.
        if (associated(list % head % next)) then
            node => list % head % next
            nodePrev => list % head
        else
            nullify(node)
            nullify(nodePrev)
        end if

        ! Loop through all streams until we find the one with StreamName
        do while(associated(node))
            if (trim(node % name) == trim(streamName)) then
                nodePrev % next => node % next
                call mpas_closeStream(node % stream)
                deallocate(node)
                list % nStreams = list % nStreams - 1
                return
            end if

            nodePrev => node
            node => node % next
        end do

        ! If the routine hasn't returned yet, the stream was not found. Return an error.
        if (present(ierr)) ierr = MPAS_STREAM_LIST_FATAL

    end subroutine mpas_stream_list_remove_stream!}}}

    !-----------------------------------------------------------------------
    !  routine mpas_stream_list_get_stream
    !
    !> \brief Get a stream from a stream list
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/06/2014
    !> \details
    !>  Searches through a stream list, and returns a pointer for the stream with a matching name.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_stream_list_get_stream(list, streamName, stream, ierr) !{{{
        type (mpas_stream_list_type), pointer :: list
        character (len=*), intent(in) :: streamName
        type (mpas_stream_type), pointer :: stream
        integer, intent(out), optional :: ierr

        type (mpas_stream_list_node_type), pointer :: node

        if (present(ierr)) ierr = MPAS_STREAM_LIST_NOERR

        ! Return if no streams exist in stream list
        if ( list % nStreams == 0 ) then
            write(stderrUnit, *) 'ERROR: mpas_stream_list_get_stream -- no streams exist in stream list yet.'
            if (present(ierr)) ierr = MPAS_STREAM_LIST_FATAL
            return
        end if

        ! Search through existing streams, returning a pointer to one that has a matching name.
        node => list % head
        do while(associated(node))
            if (trim(node % name) == trim(streamName)) then
                stream => node % stream
                return
            end if
            node => node % next
        end do

        ! If the routine hasn't returned yet, stream is not found.
        write(stderrUnit, *) 'ERROR: mpas_stream_list_get_stream -- no streams exist in stream list with the name ', trim(streamName), '.'
    end subroutine mpas_stream_list_get_stream!}}}


end module mpas_stream_list

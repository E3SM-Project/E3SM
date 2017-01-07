      threadNum = mpas_threading_get_thread_num()

      if ( threadNum == 0 ) then
         if (present(copy_array_only)) then
            local_copy_only = copy_array_only
         else
            local_copy_only = .false.
         end if
   
         
         src_cursor => src
         if (.not. local_copy_only) then
            nullify(dst_cursor)
         else
            dst_cursor => dst
         end if
   
!        do while (associated(src_cursor))
   
            if (.not. local_copy_only) then
               if (associated(dst_cursor)) then
                  allocate(dst_cursor % next)
                  dst_cursor % next % prev => dst_cursor
                  dst_cursor => dst_cursor % next
               else
                  allocate(dst)
                  nullify(dst % prev)
                  dst_cursor => dst
               end if
               nullify(dst_cursor % next)
            end if
   
   
            !
            ! Fill in members of dst_cursor from src_cursor
            !
            if (.not. local_copy_only) then
               dst_cursor % block => src_cursor % block
               dst_cursor % fieldName = src_cursor % fieldName
               dst_cursor % isVarArray = src_cursor % isVarArray
               dst_cursor % isActive = src_cursor % isActive
               dst_cursor % isDecomposed = src_cursor % isDecomposed
               dst_cursor % hasTimeDimension = src_cursor % hasTimeDimension
               dst_cursor % sendList => src_cursor % sendList
               dst_cursor % recvList => src_cursor % recvList
               dst_cursor % copyList => src_cursor % copyList

               if ( associated(dst_cursor % attLists) ) then
                  deallocate(dst_cursor % attLists)
               end if

               allocate(dst_cursor % attLists( size(src_cursor % attLists) ) )

               do iConstituent = 1, size(src_cursor % attLists)
                  call mpas_duplicate_attlist(src_cursor % attLists(iConstituent) % attList, dst_cursor % attLists(iConstituent) % attList)
               end do
            end if
            if ( dst_cursor % isActive .and. src_cursor % isActive ) then
               dst_cursor % scalar = src_cursor % scalar
            end if

!           src_cursor => src_cursor % next
!           if (.not. local_copy_only) then
!              dst_cursor => dst_cursor % next
!           end if
     
!        end do
      end if

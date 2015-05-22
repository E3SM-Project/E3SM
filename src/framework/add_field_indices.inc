      idim = ndims
      totalDimSize = 0
      field_ptr => field
      if ( field % isDecomposed ) then
         if (trim(field % dimNames(idim)) == 'nCells') then
!write(0,*) '... outer dimension is nCells'
            allocate(indices(0))
            do while (associated(field_ptr))
               call mpas_pool_get_array(field_ptr % block % allFields, 'indexToCellID', indexArray)
               call mpas_pool_get_dimension(field_ptr % block % dimensions, 'nCellsSolve', indexDimension)

               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension

               field_ptr => field_ptr % next
            end do
            call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         else if (trim(field % dimNames(idim)) == 'nEdges') then
!write(0,*) '... outer dimension is nEdges'
            allocate(indices(0))
            do while (associated(field_ptr))
               call mpas_pool_get_array(field_ptr % block % allFields, 'indexToEdgeID', indexArray)
               call mpas_pool_get_dimension(field_ptr % block % dimensions, 'nEdgesSolve', indexDimension)

               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension

               field_ptr => field_ptr % next
            end do
            call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         else if (trim(field % dimNames(idim)) == 'nVertices') then
!write(0,*) '... outer dimension is nVertices'
            allocate(indices(0))
            do while (associated(field_ptr))
               call mpas_pool_get_array(field_ptr % block % allFields, 'indexToVertexID', indexArray)
               call mpas_pool_get_dimension(field_ptr % block % dimensions, 'nVerticesSolve', indexDimension)

               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension

               field_ptr => field_ptr % next
            end do
            call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         else ! Use defined decomposition
            allocate(indices(0))
            do while (associated(field_ptr))
               call mpas_pool_get_array(field_ptr % block % allFields, trim(field_ptr % dimNames(idim))//'OwnedIndices', indexArray)

               call mpas_pool_get_dimension(field_ptr % block % dimensions, trim(field_ptr % dimNames(idim)), indexDimension)

               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension


               field_ptr => field_ptr % next
            end do
            call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         end if
      else
         globalDimSize = field % dimSizes(idim)
         totalDimSize = globalDimSize

         if (field % block % domain % dminfo % my_proc_id == IO_NODE) then
            allocate(indices(field % dimSizes(ndims)))
            do i=1,field % dimSizes(ndims)
               indices(i) = i
            end do
         else
            allocate(indices(0))
         end if
      end if

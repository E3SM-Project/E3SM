      idim = ndims
      totalDimSize = 0
      field_ptr => field
      if ( field % isDecomposed ) then
         meshFieldDim = .false.
         if (trim(field % dimNames(idim)) == 'nCells') then
            elementName = 'Cell'
            elementNamePlural = 'Cells'
            meshFieldDim = .true.
         else if (trim(field % dimNames(idim)) == 'nEdges') then
            elementName = 'Edge'
            elementNamePlural = 'Edges'
            meshFieldDim = .true.
         else if (trim(field % dimNames(idim)) == 'nVertices') then
            elementName = 'Vertex'
            elementNamePlural = 'Vertices'
            meshFieldDim = .true.
         end if

         if ( meshFieldDim ) then
            allocate(indices(0))
            if ( .not. stream % blockWrite ) then
               do while (associated(field_ptr))
                  call mpas_pool_get_array(field_ptr % block % allFields, 'indexTo' // trim(elementName) // 'ID', indexArray)
                  call mpas_pool_get_dimension(field_ptr % block % dimensions, 'n' // trim(elementNamePlural) // 'Solve', &
                                               indexDimension)

                  call mergeArrays(indices, indexArray(1:indexDimension))
                  totalDimSize = totalDimSize + indexDimension

                  field_ptr => field_ptr % next
               end do
               call mpas_dmpar_sum_int(stream % fileHandle % ioContext % dminfo, totalDimSize, globalDimSize)
            else
               call mpas_pool_get_dimension(field_ptr % block % dimensions, 'n' // trim(elementNamePlural), &
                                            indexDimension)

               allocate(indexArray(indexDimension))
               do i = 1, indexDimension
                  indexArray(i) = i
               end do
               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension
               globalDimSize = totalDimSize

               deallocate(indexArray)
            end if
         else ! Use defined decomposition
            allocate(indices(0))
            do while (associated(field_ptr))
               call mpas_pool_get_array(field_ptr % block % allFields, trim(field_ptr % dimNames(idim))//'OwnedIndices', indexArray)

               call mpas_pool_get_dimension(field_ptr % block % dimensions, trim(field_ptr % dimNames(idim)), indexDimension)

               call mergeArrays(indices, indexArray(1:indexDimension))
               totalDimSize = totalDimSize + indexDimension

               field_ptr => field_ptr % next
            end do
            call mpas_dmpar_sum_int(stream % fileHandle % ioContext % dminfo, totalDimSize, globalDimSize)
         end if
      else
         globalDimSize = field % dimSizes(idim)
         totalDimSize = globalDimSize

         if (stream % fileHandle % ioContext % dminfo % my_proc_id == IO_NODE) then
            allocate(indices(field % dimSizes(ndims)))
            do i=1,field % dimSizes(ndims)
               indices(i) = i
            end do
         else
            allocate(indices(0))
         end if
      end if

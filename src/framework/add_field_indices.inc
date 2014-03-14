      idim = ndims
      totalDimSize = 0
      field_ptr => field
      if (trim(field % dimNames(idim)) == 'nCells') then
!write(0,*) '... outer dimension is nCells'
         allocate(indices(0))
         do while (associated(field_ptr))
            call mergeArrays(indices, field_ptr % block % mesh % indexToCellID % array(1:field_ptr % block % mesh % nCellsSolve))
            totalDimSize = totalDimSize + field_ptr % block % mesh % nCellsSolve
            field_ptr => field_ptr % next
         end do
         call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         isDecomposed = .true.
      else if (trim(field % dimNames(idim)) == 'nEdges') then
!write(0,*) '... outer dimension is nEdges'
         allocate(indices(0))
         do while (associated(field_ptr))
            call mergeArrays(indices, field_ptr % block % mesh % indexToEdgeID % array(1:field_ptr % block % mesh % nEdgesSolve))
            totalDimSize = totalDimSize + field_ptr % block % mesh % nEdgesSolve
            field_ptr => field_ptr % next
         end do
         call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         isDecomposed = .true.
      else if (trim(field % dimNames(idim)) == 'nVertices') then
!write(0,*) '... outer dimension is nVertices'
         allocate(indices(0))
         do while (associated(field_ptr))
            call mergeArrays(indices, field_ptr % block % mesh % indexToVertexID % array(1:field_ptr % block % mesh % nVerticesSolve))
            totalDimSize = totalDimSize + field_ptr % block % mesh % nVerticesSolve
            field_ptr => field_ptr % next
         end do
         call mpas_dmpar_sum_int(field % block % domain % dminfo, totalDimSize, globalDimSize)
         isDecomposed = .true.
      else
         isDecomposed = .false.
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

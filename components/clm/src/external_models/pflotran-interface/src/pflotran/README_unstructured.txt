1.  The unstructured grid is read in from a file:

  A.  ASCII txt file read in UGridRead()

    Format of unstructured grid file
    -----------------------------------------------------------------
    num_cells num_vertices  ! (integers)
    element_type1 vert1 vert2 vert3 ... vert8  ! for cell 1 (char, integers)
    element_type1 vert1 vert2 vert3 ... vert8  ! for cell 2
    ...
    ...
    element_type1 vert1 vert2 vert3 ... vert8  ! for cell num_cells
    xcoord ycoord zcoord ! coordinates of vertex 1 (real)
    xcoord ycoord zcoord ! coordinates of vertex 2 (real)
    ...
    xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
    -----------------------------------------------------------------

    element types:
      H = hexahedron with 8 vertices
      T = tetrahedron with 4 vertices
      W = wedge with 6 vertices
      P = pyramid with 5 vertices (not yet supported)

    Cells are read through p0 and distributed as evenly as possible across all
    processors with the remaining cells spread across remainder (the #) cores 
    in the communicator

    Same goes for vertices.

2.  The unstructured grid is decomposed in UGridDecompose()

  A.  Package cells/vertices/duals into a PETSc vec split into 3 sections for 
      each cell (cell id, vertex ids, duals) as follows:

      ! Information for each cell is packed in a strided petsc vec
      ! The information is ordered within each stride as follows:
      ! -cell_N   ! global cell id (negative indicates 1-based)
      ! -777      ! separator between cell id and vertex ids for cell_N
      ! vertex1   ! in cell_N
      ! vertex2
      ! ...
      ! vertexN   
      ! -888      ! separator between vertex and dual ids
      ! dual1     ! dual ids between cell_N and others
      ! dual2
      ! ...
      ! dualN     
      ! -999    ! separator indicating end of information for cell_N
      
      ! the purpose of -777, -888, and -999 is to allow one to use cells of 
      ! various geometry.  Currently, the max # vertices = 8 and max # duals 
      ! = 6. But this will be generalized in the future.

      The variable "stride" is based on the block size of each cell.

    1.  Load vertex ids in local arrays.
    2.  MPI_Exscan to get global offset for first cell on processor
    3.  Create and fill MPIAdj (adjacency matrix)
    4.  Call MatMeshToCellGraph() to produce the dual graph
    5.  Create and set matrix partitioning based on dual produced in 4.
      a.  MatPartitioningCreate()
      b.  MatPartitioningSetAdjacency()
    6.  Call MatPartitioningApply() to produce an IS (is_new) that indicates the
        processor id of each cell.
    7.  Calculate number of elements on each processor with 
        ISPartitioningCount().  Should be able to use ISGetLocalSize() to 
        accomplish the same....
    8.  Create a vec (elements_natural) of size num_cells_local_new*stride
    9.  Create a new IS with the new global number of each index (cell) on
        each processor (is_num)
    10. Create a Blocked IS (is_scatter geh:name not clear)) for a new PetscVec
        geh: The current code passes a non-strided vector of cell ids
             into ISCreateBlock(), but this conflicts with PETSc docs...?
    11. Create a vec (elements_old) of size num_cells_local_old*stride
    12. Using MatGetRowIJF90(), access the information in the Dual matrix and
        fill the old element vector
    13. Scatter the old element vector to the natural...is it really natural
        or new petsc global?
      a.  Create scatter context from is_scatter created above
      b.  Scatter elements_old -> elements_natural (geh: name appropriate?)
    14. Destroy elements_old vec.
    15. MPI_ExScan to get new global offset
    16. Reallocate array holding cell vertices (%cell_vertices_0) 
    17. Duplicate elements_natural to elements_petsc
    18. Store natural cell ids
      a.  Allocate array (%cell_ids_natural) holding cell natural ids
      b.  Extract natural cell ids from elements_natural (elements_natural)
    19. Create an AO (ao_natural_to_petsc)
      a.  Make an array of local ids
      b.  pass along with natural ids (both permuted to 0-based) to 
          AOCreateBasic()
      c.  clean up afterwards
    20. Take duals and convert to petsc ordering
      a.  Count up number of duals and cells and allocate temporary array
          geh: why include cell ids too?
      b.  Set values in array based on elements_natural vec
    21. Use AOApplicationToPetsc() with ao_natural_to_petsc to map natural to 
        petsc ordering
    22. Load into elements_petsc
    23. Count up the number of ghost cells and place in array
    24. Scott added code to avoid duplicates, but this will be slow for large
        grids.  Need to leave duplicates and remove at end.
    25. Sort ghost cells and assign back to duals based on local ids
    26. Sort vertex ids and assign to vertex ids natural
    27. Create strided vector and scatter vertices to new processors
    28. Create scatter/gather ISs and perform scatter/gather

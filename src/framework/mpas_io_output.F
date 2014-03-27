! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io_output

   use mpas_grid_types
   use mpas_dmpar
   use mpas_sort
   use mpas_configure
   use mpas_io_streams

   integer, parameter :: OUTPUT = 1
   integer, parameter :: RESTART = 2
   integer, parameter :: SFC = 3
 
   type io_output_object
      character (len=StrKIND) :: filename
      integer :: stream

      integer :: time

      type (MPAS_Stream_type) :: io_stream
   end type io_output_object

   private :: mpas_insert_string_suffix

   contains

   subroutine mpas_output_state_init(output_obj, domain, stream, outputSuffix)!{{{

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(inout) :: domain
      character (len=*) :: stream
      character (len=*), optional :: outputSuffix

      character (len=StrKIND) :: tempfilename

      type (block_type), pointer :: block_ptr

      character (len=StrKIND), pointer :: config_output_name, config_restart_name, config_restart_timestamp_name

      call mpas_pool_get_config(domain % configs, 'config_output_name', config_output_name)
      call mpas_pool_get_config(domain % configs, 'config_restart_name', config_restart_name)
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      block_ptr => domain % blocklist

      if (trim(stream) == 'OUTPUT') then
         if(present(outputSuffix)) then
            call mpas_insert_string_suffix(config_output_name, outputSuffix, tempfilename)
         else
            tempfilename = config_output_name
         end if
         output_obj % filename = trim(tempfilename)
         output_obj % stream = OUTPUT
      else if (trim(stream) == 'RESTART') then
         if(present(outputSuffix)) then
            call mpas_insert_string_suffix(config_restart_name, outputSuffix, tempfilename)
            open(22,file=trim(config_restart_timestamp_name),form='formatted',status='replace')
            write(22,*) outputSuffix
            close(22)
         else
            tempfilename = config_restart_name
         end if
         output_obj % filename = trim(tempfilename)
         output_obj % stream = RESTART
      else if (trim(stream) == 'SFC') then
         ! Keep filename as whatever was set by the user
         output_obj % stream = SFC
      end if

      ! For now, we assume that a domain consists only of one block,
      !   although in future, work needs to be done to write model state
      !   from many distributed blocks
      call mpas_io_output_init(domain, output_obj, domain % dminfo)

   end subroutine mpas_output_state_init!}}}

   subroutine mpas_insert_string_suffix(stream, suffix, filename)!{{{

      implicit none

      character (len=*), intent(in) :: stream
      character (len=*), intent(in) :: suffix
      character (len=*), intent(out) :: filename
      integer :: length, i

      filename = trim(stream) // '.' // trim(suffix)

      length = len_trim(stream)
      do i=length-1,1,-1
         if(stream(i:i) == '.') then
            filename = trim(stream(:i)) // trim(suffix) // trim(stream(i:))
            exit
         end if
      end do

      do i=1,len_trim(filename)
         if (filename(i:i) == ':') filename(i:i) = '.'
      end do

   end subroutine mpas_insert_string_suffix!}}}

   subroutine mpas_output_state_for_domain(output_obj, domain, itime)!{{{
   
      implicit none
   
      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itime

      type (block_type), pointer :: block_ptr

      integer, pointer :: nCells, nEdges, nVertices, vertexDegree
      integer, pointer :: maxEdges, maxEdges2, nEdgesSolve, nCellsSolve, nVerticesSolve
      integer :: ierr
      integer :: i, j
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, indexToCellID, indexToEdgeID, indexToVertexID
      type (field2dInteger), pointer :: cellsOnCell_save, edgesOnCell_save, verticesOnCell_save, &
                               cellsOnEdge_save, verticesOnEdge_save, edgesOnEdge_save, &
                               cellsOnVertex_save, edgesOnVertex_save

      type (field2dInteger), pointer :: cellsOnCell_ptr, edgesOnCell_ptr, verticesOnCell_ptr, &
                               cellsOnEdge_ptr, verticesOnEdge_ptr, edgesOnEdge_ptr, &
                               cellsOnVertex_ptr, edgesOnVertex_ptr

      type (field2dInteger), pointer :: cellsOnCell, edgesOnCell, verticesOnCell, &
                               cellsOnEdge, verticesOnEdge, edgesOnEdge, &
                               cellsOnVertex, edgesOnVertex

      output_obj % time = itime

      !
      ! Convert connectivity information from local to global indices
      ! Needs to be done block by block
      !
      ! Also, backup local indices to be copied back into blocks after output is complete.
      !
      allocate(cellsOnCell_save)
      allocate(edgesOnCell_save) 
      allocate(verticesOnCell_save)
      allocate(cellsOnEdge_save)
      allocate(verticesOnEdge_save)
      allocate(edgesOnEdge_save)
      allocate(cellsOnVertex_save)
      allocate(edgesOnVertex_save)

      cellsOnCell_ptr => cellsOnCell_save
      edgesOnCell_ptr => edgesOnCell_save 
      verticesOnCell_ptr => verticesOnCell_save
      cellsOnEdge_ptr => cellsOnEdge_save 
      verticesOnEdge_ptr => verticesOnEdge_save 
      edgesOnEdge_ptr => edgesOnEdge_save
      cellsOnVertex_ptr => cellsOnVertex_save 
      edgesOnVertex_ptr => edgesOnVertex_save

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_dimension(block_ptr % dimensions, 'maxEdges', maxEdges)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'maxEdges2', maxEdges2)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'vertexDegree', vertexDegree)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdges', nEdges)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertices', nVertices)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_dimension(block_ptr % dimensions, 'nVerticesSolve', nVerticesSolve)

        call mpas_pool_get_array(block_ptr % allFields, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(block_ptr % allFields, 'nEdgesOnEdge', nEdgesOnEdge)
        call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID', indexToCellID)
        call mpas_pool_get_array(block_ptr % allFields, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_array(block_ptr % allFields, 'indexToVertexID', indexToVertexID)

        nullify(cellsOncell_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnCell', cellsOnCell)
        cellsOnCell_ptr % array => cellsOnCell % array
        allocate(cellsOnCell % array(maxEdges, nCells+1))

        nullify(edgesOnCell_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnCell', edgesOnCell)
        edgesOnCell_ptr % array => edgesOnCell % array
        allocate(edgesOnCell % array(maxEdges, nCells+1))

        nullify(verticesOnCell_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'verticesOnCell', verticesOnCell)
        verticesOnCell_ptr % array => verticesOnCell % array
        allocate(verticesOnCell % array(maxEdges, nCells+1))

        nullify(cellsOnEdge_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnEdge', cellsOnEdge)
        cellsOnEdge_ptr % array => cellsOnEdge % array
        allocate(cellsOnEdge % array(2, nEdges+1))

        nullify(verticesOnEdge_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'verticesOnEdge', verticesOnEdge)
        verticesOnEdge_ptr % array => verticesOnEdge % array
        allocate(verticesOnEdge % array(2, nEdges+1))

        nullify(edgesOnEdge_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnEdge', edgesOnEdge)
        edgesOnEdge_ptr % array => edgesOnEdge % array
        allocate(edgesOnEdge % array(maxEdges2, nEdges+1))

        nullify(cellsOnVertex_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnVertex', cellsOnVertex)
        cellsOnVertex_ptr % array => cellsOnVertex % array
        allocate(cellsOnVertex % array(vertexDegree, nVertices+1))

        nullify(edgesOnVertex_ptr % ioinfo)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnVertex', edgesOnVertex)
        edgesOnVertex_ptr % array => edgesOnVertex % array
        allocate(edgesOnVertex % array(vertexDegree, nVertices+1))

        do i = 1, nCellsSolve
          do j = 1, nEdgesOnCell(i)
            cellsOnCell % array(j, i) = indexToCellID(cellsOnCell_ptr % array(j, i))
            edgesOnCell % array(j, i) = indexToEdgeID(edgesOnCell_ptr % array(j, i))
            verticesOnCell % array(j, i) = indexToVertexID(verticesOnCell_ptr % array(j, i))
          end do

          cellsOnCell % array(nEdgesOnCell(i) + 1:maxEdges, i) = nCells+1
          edgesOnCell % array(nEdgesOnCell(i) + 1:maxEdges, i) = nEdges+1
          verticesOnCell % array(nEdgesOnCell(i) + 1:maxEdges, i) = nVertices+1
        end do

        do i = 1, nEdgesSolve
          cellsOnEdge % array(1, i) = indexToCellID(cellsOnEdge_ptr % array(1, i))
          cellsOnEdge % array(2, i) = indexToCellID(cellsOnEdge_ptr % array(2, i))

          verticesOnedge % array(1, i) = indexToVertexID(verticesOnEdge_ptr % array(1,i))
          verticesOnedge % array(2, i) = indexToVertexID(verticesOnEdge_ptr % array(2,i))

          do j = 1, nEdgesOnEdge(i)
            edgesOnEdge % array(j, i) = indexToEdgeID(edgesOnEdge_ptr % array(j, i))
          end do

          edgesOnEdge % array(nEdgesOnEdge(i)+1:maxEdges2, i) = nEdges+1
        end do

        do i = 1, nVerticesSolve
          do j = 1, vertexDegree
            cellsOnVertex % array(j, i) = indexToCellID(cellsOnVertex_ptr % array(j, i))
            edgesOnVertex % array(j, i) = indexToEdgeID(edgesOnVertex_ptr % array(j, i))
          end do
        end do

        block_ptr => block_ptr % next
        if(associated(block_ptr)) then
          allocate(cellsOnCell_ptr % next)
          allocate(edgesOnCell_ptr % next)
          allocate(verticesOnCell_ptr % next)
          allocate(cellsOnEdge_ptr % next)
          allocate(verticesOnEdge_ptr % next)
          allocate(edgesOnEdge_ptr % next)
          allocate(cellsOnVertex_ptr % next)
          allocate(edgesOnVertex_ptr % next)

          cellsOnCell_ptr => cellsOnCell_ptr % next
          edgesOnCell_ptr => edgesOnCell_ptr % next
          verticesOnCell_ptr => verticesOnCell_ptr % next
          cellsOnEdge_ptr => cellsOnEdge_ptr % next
          verticesOnEdge_ptr => verticesOnEdge_ptr % next
          edgesOnEdge_ptr => edgesOnEdge_ptr % next
          cellsOnVertex_ptr => cellsOnVertex_ptr % next
          edgesOnVertex_ptr => edgesOnVertex_ptr % next
        end if

        nullify(cellsOnCell_ptr % next)
        nullify(edgesOnCell_ptr % next)
        nullify(verticesOnCell_ptr % next)
        nullify(cellsOnEdge_ptr % next)
        nullify(verticesOnEdge_ptr % next)
        nullify(edgesOnEdge_ptr % next)
        nullify(cellsOnVertex_ptr % next)
        nullify(edgesOnVertex_ptr % next)
      end do

      ! Write output file
      call MPAS_writeStream(output_obj % io_stream, output_obj % time, ierr)

      ! Converge indices back to local indices, and deallocate all temporary arrays.
      cellsOnCell_ptr => cellsOnCell_save
      edgesOnCell_ptr => edgesOnCell_save 
      verticesOnCell_ptr => verticesOnCell_save
      cellsOnEdge_ptr => cellsOnEdge_save 
      verticesOnEdge_ptr => verticesOnEdge_save 
      edgesOnEdge_ptr => edgesOnEdge_save
      cellsOnVertex_ptr => cellsOnVertex_save 
      edgesOnVertex_ptr => edgesOnVertex_save

      block_ptr => domain % blocklist
      do while(associated(block_ptr))

        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnCell', edgesOnCell)
        call mpas_pool_get_field(block_ptr % allFields, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_field(block_ptr % allFields, 'verticesOnEdge', verticesOnEdge)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnEdge', edgesOnEdge)
        call mpas_pool_get_field(block_ptr % allFields, 'cellsOnVertex', cellsOnVertex)
        call mpas_pool_get_field(block_ptr % allFields, 'edgesOnVertex', edgesOnVertex)

        deallocate(cellsOnCell % array)
        deallocate(edgesOnCell % array)
        deallocate(verticesOnCell % array)
        deallocate(cellsOnEdge % array)
        deallocate(verticesOnEdge % array)
        deallocate(edgesOnEdge % array)
        deallocate(cellsOnVertex % array)
        deallocate(edgesOnVertex % array)

        cellsOncell % array => cellsOnCell_ptr % array
        edgesOnCell % array => edgesOnCell_ptr % array
        verticesOnCell % array => verticesOnCell_ptr % array
        cellsOnEdge % array => cellsOnEdge_ptr % array
        verticesOnEdge % array => verticesOnEdge_ptr % array
        edgesOnEdge % array => edgesOnEdge_ptr % array
        cellsOnVertex % array => cellsOnVertex_ptr % array
        edgesOnVertex % array => edgesOnVertex_ptr % array

        nullify(cellsOnCell_ptr % array)
        nullify(edgesOnCell_ptr % array)
        nullify(verticesOnCell_ptr % array)
        nullify(cellsOnEdge_ptr % array)
        nullify(verticesOnEdge_ptr % array)
        nullify(edgesOnEdge_ptr % array)
        nullify(cellsOnVertex_ptr % array)
        nullify(edgesOnVertex_ptr % array)

        block_ptr => block_ptr % next
        cellsOnCell_ptr => cellsOnCell_ptr % next
        edgesOnCell_ptr => edgesOnCell_ptr % next
        verticesOnCell_ptr => verticesOnCell_ptr % next
        cellsOnEdge_ptr => cellsOnEdge_ptr % next
        verticesOnEdge_ptr => verticesOnEdge_ptr % next
        edgesOnEdge_ptr => edgesOnEdge_ptr % next
        cellsOnVertex_ptr => cellsOnVertex_ptr % next
        edgesOnVertex_ptr => edgesOnVertex_ptr % next
      end do

      call mpas_deallocate_field(cellsOnCell_save)
      call mpas_deallocate_field(edgesOnCell_save) 
      call mpas_deallocate_field(verticesOnCell_save)
      call mpas_deallocate_field(cellsOnEdge_save)
      call mpas_deallocate_field(verticesOnEdge_save)
      call mpas_deallocate_field(edgesOnEdge_save)
      call mpas_deallocate_field(cellsOnVertex_save)
      call mpas_deallocate_field(edgesOnVertex_save)

   end subroutine mpas_output_state_for_domain!}}}

   subroutine mpas_output_state_finalize(output_obj, dminfo)!{{{

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo

      call mpas_io_output_finalize(output_obj, dminfo)

   end subroutine mpas_output_state_finalize!}}}

   subroutine mpas_io_output_init( domain, output_obj, dminfo)!{{{
      implicit none
 
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr, ierr
      integer, dimension(10) :: dimlist
      character (len=StrKIND*4) :: runCmd

      if(len(trim(domain % history)) > 0) then
          write(runCmd,'(a,i0,a,a,a,a)')  ' mpirun -n ', domain % dminfo % nProcs, ' ', trim(domain % executableName), new_line('A'), trim(domain % history)
      else
          write(runCmd,'(a,i0,a,a)') 'mpirun -n ',domain % dminfo % nProcs, ' ', trim(domain % executableName)
      end if
 
      call MPAS_createStream(output_obj % io_stream, trim(output_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_WRITE, 1, nferr)

      if (domain % on_a_sphere) then
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'YES             ', nferr)
      else
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'NO              ', nferr)
      end if
      call MPAS_writeStreamAtt(output_obj % io_stream, 'sphere_radius', domain % sphere_radius, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'model_name', domain % modelName, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'core_name', domain % coreName, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'model_version', domain % modelVersion, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'history', runCmd, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'source', domain % source, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'Conventions', domain % Conventions, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'mesh_id', domain % mesh_id, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'parent_id', domain % parent_id, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'mesh_spec', domain % mesh_spec, nferr)
      call MPAS_writeStreamAtt(output_obj % io_stream, 'git_version', domain % git_version, nferr)

      call mpas_add_output_atts(domain % configs, output_obj)

      call mpas_add_output_fields(domain, output_obj)
 
   end subroutine mpas_io_output_init!}}}

   subroutine mpas_io_output_finalize(output_obj, dminfo)!{{{
 
      implicit none
 
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr
 
      call MPAS_closeStream(output_obj % io_stream, nferr)
 
   end subroutine mpas_io_output_finalize!}}}


#include "add_output_atts.inc"
#include "add_output_fields.inc"
 
end module mpas_io_output

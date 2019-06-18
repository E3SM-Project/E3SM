module boundary_exchange_ut

  use edgetype_mod,           only: EdgeBuffer_t
  use kinds,                  only: real_kind

  implicit none

  type (EdgeBuffer_t) :: edge
  type (EdgeBuffer_t) :: edgeMinMax
contains

  subroutine init_edges_structs_f90 (num_min_max_fields_1d,num_scalar_fields_2d,num_scalar_fields_3d,num_scalar_fields_3d_int,num_vector_fields_3d,vector_dim) bind(c)
    use iso_c_binding,  only : c_int
    use dimensions_mod, only : nlev, nlevp, qsize
    use edge_mod_base,  only : initEdgeBuffer, initEdgeSBuffer
    use geometry_interface_mod, only: par, elem
    !
    ! Inputs
    !
    integer(kind=c_int), intent(in) :: num_min_max_fields_1d, num_scalar_fields_2d, num_scalar_fields_3d, num_scalar_fields_3d_int, num_vector_fields_3d, vector_dim

    call initEdgeBuffer(par,edge,elem,num_scalar_fields_2d + num_scalar_fields_3d*nlev + num_scalar_fields_3d_int*nlevp + vector_dim*num_vector_fields_3d*nlev)
    call initEdgeSBuffer(par,edgeMinMax,elem,num_min_max_fields_1d*nlev*2)

    qsize = num_min_max_fields_1d
  end subroutine init_edges_structs_f90

  subroutine boundary_exchange_test_f90 (field_min_1d_ptr, field_max_1d_ptr,      &
                                         field_2d_ptr, field_3d_ptr,              &
                                         field_3d_int_ptr, field_4d_ptr,          &
                                         inner_dim_4d, num_time_levels,           &
                                         idim_2d, idim_3d, idim_4d, minmax_split) bind(c)
    use iso_c_binding,      only : c_ptr, c_f_pointer, c_int
    use dimensions_mod,     only : np, nlev, nlevp, nelemd, qsize
    use edge_mod_base,      only : edgevpack, edgevunpack
    use bndry_mod,          only : bndry_exchangev
    use geometry_interface_mod, only: hybrid
    use viscosity_base,     only : neighbor_minmax, neighbor_minmax_start, neighbor_minmax_finish
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: field_min_1d_ptr, field_max_1d_ptr
    type (c_ptr), intent(in) :: field_2d_ptr, field_3d_ptr, field_3d_int_ptr, field_4d_ptr
    integer (kind=c_int), intent(in) :: inner_dim_4d, num_time_levels
    integer (kind=c_int), intent(in) :: idim_2d, idim_3d, idim_4d
    integer (kind=c_int), intent(in) :: minmax_split
    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:),       pointer :: field_min_1d
    real (kind=real_kind), dimension(:,:,:),       pointer :: field_max_1d
    real (kind=real_kind), dimension(:,:,:,:),     pointer :: field_2d
    real (kind=real_kind), dimension(:,:,:,:,:),   pointer :: field_3d
    real (kind=real_kind), dimension(:,:,:,:,:),   pointer :: field_3d_int
    real (kind=real_kind), dimension(:,:,:,:,:,:), pointer :: field_4d
    integer :: ie, kptr

    call c_f_pointer(field_min_1d_ptr, field_min_1d, [nlev,qsize,nelemd])
    call c_f_pointer(field_max_1d_ptr, field_max_1d, [nlev,qsize,nelemd])
    call c_f_pointer(field_2d_ptr,     field_2d,     [np,np,num_time_levels,nelemd])
    call c_f_pointer(field_3d_ptr,     field_3d,     [np,np,nlev,num_time_levels,nelemd])
    call c_f_pointer(field_3d_int_ptr, field_3d_int, [np,np,nlevp,num_time_levels,nelemd])
    call c_f_pointer(field_4d_ptr,     field_4d,     [np,np,nlev,inner_dim_4d,num_time_levels,nelemd])

    ! Perform 'standard' 2d/3d boundary exchange
    do ie=1,nelemd
      kptr = 0
      call edgeVpack(edge,field_2d(:,:,idim_2d,ie),1,kptr,ie)

      kptr = 1
      call edgeVpack(edge,field_3d(:,:,:,idim_3d,ie),nlev,kptr,ie)

      kptr = 1 + nlev
      call edgeVpack(edge,field_3d_int(:,:,:,idim_3d,ie),nlevp,kptr,ie)

      kptr = 1 + nlev + nlevp
      call edgeVpack(edge,field_4d(:,:,:,:,idim_4d,ie),inner_dim_4d*nlev,kptr,ie)
    enddo

    call bndry_exchangev(hybrid,edge)

    do ie=1,nelemd
      kptr = 0
      call edgeVunpack(edge,field_2d(:,:,idim_2d,ie),1,kptr,ie)

      kptr = 1
      call edgeVunpack(edge,field_3d(:,:,:,idim_3d,ie),nlev,kptr,ie)

      kptr = 1 + nlev
      call edgeVunpack(edge,field_3d_int(:,:,:,idim_3d,ie),nlevp,kptr,ie)

      kptr = 1 + nlev + nlevp
      call edgeVunpack(edge,field_4d(:,:,:,:,idim_4d,ie),inner_dim_4d*nlev,kptr,ie)
    enddo

    ! Perform min/max boundary exchange
    if (minmax_split .eq. 0) then
      call neighbor_minmax (hybrid, edgeMinMax, 1, nelemd, field_min_1d, field_max_1d)
    else
      call neighbor_minmax_start  (hybrid, edgeMinMax, 1, nelemd, field_min_1d, field_max_1d)
      call neighbor_minmax_finish (hybrid, edgeMinMax, 1, nelemd, field_min_1d, field_max_1d)
    endif
  end subroutine boundary_exchange_test_f90

  subroutine cleanup_f90 () bind(c)
    use edge_mod_base, only : FreeEdgeBuffer
    use geometry_interface_mod, only: cleanup_geometry_f90

    call FreeEdgeBuffer(edge)
    call cleanup_geometry_f90()

  end subroutine cleanup_f90

end module boundary_exchange_ut

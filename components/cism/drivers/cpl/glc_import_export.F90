module glc_import_export

  use shr_sys_mod
  use shr_kind_mod,        only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod,        only: CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use glc_constants,       only: verbose, stdout, stderr, tkfrz
  use glc_communicate,     only: my_task, master_task
  use glc_cpl_indices

  implicit none
  save
  public

  ! Public interfaces
  public :: glc_import

!=================================================================================
contains
!=================================================================================

   subroutine glc_import(x2g)

    !-------------------------------------------------------------------
     use glc_indexing_info, only : nx, ny, local_indices
     use glc_fields, only: tsfc, qsmb 

    real(r8)   , intent(in) :: x2g(:,:)

    integer(IN) :: i,j,n
    character(*), parameter :: subName = "(glc_import) "
    !-------------------------------------------------------------------

    do j = 1, ny
       do i = 1, nx
          n = local_indices(i,j)
          tsfc(i,j) = x2g(index_x2g_Sl_tsrf, n) - tkfrz
          qsmb(i,j) = x2g(index_x2g_Flgl_qice, n)
       enddo
    enddo

    !Jer hack fix: 
    !For some land points where CLM sees ocean, and all ocean points, CLM doesn't provide a temperature,
    !and so the incoming temperature is 0.d0.  This gets dropped to -273.15, in the above code.  So,
    !manually reverse this, below, to set to 0C.
    where (tsfc < -250.d0) tsfc=0.d0 

  end subroutine glc_import

!=================================================================================

  subroutine glc_export(g2x)

    !-------------------------------------------------------------------
    use glc_indexing_info, only : nx, ny, local_indices
    use glc_fields   , only: ice_covered, topo, rofi, rofl, hflx, &
                             ice_sheet_grid_mask, icemask_coupled_fluxes   ! to coupler
    use glc_route_ice_runoff, only: route_ice_runoff    
    use glc_override_frac   , only: frac_overrides_enabled, do_frac_overrides
    
    real(r8)    ,intent(inout) :: g2x(:,:)

    ! if doing frac overrides, these are the modified versions sent to the coupler;
    ! otherwise they point to the real fields
    real(r8), pointer :: ice_covered_to_cpl(:,:)
    real(r8), pointer :: topo_to_cpl(:,:)
    logical :: fields_to_cpl_allocated  ! whether we allocated the above fields

    integer(IN) :: i,j,n
    character(*), parameter :: subName = "(glc_export) "
    !-------------------------------------------------------------------

    ! If overrides of glc fraction are enabled (for testing purposes), then apply
    ! these overrides, otherwise use the real version of ice_covered and topo
    if (frac_overrides_enabled()) then
       allocate(ice_covered_to_cpl(lbound(ice_covered,1):ubound(ice_covered,1), &
                                   lbound(ice_covered,2):ubound(ice_covered,2)))
       allocate(topo_to_cpl(lbound(topo,1):ubound(topo,1), &
                            lbound(topo,2):ubound(topo,2)))
            
       ice_covered_to_cpl = ice_covered
       topo_to_cpl = topo
       call do_frac_overrides(ice_covered_to_cpl, topo_to_cpl, ice_sheet_grid_mask)
       fields_to_cpl_allocated = .true.
    else
       ice_covered_to_cpl => ice_covered
       topo_to_cpl => topo
       fields_to_cpl_allocated = .false.
    end if

    do j = 1, ny
       do i = 1, nx
          n = local_indices(i,j)

          call route_ice_runoff(rofi(i,j), &
               rofi_to_ocn=g2x(index_g2x_Fogg_rofi, n), &
               rofi_to_ice=g2x(index_g2x_Figg_rofi, n))
          
          g2x(index_g2x_Fogg_rofl, n) = rofl(i,j)

          g2x(index_g2x_Sg_ice_covered, n) = ice_covered_to_cpl(i,j)
          g2x(index_g2x_Sg_topo, n) = topo_to_cpl(i,j)
          g2x(index_g2x_Flgg_hflx, n) = hflx(i,j)
	  
          g2x(index_g2x_Sg_icemask, n) = ice_sheet_grid_mask(i,j)
          g2x(index_g2x_Sg_icemask_coupled_fluxes, n) = icemask_coupled_fluxes(i,j)
          
       enddo
    enddo

    if (fields_to_cpl_allocated) then
       deallocate(ice_covered_to_cpl)
       deallocate(topo_to_cpl)
    end if

  end subroutine glc_export

end module glc_import_export

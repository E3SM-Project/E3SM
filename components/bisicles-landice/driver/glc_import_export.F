module glc_import_export

  use shr_sys_mod
  use shr_kind_mod,        only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod,        only: CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use glc_constants,       only: verbose, stdout, stderr, tkfrz, zero_gcm_fluxes
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
     use glc_indexing, only : vector_to_spatial
     use glc_fields, only: tsfc, qsmb 

    real(r8)   , intent(in) :: x2g(:,:)

    character(*), parameter :: subName = "(glc_import) "
    !-------------------------------------------------------------------

    call vector_to_spatial(x2g(index_x2g_Sl_tsrf,:), tsfc)
    call vector_to_spatial(x2g(index_x2g_Flgl_qice,:), qsmb)

    tsfc = tsfc - tkfrz

    !Jer hack fix: 
    !For some land points where CLM sees ocean, and all ocean points, CLM doesn't provide a temperature,
    !and so the incoming temperature is 0.d0.  This gets dropped to -273.15, in the above code.  So,
    !manually reverse this, below, to set to 0C.
    where (tsfc < -250.d0) tsfc=0.d0 

  end subroutine glc_import

!=================================================================================

  subroutine glc_export(g2x)

    !-------------------------------------------------------------------
    use glc_indexing, only : nx, ny, spatial_to_vector
    use glc_fields   , only: ice_covered, topo, rofi, rofl, hflx, &
                             ice_sheet_grid_mask
    use glc_route_ice_runoff, only: route_ice_runoff    
    use glc_override_frac   , only: frac_overrides_enabled, do_frac_overrides
    
    real(r8)    ,intent(inout) :: g2x(:,:)

    ! if doing frac overrides, these are the modified versions sent to the coupler;
    ! otherwise they point to the real fields
    real(r8), pointer :: ice_covered_to_cpl(:,:)
    real(r8), pointer :: topo_to_cpl(:,:)
    logical :: fields_to_cpl_allocated  ! whether we allocated the above fields

    ! Note that there are two separate mask fields. Both of them provide information about
    ! where CISM is running. The difference is that ice_sheet_grid_mask includes icesheet
    ! areas that are diagtostic-only, whereas icemask_coupled_fluxes excludes icesheet
    ! areas where we are zeroing the fluxes sent to the coupler (thus, icesheets that are
    ! "diagnostic" in some sense). We need two separate maps, as opposed to a single map
    ! plus a scalar logical variable, in case we're running with multiple icesheet
    ! instances (e.g., Greenland & Antarctica), one of which is fully prognostic and one
    ! of which is diagnostic-only: in that case, ice_sheet_grid_mask would be non-zero
    ! over both Greenland and Antarctica, whereas icemask_coupled_fluxes would be non-zero
    ! over (e.g.)  Greenland, but 0 over Antarctica.
    real(r8), allocatable :: icemask_coupled_fluxes(:,:)  ! mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

    real(r8), allocatable :: hflx_to_cpl(:,:)
    real(r8), allocatable :: rofl_to_cpl(:,:)
    real(r8), allocatable :: rofi_to_ocn(:,:)
    real(r8), allocatable :: rofi_to_ice(:,:)

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

    allocate(icemask_coupled_fluxes(nx, ny))
    allocate(hflx_to_cpl(nx, ny))
    allocate(rofl_to_cpl(nx, ny))
    allocate(rofi_to_ocn(nx, ny))
    allocate(rofi_to_ice(nx, ny))

    if (zero_gcm_fluxes) then
       icemask_coupled_fluxes = 0._r8
       hflx_to_cpl = 0._r8
       rofl_to_cpl = 0._r8
       rofi_to_ocn = 0._r8
       rofi_to_ice = 0._r8
    else
       icemask_coupled_fluxes = ice_sheet_grid_mask
       hflx_to_cpl = hflx
       rofl_to_cpl = rofl
       call route_ice_runoff(rofi, rofi_to_ocn, rofi_to_ice)
    end if

    call spatial_to_vector(rofi_to_ocn, g2x(index_g2x_Fogg_rofi,:))
    call spatial_to_vector(rofi_to_ice, g2x(index_g2x_Figg_rofi,:))
    call spatial_to_vector(rofl_to_cpl, g2x(index_g2x_Fogg_rofl,:))

    call spatial_to_vector(ice_covered_to_cpl, g2x(index_g2x_Sg_ice_covered,:))
    call spatial_to_vector(topo_to_cpl, g2x(index_g2x_Sg_topo,:))
    call spatial_to_vector(hflx_to_cpl, g2x(index_g2x_Flgg_hflx,:))

    call spatial_to_vector(ice_sheet_grid_mask, g2x(index_g2x_Sg_icemask,:))
    call spatial_to_vector(icemask_coupled_fluxes, g2x(index_g2x_Sg_icemask_coupled_fluxes,:))

    deallocate(icemask_coupled_fluxes)
    deallocate(hflx_to_cpl)
    deallocate(rofl_to_cpl)
    deallocate(rofi_to_ocn)
    deallocate(rofi_to_ice)
    if (fields_to_cpl_allocated) then
       deallocate(ice_covered_to_cpl)
       deallocate(topo_to_cpl)
    end if

  end subroutine glc_export

end module glc_import_export

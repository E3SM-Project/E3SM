module glc_import_export

  use shr_sys_mod
  use shr_kind_mod,        only: IN=>SHR_KIND_IN
!  use shr_kind_mod,        only: CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
!  use glc_constants,       only: verbose, stdout, stderr, tkfrz, glc_nec
!  use glc_global_grid,     only: glc_grid
!  use glc_cpl_indices
  use mpas_kind_types

  implicit none
  save

  ! Public interfaces
  public :: glc_import_mct, glc_export_mct

!=================================================================================
contains
!=================================================================================

   subroutine glc_import_mct(x2g, errorCode)

    !-------------------------------------------------------------------
!    use glc_global_fields, only: tsfc, topo, qsmb 

    real(RKIND)   , intent(in) :: x2g(:,:)
    integer, intent(out) :: errorCode

    integer(IN) :: j,jj,i,g,nxg,nyg,n,elev_class
    character(*), parameter :: subName = "(glc_import) "
    !-------------------------------------------------------------------

    errorCode = 0

!    nxg = glc_grid%nx
!    nyg = glc_grid%ny
!    do j = 1, nyg           ! S to N
!       jj = nyg - j + 1     ! reverse j index for glint grid (N to S)
!       do i = 1, nxg
!          g = (j-1)*nxg + i   ! global index (W to E, S to N)
!          do elev_class = 0, glc_nec
!             tsfc(i,jj,elev_class) = x2g(index_x2g_Sl_tsrf(elev_class), g) - tkfrz
!             topo(i,jj,elev_class) = x2g(index_x2g_Sl_topo(elev_class), g)
!             qsmb(i,jj,elev_class) = x2g(index_x2g_Flgl_qice(elev_class), g)
!          enddo
!       enddo
!    enddo
!    
!    !Jer hack fix: 
!    !For some land points where CLM sees ocean, and all ocean points, CLM doesn't provide a temperature,
!    !and so the incoming temperature is 0.d0.  This gets dropped to -273.15, in the above code.  So,
!    !manually reverse this, below, to set to 0C.
!    where (tsfc < -250.d0) tsfc=0.d0 

  end subroutine glc_import_mct

!=================================================================================

  subroutine glc_export_mct(g2x, errorCode)

    !-------------------------------------------------------------------
 !   use glc_global_fields   , only: gfrac, gtopo, grofi, grofl, ghflx, &
 !                                   ice_sheet_grid_mask, icemask_coupled_fluxes   ! to coupler
 !   use glc_route_ice_runoff, only: route_ice_runoff    
 !   use glc_override_frac   , only: frac_overrides_enabled, do_frac_overrides
    
    real(RKIND)    ,intent(inout) :: g2x(:,:)
    integer, intent(out) :: errorCode

    real(RKIND), pointer :: gfrac_to_cpl(:,:,:)   ! if overriding gfrac, this is the modified version, 
                                               ! sent to the coupler; otherwise it points to the real gfrac
    logical :: gfrac_to_cpl_allocated          ! whether we allocated gfrac_to_cpl
    integer(IN) :: j,jj,i,g,nxg,nyg,n,elev_class
    character(*), parameter :: subName = "(glc_export) "
    !-------------------------------------------------------------------

    errorCode = 0

    ! If overrides of glc fraction are enabled (for testing purposes), then apply
    ! these overrides, otherwise use the real version of gfrac
 !   if (frac_overrides_enabled()) then
 !      allocate(gfrac_to_cpl(lbound(gfrac,1):ubound(gfrac,1), &
 !                            lbound(gfrac,2):ubound(gfrac,2), &
 !                            lbound(gfrac,3):ubound(gfrac,3)))
 !      gfrac_to_cpl = gfrac
 !      call do_frac_overrides(gfrac_to_cpl, ice_sheet_grid_mask)
 !      gfrac_to_cpl_allocated = .true.
 !   else
 !      gfrac_to_cpl => gfrac
 !      gfrac_to_cpl_allocated = .false.
 !   end if

 !   nxg = glc_grid%nx
 !   nyg = glc_grid%ny
 !   do j = 1, nyg           ! S to N
 !      jj = nyg - j + 1     ! reverse j index for glint grid (N to S)
 !      do i = 1, nxg
 !         g = (j-1)*nxg + i ! global index (W to E, S to N)
!
!          call route_ice_runoff(grofi(i,jj), &
!               rofi_to_ocn=g2x(index_g2x_Fogg_rofi, g), &
!               rofi_to_ice=g2x(index_g2x_Figg_rofi, g))
!          
!          g2x(index_g2x_Fogg_rofl, g) = grofl(i,jj)
!
!          do elev_class = 0, glc_nec !Jer
!             g2x(index_g2x_Sg_frac(elev_class), g) = gfrac_to_cpl(i,jj,elev_class)
!             g2x(index_g2x_Sg_topo(elev_class), g) = gtopo(i,jj,elev_class)
!             g2x(index_g2x_Flgg_hflx(elev_class), g) = ghflx(i,jj,elev_class)
!          enddo
	  
!          g2x(index_g2x_Sg_icemask, g) = ice_sheet_grid_mask(i,jj)
!          g2x(index_g2x_Sg_icemask_coupled_fluxes, g) = icemask_coupled_fluxes(i,jj)
          
!       enddo
!    enddo

!    if (gfrac_to_cpl_allocated) then
!       deallocate(gfrac_to_cpl)
!    end if

  end subroutine glc_export_mct

end module glc_import_export

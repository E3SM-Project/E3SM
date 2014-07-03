  ! Read the average temperature profile from the initial condition file.
  !
  ! NOTE: This needs to be in its own file to avoid circular references.
  subroutine carma_getT(T)
    use shr_kind_mod, only: r8 => shr_kind_r8
    use cam_initfiles,only: initial_file_get_id
    use pio,          only: file_desc_t
    use ncdio_atm,    only: infld
    use pmgrid,       only: plat, plev, plevp, plon
    use ppgrid,       only: pcols, pver, pverp
    use abortutils,   only: endrun

    real(r8), intent(out)   :: T(pver)      ! midpoint temperature (Pa)

    integer                    :: iz           ! vertical index
    type(file_desc_t), pointer :: ncid_ini
    logical                    :: found
    real(r8), pointer          :: init_t(:,:,:)

    ! For an initial run, if the file is missing, then create one using the average
    ! temperature from the initial condition file.
    ncid_ini  => initial_file_get_id()
    nullify(init_t)

    allocate(init_t(plon,pver,plat))
    call infld('T', ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, pver, 1, &
         plat, init_t, found, grid_map='GLOBAL', array_order_in='xyz')
            
    if (.not. found) then
      call endrun('carma_init::infld failed to find field T.')
    end if

    ! Just do a simple average. Could get gw and do a weighted average.
    do iz = 1, pver
      T(iz) = sum(init_t(:, iz, :)) / plat / plon
    end do

    deallocate(init_t)
    
    return
  end

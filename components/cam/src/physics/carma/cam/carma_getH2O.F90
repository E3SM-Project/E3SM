  ! Read the average water vapor profile from the initial condition file.
  !
  ! NOTE: This needs to be in its own file to avoid circular references.
  subroutine carma_getH2O(h2o)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use cam_initfiles,  only: initial_file_get_id
    use pio,            only: file_desc_t
    use cam_pio_utils,  only: cam_pio_get_var
    use pmgrid,         only: plat, plev, plevp, plon
    use ppgrid,         only: pcols, pver, pverp
    use cam_abortutils, only: endrun

    real(r8), intent(out)   :: h2o(pver)      ! midpoint h2o mmr (kg/kg)

    integer                    :: iz           ! vertical index
    type(file_desc_t), pointer :: ncid_ini
    logical                    :: found
    real(r8), pointer          :: init_h2o(:,:,:)

    ! For an initial run, if the file is missing, then create one using the
    ! average concentration from the initial condition file.
    ncid_ini  => initial_file_get_id()
    nullify(init_h2o)

    allocate(init_h2o(plon,pver,plat))
    call cam_pio_get_var('Q', ncid_ini, init_h2o, found=found)
            
    if (.not. found) then
      call endrun('carma_init::cam_pio_get_var failed to find field Q.')
    end if

    ! Just do a simple average. Could get gw and do a weighted average.
    do iz = 1, pver
      h2o(iz) = sum(init_h2o(:, iz, :)) / plat / plon
    end do

    deallocate(init_h2o)
    
    return
  end

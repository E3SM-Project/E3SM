! Test-only C bindings to individual routines of the Fortran COSP library, used to check
! the C++/Kokkos port BFB against the Fortran. Arrays are column-major, so the C++ side
! passes LayoutLeft data.
module cosp_f90_shims
  use iso_c_binding
  use cosp_kinds,                only: wp
  use mod_rng,                   only: rng_state, init_rng, get_rng
  use cosp_optics,               only: cosp_simulator_optics, modis_optics_partition, modis_optics
  use mod_modis_sim,             only: get_g_nir, get_ssa_nir
  use mod_cosp_modis_interface,  only: cosp_modis_init
  use mod_cosp_stats,            only: hist2D
  use mod_cosp_config,           only: R_UNDEF, tau_binBounds, pres_binBounds, hgt_binBounds, &
                                       ntau, npres, nhgt, numMISRHgtBins, numMISRTauBins
  use mod_misr_simulator,        only: misr_subcolumn, misr_column, misr_tauchk => tauchk
  use mod_icarus,                only: icarus_subcolumn, icarus_column, isccp_tauchk => tauchk
  use mod_cosp_isccp_interface,  only: cosp_isccp_init

  implicit none

contains

  ! Draw ndraw numbers from each of n streams, seeded with seeds(:)
  subroutine cosp_rng_f(n, ndraw, seeds, rnd) bind(C, name='cosp_rng_f')
    integer(c_int), value, intent(in) :: n, ndraw
    integer(c_int), intent(in)        :: seeds(n)
    real(c_double), intent(out)       :: rnd(n,ndraw)
    type(rng_state) :: rngs(n)
    integer :: seed(n), i
    seed = seeds
    call init_rng(rngs, seed)
    do i = 1, ndraw
       rnd(:,i) = get_rng(rngs)
    enddo
  end subroutine cosp_rng_f

  subroutine cosp_simulator_optics_f(dim1, dim2, dim3, flag, var1, var2, varout) &
       bind(C, name='cosp_simulator_optics_f')
    integer(c_int), value, intent(in) :: dim1, dim2, dim3
    real(c_double), intent(in)        :: flag(dim1,dim2,dim3), var1(dim1,dim3), var2(dim1,dim3)
    real(c_double), intent(out)       :: varout(dim1,dim2,dim3)
    call cosp_simulator_optics(dim1, dim2, dim3, flag, var1, var2, varout)
  end subroutine cosp_simulator_optics_f

  subroutine cosp_modis_optics_partition_f(npoints, ncolumns, nlev, cloud_water, cloud_ice, &
       water_size, ice_size, tau, tau_liq, tau_ice) bind(C, name='cosp_modis_optics_partition_f')
    integer(c_int), value, intent(in) :: npoints, ncolumns, nlev
    real(c_double), intent(in),  dimension(npoints,ncolumns,nlev) :: cloud_water, cloud_ice, &
         water_size, ice_size, tau
    real(c_double), intent(out), dimension(npoints,ncolumns,nlev) :: tau_liq, tau_ice
    ! Same argument order (and dimension swap) as the call in cosp_c2f.F90
    call modis_optics_partition(npoints, nlev, ncolumns, cloud_water, cloud_ice, water_size, &
         ice_size, tau, tau_liq, tau_ice)
  end subroutine cosp_modis_optics_partition_f

  subroutine cosp_modis_optics_f(npoints, ncolumns, nlev, tau_liq, size_liq, tau_ice, size_ice, &
       frac_liq, g, w0) bind(C, name='cosp_modis_optics_f')
    integer(c_int), value, intent(in) :: npoints, ncolumns, nlev
    real(c_double), intent(in),  dimension(npoints,ncolumns,nlev) :: tau_liq, size_liq, &
         tau_ice, size_ice
    real(c_double), intent(out), dimension(npoints,ncolumns,nlev) :: frac_liq, g, w0
    call cosp_modis_init()
    call modis_optics(npoints, nlev, ncolumns, tau_liq, size_liq, tau_ice, size_ice, &
         frac_liq, g, w0)
  end subroutine cosp_modis_optics_f

  subroutine cosp_nir_fits_f(n, phase, re, g, w0) bind(C, name='cosp_nir_fits_f')
    integer(c_int), value, intent(in) :: n, phase
    real(c_double), intent(in)        :: re(n)
    real(c_double), intent(out)       :: g(n), w0(n)
    call cosp_modis_init()
    g  = get_g_nir(phase, re)
    w0 = get_ssa_nir(phase, re)
  end subroutine cosp_nir_fits_f

  subroutine cosp_hist2d_f(npts, var1, var2, nbin1, bin1, nbin2, bin2, hist) &
       bind(C, name='cosp_hist2d_f')
    integer(c_int), value, intent(in) :: npts, nbin1, nbin2
    real(c_double), intent(in)        :: var1(npts), var2(npts), bin1(nbin1+1), bin2(nbin2+1)
    real(c_double), intent(out)       :: hist(nbin1,nbin2)
    call hist2D(var1, var2, npts, bin1, nbin1, bin2, nbin2, hist)
  end subroutine cosp_hist2d_f

  ! Constants and bin tables, to check the C++ copies of them
  subroutine cosp_constants_f(r_undef_out, tau_bnds, pres_bnds, hgt_bnds, misr_tauchk_out) &
       bind(C, name='cosp_constants_f')
    real(c_double), intent(out) :: r_undef_out, misr_tauchk_out
    real(c_double), intent(out) :: tau_bnds(ntau+1), pres_bnds(npres+1), hgt_bnds(nhgt+1)
    r_undef_out     = R_UNDEF
    tau_bnds        = tau_binBounds
    pres_bnds       = pres_binBounds
    hgt_bnds        = hgt_binBounds
    misr_tauchk_out = misr_tauchk
  end subroutine cosp_constants_f

  subroutine cosp_misr_subcolumn_f(npoints, ncol, nlev, dtau, zfull, at, sunlit, tau_out, ztop) &
       bind(C, name='cosp_misr_subcolumn_f')
    integer(c_int), value, intent(in) :: npoints, ncol, nlev
    real(c_double), intent(in)        :: dtau(npoints,ncol,nlev), zfull(npoints,nlev), at(npoints,nlev)
    integer(c_int), intent(in)        :: sunlit(npoints)
    real(c_double), intent(out)       :: tau_out(npoints,ncol), ztop(npoints,ncol)
    real(wp) :: dist_model_layertops(npoints,numMISRHgtBins)
    integer  :: sunlit_i(npoints)
    sunlit_i = sunlit
    call misr_subcolumn(npoints, ncol, nlev, dtau, zfull, at, sunlit_i, tau_out, &
         dist_model_layertops, ztop)
  end subroutine cosp_misr_subcolumn_f

  subroutine cosp_misr_column_f(npoints, ncol, ztop, sunlit, tau, fq) bind(C, name='cosp_misr_column_f')
    integer(c_int), value, intent(in) :: npoints, ncol
    real(c_double), intent(in)        :: ztop(npoints,ncol), tau(npoints,ncol)
    integer(c_int), intent(in)        :: sunlit(npoints)
    real(c_double), intent(out)       :: fq(npoints,numMISRTauBins,numMISRHgtBins)
    real(wp) :: cldarea(npoints), mean_ztop(npoints)
    integer  :: sunlit_i(npoints)
    sunlit_i = sunlit
    call misr_column(npoints, ncol, ztop, sunlit_i, tau, cldarea, mean_ztop, fq)
  end subroutine cosp_misr_column_f

  subroutine cosp_isccp_constants_f(tauchk_out) bind(C, name='cosp_isccp_constants_f')
    real(c_double), intent(out) :: tauchk_out
    tauchk_out = isccp_tauchk
  end subroutine cosp_isccp_constants_f

  ! ICARUS subcolumn and column parts, with the settings EAMxx uses (cosp_c2f.F90)
  subroutine cosp_icarus_f(npoints, ncol, nlev, sunlit, dtau, dem, skt, emsfc_lw, qv, at, &
       pfull, phalf, boxtau, boxptop, boxttop, fq, totalcldarea) bind(C, name='cosp_icarus_f')
    integer(c_int), value, intent(in) :: npoints, ncol, nlev
    integer(c_int), intent(in)        :: sunlit(npoints)
    real(c_double), value, intent(in) :: emsfc_lw
    real(c_double), intent(in)        :: dtau(npoints,ncol,nlev), dem(npoints,ncol,nlev), skt(npoints)
    real(c_double), intent(in), dimension(npoints,nlev) :: qv, at, pfull
    real(c_double), intent(in)        :: phalf(npoints,nlev+1)
    real(c_double), intent(out)       :: boxtau(npoints,ncol), boxptop(npoints,ncol), boxttop(npoints,ncol)
    real(c_double), intent(out)       :: fq(npoints,7,7), totalcldarea(npoints)
    real(wp) :: frac_out(npoints,ncol,nlev), meantbclr(npoints)
    real(wp), dimension(npoints) :: meanalbedocld, meanptop, meantaucld, meantb
    integer  :: levmatch(npoints,ncol), sunlit_i(npoints)
    sunlit_i = sunlit
    frac_out = 1._wp
    call cosp_isccp_init(1, 2)
    call icarus_subcolumn(npoints, ncol, nlev, sunlit_i, dtau, dem, skt, emsfc_lw, qv, at, &
         pfull, phalf, frac_out, levmatch, boxtau, boxptop, boxttop, meantbclr)
    call icarus_column(npoints, ncol, boxtau, boxptop/100._wp, sunlit_i, boxttop, fq, &
         meanalbedocld, meanptop, meantaucld, totalcldarea, meantb)
  end subroutine cosp_icarus_f

end module cosp_f90_shims

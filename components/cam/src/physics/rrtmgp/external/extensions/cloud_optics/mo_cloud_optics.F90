! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Provides cloud optical properties as a function of effective radius for the RRTMGP bands
!   Based on Mie calculations for liquid
!     and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
!   Can use either look-up tables or Pade approximates according to which data has been loaded
!   Mike Iacono (AER) is the original author
!
! The class can be used as-is but is also intended as an example of how to extend the RTE framework
! -------------------------------------------------------------------------------------------------

module mo_cloud_optics
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private
  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_cloud_optics
    private
    !
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    !
    integer            :: icergh = 0  ! (1 = none, 2 = medium, 3 = high)
    !
    ! Lookup table information
    !
    ! Upper and lower limits of the tables
    real(wp) :: radliq_lwr = 0._wp, radliq_upr = 0._wp
    real(wp) :: radice_lwr = 0._wp, radice_upr = 0._wp
    ! How many steps in the table? (for convenience)
    integer  :: liq_nsteps = 0,        ice_nsteps = 0
    ! How big is each step in the table?
    real(wp) :: liq_step_size = 0._wp, ice_step_size = 0._wp
    !
    ! The tables themselves.
    !
    real(wp), dimension(:,:    ), allocatable :: lut_extliq, lut_ssaliq, lut_asyliq ! (nsize_liq, nbnd)
    real(wp), dimension(:,:,:  ), allocatable :: lut_extice, lut_ssaice, lut_asyice ! (nsize_ice, nbnd, nrghice)

    !
    ! Pade approximant coefficients
    !
    real(wp), dimension(:,:,:  ), allocatable :: pade_extliq                 ! (nbnd, nsizereg, ncoeff_ext)
    real(wp), dimension(:,:,:  ), allocatable :: pade_ssaliq,  pade_asyliq   ! (nbnd, nsizereg, ncoeff_ssa_g)
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice                 ! (nbnd, nsizereg, ncoeff_ext, nrghice)
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice, pade_asyice    ! (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
    ! Particle size regimes for Pade formulations
    real(wp), dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq  ! (nbound)
    real(wp), dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice  ! (nbound)
    ! -----
  contains
    generic,   public :: load  => load_lut, load_pade
    procedure, public :: finalize
    procedure, public :: cloud_optics
    procedure, public :: get_min_radius_liq
    procedure, public :: get_min_radius_ice
    procedure, public :: get_max_radius_liq
    procedure, public :: get_max_radius_ice
    procedure, public :: get_num_ice_roughness_types
    procedure, public :: set_ice_roughness
    ! Internal procedures
    procedure, private :: load_lut
    procedure, private :: load_pade
  end type ty_cloud_optics

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routines to load data needed for cloud optics calculations. Two routines: one to load
  !    lookup-tables and one for coefficients for Pade approximates
  !
  ! ------------------------------------------------------------------------------
  function load_lut(this, band_lims_wvn, &
                    radliq_lwr, radliq_upr, radliq_fac, &
                    radice_lwr, radice_upr, radice_fac, &
                    lut_extliq, lut_ssaliq, lut_asyliq, &
                    lut_extice, lut_ssaice, lut_asyice) result(error_msg)
    class(ty_cloud_optics),     intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! Spectral discretization
    ! Lookup table interpolation constants
    ! Lower and upper bounds of the tables; also the constant for calculating interpolation indices for liquid
    real(wp),                   intent(in   ) :: radliq_lwr, radliq_upr, radliq_fac
    real(wp),                   intent(in   ) :: radice_lwr, radice_upr, radice_fac
    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for liquid and ice respectively
    real(wp), dimension(:,:),   intent(in)    :: lut_extliq, lut_ssaliq, lut_asyliq
    real(wp), dimension(:,:,:), intent(in)    :: lut_extice, lut_ssaice, lut_asyice
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer               :: nbnd, nrghice, nsize_liq, nsize_ice

    error_msg = this%init(band_lims_wvn, name="RRTMGP cloud optics")
    !
    ! LUT coefficient dimensions
    !
    nsize_liq = size(lut_extliq,dim=1)
    nsize_ice = size(lut_extice,dim=1)
    nbnd     = size(lut_extliq,dim=2)
    nrghice   = size(lut_extice,dim=3)
    !
    ! Error checking
    !   Can we check for consistency between table bounds and _fac?
    !
    if(nbnd /= this%get_nband()) &
      error_msg = "cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
    if(size(lut_extice, 2) /= nbnd) &
      error_msg = "cloud_optics%init(): array lut_extice has the wrong number of bands"
    if(any([size(lut_ssaliq, 1), size(lut_ssaliq, 2)] /= [nsize_liq, nbnd])) &
      error_msg = "cloud_optics%init(): array lut_ssaliq isn't consistently sized"
    if(any([size(lut_asyliq, 1), size(lut_asyliq, 2)] /= [nsize_liq, nbnd])) &
      error_msg = "cloud_optics%init(): array lut_asyliq isn't consistently sized"
    if(any([size(lut_ssaice, 1), size(lut_ssaice, 2), size(lut_ssaice, 3)] /= [nsize_ice, nbnd, nrghice])) &
      error_msg = "cloud_optics%init(): array lut_ssaice  isn't consistently sized"
    if(any([size(lut_asyice, 1), size(lut_asyice, 2), size(lut_asyice, 3)] /= [nsize_ice, nbnd, nrghice])) &
      error_msg = "cloud_optics%init(): array lut_asyice  isn't consistently sized"
    if(error_msg /= "") return

    this%liq_nsteps = nsize_liq
    this%ice_nsteps = nsize_ice
    this%liq_step_size = (radliq_upr - radliq_lwr)/real(nsize_liq-1,wp)
    this%ice_step_size = (radice_upr - radice_lwr)/real(nsize_ice-1,wp)
    ! Allocate LUT coefficients
    allocate(this%lut_extliq(nsize_liq, nbnd), &
             this%lut_ssaliq(nsize_liq, nbnd), &
             this%lut_asyliq(nsize_liq, nbnd), &
             this%lut_extice(nsize_ice, nbnd, nrghice), &
             this%lut_ssaice(nsize_ice, nbnd, nrghice), &
             this%lut_asyice(nsize_ice, nbnd, nrghice))

    ! Load LUT constants
    this%radliq_lwr = radliq_lwr
    this%radliq_upr = radliq_upr
    this%radice_lwr = radice_lwr
    this%radice_upr = radice_upr

    ! Load LUT coefficients
    this%lut_extliq = lut_extliq
    this%lut_ssaliq = lut_ssaliq
    this%lut_asyliq = lut_asyliq
    this%lut_extice = lut_extice
    this%lut_ssaice = lut_ssaice
    this%lut_asyice = lut_asyice
  end function load_lut
  ! ------------------------------------------------------------------------------
  !
  ! Cloud optics initialization function - Pade
  !
  ! ------------------------------------------------------------------------------
  function load_pade(this, band_lims_wvn, &
                     pade_extliq, pade_ssaliq, pade_asyliq, &
                     pade_extice, pade_ssaice, pade_asyice, &
                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice) &
                     result(error_msg)
    class(ty_cloud_optics),       intent(inout) :: this          ! cloud specification data
    real(wp), dimension(:,:),     intent(in   ) :: band_lims_wvn ! Spectral discretization
    !
    ! Pade coefficients: extinction, single-scattering albedo, and asymmetry factor for liquid and ice
    !
    real(wp), dimension(:,:,:),   intent(in)    :: pade_extliq, pade_ssaliq, pade_asyliq
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_extice, pade_ssaice, pade_asyice
    !
    ! Boundaries of size regimes. Liquid and ice are separate;
    !   extinction is fit to different numbers of size bins than single-scattering albedo and asymmetry factor
    !
    real(wp),  dimension(:),       intent(in)    :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq
    real(wp),  dimension(:),       intent(in)    :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice
    character(len=128)    :: error_msg

! ------- Local -------

    integer               :: nbnd, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

! ------- Definitions -------

    ! Pade coefficient dimensions
    nbnd        = size(pade_extliq,dim=1)
    nsizereg     = size(pade_extliq,dim=2)
    ncoeff_ext   = size(pade_extliq,dim=3)
    ncoeff_ssa_g = size(pade_ssaliq,dim=3)
    nrghice      = size(pade_extice,dim=4)
    nbound       = size(pade_sizreg_extliq)
    error_msg = this%init(band_lims_wvn, name="RRTMGP cloud optics")
    !
    ! Error checking
    !
    if(nbnd /= this%get_nband()) &
      error_msg = "cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
    if(any([size(pade_ssaliq, 1), size(pade_ssaliq, 2), size(pade_ssaliq, 3)] /= [nbnd, nsizereg, ncoeff_ssa_g])) &
      error_msg = "cloud_optics%init(): array pade_ssaliq isn't consistently sized"
    if(any([size(pade_asyliq, 1), size(pade_asyliq, 2), size(pade_asyliq, 3)] /= [nbnd, nsizereg, ncoeff_ssa_g])) &
      error_msg = "cloud_optics%init(): array pade_asyliq isn't consistently sized"
    if(any([size(pade_extice, 1), size(pade_extice, 2), size(pade_extice, 3), size(pade_extice, 4)] /= &
           [nbnd,                nsizereg,             ncoeff_ext,           nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_extice isn't consistently sized"
    if(any([size(pade_ssaice, 1), size(pade_ssaice, 2), size(pade_ssaice, 3), size(pade_ssaice, 4)] /= &
           [nbnd,                nsizereg,             ncoeff_ssa_g,         nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_ssaice isn't consistently sized"
    if(any([size(pade_asyice, 1), size(pade_asyice, 2), size(pade_asyice, 3), size(pade_asyice, 4)] /= &
           [nbnd,                nsizereg,             ncoeff_ssa_g,         nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_asyice isn't consistently sized"
    if(any([                          size(pade_sizreg_ssaliq), size(pade_sizreg_asyliq),               &
            size(pade_sizreg_extice), size(pade_sizreg_ssaice), size(pade_sizreg_asyice)] /= nbound))   &
      error_msg = "cloud_optics%init(): one or more Pade size regime arrays are inconsistently sized"
    if(nsizereg /= 3) &
        error_msg = "cloud_optics%init(): Expecting precisely three size regimes for Pade approximants"
    if(error_msg /= "") return
    !
    ! Consistency among size regimes
    !
    this%radliq_lwr = pade_sizreg_extliq(1)
    this%radliq_upr = pade_sizreg_extliq(nbound)
    this%radice_lwr = pade_sizreg_extice(1)
    this%radice_upr = pade_sizreg_extice(nbound)

    if(any([pade_sizreg_ssaliq(1), pade_sizreg_asyliq(1)] < this%radliq_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lowest values"
    if(any([pade_sizreg_ssaice(1), pade_sizreg_asyice(1)] < this%radice_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lower values"

    if(any([pade_sizreg_ssaliq(nbound), pade_sizreg_asyliq(nbound)] > this%radliq_upr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr"
    if(any([pade_sizreg_ssaice(nbound), pade_sizreg_asyice(nbound)] > this%radice_upr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr"
    if(error_msg /= "") return
    !
    ! Allocate Pade coefficients
    !
    allocate(this%pade_extliq(nbnd, nsizereg, ncoeff_ext),   &
             this%pade_ssaliq(nbnd, nsizereg, ncoeff_ssa_g), &
             this%pade_asyliq(nbnd, nsizereg, ncoeff_ssa_g), &
             this%pade_extice(nbnd, nsizereg, ncoeff_ext,   nrghice), &
             this%pade_ssaice(nbnd, nsizereg, ncoeff_ssa_g, nrghice), &
             this%pade_asyice(nbnd, nsizereg, ncoeff_ssa_g, nrghice))
    !
    ! Allocate Pade coefficient particle size regime boundaries
    !
    allocate(this%pade_sizreg_extliq(nbound), &
             this%pade_sizreg_ssaliq(nbound), &
             this%pade_sizreg_asyliq(nbound), &
             this%pade_sizreg_extice(nbound), &
             this%pade_sizreg_ssaice(nbound), &
             this%pade_sizreg_asyice(nbound))
    !
    ! Load data
    !
    this%pade_extliq = pade_extliq
    this%pade_ssaliq = pade_ssaliq
    this%pade_asyliq = pade_asyliq
    this%pade_extice = pade_extice
    this%pade_ssaice = pade_ssaice
    this%pade_asyice = pade_asyice
    this%pade_sizreg_extliq = pade_sizreg_extliq
    this%pade_sizreg_ssaliq = pade_sizreg_ssaliq
    this%pade_sizreg_asyliq = pade_sizreg_asyliq
    this%pade_sizreg_extice = pade_sizreg_extice
    this%pade_sizreg_ssaice = pade_sizreg_ssaice
    this%pade_sizreg_asyice = pade_sizreg_asyice
  end function load_pade
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Finalize
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize(this)
    class(ty_cloud_optics), intent(inout) :: this

    this%radliq_lwr = 0._wp
    this%radliq_upr = 0._wp
    this%radice_lwr = 0._wp
    this%radice_upr = 0._wp

    ! Lookup table cloud optics coefficients
    if(allocated(this%lut_extliq)) then
      deallocate(this%lut_extliq, this%lut_ssaliq, this%lut_asyliq, &
                 this%lut_extice, this%lut_ssaice, this%lut_asyice)
      this%liq_nsteps = 0
      this%ice_nsteps = 0
      this%liq_step_size = 0._wp
      this%ice_step_size = 0._wp
    end if

    ! Pade cloud optics coefficients
    if(allocated(this%pade_extliq)) then
      deallocate(this%pade_extliq, this%pade_ssaliq, this%pade_asyliq, &
                 this%pade_extice, this%pade_ssaice, this%pade_asyice, &
                 this%pade_sizreg_extliq, this%pade_sizreg_ssaliq, this%pade_sizreg_asyliq, &
                 this%pade_sizreg_extice, this%pade_sizreg_ssaice, this%pade_sizreg_asyice)
    end if
  end subroutine finalize
  ! ------------------------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud physical properties
  !
  ! ------------------------------------------------------------------------------
  !
  ! Compute single-scattering properties
  !
  function cloud_optics(this, &
                        ncol, nlay, nbnd, nrghice,                &
                        liqmsk, icemsk, clwp, ciwp, reliq, reice, &
                        optical_props) result(error_msg)
    class(ty_cloud_optics), &
              intent(in   ) :: this
    integer,  intent(in   ) :: ncol, nlay, nbnd
    integer,  intent(in   ) :: nrghice              ! number of ice roughness categories
    logical,  intent(in   ) :: liqmsk(ncol,nlay), & ! Cloud mask for liquid and ice clouds respectively
                               icemsk(ncol,nlay)
    real(wp), intent(in   ) :: ciwp  (ncol,nlay), &     ! cloud ice water path
                               clwp  (ncol,nlay), &     ! cloud liquid water path
                               reice (ncol,nlay), &      ! cloud ice particle effective size (microns)
                               reliq (ncol,nlay)      ! cloud liquid particle effective radius (microns)
    class(ty_optical_props_arry), &
              intent(inout) :: optical_props
                                               ! Dimensions: (ncol,nlay,nbnd)

    character(len=128)      :: error_msg
    ! ------- Local -------
    type(ty_optical_props_2str) :: clouds_liq, clouds_ice
    integer  :: nsizereg, ibnd, imom
    ! ----------------------------------------
    !
    ! Error checking
    !
    ! ----------------------------------------
    error_msg = ''
    if(.not.(allocated(this%lut_extliq) .or. allocated(this%pade_extliq))) then
      error_msg = 'cloud optics: no data has been initialized'
      return
    end if

    if(.not. this%bands_are_equal(optical_props)) &
      error_msg = "cloud optics: optical properties don't have the same band structure"

    if(optical_props%get_nband() /= optical_props%get_ngpt() ) &
      error_msg = "cloud optics: optical properties must be requested by band not g-points"

    if (this%icergh < 1 .or. this%icergh > this%get_num_ice_roughness_types()) &
       error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'

    if(any(liqmsk .and. (reliq < this%radliq_lwr .or. reliq > this%radliq_upr))) &
      error_msg = 'cloud optics: liquid effective radius is out of bounds'

    if(any(icemsk .and. (reice < this%radice_lwr .or. reice > this%radice_upr))) &
      error_msg = 'cloud optics: ice effective radius is out of bounds'

    if(any((liqmsk .and.  clwp < 0._wp) .or. (icemsk .and.  ciwp < 0._wp))) &
      error_msg = 'cloud optics: negative clwp or ciwp where clouds are supposed to be'

    if(error_msg /= "") return

    ! ----------------------------------------
    !
    ! Compute cloud optical properties. Use lookup tables if available, Pade approximants if not.
    !
    error_msg = clouds_liq%alloc_2str(ncol,nlay, this)
    if(error_msg /= "") return
    error_msg = clouds_ice%alloc_2str(ncol,nlay, this)
    if (allocated(this%lut_extliq)) then
      !
      ! Liquid
      !
      clouds_liq%tau = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this%liq_nsteps,this%liq_step_size,this%radliq_lwr, &
                                  this%lut_extliq)
      clouds_liq%ssa = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this%liq_nsteps,this%liq_step_size,this%radliq_lwr, &
                                  this%lut_ssaliq)
      clouds_liq%g   = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this%liq_nsteps,this%liq_step_size,this%radliq_lwr, &
                                  this%lut_asyliq  )
      do ibnd = 1,nbnd
        where(liqmsk) clouds_liq%tau(1:ncol,1:nlay,ibnd) = clouds_liq%tau(1:ncol,1:nlay,ibnd) * clwp(1:ncol,1:nlay)
      end do
      !
      ! Ice
      !
      clouds_ice%tau = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this%ice_nsteps,this%ice_step_size,this%radice_lwr, &
                                  this%lut_extice(:,:,this%icergh))
      clouds_ice%ssa = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this%ice_nsteps,this%ice_step_size,this%radice_lwr, &
                                  this%lut_ssaice (:,:,this%icergh))
      clouds_ice%g   = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this%ice_nsteps,this%ice_step_size,this%radice_lwr, &
                                  this%lut_asyice  (:,:,this%icergh))
      do ibnd = 1, nbnd
        where(icemsk) clouds_ice%tau(1:ncol,1:nlay,ibnd) = clouds_ice%tau(1:ncol,1:nlay,ibnd) * ciwp(1:ncol,1:nlay)
      end do
    else
      !
      ! Cloud optical properties from Pade coefficient method
      !
      ! This assumes that all the Pade treaments have the same number of size regimes
      nsizereg = size(this%pade_extliq,2)
      !
      ! Liquid
      !
      clouds_liq%tau = compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this%pade_sizreg_extliq,2,3,this%pade_extliq)
      clouds_liq%ssa = 1._wp - max(0._wp,                                             &
                       compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this%pade_sizreg_ssaliq,2,2,this%pade_ssaliq) &
                                  )
      clouds_liq%g   = compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this%pade_sizreg_asyliq,2,2,this%pade_asyliq)
      do ibnd = 1, nbnd
        where(liqmsk) clouds_liq%tau(1:ncol,1:nlay,ibnd) = clouds_liq%tau(1:ncol,1:nlay,ibnd) * clwp(1:ncol,1:nlay)
      end do
      !
      ! Ice
      !
      clouds_ice%tau = compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this%pade_sizreg_extice,2,3, &
                                         this%pade_extice(:,:,:,this%icergh))
      clouds_ice%ssa = 1._wp - max(0._wp,                                                                  &
                       compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this%pade_sizreg_ssaice,2,2, &
                                         this%pade_ssaice(:,:,:,this%icergh)))
      clouds_ice%g   = compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this%pade_sizreg_asyice,2,2, &
                                         this%pade_asyice(:,:,:,this%icergh))
      do ibnd = 1, nbnd
        where(icemsk) clouds_ice%tau(1:ncol,1:nlay,ibnd) = clouds_ice%tau(1:ncol,1:nlay,ibnd) * ciwp(1:ncol,1:nlay)
      end do
    endif

    !
    ! Combine liquid and ice contributions into total cloud optical properties
    !
   error_msg = clouds_ice%increment(clouds_liq)
   !
   ! Copy total cloud properties onto outputs
   !
   select type(optical_props)
   type is (ty_optical_props_1scl)
     optical_props%tau(1:ncol,1:nlay,1:nbnd) = clouds_liq%tau(1:ncol,1:nlay,1:nbnd) * &
                                      (1._wp - clouds_liq%ssa(1:ncol,1:nlay,1:nbnd))
   type is (ty_optical_props_2str)
     optical_props%tau(1:ncol,1:nlay,1:nbnd) = clouds_liq%tau(1:ncol,1:nlay,1:nbnd)
     optical_props%ssa(1:ncol,1:nlay,1:nbnd) = clouds_liq%ssa(1:ncol,1:nlay,1:nbnd)
     optical_props%g  (1:ncol,1:nlay,1:nbnd) = clouds_liq%g  (1:ncol,1:nlay,1:nbnd)
   type is (ty_optical_props_nstr)
     optical_props%tau(  1:ncol,1:nlay,1:nbnd) = clouds_liq%tau(1:ncol,1:nlay,1:nbnd)
     optical_props%ssa(  1:ncol,1:nlay,1:nbnd) = clouds_liq%ssa(1:ncol,1:nlay,1:nbnd)
     optical_props%p  (1,1:ncol,1:nlay,1:nbnd) = clouds_liq%g  (1:ncol,1:nlay,1:nbnd)
     do imom = 2, optical_props%get_nmom()
       optical_props%p(imom,1:ncol,1:nlay,1:nbnd) = clouds_liq%g   (       1:ncol,1:nlay,1:nbnd) * &
                                                    optical_props%p(imom-1,1:ncol,1:nlay,1:nbnd)
     end do
   end select

  end function cloud_optics
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function set_ice_roughness(this, icergh) result(error_msg)
    class(ty_cloud_optics), intent(inout) :: this
    integer,                intent(in   ) :: icergh
    character(len=128)                    :: error_msg

    error_msg = ""
    if(icergh < 1) &
      error_msg = "cloud_optics%set_ice_roughness(): must be > 0"
    if(error_msg /= "") return

    this%icergh = icergh
  end function set_ice_roughness
  !-----------------------------------------------
  function get_num_ice_roughness_types(this) result(i)
    class(ty_cloud_optics), intent(in   ) :: this
    integer                               :: i

    i = 0
    if(allocated(this%pade_extice)) i = size(this%pade_extice, dim=4)
    if(allocated(this%lut_extice )) i = size(this%lut_extice,  dim=3)
  end function get_num_ice_roughness_types
  !-----------------------------------------------
  function get_min_radius_liq(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_lwr
  end function get_min_radius_liq
  !-----------------------------------------------
  function get_max_radius_liq(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_upr
  end function get_max_radius_liq
  !-----------------------------------------------
  function get_min_radius_ice(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radice_lwr
  end function get_min_radius_ice
  !-----------------------------------------------
  function get_max_radius_ice(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radice_upr
  end function get_max_radius_ice
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
  !   elements starting at "offset." The table's second dimension is band.
  ! Returns 0 where the mask is false.
  ! We could also try gather/scatter for efficiency
  !
  function compute_from_table(ncol, nlay, nbnd, mask, size, nsteps, step_size, offset, table)
    integer,                          intent(in) :: ncol, nlay, nbnd, nsteps
    logical,  dimension(ncol,  nlay), intent(in) :: mask
    real(wp), dimension(ncol,  nlay), intent(in) :: size
    real(wp),                         intent(in) :: step_size, offset
    real(wp), dimension(nsteps,nbnd), intent(in) :: table
    real(wp), dimension(ncol,nlay,nbnd)          :: compute_from_table
    ! ---------------------------
    integer  :: icol, ilay, ibnd
    integer  :: index
    real(wp) :: fint
    ! ---------------------------
    do ilay = 1,nlay
      do icol = 1, ncol
        if(mask(icol,ilay)) then
          index = min(floor((size(icol,ilay) - offset)/step_size)+1, nsteps-1)
          fint = (size(icol,ilay) - offset)/step_size - (index-1)
          do ibnd = 1, nbnd
            compute_from_table(icol,ilay,ibnd) = table(index,  ibnd) + &
                                         fint * (table(index+1,ibnd) - table(index,ibnd))
          end do
        else
          do ibnd = 1, nbnd
            compute_from_table(icol,ilay,ibnd) = 0._wp
          end do
        end if
      end do
    end do
  end function compute_from_table
  !---------------------------------------------------------------------------
  !
  ! Pade functions
  !
  !---------------------------------------------------------------------------
  function compute_from_pade(ncol, nlay, nbnd, mask, size, nsizes, size_bounds, m, n, pade_coeffs)
    integer,                        intent(in) :: ncol, nlay, nbnd, nsizes
    logical,  dimension(ncol,nlay), intent(in) :: mask
    real(wp), dimension(ncol,nlay), intent(in) :: size
    real(wp), dimension(nsizes+1),  intent(in) :: size_bounds
    integer,                        intent(in) :: m, n
    real(wp), dimension(nbnd,nsizes,0:m+n), &
                                    intent(in) :: pade_coeffs
    real(wp), dimension(ncol,nlay,nbnd)        :: compute_from_pade
    ! ---------------------------
    integer  :: icol, ilay, ibnd, irad

    do ilay = 1,nlay
      do icol = 1, ncol
        if(mask(icol,ilay)) then
          !
          ! Finds index into size regime table
          ! This works only if there are precisely three size regimes (four bounds) and it's
          !   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
          !
          irad = min(floor((size(icol,ilay) - size_bounds(2))/size_bounds(3))+2, 3)
          compute_from_pade(icol,ilay,1:nbnd) = &
               pade_eval(nbnd, nsizes, m, n, irad, size(icol,ilay), pade_coeffs)
        else
          do ibnd = 1, nbnd
            compute_from_pade(icol,ilay,ibnd) = 0._wp
          end do
        end if
      end do
    end do

  end function compute_from_pade
  !---------------------------------------------------------------------------
  !
  ! Evaluate Pade approximant of order [m/n]
  !
  function pade_eval(nbnd, nrads, m, n, irad, re, pade_coeffs)
    integer,                intent(in) :: nbnd, nrads, m, n, irad
    real(wp), dimension(nbnd, nrads, 0:m+n), &
                            intent(in) :: pade_coeffs
    real(wp),               intent(in) :: re
    real(wp), dimension(nbnd)         :: pade_eval

    integer :: iband
    real(wp) :: numer, denom
    integer  :: i

    do iband = 1, nbnd
      denom = pade_coeffs(iband,irad,n+m)
      do i = n-1+m, 1+m, -1
        denom = pade_coeffs(iband,irad,i)+re*denom
      end do
      denom =  1._wp                     +re*denom

      numer = pade_coeffs(iband,irad,m)
      do i = m-1, 1, -1
        numer = pade_coeffs(iband,irad,i)+re*numer
      end do
      numer = pade_coeffs(iband,irad,0)  +re*numer

      pade_eval(iband) = numer/denom
    end do
  end function pade_eval
end module mo_cloud_optics

module LakeBGCType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Hold lake biogeochemical state and flux variables 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use elm_varcon        , only : zisoi, spval
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use ColumnType        , only : col_pp
  use LandunitType      , only : lun_pp
  !
  implicit none
  save
  private
  ! gas identifier
  integer, parameter, public :: gn2lak = 1
  integer, parameter, public :: go2lak = 2
  integer, parameter, public :: gco2lak = 3
  integer, parameter, public :: gch4lak = 4
  ! phytoplankton group identifier
  integer, parameter, public :: small_phyto = 1
  integer, parameter, public :: large_phyto = 2
  ! sediment C pool identifier
  integer, parameter, public :: pasC = 1
  integer, parameter, public :: actC = 2
  ! lake type identifier
  integer, parameter, public :: yedoma_lake = 1
  integer, parameter, public :: thaw_lake = 2
  ! sub-cycle time step
  integer, parameter, public :: nsubstep = 10
  !
  ! !PUBLIC TYPES:
  type, public :: lakebgc_type

     ! Lake BGC input variables
     real(r8), pointer :: tp_col(:)                   ! col epilimnion average total phosphorus  (gP/m3)
     real(r8), pointer :: ph_col(:)                   ! col water average pH
     real(r8), pointer :: sdep_col(:)                 ! col sediment deposition rate (g/m2/yr)
     integer,  pointer :: ltype_col(:)                ! col lake type (regular lake = 0) 

     ! Lake BGC state variables
     real(r8), pointer :: conc_wat_col(:,:,:)         ! col water-column depth-resolved dissolved gas conc (mol/m3)
     real(r8), pointer :: conc_sed_col(:,:,:)         ! col sed-column depth-resolved dissloved gas conc (mol/m3)
     real(r8), pointer :: conc_bubl_col(:,:,:)        ! col depth-resolved bubble gas conc (mol/m3)
     real(r8), pointer :: conc_iceb_col(:,:)          ! col total bubble gas trapped in ice layers (mol/m2)
     real(r8), pointer :: biomas_phyto_col(:,:,:)     ! col phytoplankton biomass (gC/m3)
     real(r8), pointer :: chla_col(:,:)               ! col chlorophyll-a conc (g/m3)
     real(r8), pointer :: soilc_col(:,:,:)            ! col sediment C pools (gC/m3)

     ! Lake BGC flux variables
     real(r8), pointer :: ch4_sed_diff_col(:)         ! col CH4 diffusion at the water-sediment interface (mol/m2/s)
     real(r8), pointer :: ch4_surf_diff_col(:)        ! col CH4 diffusion at the air-water interface (mol/m2/s)
     real(r8), pointer :: ch4_sed_ebul_col(:)         ! col CH4 ebullition at the water-sediment interface (mol/m2/s)
     real(r8), pointer :: ch4_surf_ebul_col(:)        ! col CH4 ebullition at the air-water interface (mol/m2/s)
     real(r8), pointer :: ch4_surf_flux_tot_col(:)    ! col total CH4 flux at the air-water interface (kg C/m**2/s)
     real(r8), pointer :: ch4_prod_wat_col(:,:)       ! col water-column depth-resolved CH4 production (mol/m3/s)
     real(r8), pointer :: ch4_prod_sed_col(:,:)       ! col sed-column depth-resolved CH4 production (mol/m3/s)
     real(r8), pointer :: ch4_prod_tot_col(:)         ! col integrated CH4 production (gC/m2/s)
     real(r8), pointer :: ch4_oxid_tot_col(:)         ! col integrated CH4 oxidation (gC/m2/s)
     real(r8), pointer :: nem_col(:)                  ! net adjustment to atm. C flux from methane production (g C/m**2/s)
     real(r8), pointer :: ch4_oxid_wat_col(:,:)       ! col water-column depth-resolved CH4 oxidation (mol/m3/s)
     real(r8), pointer :: ch4_oxid_sed_col(:,:)       ! col sed-column depth-resolved CH4 oxidation (mol/m3/s)
     real(r8), pointer :: gpp_col(:)                  ! col depth-integrated gross primary production (gC/m2/s)
     real(r8), pointer :: npp_col(:)                  ! col depth-integrated net primary production (gC/m2/s)

  contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: Restart

  end type lakebgc_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    
    class(lakebgc_type)            :: this
    type(bounds_type) , intent(in) :: bounds

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds )

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use elm_varpar     , only : nlevlak, nlevgrnd, nlevsoi
    use elm_varpar     , only : ngaslak, nphytolak, nsoilclak
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(lakebgc_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate( this%tp_col              (begc:endc))                                  ; this%tp_col             (:)     = nan
    allocate( this%ph_col              (begc:endc))                                  ; this%ph_col             (:)     = nan
    allocate( this%sdep_col            (begc:endc))                                  ; this%sdep_col           (:)     = nan
    allocate( this%ltype_col           (begc:endc))                                  ; this%ltype_col          (:)     = 0

    allocate( this%conc_wat_col        (begc:endc,1:nlevlak,1:ngaslak))              ; this%conc_wat_col       (:,:,:) = nan
    allocate( this%conc_sed_col        (begc:endc,1:nlevgrnd,1:ngaslak))             ; this%conc_sed_col       (:,:,:) = nan
    allocate( this%conc_bubl_col       (begc:endc,1:nlevlak,1:ngaslak))              ; this%conc_bubl_col      (:,:,:) = nan 
    allocate( this%biomas_phyto_col    (begc:endc,1:nlevlak,1:nphytolak))            ; this%biomas_phyto_col   (:,:,:) = nan
    allocate( this%soilc_col           (begc:endc,1:nlevgrnd,1:nsoilclak))           ; this%soilc_col          (:,:,:) = nan
    allocate( this%chla_col            (begc:endc,1:nlevlak))                        ; this%chla_col           (:,:)   = nan
    allocate( this%conc_iceb_col       (begc:endc,1:ngaslak))                        ; this%conc_iceb_col      (:,:)   = nan
   
    allocate( this%ch4_sed_diff_col    (begc:endc))                                  ; this%ch4_sed_diff_col   (:)     = nan
    allocate( this%ch4_surf_diff_col   (begc:endc))                                  ; this%ch4_surf_diff_col  (:)     = nan
    allocate( this%ch4_sed_ebul_col    (begc:endc))                                  ; this%ch4_sed_ebul_col   (:)     = nan
    allocate( this%ch4_surf_ebul_col   (begc:endc))                                  ; this%ch4_surf_ebul_col  (:)     = nan
    allocate( this%ch4_surf_flux_tot_col (begc:endc))                                ; this%ch4_surf_flux_tot_col (:)  = nan
    allocate( this%gpp_col             (begc:endc))                                  ; this%gpp_col            (:)     = nan
    allocate( this%npp_col             (begc:endc))                                  ; this%npp_col            (:)     = nan
    allocate( this%ch4_prod_wat_col    (begc:endc,1:nlevlak))                        ; this%ch4_prod_wat_col   (:,:)   = nan
    allocate( this%ch4_oxid_wat_col    (begc:endc,1:nlevlak))                        ; this%ch4_oxid_wat_col   (:,:)   = nan
    allocate( this%ch4_prod_sed_col    (begc:endc,1:nlevgrnd))                       ; this%ch4_prod_sed_col   (:,:)   = nan
    allocate( this%ch4_oxid_sed_col    (begc:endc,1:nlevgrnd))                       ; this%ch4_oxid_sed_col   (:,:)   = nan
    allocate( this%ch4_prod_tot_col    (begc:endc))                                  ; this%ch4_prod_tot_col   (:)     = nan
    allocate( this%ch4_oxid_tot_col    (begc:endc))                                  ; this%ch4_oxid_tot_col   (:)     = nan
    allocate( this%nem_col             (begc:endc))                                  ; this%nem_col            (:)     = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use elm_varpar     , only : nlevlak, nlevgrnd
    use elm_varpar     , only : ngaslak, nphytolak, nsoilclak
    use histFileMod    , only : hist_addfld1d, hist_addfld2d
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(lakebgc_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc, k
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays 
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    this%ch4_sed_diff_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SED_DIFF_LAKE',  units='mol/m^2/s',  &
         avgflag='A', long_name='diffusive CH4 flux at the water-sediment interface', &
         ptr_col=this%ch4_sed_diff_col, default='inactive')

    this%ch4_surf_diff_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_DIFF_LAKE',  units='mol/m^2/s',  &
         avgflag='A', long_name='diffusive CH4 flux at the lake surface', &
         ptr_col=this%ch4_surf_diff_col, default='inactive')

    this%ch4_sed_ebul_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SED_EBUL_LAKE',  units='mol/m^2/s',  &
         avgflag='A', long_name='ebullition CH4 flux at the water-sediment interface', &
         ptr_col=this%ch4_sed_ebul_col, default='inactive')

    this%ch4_surf_ebul_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_EBUL_LAKE',  units='mol/m^2/s',  &
         avgflag='A', long_name='ebullition CH4 flux at the lake surface', &
         ptr_col=this%ch4_surf_ebul_col, default='inactive')

    this%gpp_col(begc:endc) = spval
    call hist_addfld1d (fname='GPP_LAKE',  units='gC/m^2/s',  &
         avgflag='A', long_name='lake gross primary production', &
         ptr_col=this%gpp_col, default='inactive')

    this%npp_col(begc:endc) = spval
    call hist_addfld1d (fname='NPP_LAKE',  units='gC/m^2/s',  &
         avgflag='A', long_name='lake net primary production', &
         ptr_col=this%npp_col, default='inactive')

    this%conc_wat_col(begc:endc,1:nlevlak,1:ngaslak) = spval
    data2dptr => this%conc_wat_col(:,:,gch4lak)
    call hist_addfld2d (fname='CONC_CH4_LAKE',  units='mol/m^3', type2d='levlak', & 
         avgflag='A', long_name='CH4 concentration in the lake water', &
         ptr_col=data2dptr, default='inactive')

    data2dptr => this%conc_wat_col(:,:,go2lak)
    call hist_addfld2d (fname='CONC_O2_LAKE',  units='mol/m^3', type2d='levlak', &
         avgflag='A', long_name='O2 concentration in the lake water', &
         ptr_col=data2dptr, default='inactive')

    this%conc_sed_col(begc:endc,1:nlevgrnd,1:ngaslak) = spval
    data2dptr => this%conc_sed_col(:,:,gch4lak)
    call hist_addfld2d (fname='CONC_CH4_SED',  units='mol/m^3', type2d='levgrnd', & 
         avgflag='A', long_name='CH4 concentration in the lake sediment', &
         ptr_col=data2dptr, default='inactive')

    this%conc_bubl_col(begc:endc,1:nlevlak,1:ngaslak) = spval
    data2dptr => this%conc_bubl_col(:,:,gch4lak)
    call hist_addfld2d (fname='BUBL_CH4_LAKE',  units='mol/m^3', type2d='levlak', &
         avgflag='A', long_name='Bubble CH4 concentration in the lake water', &
         ptr_col=data2dptr, default='inactive')

    this%conc_iceb_col(begc:endc,1:ngaslak) = spval
    data1dptr => this%conc_iceb_col(:,gch4lak)
    call hist_addfld1d (fname='BUBL_ICE_CH4',  units='mol/m^2',  & 
         avgflag='A', long_name='lake CH4 trapped in ice layers', &
         ptr_col=data1dptr, default='inactive')

    this%ch4_prod_wat_col(begc:endc,1:nlevlak) = spval
    call hist_addfld2d (fname='CH4_PROD_LAKE',  units='mol/m^3/s', type2d='levlak', &
         avgflag='A', long_name='CH4 production in the lake water', &
         ptr_col=this%ch4_prod_wat_col, default='inactive')

    this%ch4_prod_sed_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CH4_PROD_SED',  units='mol/m^3/s', type2d='levgrnd', &
         avgflag='A', long_name='CH4 production in the lake sediment', &
         ptr_col=this%ch4_prod_sed_col, default='inactive')

    this%ch4_oxid_wat_col(begc:endc,1:nlevlak) = spval
    call hist_addfld2d (fname='CH4_OXID_LAKE',  units='mol/m^3/s', type2d='levlak', &
         avgflag='A', long_name='CH4 oxidation in the lake water', &
         ptr_col=this%ch4_oxid_wat_col, default='inactive')

    this%ch4_oxid_sed_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CH4_OXID_SED',  units='mol/m^3/s', type2d='levgrnd', &
         avgflag='A', long_name='CH4 oxidation in the lake sediment', &
         ptr_col=this%ch4_oxid_sed_col, default='inactive')

    this%soilc_col(begc:endc,1:nlevgrnd,1:nsoilclak) = spval
    do k = 1, nsoilclak
       data2dptr => this%soilc_col(:,:,k)
       write(fieldname, "(A,I0,A)") 'LAKE_SED', k, 'C'
       write(longname, "(A,I0,A)") 'lake sediment ', k, ' C'
       call hist_addfld2d (fname=fieldname,  units='gC/m^3', type2d='levgrnd', &
            avgflag='A', long_name=longname, ptr_col=data2dptr, default='inactive')
    end do 

    this%biomas_phyto_col(begc:endc,1:nlevlak,1:nphytolak) = spval
    do k = 1, nphytolak
       data2dptr => this%biomas_phyto_col(:,:,k)
       write(fieldname, "(A,I0)") 'LAKE_PHYTO', k
       write(longname, "(A,I0)") 'biomass of lake phytoplankton group', k
       call hist_addfld2d (fname=fieldname,  units='gC/m^3', type2d='levlak', &
            avgflag='A', long_name=longname, ptr_col=data2dptr, default='inactive')
    end do

    this%chla_col(begc:endc,1:nlevlak) = spval
    call hist_addfld2d (fname='CHLA_LAKE',  units='g/m^3', type2d='levlak', &
         avgflag='A', long_name='chlorophyll-a concentration in the lake water', &
         ptr_col=this%chla_col, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use elm_varcon     , only : grlnd
    use elm_varctl     , only : fsurdat
    use elm_varpar     , only : nlevlak, nlevgrnd, nlevsoi, nlevsoifl
    use elm_varpar     , only : ngaslak, nphytolak, nsoilclak
    use elm_varpar     , only : more_vertlayers
    use landunit_varcon, only : istdlak
    use fileutils      , only : getfil
    use ncdio_pio      , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    use ncdio_pio      , only : ncd_inqdlen
    !
    ! !ARGUMENTS:
    class(lakebgc_type) :: this
    type(bounds_type) , intent(in) :: bounds
    ! !CONSTANTS
    real(r8), parameter :: carb0 = 8.9386e3_r8  ! gC/m2
    !
    ! !LOCAL VARIABLES:
    integer            :: c, g, l
    integer            :: j, lev
    integer            :: dimid                         ! dimension id
    logical            :: readvar
    type(file_desc_t)  :: ncid
    character(len=256) :: locfn
    real(r8) ,pointer  :: zsoifl     (:)      ! Output: [real(r8) (:)]  original soil midpoint 
    real(r8) ,pointer  :: zisoifl    (:)      ! Output: [real(r8) (:)]  original soil interface depth 
    real(r8) ,pointer  :: dzsoifl    (:)      ! Output: [real(r8) (:)]  original soil thickness 
    real(r8) ,pointer  :: tp2d       (:)      ! read in - TP 
    real(r8) ,pointer  :: ph2d       (:)      ! read in - pH
    real(r8) ,pointer  :: sdep2d     (:)      ! read in - sediment deposition
    real(r8) ,pointer  :: type2d     (:)      ! read in - lake type 
    real(r8) :: carbon
    !-----------------------------------------------------------------------

    allocate(tp2d(bounds%begg:bounds%endg))
    allocate(ph2d(bounds%begg:bounds%endg))
    allocate(sdep2d(bounds%begg:bounds%endg))
    allocate(type2d(bounds%begg:bounds%endg))

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid,dimid,nlevsoifl,name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call endrun(msg=' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if

    call ncd_io(ncid=ncid, varname='LAKE_TP', flag='read', data=tp2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) tp2d(:) = 0._r8
    call ncd_io(ncid=ncid, varname='LAKE_PH', flag='read', data=ph2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) ph2d(:) = 7._r8
    call ncd_io(ncid=ncid, varname='LAKE_SDEP', flag='read', data=sdep2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) sdep2d(:) = 0._r8
    call ncd_io(ncid=ncid, varname='LAKE_TYPE', flag='read', data=type2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) type2d(:) = 0
 
    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! get original soil depths to be used in interpolation of sediment C 
    ! --------------------------------------------------------------------

    allocate(zsoifl(1:nlevsoifl), zisoifl(1:nlevsoifl+1), dzsoifl(1:nlevsoifl))
    do j = 1, nlevsoifl
       zsoifl(j) = 0.025*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
    enddo

    dzsoifl(1) = 0.5_r8*(zsoifl(1)+zsoifl(2))              !thickness b/n two interfaces
    do j = 2,nlevsoifl-1
       dzsoifl(j)= 0.5_r8*(zsoifl(j+1)-zsoifl(j-1))
    enddo
    dzsoifl(nlevsoifl) = zsoifl(nlevsoifl)-zsoifl(nlevsoifl-1)

    zisoifl(1) = 0._r8
    do j = 2, nlevsoifl
       zisoifl(j) = 0.5_r8*(zsoifl(j-1)+zsoifl(j))         !interface depths
    enddo
    zisoifl(nlevsoifl+1) = zsoifl(nlevsoifl) + 0.5_r8*dzsoifl(nlevsoifl)

    carbon = carb0 / zisoifl(nlevsoi+1) / 0.6321_r8

    do c = bounds%begc, bounds%endc

       this%ch4_sed_diff_col(c)                       = 0._r8
       this%ch4_surf_diff_col(c)                      = 0._r8
       this%ch4_sed_ebul_col(c)                       = 0._r8
       this%ch4_surf_ebul_col(c)                      = 0._r8
       this%ch4_surf_flux_tot_col(c)                  = 0._r8
       this%gpp_col(c)                                = 0._r8
       this%npp_col(c)                                = 0._r8
       this%ch4_prod_wat_col(c,1:nlevlak)             = 0._r8
       this%ch4_oxid_wat_col(c,1:nlevlak)             = 0._r8
       this%ch4_prod_sed_col(c,1:nlevgrnd)            = 0._r8
       this%ch4_oxid_sed_col(c,1:nlevgrnd)            = 0._r8
       this%ch4_prod_tot_col(c)                       = 0._r8
       this%ch4_oxid_tot_col(c)                       = 0._r8
       this%nem_col(c)                                = 0._r8
       this%conc_wat_col(c,1:nlevlak,1:ngaslak)       = 0._r8
       this%conc_sed_col(c,1:nlevgrnd,1:ngaslak)      = 0._r8
       this%conc_bubl_col(c,1:nlevlak,1:ngaslak)      = 0._r8
       this%conc_iceb_col(c,1:ngaslak)                = 0._r8
       this%biomas_phyto_col(c,1:nlevlak,1:nphytolak) = 0._r8
       this%chla_col(c,1:nlevlak)                     = 0._r8
       this%soilc_col(c,1:nlevgrnd,1:nsoilclak)       = 0._r8

       g = col_pp%gridcell(c)
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) /= istdlak) then
          this%tp_col(c)                              = spval
          this%ph_col(c)                              = spval          
          this%sdep_col(c)                            = spval
          this%ltype_col(c)                           = 0
          this%ch4_sed_diff_col(c)                    = spval
          this%ch4_surf_diff_col(c)                   = spval
          this%ch4_sed_ebul_col(c)                    = spval
          this%ch4_surf_ebul_col(c)                   = spval
          this%ch4_surf_flux_tot_col(c)               = spval
          this%gpp_col(c)                             = spval
          this%npp_col(c)                             = spval
          this%ch4_prod_wat_col(c,:)                  = spval
          this%ch4_oxid_wat_col(c,:)                  = spval
          this%ch4_prod_sed_col(c,:)                  = spval
          this%ch4_oxid_sed_col(c,:)                  = spval
          this%ch4_prod_tot_col(c)                    = spval
          this%ch4_oxid_tot_col(c)                    = spval
          this%conc_wat_col(c,:,:)                    = spval
          this%conc_sed_col(c,:,:)                    = spval
          this%conc_bubl_col(c,:,:)                   = spval
          this%conc_iceb_col(c,:)                     = spval
          this%biomas_phyto_col(c,:,:)                = spval
          this%chla_col(c,:)                          = spval
          this%soilc_col(c,:,:)                       = spval
       else
          this%tp_col(c)                              = tp2d(g)
          this%ph_col(c)                              = ph2d(g)
          this%sdep_col(c)                            = sdep2d(g)
          this%ltype_col(c)                           = type2d(g)
          this%conc_wat_col(c,:,gn2lak)               = 0.347_r8
          this%conc_wat_col(c,:,go2lak)               = 0.425_r8 
          this%conc_wat_col(c,:,gco2lak)              = 0.032_r8
          this%conc_wat_col(c,:,gch4lak)              = 1.544e-6_r8
          this%biomas_phyto_col(c,:,small_phyto)      = 0.025_r8
          this%biomas_phyto_col(c,:,large_phyto)      = 0.025_r8
          this%chla_col(c,:)                          = 5.e-4_r8 
          this%conc_sed_col(c,:,gn2lak)               = 0.139_r8
          this%conc_sed_col(c,:,go2lak)               = 0._r8
          this%conc_sed_col(c,:,gco2lak)              = 0.013_r8
          this%conc_sed_col(c,:,gch4lak)              = 0.617e-6_r8
          
          do j = 1, nlevgrnd
             if (j<=nlevsoi) then
                this%soilc_col(c,j,pasC)              = 0.95e3 * carbon * exp(-zisoifl(j)) 
                this%soilc_col(c,j,actC)              = 0.05e3 * carbon * exp(-zisoifl(j))
             else
                this%soilc_col(c,j,pasC)              = 0._r8
                this%soilc_col(c,j,actC)              = 0._r8
             end if
          end do
       end if
    end do

    deallocate(tp2d, ph2d, sdep2d, type2d)

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t, ncd_double
    use restUtilMod 
    !
    ! !ARGUMENTS:
    class(lakebgc_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    logical           :: readvar
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------

    data2dptr => this%conc_wat_col(:,:,gn2lak) 
    call restartvar(ncid=ncid, flag=flag, varname='CONC_N2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='nitrogen lake concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr) 

    data2dptr => this%conc_wat_col(:,:,go2lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='oxygen lake concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_wat_col(:,:,gco2lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CO2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='carbon dioxide lake concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_wat_col(:,:,gch4lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='methane lake concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_sed_col(:,:,gn2lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_N2_SED', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='nitrogen lake sediment concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_sed_col(:,:,go2lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_SED', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen lake sediment concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_sed_col(:,:,gco2lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CO2_SED', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='carbon dioxide lake sediment concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_sed_col(:,:,gch4lak)
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_SED', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='methane lake sediment concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)

    data2dptr => this%conc_bubl_col(:,:,gn2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_N2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='bubble nitrogen concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)           

    data2dptr => this%conc_bubl_col(:,:,go2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_O2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='bubble oxygen concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_bubl_col(:,:,gco2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_CO2_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='bubble carbon dioxide concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  

    data2dptr => this%conc_bubl_col(:,:,gch4lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_CH4_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='bubble methane concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)  
 
    data1dptr => this%conc_iceb_col(:,gn2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_ICE_N2', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ice-trapped bubble nitrogen concentration', units='mol/m^2', &
         interpinic_flag='interp', readvar=readvar, data=data1dptr)

    data1dptr => this%conc_iceb_col(:,go2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_ICE_O2', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ice-trapped bubble oxygen concentration', units='mol/m^2', &
         interpinic_flag='interp', readvar=readvar, data=data1dptr)

    data1dptr => this%conc_iceb_col(:,gco2lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_ICE_CO2', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ice-trapped bubble carbon dioxide concentration', units='mol/m^2', &
         interpinic_flag='interp', readvar=readvar, data=data1dptr)

    data1dptr => this%conc_iceb_col(:,gch4lak)
    call restartvar(ncid=ncid, flag=flag, varname='BUBL_ICE_CH4', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ice-trapped bubble methane concentration', units='mol/m^2', &
         interpinic_flag='interp', readvar=readvar, data=data1dptr)

    call restartvar(ncid=ncid, flag=flag, varname='CHLA_LAKE', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='chlorophyll-a lake concentration', units='g/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%chla_col)

    data2dptr => this%biomas_phyto_col(:,:,small_phyto)
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_PHYTO1', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='small phytoplankton lake biomass', units='gC/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)    

    data2dptr => this%biomas_phyto_col(:,:,large_phyto)
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_PHYTO2', xtype=ncd_double, &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='large phytoplankton lake biomass', units='gC/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)

    data2dptr => this%soilc_col(:,:,pasC)
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_SED1C', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='passive sediment C', units='gC/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)

    data2dptr => this%soilc_col(:,:,actC)
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_SED2C', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='active sediment C', units='gC/m^3', &
         readvar=readvar, interpinic_flag='interp', data=data2dptr)

  end subroutine Restart

end module LakeBGCType

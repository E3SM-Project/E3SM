module EcophysConType

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use abortutils     , only : endrun
  !
  implicit none
  save
  public 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ecophysconInit
  !
  ! !PUBLIC TYPES:
  type, public :: ecophyscon_type
     integer , allocatable :: noveg         (:)   ! value for not vegetated
     integer , allocatable :: tree          (:)   ! tree or not?
     real(r8), allocatable :: smpso         (:)   ! soil water potential at full stomatal opening (mm)
     real(r8), allocatable :: smpsc         (:)   ! soil water potential at full stomatal closure (mm)
     real(r8), allocatable :: fnitr         (:)   ! foliage nitrogen limitation factor (-)
     real(r8), allocatable :: foln          (:)   ! foliage nitrogen (%)
     real(r8), allocatable :: dleaf         (:)   ! characteristic leaf dimension (m)
     real(r8), allocatable :: c3psn         (:)   ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), allocatable :: xl            (:)   ! leaf/stem orientation index
     real(r8), allocatable :: rhol          (:,:) ! leaf reflectance: 1=vis, 2=nir   (numrad)
     real(r8), allocatable :: rhos          (:,:) ! stem reflectance: 1=vis, 2=nir   (numrad)
     real(r8), allocatable :: taul          (:,:) ! leaf transmittance: 1=vis, 2=nir (numrad)
     real(r8), allocatable :: taus          (:,:) ! stem transmittance: 1=vis, 2=nir (numrad)
     real(r8), allocatable :: z0mr          (:)   ! ratio of momentum roughness length to canopy top height (-)
     real(r8), allocatable :: displar       (:)   ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: roota_par     (:)   ! rooting distribution parameter [1/m]
     real(r8), allocatable :: rootb_par     (:)   ! rooting distribution parameter [1/m]
     real(r8), allocatable :: rootprof_beta (:)   ! rooting distribution parameter for C and N inputs [unitless]
     real(r8), allocatable :: dwood         (:)   ! wood density (gC/m3)
     real(r8), allocatable :: slatop        (:)   ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real(r8), allocatable :: dsladlai      (:)   ! dSLA/dLAI, projected area basis [m^2/gC]
     real(r8), allocatable :: leafcn        (:)   ! leaf C:N (gC/gN)
     real(r8), allocatable :: flnr          (:)   ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real(r8), allocatable :: woody         (:)   ! binary flag for woody lifeform (1=woody, 0=not woody)
     real(r8), allocatable :: lflitcn       (:)   ! leaf litter C:N (gC/gN)
     real(r8), allocatable :: frootcn       (:)   ! fine root C:N (gC/gN)
     real(r8), allocatable :: livewdcn      (:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), allocatable :: deadwdcn      (:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), allocatable :: graincn       (:)   ! grain C:N (gC/gN) for prognostic crop model
     real(r8), allocatable :: froot_leaf    (:)   ! allocation parameter: new fine root C per new leaf C (gC/gC)
     real(r8), allocatable :: stem_leaf     (:)   ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), allocatable :: croot_stem    (:)   ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), allocatable :: flivewd       (:)   ! allocation parameter: fraction of new wood that is live 
                                                  ! (phloem and ray parenchyma) (no units)
     real(r8), allocatable :: fcur          (:)   ! allocation parameter: fraction of allocation that goes 
                                                  ! to currently displayed growth, remainder to storage
     real(r8), allocatable :: lf_flab       (:)   ! leaf litter labile fraction
     real(r8), allocatable :: lf_fcel       (:)   ! leaf litter cellulose fraction
     real(r8), allocatable :: lf_flig       (:)   ! leaf litter lignin fraction
     real(r8), allocatable :: fr_flab       (:)   ! fine root litter labile fraction
     real(r8), allocatable :: fr_fcel       (:)   ! fine root litter cellulose fraction
     real(r8), allocatable :: fr_flig       (:)   ! fine root litter lignin fraction
     real(r8), allocatable :: leaf_long     (:)   ! leaf longevity (yrs)
     real(r8), allocatable :: evergreen     (:)   ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), allocatable :: stress_decid  (:)   ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: season_decid  (:)   ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: cc_leaf       (:)   ! combustion completeness factor for leaf (0 to 1)
     real(r8), allocatable :: cc_lstem      (:)   ! combustion completeness factor for live stem (0 to 1)
     real(r8), allocatable :: cc_dstem      (:)   ! combustion completeness factor for dead stem (0 to 1)
     real(r8), allocatable :: cc_other      (:)   ! combustion completeness factor for other plant tissues (0 to 1)
     real(r8), allocatable :: fm_leaf       (:)   ! fire-related mortality factor for leaf (0 to 1)
     real(r8), allocatable :: fm_lstem      (:)   ! fire-related mortality factor for live stem (0 to 1)
     real(r8), allocatable :: fm_dstem      (:)   ! fire-related mortality factor for dead stem (0 to 1)
     real(r8), allocatable :: fm_other      (:)   ! fire-related mortality factor for other plant tissues (0 to 1)
     real(r8), allocatable :: fm_root       (:)   ! fire-related mortality factor for fine roots (0 to 1)
     real(r8), allocatable :: fm_lroot      (:)   ! fire-related mortality factor for live roots (0 to 1)
     real(r8), allocatable :: fm_droot      (:)   ! fire-related mortality factor for dead roots (0 to 1)
     real(r8), allocatable :: fertnitro     (:)   ! fertilizer applied (crop)
     real(r8), allocatable :: fleafcn       (:)   ! C:N during grain fill; leaf (crop)
     real(r8), allocatable :: ffrootcn      (:)   ! C:N during grain fill; froot (crop)
     real(r8), allocatable :: fstemcn       (:)   ! C:N during grain fill; stem (crop)
  end type Ecophyscon_type

  type(ecophyscon_type), public :: ecophyscon ! patch ecophysiological constants structure

contains

  !-----------------------------------------------------------------------
  subroutine ecophysconInit()
    !
    ! !USES:
    use clm_varpar, only : numrad, numpft 
    use pftvarcon , only : ntree, smpso, smpsc, fnitr
    use pftvarcon , only : z0mr, displar, dleaf, rhol, rhos, taul, taus, xl
    use pftvarcon , only : c3psn, slatop, dsladlai, leafcn, flnr, woody
    use pftvarcon , only : lflitcn, frootcn, livewdcn, deadwdcn, froot_leaf, stem_leaf, croot_stem
    use pftvarcon , only : flivewd, fcur, lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig
    use pftvarcon , only : leaf_long, evergreen, stress_decid, season_decid
    use pftvarcon , only : fertnitro, graincn, fleafcn, ffrootcn, fstemcn, dwood
    !
    ! !LOCAL VARIABLES:
    integer :: m, ib
    !------------------------------------------------------------------------

    allocate(ecophyscon%noveg         (0:numpft))        ; ecophyscon%noveg        (:)   =huge(1)
    allocate(ecophyscon%tree          (0:numpft))        ; ecophyscon%tree         (:)   =huge(1)
    allocate(ecophyscon%smpso         (0:numpft))        ; ecophyscon%smpso        (:)   =nan
    allocate(ecophyscon%smpsc         (0:numpft))        ; ecophyscon%smpsc        (:)   =nan
    allocate(ecophyscon%fnitr         (0:numpft))        ; ecophyscon%fnitr        (:)   =nan
    allocate(ecophyscon%foln          (0:numpft))        ; ecophyscon%foln         (:)   =nan
    allocate(ecophyscon%dleaf         (0:numpft))        ; ecophyscon%dleaf        (:)   =nan
    allocate(ecophyscon%c3psn         (0:numpft))        ; ecophyscon%c3psn        (:)   =nan
    allocate(ecophyscon%xl            (0:numpft))        ; ecophyscon%xl           (:)   =nan
    allocate(ecophyscon%rhol          (0:numpft,numrad)) ; ecophyscon%rhol         (:,:) =nan
    allocate(ecophyscon%rhos          (0:numpft,numrad)) ; ecophyscon%rhos         (:,:) =nan
    allocate(ecophyscon%taul          (0:numpft,numrad)) ; ecophyscon%taul         (:,:) =nan
    allocate(ecophyscon%taus          (0:numpft,numrad)) ; ecophyscon%taus         (:,:) =nan
    allocate(ecophyscon%z0mr          (0:numpft))        ; ecophyscon%z0mr         (:)   =nan
    allocate(ecophyscon%displar       (0:numpft))        ; ecophyscon%displar      (:)   =nan
    allocate(ecophyscon%roota_par     (0:numpft))        ; ecophyscon%roota_par    (:)   =nan
    allocate(ecophyscon%rootb_par     (0:numpft))        ; ecophyscon%rootb_par    (:)   =nan
    allocate(ecophyscon%slatop        (0:numpft))        ; ecophyscon%slatop       (:)   =nan
    allocate(ecophyscon%dsladlai      (0:numpft))        ; ecophyscon%dsladlai     (:)   =nan
    allocate(ecophyscon%leafcn        (0:numpft))        ; ecophyscon%leafcn       (:)   =nan
    allocate(ecophyscon%flnr          (0:numpft))        ; ecophyscon%flnr         (:)   =nan
    allocate(ecophyscon%woody         (0:numpft))        ; ecophyscon%woody        (:)   =nan
    allocate(ecophyscon%lflitcn       (0:numpft))        ; ecophyscon%lflitcn      (:)   =nan
    allocate(ecophyscon%frootcn       (0:numpft))        ; ecophyscon%frootcn      (:)   =nan
    allocate(ecophyscon%livewdcn      (0:numpft))        ; ecophyscon%livewdcn     (:)   =nan
    allocate(ecophyscon%deadwdcn      (0:numpft))        ; ecophyscon%deadwdcn     (:)   =nan
    allocate(ecophyscon%graincn       (0:numpft))        ; ecophyscon%graincn      (:)   =nan
    allocate(ecophyscon%froot_leaf    (0:numpft))        ; ecophyscon%froot_leaf   (:)   =nan
    allocate(ecophyscon%stem_leaf     (0:numpft))        ; ecophyscon%stem_leaf    (:)   =nan
    allocate(ecophyscon%croot_stem    (0:numpft))        ; ecophyscon%croot_stem   (:)   =nan
    allocate(ecophyscon%flivewd       (0:numpft))        ; ecophyscon%flivewd      (:)   =nan
    allocate(ecophyscon%fcur          (0:numpft))        ; ecophyscon%fcur         (:)   =nan
    allocate(ecophyscon%lf_flab       (0:numpft))        ; ecophyscon%lf_flab      (:)   =nan
    allocate(ecophyscon%lf_fcel       (0:numpft))        ; ecophyscon%lf_fcel      (:)   =nan
    allocate(ecophyscon%lf_flig       (0:numpft))        ; ecophyscon%lf_flig      (:)   =nan
    allocate(ecophyscon%fr_flab       (0:numpft))        ; ecophyscon%fr_flab      (:)   =nan
    allocate(ecophyscon%fr_fcel       (0:numpft))        ; ecophyscon%fr_fcel      (:)   =nan
    allocate(ecophyscon%fr_flig       (0:numpft))        ; ecophyscon%fr_flig      (:)   =nan
    allocate(ecophyscon%leaf_long     (0:numpft))        ; ecophyscon%leaf_long    (:)   =nan
    allocate(ecophyscon%evergreen     (0:numpft))        ; ecophyscon%evergreen    (:)   =nan
    allocate(ecophyscon%stress_decid  (0:numpft))        ; ecophyscon%stress_decid (:)   =nan
    allocate(ecophyscon%season_decid  (0:numpft))        ; ecophyscon%season_decid (:)   =nan
    allocate(ecophyscon%dwood         (0:numpft))        ; ecophyscon%dwood        (:)   =nan
    allocate(ecophyscon%rootprof_beta (0:numpft))        ; ecophyscon%rootprof_beta(:)   =nan
    allocate(ecophyscon%fertnitro     (0:numpft))        ; ecophyscon%fertnitro    (:)   =nan
    allocate(ecophyscon%fleafcn       (0:numpft))        ; ecophyscon%fleafcn      (:)   =nan
    allocate(ecophyscon%ffrootcn      (0:numpft))        ; ecophyscon%ffrootcn     (:)   =nan
    allocate(ecophyscon%fstemcn       (0:numpft))        ; ecophyscon%fstemcn      (:)   =nan

    do m = 0,numpft

       if (m <= ntree) then
          ecophyscon%tree(m) = 1
       else
          ecophyscon%tree(m) = 0
       end if

       do ib = 1,numrad
          ecophyscon%rhol(m,ib)   = rhol(m,ib)
          ecophyscon%rhos(m,ib)   = rhos(m,ib)
          ecophyscon%taul(m,ib)   = taul(m,ib)
          ecophyscon%taus(m,ib)   = taus(m,ib)
       end do

       ecophyscon%z0mr(m)         = z0mr(m)
       ecophyscon%displar(m)      = displar(m)
       ecophyscon%dleaf(m)        = dleaf(m)
       ecophyscon%xl(m)           = xl(m)
       ecophyscon%c3psn(m)        = c3psn(m)
       ecophyscon%slatop(m)       = slatop(m)
       ecophyscon%dsladlai(m)     = dsladlai(m)
       ecophyscon%leafcn(m)       = leafcn(m)
       ecophyscon%flnr(m)         = flnr(m)
       ecophyscon%smpso(m)        = smpso(m)
       ecophyscon%smpsc(m)        = smpsc(m)
       ecophyscon%fnitr(m)        = fnitr(m)
       ecophyscon%woody(m)        = woody(m)
       ecophyscon%lflitcn(m)      = lflitcn(m)
       ecophyscon%frootcn(m)      = frootcn(m)
       ecophyscon%livewdcn(m)     = livewdcn(m)
       ecophyscon%deadwdcn(m)     = deadwdcn(m)
       ecophyscon%graincn(m)      = graincn(m)
       ecophyscon%froot_leaf(m)   = froot_leaf(m)
       ecophyscon%stem_leaf(m)    = stem_leaf(m)
       ecophyscon%croot_stem(m)   = croot_stem(m)
       ecophyscon%flivewd(m)      = flivewd(m)
       ecophyscon%fcur(m)         = fcur(m)
       ecophyscon%lf_flab(m)      = lf_flab(m)
       ecophyscon%lf_fcel(m)      = lf_fcel(m)
       ecophyscon%lf_flig(m)      = lf_flig(m)
       ecophyscon%fr_flab(m)      = fr_flab(m)
       ecophyscon%fr_fcel(m)      = fr_fcel(m)
       ecophyscon%fr_flig(m)      = fr_flig(m)
       ecophyscon%leaf_long(m)    = leaf_long(m)
       ecophyscon%evergreen(m)    = evergreen(m)
       ecophyscon%stress_decid(m) = stress_decid(m)
       ecophyscon%season_decid(m) = season_decid(m)
       ecophyscon%dwood(m)        = dwood
       ecophyscon%fertnitro(m)    = fertnitro(m)
       ecophyscon%fleafcn(m)      = fleafcn(m)
       ecophyscon%ffrootcn(m)     = ffrootcn(m)
       ecophyscon%fstemcn(m)      = fstemcn(m)
    end do
  end subroutine ecophysconInit

end module EcophysConType

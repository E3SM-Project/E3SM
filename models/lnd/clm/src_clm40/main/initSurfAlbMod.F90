module initSurfalbMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initSurfalbMod
!
! !DESCRIPTION:
! Computes initial surface albedo calculation - 
! Initialization of ecosystem dynamics is needed for this
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use clm_varctl,   only : iulog, use_cn, use_cndv
!
! !PUBLIC TYPES:
  implicit none
  logical, public :: do_initsurfalb
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: InitSurfAlb
!
! !REVISION HISTORY:
! 2005-06-12: Created by Mariana Vertenstein
! 2008-02-29: Revised snow cover fraction from Niu and Yang, 2007
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: initSurfalb
!
! !INTERFACE:
  subroutine initSurfalb( calday, declin, declinm1 )
!
! !DESCRIPTION:
! The variable, h2osoi_vol, is needed by the soil albedo routine - this is not needed
! on restart since it is computed before the soil albedo computation is called.
! The remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines.
!
! !USES:
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use shr_orb_mod         , only : shr_orb_decl
    use shr_const_mod       , only : SHR_CONST_PI
    use clmtype
    use spmdMod             , only : masterproc,iam
    use decompMod           , only : get_proc_clumps, get_clump_bounds
    use filterMod           , only : filter
    use clm_varpar          , only : nlevsoi, nlevsno, nlevlak, nlevgrnd
    use clm_varcon          , only : zlnd, istsoil, isturb, denice, denh2o, &
                                     icol_roof, icol_road_imperv, &
                                     icol_road_perv
    use clm_varcon          , only : istcrop
    use clm_time_manager        , only : get_step_size
    use FracWetMod          , only : FracWet
    use SurfaceAlbedoMod    , only : SurfaceAlbedo
    use CNEcosystemDynMod   , only : CNEcosystemDyn
    use CNVegStructUpdateMod, only : CNVegStructUpdate
    use STATICEcosysDynMod  , only : EcosystemDyn, interpMonthlyVeg
    use UrbanMod            , only : UrbanAlbedo
    use abortutils          , only : endrun
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: calday               ! calendar day for declin
    real(r8), intent(in) :: declin               ! declination angle (radians) for calday
    real(r8), intent(in), optional :: declinm1   ! declination angle (radians) for caldaym1
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: plandunit(:)      ! landunit index associated with each pft
    integer , pointer :: ctype(:)          ! column type
    integer , pointer :: clandunit(:)      ! landunit index associated with each column
    integer,  pointer :: pgridcell(:)      ! gridcell associated with each pft
    integer , pointer :: itypelun(:) 	   ! landunit type
    logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
    real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
    integer , pointer :: frac_veg_nosno_alb(:) ! fraction of vegetation not covered by snow (0 OR 1) [-] 
    real(r8), pointer :: dayl(:)           ! daylength (seconds)
    real(r8), pointer :: latdeg(:)         ! latitude (degrees)
    integer , pointer :: pcolumn(:)        ! index into column level quantities
    real(r8), pointer :: soilpsi(:,:)      ! soil water potential in each soil layer (MPa)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: snowdp(:)         ! snow height (m)
    real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(r8), pointer :: fwet(:)           ! fraction of canopy that is wet (0 to 1) (pft-level)
    real(r8), pointer :: decl(:)           ! solar declination angle (radians)
!
! local pointers to implicit out arguments (lake points only)
!
    real(r8), pointer :: fdry(:)     ! fraction of foliage that is green and dry [-] (new)
    real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)     ! one-sided stem area index, no burying by snow
    real(r8), pointer :: htop(:)     ! canopy top (m)
    real(r8), pointer :: hbot(:)     ! canopy bottom (m)
    real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: nc,j,l,c,p,fc ! indices
    integer :: nclumps       ! number of clumps on this processor
    integer :: begp, endp    ! per-clump beginning and ending pft indices
    integer :: begc, endc    ! per-clump beginning and ending column indices
    integer :: begl, endl    ! per-clump beginning and ending landunit indices
    integer :: begg, endg    ! per-clump gridcell ending gridcell indices
    integer :: ier           ! MPI return code
    real(r8):: lat           ! latitude (radians) for daylength calculation
    real(r8):: temp          ! temporary variable for daylength
    real(r8):: snowbd        ! temporary calculation of snow bulk density (kg/m3)
    real(r8):: fmelt         ! snowbd/100
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (landunit-level)
    
    lakpoi              => lun%lakpoi
    itypelun            => lun%itype

    ! Assign local pointers to derived subtypes components (column-level)

    dz                  => cps%dz
    h2osoi_ice          => cws%h2osoi_ice
    h2osoi_liq          => cws%h2osoi_liq
    h2osoi_vol          => cws%h2osoi_vol
    snowdp              => cps%snowdp
    h2osno              => cws%h2osno
    frac_sno            => cps%frac_sno 
    ctype               => col%itype
    clandunit           => col%landunit
    soilpsi             => cps%soilpsi

    ! Assign local pointers to derived subtypes components (pft-level)

    plandunit          => pft%landunit
    frac_veg_nosno_alb => pps%frac_veg_nosno_alb
    frac_veg_nosno     => pps%frac_veg_nosno
    fwet               => pps%fwet

    ! Assign local pointers to derived subtypes components (pft-level)
    ! The folowing pointers will only be used for lake points in this routine

    htop               => pps%htop
    hbot               => pps%hbot
    tlai               => pps%tlai
    tsai               => pps%tsai
    elai               => pps%elai
    esai               => pps%esai
    fdry               => pps%fdry

    decl      => cps%decl
    dayl      => pepv%dayl
    pcolumn   => pft%column
    pgridcell => pft%gridcell
    latdeg    => grc%latdeg 

    ! ========================================================================
    ! Determine surface albedo - initialized by calls to ecosystem dynamics and
    ! albedo subroutines. Note: elai, esai, frac_veg_nosno_alb are computed in
    ! Ecosysdyn and needed by routines FracWet and SurfaceAlbedo and 
    ! frac_veg_nosno is needed by FracWet
    ! fwet is needed in routine TwoStream (called by SurfaceAlbedo)
    ! frac_sno is needed by SoilAlbedo (called by SurfaceAlbedo)
    ! ========================================================================

    if (.not. use_cn) then
       ! the default mode uses prescribed vegetation structure
       ! Read monthly vegetation data for interpolation to daily values
       
       call interpMonthlyVeg()
    end if

    ! Determine clump bounds for this processor

    nclumps = get_proc_clumps()

    ! Loop over clumps on this processor
!$OMP PARALLEL DO PRIVATE (nc,p,j,l,c,fc,begg,endg,begl,endl,begc,endc,begp,endp,lat,temp,snowbd,fmelt)
    do nc = 1,nclumps

       ! Determine clump bounds

       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

       ! Determine variables needed by SurfaceAlbedo for lake points

       do p = begp,endp
          l = plandunit(p)
          if (lakpoi(l)) then
             fwet(p) = 0._r8
             fdry(p) = 0._r8
             elai(p) = 0._r8
             esai(p) = 0._r8
             htop(p) = 0._r8
             hbot(p) = 0._r8
             tlai(p) = 0._r8
             tsai(p) = 0._r8
             frac_veg_nosno_alb(p) = 0._r8
             frac_veg_nosno(p) = 0._r8
          end if
       end do

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, or static parameterizations
       ! ============================================================================

       if (use_cn) then
          do j = 1, nlevgrnd
             do fc = 1, filter(nc)%num_soilc
                c = filter(nc)%soilc(fc)
                soilpsi(c,j) = -15.0_r8
             end do
          end do
       end if

       ! Determine variables needed for SurfaceAlbedo for non-lake points

       if (use_cn) then
          ! CN initialization is done only on the soil landunits.

          if (.not. present(declinm1)) then
             write(iulog,*)'declination for the previous timestep (declinm1) must be ',&
                  ' present as argument in CN mode'
             call endrun()
          end if
          
          ! it is necessary to initialize the solar declination for the previous
          ! timestep (caldaym1) so that the CNphenology routines know if this is 
          ! before or after the summer solstice.
          
          ! declination for previous timestep
          do c = begc, endc
             l = clandunit(c)
             if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then
                decl(c) = declinm1
             end if
          end do
          
          ! daylength for previous timestep
          do p = begp, endp
             c = pcolumn(p)
             l = plandunit(p)
             if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then
                lat = latdeg(pgridcell(p)) * SHR_CONST_PI / 180._r8
                temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
                temp = min(1._r8,max(-1._r8,temp))
                dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp) 
             end if
          end do
          
          ! declination for current timestep
          do c = begc, endc
             l = clandunit(c)
             if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then
                decl(c) = declin
             end if
          end do
          
          call CNEcosystemDyn(begc, endc, begp, endp, filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb=.true.)

       else

          ! this is the default call if CN not set

          call EcosystemDyn(begp, endp, filter(nc)%num_nolakep, filter(nc)%nolakep, &
               doalb=.true.)
       end if

       do p = begp, endp
          l = plandunit(p)
          if (.not. lakpoi(l)) then
             frac_veg_nosno(p) = frac_veg_nosno_alb(p)
             fwet(p) = 0._r8
          end if
       end do
       
       call FracWet(filter(nc)%num_nolakep, filter(nc)%nolakep)
       
       ! Compute Surface Albedo - all land points (including lake) other than urban
       ! Needs as input fracion of soil covered by snow (Z.-L. Yang U. Texas)

       do c = begc, endc
          l = clandunit(c)
          if (itypelun(l) == isturb) then
             ! From Bonan 1996 (LSM technical note)
             frac_sno(c) = min( snowdp(c)/0.05_r8, 1._r8)
          else
             frac_sno(c) = 0._r8
             ! snow cover fraction as in Niu and Yang 2007
             if(snowdp(c) .gt. 0.0)  then
                snowbd   = min(800._r8,h2osno(c)/snowdp(c)) !bulk density of snow (kg/m3)
                fmelt    = (snowbd/100.)**1.
                ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
                ! reconsidered, optimal value of 1.5 in Niu et al., 2007
                frac_sno(c) = tanh( snowdp(c) /(2.5 * zlnd * fmelt) )
             endif
          end if
       end do
       call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                          filter(nc)%num_nourbanc, filter(nc)%nourbanc, &
                          filter(nc)%num_nourbanp, filter(nc)%nourbanp, &
                          calday, declin)
       

       ! Determine albedos for urban landunits

       if (filter(nc)%num_urbanl > 0) then
          call UrbanAlbedo(nc, begl, endl, begc, endc, begp, endp, &
                           filter(nc)%num_urbanl, filter(nc)%urbanl, &
                           filter(nc)%num_urbanc, filter(nc)%urbanc, &
                           filter(nc)%num_urbanp, filter(nc)%urbanp )

       end if

    end do   ! end of loop over clumps
!$OMP END PARALLEL DO

  end subroutine initSurfalb

end module initSurfalbMod

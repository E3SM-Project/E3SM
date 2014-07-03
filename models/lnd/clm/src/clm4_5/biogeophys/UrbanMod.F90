module UrbanMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanMod
! 
! !DESCRIPTION: 
! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : numrad
  use clm_varcon  , only : isecspday, degpsec
  use clm_varctl  , only : iulog
  use abortutils  , only : endrun  
  use shr_sys_mod , only : shr_sys_flush 
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanClumpInit    ! Initialization of urban clump data structure
  public :: UrbanRadiation    ! Urban radiative fluxes
  public :: UrbanAlbedo       ! Urban albedos  
  public :: UrbanSnowAlbedo   ! Urban snow albedos
  public :: UrbanFluxes       ! Urban turbulent fluxes

! !Urban control variables
  character(len= *), parameter, public :: urban_hac_off = 'OFF'               ! 
  character(len= *), parameter, public :: urban_hac_on =  'ON'                ! 
  character(len= *), parameter, public :: urban_wasteheat_on = 'ON_WASTEHEAT' ! 
  character(len= 16), public :: urban_hac = urban_hac_off
  logical, public :: urban_traffic = .false.        ! urban traffic fluxes
!
! !REVISION HISTORY:
! Created by Gordon Bonan and Mariana Vertenstein and Keith Oleson 04/2003
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS
  private :: view_factor      ! View factors for road and one wall
  private :: incident_direct  ! Direct beam solar rad incident on walls and road in urban canyon 
  private :: incident_diffuse ! Diffuse solar rad incident on walls and road in urban canyon
  private :: net_solar        ! Solar radiation absorbed by road and both walls in urban canyon 
  private :: net_longwave     ! Net longwave radiation for road and both walls in urban canyon 

! PRIVATE TYPES
  private
  type urban_clump_t
     real(r8), pointer :: canyon_hwr(:)            ! ratio of building height to street width 
     real(r8), pointer :: wtroad_perv(:)           ! weight of pervious road wrt total road
     real(r8), pointer :: ht_roof(:)               ! height of urban roof (m)
     real(r8), pointer :: wtlunit_roof(:)          ! weight of roof with respect to landunit
     real(r8), pointer :: wind_hgt_canyon(:)       ! height above road at which wind in canyon is to be computed (m)
     real(r8), pointer :: em_roof(:)               ! roof emissivity
     real(r8), pointer :: em_improad(:)            ! impervious road emissivity
     real(r8), pointer :: em_perroad(:)            ! pervious road emissivity
     real(r8), pointer :: em_wall(:)               ! wall emissivity
     real(r8), pointer :: alb_roof_dir(:,:)        ! direct  roof albedo
     real(r8), pointer :: alb_roof_dif(:,:)        ! diffuse roof albedo
     real(r8), pointer :: alb_improad_dir(:,:)     ! direct  impervious road albedo
     real(r8), pointer :: alb_improad_dif(:,:)     ! diffuse impervious road albedo
     real(r8), pointer :: alb_perroad_dir(:,:)     ! direct  pervious road albedo
     real(r8), pointer :: alb_perroad_dif(:,:)     ! diffuse pervious road albedo
     real(r8), pointer :: alb_wall_dir(:,:)        ! direct  wall albedo
     real(r8), pointer :: alb_wall_dif(:,:)        ! diffuse wall albedo
  end type urban_clump_t

  type (urban_clump_t), private, pointer :: urban_clump(:)  ! array of urban clumps for this processor

  integer,  private, parameter :: noonsec   = isecspday / 2 ! seconds at local noon
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanAlbedo
!
! !INTERFACE:
  subroutine UrbanAlbedo (nc, lbl, ubl, lbc, ubc, lbp, ubp, &
                          num_urbanl, filter_urbanl, &
                          num_urbanc, filter_urbanc, &
                          num_urbanp, filter_urbanp)
!
! !DESCRIPTION: 
! Determine urban landunit component albedos
!
! !USES:
    use clmtype
    use shr_orb_mod  , only : shr_orb_decl, shr_orb_cosz
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv, &
                              sb
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nc                        ! clump index
    integer,  intent(in) :: lbl, ubl                  ! landunit-index bounds
    integer,  intent(in) :: lbc, ubc                  ! column-index bounds
    integer,  intent(in) :: lbp, ubp                  ! pft-index bounds
    integer , intent(in) :: num_urbanl                ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(ubl-lbl+1) ! urban landunit filter
    integer , intent(in) :: num_urbanc                ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(ubc-lbc+1) ! urban column filter
    integer , intent(in) :: num_urbanp                ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(ubp-lbp+1) ! urban pft filter
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 03/2003, Mariana Vertenstein: Migrated to clm2.2 
! 01/2008, Erik Kluzek:         Migrated to clm3.5.15
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: pgridcell(:) ! gridcell of corresponding pft
    integer , pointer :: lgridcell(:) ! gridcell of corresponding landunit
    integer , pointer :: clandunit(:) ! column's landunit
    integer , pointer :: cgridcell(:) ! gridcell of corresponding column
    integer , pointer :: coli(:)      ! beginning column index for landunit 
    integer , pointer :: colf(:)      ! ending column index for landunit
    integer , pointer :: ctype(:)     ! column type
    integer , pointer :: pcolumn(:)   ! column of corresponding pft
    real(r8), pointer :: czen(:)      ! cosine of solar zenith angle for each column
    real(r8), pointer :: lat(:)       ! latitude (radians)
    real(r8), pointer :: lon(:)       ! longitude (radians)
    real(r8), pointer :: frac_sno(:)  ! fraction of ground covered by snow (0 to 1)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: albgrd(:,:)  ! ground albedo (direct)
    real(r8), pointer :: albgri(:,:)  ! ground albedo (diffuse)
    real(r8), pointer :: albd(:,:)    ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)    ! surface albedo (diffuse)
    real(r8), pointer :: fabd(:,:)    ! flux absorbed by veg per unit direct  flux
    real(r8), pointer :: fabd_sun(:,:)! flux absorbed by sunlit leaf per unit direct flux
    real(r8), pointer :: fabd_sha(:,:)! flux absorbed by shaded leaf per unit direct flux
    real(r8), pointer :: fabi(:,:)    ! flux absorbed by veg per unit diffuse flux
    real(r8), pointer :: fabi_sun(:,:)! flux absorbed by sunlit leaf per unit diffuse flux
    real(r8), pointer :: fabi_sha(:,:)! flux absorbed by shaded leaf per unit diffuse flux
    real(r8), pointer :: ftdd(:,:)    ! down direct  flux below veg per unit dir flx
    real(r8), pointer :: ftid(:,:)    ! down diffuse flux below veg per unit dir flx
    real(r8), pointer :: ftii(:,:)    ! down diffuse flux below veg per unit dif flx
    real(r8), pointer :: fsun(:)      ! sunlit fraction of canopy
    real(r8), pointer :: vf_sr(:)     ! view factor of sky for road
    real(r8), pointer :: vf_wr(:)     ! view factor of one wall for road
    real(r8), pointer :: vf_sw(:)     ! view factor of sky for one wall
    real(r8), pointer :: vf_rw(:)     ! view factor of road for one wall
    real(r8), pointer :: vf_ww(:)     ! view factor of opposing wall for one wall
    real(r8), pointer :: sabs_roof_dir(:,:)       ! direct  solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_roof_dif(:,:)       ! diffuse solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dir(:,:)    ! direct  solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dif(:,:)    ! diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dir(:,:)  ! direct  solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dif(:,:)  ! diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_improad_dir(:,:)    ! direct  solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_improad_dif(:,:)    ! diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dir(:,:)    ! direct  solar absorbed  by pervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dif(:,:)    ! diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
!
!
! !OTHER LOCAL VARIABLES
!EOP
!
    real(r8) :: coszen(num_urbanl)                 ! cosine solar zenith angle
    real(r8) :: coszen_pft(num_urbanp)             ! cosine solar zenith angle for next time step (pft level)
    real(r8) :: zen(num_urbanl)                    ! solar zenith angle (radians)
    real(r8) :: sdir(num_urbanl, numrad)           ! direct beam solar radiation on horizontal surface
    real(r8) :: sdif(num_urbanl, numrad)           ! diffuse solar radiation on horizontal surface

    real(r8) :: sdir_road(num_urbanl, numrad)      ! direct beam solar radiation incident on road
    real(r8) :: sdif_road(num_urbanl, numrad)      ! diffuse solar radiation incident on road
    real(r8) :: sdir_sunwall(num_urbanl, numrad)   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdif_sunwall(num_urbanl, numrad)   ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdir_shadewall(num_urbanl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_shadewall(num_urbanl, numrad) ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: albsnd_roof(num_urbanl,numrad)     ! snow albedo for roof (direct)
    real(r8) :: albsni_roof(num_urbanl,numrad)     ! snow albedo for roof (diffuse)
    real(r8) :: albsnd_improad(num_urbanl,numrad)  ! snow albedo for impervious road (direct)
    real(r8) :: albsni_improad(num_urbanl,numrad)  ! snow albedo for impervious road (diffuse)
    real(r8) :: albsnd_perroad(num_urbanl,numrad)  ! snow albedo for pervious road (direct)
    real(r8) :: albsni_perroad(num_urbanl,numrad)  ! snow albedo for pervious road (diffuse)
                                       
    integer  :: fl,fp,fc,g,l,p,c,ib                ! indices
    integer  :: ic                                 ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: num_solar                          ! counter
    real(r8) :: alb_roof_dir_s(num_urbanl,numrad)    ! direct roof albedo with snow effects
    real(r8) :: alb_roof_dif_s(num_urbanl,numrad)    ! diffuse roof albedo with snow effects
    real(r8) :: alb_improad_dir_s(num_urbanl,numrad) ! direct impervious road albedo with snow effects
    real(r8) :: alb_perroad_dir_s(num_urbanl,numrad) ! direct pervious road albedo with snow effects
    real(r8) :: alb_improad_dif_s(num_urbanl,numrad) ! diffuse impervious road albedo with snow effects
    real(r8) :: alb_perroad_dif_s(num_urbanl,numrad) ! diffuse pervious road albedo with snow effects
    real(r8) :: sref_roof_dir(num_urbanl,numrad)      ! direct  solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_roof_dif(num_urbanl,numrad)      ! diffuse solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_sunwall_dir(num_urbanl,numrad)   ! direct  solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_sunwall_dif(num_urbanl,numrad)   ! diffuse solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dir(num_urbanl,numrad) ! direct  solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dif(num_urbanl,numrad) ! diffuse solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_improad_dir(num_urbanl,numrad)   ! direct  solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_improad_dif(num_urbanl,numrad)   ! diffuse solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dir(num_urbanl,numrad)   ! direct  solar reflected by pervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dif(num_urbanl,numrad)   ! diffuse solar reflected by pervious road per unit ground area per unit incident flux   
    real(r8), pointer :: canyon_hwr(:)             ! ratio of building height to street width
    real(r8), pointer :: wtroad_perv(:)            ! weight of pervious road wrt total road
    real(r8), pointer :: alb_roof_dir(:,:)         ! direct roof albedo
    real(r8), pointer :: alb_roof_dif(:,:)         ! diffuse roof albedo
    real(r8), pointer :: alb_improad_dir(:,:)      ! direct impervious road albedo
    real(r8), pointer :: alb_perroad_dir(:,:)      ! direct pervious road albedo
    real(r8), pointer :: alb_improad_dif(:,:)      ! diffuse imprevious road albedo
    real(r8), pointer :: alb_perroad_dif(:,:)      ! diffuse pervious road albedo
    real(r8), pointer :: alb_wall_dir(:,:)         ! direct wall albedo
    real(r8), pointer :: alb_wall_dif(:,:)         ! diffuse wall albedo
!-----------------------------------------------------------------------

    ! Assign pointers into module urban clumps

    canyon_hwr         => urban_clump(nc)%canyon_hwr
    wtroad_perv        => urban_clump(nc)%wtroad_perv
    alb_roof_dir       => urban_clump(nc)%alb_roof_dir
    alb_roof_dif       => urban_clump(nc)%alb_roof_dif  
    alb_improad_dir    => urban_clump(nc)%alb_improad_dir  
    alb_improad_dif    => urban_clump(nc)%alb_improad_dif  
    alb_perroad_dir    => urban_clump(nc)%alb_perroad_dir  
    alb_perroad_dif    => urban_clump(nc)%alb_perroad_dif  
    alb_wall_dir       => urban_clump(nc)%alb_wall_dir  
    alb_wall_dif       => urban_clump(nc)%alb_wall_dif  

    ! Assign gridcell level pointers

    lat                => clm3%g%lat
    lon                => clm3%g%lon

    ! Assign landunit level pointer

    lgridcell          => clm3%g%l%gridcell
    coli               => clm3%g%l%coli
    colf               => clm3%g%l%colf
    vf_sr              => clm3%g%l%lps%vf_sr
    vf_wr              => clm3%g%l%lps%vf_wr
    vf_sw              => clm3%g%l%lps%vf_sw
    vf_rw              => clm3%g%l%lps%vf_rw
    vf_ww              => clm3%g%l%lps%vf_ww
    sabs_roof_dir      => clm3%g%l%lps%sabs_roof_dir
    sabs_roof_dif      => clm3%g%l%lps%sabs_roof_dif
    sabs_sunwall_dir   => clm3%g%l%lps%sabs_sunwall_dir
    sabs_sunwall_dif   => clm3%g%l%lps%sabs_sunwall_dif
    sabs_shadewall_dir => clm3%g%l%lps%sabs_shadewall_dir
    sabs_shadewall_dif => clm3%g%l%lps%sabs_shadewall_dif
    sabs_improad_dir   => clm3%g%l%lps%sabs_improad_dir
    sabs_improad_dif   => clm3%g%l%lps%sabs_improad_dif
    sabs_perroad_dir   => clm3%g%l%lps%sabs_perroad_dir
    sabs_perroad_dif   => clm3%g%l%lps%sabs_perroad_dif

    ! Assign column level pointers

    ctype              => clm3%g%l%c%itype
    albgrd             => clm3%g%l%c%cps%albgrd
    albgri             => clm3%g%l%c%cps%albgri
    frac_sno           => clm3%g%l%c%cps%frac_sno
    clandunit          => clm3%g%l%c%landunit
    cgridcell          => clm3%g%l%c%gridcell
    czen               => clm3%g%l%c%cps%coszen

    ! Assign pft  level pointers

    pgridcell          => clm3%g%l%c%p%gridcell
    pcolumn            => clm3%g%l%c%p%column
    albd               => clm3%g%l%c%p%pps%albd
    albi               => clm3%g%l%c%p%pps%albi
    fabd               => clm3%g%l%c%p%pps%fabd
    fabd_sun           => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha           => clm3%g%l%c%p%pps%fabd_sha
    fabi               => clm3%g%l%c%p%pps%fabi
    fabi_sun           => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha           => clm3%g%l%c%p%pps%fabi_sha
    ftdd               => clm3%g%l%c%p%pps%ftdd
    ftid               => clm3%g%l%c%p%pps%ftid
    ftii               => clm3%g%l%c%p%pps%ftii
    fsun               => clm3%g%l%c%p%pps%fsun


    ! ----------------------------------------------------------------------------
    ! Solar declination and cosine solar zenith angle and zenith angle for 
    ! next time step
    ! ----------------------------------------------------------------------------

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)
       coszen(fl) = czen(coli(l))               ! Assumes coszen for each column are the same
       zen(fl)    = acos(coszen(fl))
    end do

    do fp = 1,num_urbanp
       p = filter_urbanp(fp)
       g = pgridcell(p)
       c = pcolumn(p)
       coszen_pft(fp) = czen(c)
    end do

    ! ----------------------------------------------------------------------------
    ! Initialize clmtype output since solar radiation is only done if coszen > 0
    ! ----------------------------------------------------------------------------
    
    do ib = 1,numrad
       do fc = 1,num_urbanc
          c = filter_urbanc(fc)

          albgrd(c,ib) = 0._r8
          albgri(c,ib) = 0._r8
       end do

       do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          g = pgridcell(p)
          albd(p,ib)    = 1._r8
          albi(p,ib)    = 1._r8
          fabd(p,ib)    = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib)    = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          if (coszen_pft(fp) > 0._r8) then
             ftdd(p,ib) = 1._r8
          else
             ftdd(p,ib) = 0._r8
          end if
          ftid(p,ib)    = 0._r8
          if (coszen_pft(fp) > 0._r8) then
             ftii(p,ib) = 1._r8
          else
             ftii(p,ib) = 0._r8
          end if
          if (ib == 1) then
             fsun(p)    = 0._r8
          end if
       end do
    end do

    ! ----------------------------------------------------------------------------
    ! Urban Code
    ! ----------------------------------------------------------------------------

    num_solar = 0
    do fl = 1,num_urbanl
       if (coszen(fl) > 0._r8) num_solar = num_solar + 1
    end do

    ! Initialize urban clump components

    do ib = 1,numrad
       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          sabs_roof_dir(l,ib)      = 0._r8
          sabs_roof_dif(l,ib)      = 0._r8
          sabs_sunwall_dir(l,ib)   = 0._r8
          sabs_sunwall_dif(l,ib)   = 0._r8
          sabs_shadewall_dir(l,ib) = 0._r8
          sabs_shadewall_dif(l,ib) = 0._r8
          sabs_improad_dir(l,ib)   = 0._r8
          sabs_improad_dif(l,ib)   = 0._r8
          sabs_perroad_dir(l,ib)   = 0._r8
          sabs_perroad_dif(l,ib)   = 0._r8
          sref_roof_dir(fl,ib)      = 1._r8
          sref_roof_dif(fl,ib)      = 1._r8
          sref_sunwall_dir(fl,ib)   = 1._r8
          sref_sunwall_dif(fl,ib)   = 1._r8
          sref_shadewall_dir(fl,ib) = 1._r8
          sref_shadewall_dif(fl,ib) = 1._r8
          sref_improad_dir(fl,ib)   = 1._r8
          sref_improad_dif(fl,ib)   = 1._r8
          sref_perroad_dir(fl,ib)   = 1._r8
          sref_perroad_dif(fl,ib)   = 1._r8
       end do
    end do

    ! View factors for road and one wall in urban canyon (depends only on canyon_hwr)
    
    if (num_urbanl .gt. 0) then
       call view_factor (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr)
    end if

    ! ----------------------------------------------------------------------------
    ! Only do the rest if all coszen are positive 
    ! ----------------------------------------------------------------------------

    if (num_solar > 0)then

       ! Set constants - solar fluxes are per unit incoming flux

       do ib = 1,numrad
          do fl = 1,num_urbanl
             sdir(fl,ib) = 1._r8
             sdif(fl,ib) = 1._r8
          end do
       end do

       ! Incident direct beam radiation for 
       ! (a) roof and (b) road and both walls in urban canyon
          
       if (num_urbanl .gt. 0) then
          call incident_direct (lbl, ubl, num_urbanl, canyon_hwr, coszen, zen, sdir, sdir_road, sdir_sunwall, sdir_shadewall)
       end if
       
       ! Incident diffuse radiation for 
       ! (a) roof and (b) road and both walls in urban canyon.
       
       if (num_urbanl .gt. 0) then
          call incident_diffuse (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr, sdif, sdif_road, &
                                 sdif_sunwall, sdif_shadewall)
       end if

       ! Get snow albedos for roof and impervious and pervious road
       if (num_urbanl .gt. 0) then
          ic = 0; call UrbanSnowAlbedo(lbl, ubl, num_urbanl, filter_urbanl, coszen, ic, albsnd_roof, albsnd_improad, albsnd_perroad)
          ic = 1; call UrbanSnowAlbedo(lbl, ubl, num_urbanl, filter_urbanl, coszen, ic, albsni_roof, albsni_improad, albsni_perroad)
       end if

       ! Combine snow-free and snow albedos
       do ib = 1,numrad
          do fl = 1,num_urbanl
             l = filter_urbanl(fl)
             do c = coli(l),colf(l)
                if (ctype(c) == icol_roof) then    
                   alb_roof_dir_s(fl,ib) = alb_roof_dir(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsnd_roof(fl,ib)*frac_sno(c)
                   alb_roof_dif_s(fl,ib) = alb_roof_dif(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsni_roof(fl,ib)*frac_sno(c)
                else if (ctype(c) == icol_road_imperv) then    
                   alb_improad_dir_s(fl,ib) = alb_improad_dir(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsnd_improad(fl,ib)*frac_sno(c)
                   alb_improad_dif_s(fl,ib) = alb_improad_dif(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsni_improad(fl,ib)*frac_sno(c)
                else if (ctype(c) == icol_road_perv) then    
                   alb_perroad_dir_s(fl,ib) = alb_perroad_dir(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsnd_perroad(fl,ib)*frac_sno(c)
                   alb_perroad_dif_s(fl,ib) = alb_perroad_dif(fl,ib)*(1._r8-frac_sno(c))  &
                                        + albsni_perroad(fl,ib)*frac_sno(c)
                end if
             end do
          end do
       end do
       
       ! Reflected and absorbed solar radiation per unit incident radiation 
       ! for road and both walls in urban canyon allowing for multiple reflection
       ! Reflected and absorbed solar radiation per unit incident radiation for roof
       
       if (num_urbanl .gt. 0) then
          call net_solar (lbl, ubl, num_urbanl, filter_urbanl, coszen, canyon_hwr, wtroad_perv, sdir, sdif, &
                          alb_improad_dir_s, alb_perroad_dir_s, alb_wall_dir, alb_roof_dir_s, &
                          alb_improad_dif_s, alb_perroad_dif_s, alb_wall_dif, alb_roof_dif_s, &
                          sdir_road, sdir_sunwall, sdir_shadewall,  &
                          sdif_road, sdif_sunwall, sdif_shadewall,  &
                          sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir, sref_roof_dir, &
                          sref_improad_dif, sref_perroad_dif, sref_sunwall_dif, sref_shadewall_dif, sref_roof_dif)
       end if

       ! ----------------------------------------------------------------------------
       ! Map urban output to clmtype components 
       ! ----------------------------------------------------------------------------

       !  Set albgrd and albgri (ground albedos) and albd and albi (surface albedos)

       do ib = 1,numrad
          do fl = 1,num_urbanl
             l = filter_urbanl(fl)
             do c = coli(l),colf(l)
                if (ctype(c) == icol_roof) then    
                   albgrd(c,ib) = sref_roof_dir(fl,ib) 
                   albgri(c,ib) = sref_roof_dif(fl,ib) 
                else if (ctype(c) == icol_sunwall) then   
                   albgrd(c,ib) = sref_sunwall_dir(fl,ib)
                   albgri(c,ib) = sref_sunwall_dif(fl,ib)
                else if (ctype(c) == icol_shadewall) then 
                   albgrd(c,ib) = sref_shadewall_dir(fl,ib)
                   albgri(c,ib) = sref_shadewall_dif(fl,ib)
                else if (ctype(c) == icol_road_perv) then
                   albgrd(c,ib) = sref_perroad_dir(fl,ib)
                   albgri(c,ib) = sref_perroad_dif(fl,ib)
                else if (ctype(c) == icol_road_imperv) then
                   albgrd(c,ib) = sref_improad_dir(fl,ib)
                   albgri(c,ib) = sref_improad_dif(fl,ib)
                endif
             end do
          end do
          do fp = 1,num_urbanp
             p = filter_urbanp(fp)
             c = pcolumn(p)
             albd(p,ib) = albgrd(c,ib)
             albi(p,ib) = albgri(c,ib)
          end do
       end do
    end if

  end subroutine UrbanAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanSnowAlbedo
!
! !INTERFACE:
  subroutine UrbanSnowAlbedo (lbl, ubl, num_urbanl, filter_urbanl, coszen, ind, &
                              albsn_roof, albsn_improad, albsn_perroad)
!
! !DESCRIPTION:
! Determine urban snow albedos
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon   , only : icol_roof, icol_road_perv, icol_road_imperv
!
! !ARGUMENTS:
    implicit none
    integer,  intent(in) :: lbl, ubl                    ! landunit-index bounds
    integer , intent(in) :: num_urbanl                  ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(ubl-lbl+1)    ! urban landunit filter
    integer , intent(in) :: ind                         ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(in) :: coszen(num_urbanl)          ! cosine solar zenith angle
    real(r8), intent(out):: albsn_roof(num_urbanl,2)    ! roof snow albedo by waveband (assume 2 wavebands)
    real(r8), intent(out):: albsn_improad(num_urbanl,2) ! impervious road snow albedo by waveband (assume 2 wavebands)
    real(r8), intent(out):: albsn_perroad(num_urbanl,2) ! pervious road snow albedo by waveband (assume 2 wavebands)
!
! !CALLED FROM:
! subroutine UrbanAlbedo in this module
!
! !REVISION HISTORY:
! Author: Keith Oleson 9/2005
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
    integer , pointer :: coli(:)      ! beginning column index for landunit
    integer , pointer :: colf(:)      ! ending column index for landunit
    real(r8), pointer :: h2osno(:)    ! snow water (mm H2O)
    integer , pointer :: ctype(:)     ! column type
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer  :: fl,c,l              ! indices
!
! variables and constants for snow albedo calculation
!
! These values are derived from Marshall (1989) assuming soot content of 1.5e-5 
! (three times what LSM uses globally). Note that snow age effects are ignored here.
    real(r8), parameter :: snal0 = 0.66_r8 ! vis albedo of urban snow
    real(r8), parameter :: snal1 = 0.56_r8 ! nir albedo of urban snow
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit level)

    coli       => clm3%g%l%coli
    colf       => clm3%g%l%colf

    ! Assign local pointers to derived subtypes components (column-level)

    ctype      => clm3%g%l%c%itype
    h2osno     => clm3%g%l%c%cws%h2osno

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       do c = coli(l),colf(l)
          if (coszen(fl) > 0._r8 .and. h2osno(c) > 0._r8) then
             if (ctype(c) == icol_roof) then
                albsn_roof(fl,1) = snal0
                albsn_roof(fl,2) = snal1
             else if (ctype(c) == icol_road_imperv) then
                albsn_improad(fl,1) = snal0
                albsn_improad(fl,2) = snal1
             else if (ctype(c) == icol_road_perv) then
                albsn_perroad(fl,1) = snal0
                albsn_perroad(fl,2) = snal1
             end if
          else
             if (ctype(c) == icol_roof) then
                albsn_roof(fl,1) = 0._r8
                albsn_roof(fl,2) = 0._r8
             else if (ctype(c) == icol_road_imperv) then
                albsn_improad(fl,1) = 0._r8
                albsn_improad(fl,2) = 0._r8
             else if (ctype(c) == icol_road_perv) then
                albsn_perroad(fl,1) = 0._r8
                albsn_perroad(fl,2) = 0._r8
             end if
          end if
       end do
    end do

  end subroutine UrbanSnowAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanRadiation
!
! !INTERFACE:
  subroutine UrbanRadiation (nc, lbl, ubl, lbc, ubc, lbp, ubp, &
                             num_nourbanl, filter_nourbanl, &
                             num_urbanl, filter_urbanl, &
                             num_urbanc, filter_urbanc, &
                             num_urbanp, filter_urbanp)
!
! !DESCRIPTION: 
! Solar fluxes absorbed and reflected by roof and canyon (walls, road).
! Also net and upward longwave fluxes.

! !USES:
    use clmtype
    use clm_varcon       , only : spval, icol_roof, icol_sunwall, icol_shadewall, &
                                  icol_road_perv, icol_road_imperv, sb
    use clm_varcon       , only : tfrz                ! To use new constant..
    use clm_time_manager , only : get_curr_date, get_step_size
    use clm_atmlnd       , only : clm_a2l
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nc                         ! clump index
    integer,  intent(in) :: lbl, ubl                   ! landunit-index bounds
    integer,  intent(in) :: lbc, ubc                   ! column-index bounds
    integer,  intent(in) :: lbp, ubp                   ! pft-index bounds
    integer , intent(in) :: num_nourbanl               ! number of non-urban landunits in clump
    integer , intent(in) :: filter_nourbanl(ubl-lbl+1) ! non-urban landunit filter
    integer , intent(in) :: num_urbanl                 ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(ubl-lbl+1)   ! urban landunit filter
    integer , intent(in) :: num_urbanc                 ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(ubc-lbc+1)   ! urban column filter
    integer , intent(in) :: num_urbanp                 ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(ubp-lbp+1)   ! urban pft filter
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 03/2003, Mariana Vertenstein: Migrated to clm2.2 
! 07/2004, Mariana Vertenstein: Migrated to clm3.0
! 01/2008, Erik Kluzek:         Migrated to clm3.5.15
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments (urban clump)
!
    real(r8), pointer :: canyon_hwr(:)           ! ratio of building height to street width
    real(r8), pointer :: wtroad_perv(:)          ! weight of pervious road wrt total road
    real(r8), pointer :: em_roof(:)              ! roof emissivity
    real(r8), pointer :: em_improad(:)           ! impervious road emissivity
    real(r8), pointer :: em_perroad(:)           ! pervious road emissivity
    real(r8), pointer :: em_wall(:)              ! wall emissivity
!
! local pointers to original implicit in arguments (clmtype)
!
    integer , pointer :: pgridcell(:)            ! gridcell of corresponding pft
    integer , pointer :: pcolumn(:)              ! column of corresponding pft
    integer , pointer :: lgridcell(:)            ! gridcell of corresponding landunit
    integer , pointer :: ctype(:)                ! column type
    integer , pointer :: coli(:)                 ! beginning column index for landunit 
    integer , pointer :: colf(:)                 ! ending column index for landunit
    integer , pointer :: pfti(:)                 ! beginning pfti index for landunit 
    integer , pointer :: pftf(:)                 ! ending pftf index for landunit
    real(r8), pointer :: londeg(:)               ! longitude (degrees)
    real(r8), pointer :: forc_lwrad(:)           ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: forc_solad(:,:)         ! direct beam radiation  (vis=forc_sols , nir=forc_soll ) (W/m**2)
    real(r8), pointer :: forc_solai(:,:)         ! diffuse beam radiation (vis=forc_sols , nir=forc_soll ) (W/m**2)
    real(r8), pointer :: forc_solar(:)           ! incident solar radiation (W/m**2)
    real(r8), pointer :: albd(:,:)               ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)               ! surface albedo (diffuse)
    real(r8), pointer :: t_grnd(:)               ! ground temperature (K)
    real(r8), pointer :: frac_sno(:)             ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: t_ref2m(:)              ! 2 m height surface air temperature (K)
    real(r8), pointer :: vf_sr(:)                ! view factor of sky for road
    real(r8), pointer :: vf_wr(:)                ! view factor of one wall for road
    real(r8), pointer :: vf_sw(:)                ! view factor of sky for one wall
    real(r8), pointer :: vf_rw(:)                ! view factor of road for one wall
    real(r8), pointer :: vf_ww(:)                ! view factor of opposing wall for one wall
    real(r8), pointer :: sabs_roof_dir(:,:)      ! direct  solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_roof_dif(:,:)      ! diffuse solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dir(:,:)   ! direct  solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dif(:,:)   ! diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dir(:,:) ! direct  solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dif(:,:) ! diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_improad_dir(:,:)   ! direct  solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_improad_dif(:,:)   ! diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dir(:,:)   ! direct  solar absorbed  by pervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dif(:,:)   ! diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
!
! local pointers to original implicit out arguments (clmtype)
!
    real(r8), pointer :: sabg(:)                 ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: sabv(:)                 ! solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: fsa(:)                  ! solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsa_u(:)                ! urban solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsr(:)                  ! solar radiation reflected (total) (W/m**2)
    real(r8), pointer :: fsds_vis_d(:)           ! incident direct beam vis solar radiation (W/m**2)
    real(r8), pointer :: fsds_nir_d(:)           ! incident direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsds_vis_i(:)           ! incident diffuse vis solar radiation (W/m**2)
    real(r8), pointer :: fsds_nir_i(:)           ! incident diffuse nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_vis_d(:)            ! reflected direct beam vis solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_d(:)            ! reflected direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_vis_i(:)            ! reflected diffuse vis solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_i(:)            ! reflected diffuse nir solar radiation (W/m**2)
    real(r8), pointer :: fsds_vis_d_ln(:)        ! incident direct beam vis solar rad at local noon (W/m**2)
    real(r8), pointer :: fsds_vis_i_ln(:)        ! incident diffuse beam vis solar rad at local noon (W/m**2)
    real(r8), pointer :: parveg_ln(:)            ! absorbed par by vegetation at local noon (W/m**2)
    real(r8), pointer :: fsds_nir_d_ln(:)        ! incident direct beam nir solar rad at local noon (W/m**2)
    real(r8), pointer :: fsr_vis_d_ln(:)         ! reflected direct beam vis solar rad at local noon (W/m**2)
    real(r8), pointer :: fsr_nir_d_ln(:)         ! reflected direct beam nir solar rad at local noon (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:)       ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)       ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_u(:)     ! urban net infrared (longwave) rad (W/m**2) [+ = to atm]
!
!
! !OTHER LOCAL VARIABLES
!EOP
!
    integer  :: fp,fl,p,c,l,g              ! indices
    integer  :: local_secp1                ! seconds into current date in local time
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: year,month,day             ! temporaries (not used)
    integer  :: secs                       ! seconds into current date

    real(r8), parameter :: mpe    = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: snoem  = 0.97_r8   ! snow emissivity (should use value from Biogeophysics1)
                                 
    real(r8) :: lwnet_roof(num_urbanl)     ! net (outgoing-incoming) longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwnet_improad(num_urbanl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwnet_perroad(num_urbanl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwnet_sunwall(num_urbanl)  ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwnet_shadewall(num_urbanl)! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwnet_canyon(num_urbanl)   ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: lwup_roof(num_urbanl)      ! upward longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwup_improad(num_urbanl)   ! upward longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwup_perroad(num_urbanl)   ! upward longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwup_sunwall(num_urbanl)   ! upward longwave radiation, (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwup_shadewall(num_urbanl) ! upward longwave radiation, (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwup_canyon(num_urbanl)    ! upward longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: t_roof(num_urbanl)         ! roof temperature (K)
    real(r8) :: t_improad(num_urbanl)      ! imppervious road temperature (K)
    real(r8) :: t_perroad(num_urbanl)      ! pervious road temperature (K)
    real(r8) :: t_sunwall(num_urbanl)      ! sunlit wall temperature (K)
    real(r8) :: t_shadewall(num_urbanl)    ! shaded wall temperature (K)
    real(r8) :: lwdown(num_urbanl)         ! atmospheric downward longwave radiation (W/m**2)
    real(r8) :: em_roof_s(num_urbanl)      ! roof emissivity with snow effects
    real(r8) :: em_improad_s(num_urbanl)   ! impervious road emissivity with snow effects
    real(r8) :: em_perroad_s(num_urbanl)   ! pervious road emissivity with snow effects
!-----------------------------------------------------------------------

    ! Assign pointers into module urban clumps

    if( num_urbanl > 0 )then
       canyon_hwr         => urban_clump(nc)%canyon_hwr
       wtroad_perv        => urban_clump(nc)%wtroad_perv
       em_roof            => urban_clump(nc)%em_roof
       em_improad         => urban_clump(nc)%em_improad
       em_perroad         => urban_clump(nc)%em_perroad
       em_wall            => urban_clump(nc)%em_wall
    end if

    ! Assign local pointers to multi-level derived type members (gridcell level)
    
    londeg        => clm3%g%londeg
    forc_solad    => clm_a2l%forc_solad
    forc_solai    => clm_a2l%forc_solai
    forc_solar    => clm_a2l%forc_solar
    forc_lwrad    => clm_a2l%forc_lwrad

    ! Assign local pointers to derived type members (landunit level)

    pfti           => clm3%g%l%pfti
    pftf           => clm3%g%l%pftf
    coli           => clm3%g%l%coli
    colf           => clm3%g%l%colf
    lgridcell      => clm3%g%l%gridcell
    vf_sr          => clm3%g%l%lps%vf_sr
    vf_wr          => clm3%g%l%lps%vf_wr
    vf_sw          => clm3%g%l%lps%vf_sw
    vf_rw          => clm3%g%l%lps%vf_rw
    vf_ww          => clm3%g%l%lps%vf_ww
    sabs_roof_dir      => clm3%g%l%lps%sabs_roof_dir
    sabs_roof_dif      => clm3%g%l%lps%sabs_roof_dif
    sabs_sunwall_dir   => clm3%g%l%lps%sabs_sunwall_dir
    sabs_sunwall_dif   => clm3%g%l%lps%sabs_sunwall_dif
    sabs_shadewall_dir => clm3%g%l%lps%sabs_shadewall_dir
    sabs_shadewall_dif => clm3%g%l%lps%sabs_shadewall_dif
    sabs_improad_dir   => clm3%g%l%lps%sabs_improad_dir
    sabs_improad_dif   => clm3%g%l%lps%sabs_improad_dif
    sabs_perroad_dir   => clm3%g%l%lps%sabs_perroad_dir
    sabs_perroad_dif   => clm3%g%l%lps%sabs_perroad_dif

    ! Assign local pointers to derived type members (column level)

    ctype          => clm3%g%l%c%itype
    t_grnd         => clm3%g%l%c%ces%t_grnd 
    frac_sno       => clm3%g%l%c%cps%frac_sno

    ! Assign local pointers to derived type members (pft level)

    pgridcell      => clm3%g%l%c%p%gridcell
    pcolumn        => clm3%g%l%c%p%column
    albd           => clm3%g%l%c%p%pps%albd
    albi           => clm3%g%l%c%p%pps%albi
    sabg           => clm3%g%l%c%p%pef%sabg
    sabv           => clm3%g%l%c%p%pef%sabv
    fsa            => clm3%g%l%c%p%pef%fsa
    fsa_u          => clm3%g%l%c%p%pef%fsa_u
    fsr            => clm3%g%l%c%p%pef%fsr
    fsds_vis_d     => clm3%g%l%c%p%pef%fsds_vis_d
    fsds_nir_d     => clm3%g%l%c%p%pef%fsds_nir_d
    fsds_vis_i     => clm3%g%l%c%p%pef%fsds_vis_i
    fsds_nir_i     => clm3%g%l%c%p%pef%fsds_nir_i
    fsr_vis_d      => clm3%g%l%c%p%pef%fsr_vis_d
    fsr_nir_d      => clm3%g%l%c%p%pef%fsr_nir_d
    fsr_vis_i      => clm3%g%l%c%p%pef%fsr_vis_i
    fsr_nir_i      => clm3%g%l%c%p%pef%fsr_nir_i
    fsds_vis_d_ln  => clm3%g%l%c%p%pef%fsds_vis_d_ln
    fsds_nir_d_ln  => clm3%g%l%c%p%pef%fsds_nir_d_ln
    fsds_vis_i_ln  => clm3%g%l%c%p%pef%fsds_vis_i_ln
    parveg_ln      => clm3%g%l%c%p%pef%parveg_ln
    fsr_vis_d_ln   => clm3%g%l%c%p%pef%fsr_vis_d_ln
    fsr_nir_d_ln   => clm3%g%l%c%p%pef%fsr_nir_d_ln
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out 
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net
    eflx_lwrad_net_u => clm3%g%l%c%p%pef%eflx_lwrad_net_u
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m

    ! Define fields that appear on the restart file for non-urban landunits 

    do fl = 1,num_nourbanl
       l = filter_nourbanl(fl)
       sabs_roof_dir(l,:) = spval
       sabs_roof_dif(l,:) = spval
       sabs_sunwall_dir(l,:) = spval
       sabs_sunwall_dif(l,:) = spval
       sabs_shadewall_dir(l,:) = spval
       sabs_shadewall_dif(l,:) = spval
       sabs_improad_dir(l,:) = spval
       sabs_improad_dif(l,:) = spval
       sabs_perroad_dir(l,:) = spval
       sabs_perroad_dif(l,:) = spval
       vf_sr(l) = spval
       vf_wr(l) = spval
       vf_sw(l) = spval
       vf_rw(l) = spval
       vf_ww(l) = spval
    end do

    ! Set input forcing fields
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l) 

       ! Need to set the following temperatures to some defined value even if it
       ! does not appear in the urban landunit for the net_longwave computation

       t_roof(fl)      = 19._r8 + tfrz
       t_sunwall(fl)   = 19._r8 + tfrz
       t_shadewall(fl) = 19._r8 + tfrz
       t_improad(fl)   = 19._r8 + tfrz
       t_perroad(fl)   = 19._r8 + tfrz

       ! Initial assignment of emissivity
       em_roof_s(fl) = em_roof(fl)
       em_improad_s(fl) = em_improad(fl)
       em_perroad_s(fl) = em_perroad(fl)

       ! Set urban temperatures and emissivity including snow effects.
       do c = coli(l),colf(l)
          if (ctype(c) == icol_roof       )  then
             t_roof(fl)      = t_grnd(c)
             em_roof_s(fl) = em_roof(fl)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_road_imperv) then 
             t_improad(fl)   = t_grnd(c)
             em_improad_s(fl) = em_improad(fl)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_road_perv  ) then
             t_perroad(fl)   = t_grnd(c)
             em_perroad_s(fl) = em_perroad(fl)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_sunwall    ) then
             t_sunwall(fl)   = t_grnd(c)
          else if (ctype(c) == icol_shadewall  ) then
             t_shadewall(fl) = t_grnd(c)
          end if
       end do
       lwdown(fl) = forc_lwrad(g)
    end do

    ! Net longwave radiation for road and both walls in urban canyon allowing for multiple re-emission
    
    if (num_urbanl .gt. 0) then
       call net_longwave (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr, wtroad_perv, &
                          lwdown, em_roof_s, em_improad_s, em_perroad_s, em_wall, &
                          t_roof, t_improad, t_perroad, t_sunwall, t_shadewall, &
                          lwnet_roof, lwnet_improad, lwnet_perroad, lwnet_sunwall, lwnet_shadewall, lwnet_canyon, &
                          lwup_roof, lwup_improad, lwup_perroad, lwup_sunwall, lwup_shadewall, lwup_canyon)
    end if
    
    dtime = get_step_size()
    call get_curr_date (year, month, day, secs)
    
    ! Determine clmtype variables needed for history output and communication with atm
    ! Loop over urban pfts in clump

    do fp = 1,num_urbanp
       p = filter_urbanp(fp)
       g = pgridcell(p)
       
       local_secp1 = secs + nint((londeg(g)/degpsec)/dtime)*dtime
       local_secp1 = mod(local_secp1,isecspday)

       ! Solar incident 
       
       fsds_vis_d(p) = forc_solad(g,1)
       fsds_nir_d(p) = forc_solad(g,2)    
       fsds_vis_i(p) = forc_solai(g,1)
       fsds_nir_i(p) = forc_solai(g,2)
       ! Determine local noon incident solar
       if (local_secp1 == noonsec) then
          fsds_vis_d_ln(p) = forc_solad(g,1)
          fsds_nir_d_ln(p) = forc_solad(g,2)
          fsds_vis_i_ln(p) = forc_solai(g,1)
          parveg_ln(p)     = 0._r8
       else
          fsds_vis_d_ln(p) = spval 
          fsds_nir_d_ln(p) = spval 
          fsds_vis_i_ln(p) = spval
          parveg_ln(p)     = spval
       endif
       
       ! Solar reflected 
       ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)

       fsr_vis_d(p) = albd(p,1) * forc_solad(g,1)
       fsr_nir_d(p) = albd(p,2) * forc_solad(g,2)
       fsr_vis_i(p) = albi(p,1) * forc_solai(g,1)
       fsr_nir_i(p) = albi(p,2) * forc_solai(g,2)

       ! Determine local noon reflected solar
       if (local_secp1 == noonsec) then
          fsr_vis_d_ln(p) = fsr_vis_d(p)
          fsr_nir_d_ln(p) = fsr_nir_d(p)
       else
          fsr_vis_d_ln(p) = spval 
          fsr_nir_d_ln(p) = spval 
       endif
       fsr(p) = fsr_vis_d(p) + fsr_nir_d(p) + fsr_vis_i(p) + fsr_nir_i(p)  
       
    end do

    ! Loop over urban landunits in clump

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)

       ! Solar absorbed and longwave out and net
       ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
       ! Each urban pft has its own column - this is used in the logic below
 
       do p = pfti(l), pftf(l)
          c = pcolumn(p)
          if (ctype(c) == icol_roof) then   
             eflx_lwrad_out(p) = lwup_roof(fl)
             eflx_lwrad_net(p) = lwnet_roof(fl)
             eflx_lwrad_net_u(p) = lwnet_roof(fl)
             sabg(p) = sabs_roof_dir(l,1)*forc_solad(g,1) + &
                       sabs_roof_dif(l,1)*forc_solai(g,1) + &
                       sabs_roof_dir(l,2)*forc_solad(g,2) + &
                       sabs_roof_dif(l,2)*forc_solai(g,2) 
          else if (ctype(c) == icol_sunwall) then   
             eflx_lwrad_out(p) = lwup_sunwall(fl)
             eflx_lwrad_net(p) = lwnet_sunwall(fl)
             eflx_lwrad_net_u(p) = lwnet_sunwall(fl)
             sabg(p) = sabs_sunwall_dir(l,1)*forc_solad(g,1) + &
                       sabs_sunwall_dif(l,1)*forc_solai(g,1) + &
                       sabs_sunwall_dir(l,2)*forc_solad(g,2) + &
                       sabs_sunwall_dif(l,2)*forc_solai(g,2) 
          else if (ctype(c) == icol_shadewall) then   
             eflx_lwrad_out(p) = lwup_shadewall(fl)
             eflx_lwrad_net(p) = lwnet_shadewall(fl)
             eflx_lwrad_net_u(p) = lwnet_shadewall(fl)
             sabg(p) = sabs_shadewall_dir(l,1)*forc_solad(g,1) + &
                       sabs_shadewall_dif(l,1)*forc_solai(g,1) + &
                       sabs_shadewall_dir(l,2)*forc_solad(g,2) + &
                       sabs_shadewall_dif(l,2)*forc_solai(g,2) 
          else if (ctype(c) == icol_road_perv) then       
             eflx_lwrad_out(p) = lwup_perroad(fl)
             eflx_lwrad_net(p) = lwnet_perroad(fl)
             eflx_lwrad_net_u(p) = lwnet_perroad(fl)
             sabg(p) = sabs_perroad_dir(l,1)*forc_solad(g,1) + &
                       sabs_perroad_dif(l,1)*forc_solai(g,1) + &
                       sabs_perroad_dir(l,2)*forc_solad(g,2) + &
                       sabs_perroad_dif(l,2)*forc_solai(g,2) 
          else if (ctype(c) == icol_road_imperv) then       
             eflx_lwrad_out(p) = lwup_improad(fl)
             eflx_lwrad_net(p) = lwnet_improad(fl)
             eflx_lwrad_net_u(p) = lwnet_improad(fl)
             sabg(p) = sabs_improad_dir(l,1)*forc_solad(g,1) + &
                       sabs_improad_dif(l,1)*forc_solai(g,1) + &
                       sabs_improad_dir(l,2)*forc_solad(g,2) + &
                       sabs_improad_dif(l,2)*forc_solai(g,2) 
          end if
          sabv(p)   = 0._r8
          fsa(p)    = sabv(p) + sabg(p)
          fsa_u(p)  = fsa(p)
          
       end do   ! end loop over urban pfts
       
    end do ! end loop over urban landunits

  end subroutine UrbanRadiation

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: view_factor
!
! !INTERFACE:
  subroutine view_factor (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr)

!
! !DESCRIPTION: 
! View factors for road and one wall
!                                                        WALL    |
!                  ROAD                                          |
!                                                         wall   |
!          -----\          /-----   -             -  |\----------/
!              | \  vsr   / |       |         r   |  | \  vww   /   s
!              |  \      /  |       h         o   w  |  \      /    k
!        wall  |   \    /   | wall  |         a   |  |   \    /     y
!              |vwr \  / vwr|       |         d   |  |vrw \  / vsw 
!              ------\/------       -             -  |-----\/-----
!                   road                                  wall   |
!              <----- w ---->                                    |
!                                                    <---- h --->|
!
!    vsr = view factor of sky for road          vrw = view factor of road for wall
!    vwr = view factor of one wall for road     vww = view factor of opposing wall for wall
!                                               vsw = view factor of sky for wall
!    vsr + vwr + vwr = 1                        vrw + vww + vsw = 1
!
! Source: Masson, V. (2000) A physically-based scheme for the urban energy budget in 
! atmospheric models. Boundary-Layer Meteorology 94:357-397
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer,  intent(in) :: lbl, ubl                  ! landunit-index bounds
    integer , intent(in)  :: num_urbanl               ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(ubl-lbl+1) ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr(num_urbanl)   ! ratio of building height to street width
!
! local pointers to original implicit out arguments (clmtype)
!
    real(r8), pointer :: vf_sr(:)     ! view factor of sky for road
    real(r8), pointer :: vf_wr(:)     ! view factor of one wall for road
    real(r8), pointer :: vf_sw(:)     ! view factor of sky for one wall
    real(r8), pointer :: vf_rw(:)     ! view factor of road for one wall
    real(r8), pointer :: vf_ww(:)     ! view factor of opposing wall for one wall
!
! !CALLED FROM:
! subroutine UrbanAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 03/2003, Mariana Vertenstein: Migrated to clm2.2 
! 01/2008, Erik Kluzek:         Migrated to clm3.5.15
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: l, fl   ! indices
    real(r8) :: sum    ! sum of view factors for wall or road
!-----------------------------------------------------------------------

    ! Assign landunit level pointer

    vf_sr  => clm3%g%l%lps%vf_sr
    vf_wr  => clm3%g%l%lps%vf_wr
    vf_sw  => clm3%g%l%lps%vf_sw
    vf_rw  => clm3%g%l%lps%vf_rw
    vf_ww  => clm3%g%l%lps%vf_ww

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       ! road -- sky view factor -> 1 as building height -> 0 
       ! and -> 0 as building height -> infinity

       vf_sr(l) = sqrt(canyon_hwr(fl)**2 + 1._r8) - canyon_hwr(fl)
       vf_wr(l) = 0.5_r8 * (1._r8 - vf_sr(l))

       ! one wall -- sky view factor -> 0.5 as building height -> 0 
       ! and -> 0 as building height -> infinity

       vf_sw(l) = 0.5_r8 * (canyon_hwr(fl) + 1._r8 - sqrt(canyon_hwr(fl)**2+1._r8)) / canyon_hwr(fl)
       vf_rw(l) = vf_sw(l)
       vf_ww(l) = 1._r8 - vf_sw(l) - vf_rw(l)

    end do


    ! error check -- make sure view factor sums to one for road and wall

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       sum = vf_sr(l) + 2._r8*vf_wr(l)
       if (abs(sum-1._r8) > 1.e-06_r8 ) then
          write (iulog,*) 'urban road view factor error',sum
          write (iulog,*) 'clm model is stopping'
          call endrun()
       endif
       sum = vf_sw(l) + vf_rw(l) + vf_ww(l)
       if (abs(sum-1._r8) > 1.e-06_r8 ) then
          write (iulog,*) 'urban wall view factor error',sum
          write (iulog,*) 'clm model is stopping'
          call endrun()
       endif

    end do

  end subroutine view_factor

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: incident_direct
!
! !INTERFACE:
  subroutine incident_direct (lbl, ubl, num_urbanl, canyon_hwr, coszen, zen, sdir, sdir_road, sdir_sunwall, sdir_shadewall)
!
! !DESCRIPTION: 
! Direct beam solar radiation incident on walls and road in urban canyon
!
!                           Sun
!                            /
!             roof          /
!            ------        /---            -
!                 |       / |              |
!    sunlit wall  |      /  | shaded wall  h
!                 |     /   |              |
!                 -----/-----              -
!                    road
!                 <--- w --->
!
! Method:
! Road          = Horizontal surface. Account for shading by wall. Integrate over all canyon orientations
! Wall (sunlit) = Adjust horizontal radiation for 90 degree surface. Account for shading by opposing wall.
!                 Integrate over all canyon orientations
! Wall (shaded) = 0
!
! Conservation check: Total incoming direct beam (sdir) = sdir_road + (sdir_shadewall + sdir_sunwall)*canyon_hwr
! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
!
! Source: Masson, V. (2000) A physically-based scheme for the urban energy budget in 
! atmospheric models. Boundary-Layer Meteorology 94:357-397
!
! This analytical solution from Masson (2000) agrees with the numerical solution to
! within 0.6 W/m**2 for sdir = 1000 W/m**2 and for all H/W from 0.1 to 10 by 0.1
! and all solar zenith angles from 1 to 90 deg by 1
!
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use clm_varcon    , only : rpi
    implicit none
!
! !ARGUMENTS:
    integer,  intent(in)  :: lbl, ubl                           ! landunit-index bounds
    integer , intent(in)  :: num_urbanl                         ! number of urban landunits
    real(r8), intent(in)  :: canyon_hwr(num_urbanl)             ! ratio of building height to street width
    real(r8), intent(in)  :: coszen(num_urbanl)                 ! cosine solar zenith angle
    real(r8), intent(in)  :: zen(num_urbanl)                    ! solar zenith angle (radians)
    real(r8), intent(in)  :: sdir(num_urbanl, numrad)           ! direct beam solar radiation incident on horizontal surface
    real(r8), intent(out) :: sdir_road(num_urbanl, numrad)      ! direct beam solar radiation incident on road per unit incident flux
    real(r8), intent(out) :: sdir_sunwall(num_urbanl, numrad)   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8), intent(out) :: sdir_shadewall(num_urbanl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
!
! !CALLED FROM:
! subroutine UrbanAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: l,i,ib             ! indices
!KO   logical  :: numchk = .true.     ! true => perform numerical check of analytical solution
    logical  :: numchk = .false.   ! true => perform numerical check of analytical solution
    real(r8) :: theta0(num_urbanl) ! critical canyon orientation for which road is no longer illuminated
    real(r8) :: tanzen(num_urbanl) ! tan(zenith angle)
    real(r8) :: swall_projected    ! direct beam solar radiation (per unit ground area) incident on wall
    real(r8) :: err1(num_urbanl)   ! energy conservation error
    real(r8) :: err2(num_urbanl)   ! energy conservation error
    real(r8) :: err3(num_urbanl)   ! energy conservation error
    real(r8) :: sumr               ! sum of sroad for each orientation (0 <= theta <= pi/2)
    real(r8) :: sumw               ! sum of swall for each orientation (0 <= theta <= pi/2)
    real(r8) :: num                ! number of orientations
    real(r8) :: theta              ! canyon orientation relative to sun (0 <= theta <= pi/2)
    real(r8) :: zen0               ! critical solar zenith angle for which sun begins to illuminate road
!-----------------------------------------------------------------------

    do l = 1,num_urbanl
       if (coszen(l) > 0._r8) then
          theta0(l) = asin(min( (1._r8/(canyon_hwr(l)*tan(max(zen(l),0.000001_r8)))), 1._r8 ))
          tanzen(l) = tan(zen(l))
       end if
    end do

    do ib = 1,numrad

       do l = 1,num_urbanl
          if (coszen(l) > 0._r8) then
             sdir_shadewall(l,ib) = 0._r8

             ! incident solar radiation on wall and road integrated over all canyon orientations (0 <= theta <= pi/2)
          
             sdir_road(l,ib) = sdir(l,ib) * &
                  (2._r8*theta0(l)/rpi - 2./rpi*canyon_hwr(l)*tanzen(l)*(1._r8-cos(theta0(l))))
             sdir_sunwall(l,ib) = 2._r8 * sdir(l,ib) * ((1._r8/canyon_hwr(l))* &
                  (0.5_r8-theta0(l)/rpi) + (1._r8/rpi)*tanzen(l)*(1._r8-cos(theta0(l))))

             ! conservation check for road and wall. need to use wall fluxes converted to ground area
          
             swall_projected = (sdir_shadewall(l,ib) + sdir_sunwall(l,ib)) * canyon_hwr(l)
             err1(l) = sdir(l,ib) - (sdir_road(l,ib) + swall_projected)
          else
             sdir_road(l,ib) = 0._r8
             sdir_sunwall(l,ib) = 0._r8
             sdir_shadewall(l,ib) = 0._r8
          endif
       end do

       do l = 1,num_urbanl
          if (coszen(l) > 0._r8) then
             if (abs(err1(l)) > 0.001_r8) then
                write (iulog,*) 'urban direct beam solar radiation balance error',err1(l)
                write (iulog,*) 'clm model is stopping'
                call endrun()
             endif
          endif
       end do

       ! numerical check of analytical solution
       ! sum sroad and swall over all canyon orientations (0 <= theta <= pi/2)

       if (numchk) then
          do l = 1,num_urbanl
             if (coszen(l) > 0._r8) then
             sumr = 0._r8
             sumw = 0._r8
             num  = 0._r8
             do i = 1, 9000
                theta = i/100._r8 * rpi/180._r8
                zen0 = atan(1._r8/(canyon_hwr(l)*sin(theta)))
                if (zen(l) >= zen0) then 
                   sumr = sumr + 0._r8
                   sumw = sumw + sdir(l,ib) / canyon_hwr(l)
                else
                   sumr = sumr + sdir(l,ib) * (1._r8-canyon_hwr(l)*sin(theta)*tanzen(l))
                   sumw = sumw + sdir(l,ib) * sin(theta)*tanzen(l)
                end if
                num = num + 1._r8
             end do
             err2(l) = sumr/num - sdir_road(l,ib)
             err3(l) = sumw/num - sdir_sunwall(l,ib)
             endif
          end do
          do l = 1,num_urbanl
             if (coszen(l) > 0._r8) then
             if (abs(err2(l)) > 0.0006_r8 ) then
                write (iulog,*) 'urban road incident direct beam solar radiation error',err2(l)
                write (iulog,*) 'clm model is stopping'
                call endrun
             endif
             if (abs(err3(l)) > 0.0006_r8 ) then
                write (iulog,*) 'urban wall incident direct beam solar radiation error',err3(l)
                write (iulog,*) 'clm model is stopping'
                call endrun
             end if
             end if
          end do
       end if

    end do

  end subroutine incident_direct

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: incident_diffuse
!
! !INTERFACE:
  subroutine incident_diffuse (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr, sdif, sdif_road, sdif_sunwall, sdif_shadewall)
!
! !DESCRIPTION: 
! Diffuse solar radiation incident on walls and road in urban canyon
! Conservation check: Total incoming diffuse 
! (sdif) = sdif_road + (sdif_shadewall + sdif_sunwall)*canyon_hwr
! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
!
! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: lbl, ubl                           ! landunit-index bounds
    integer , intent(in)  :: num_urbanl                         ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(ubl-lbl+1)           ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr(num_urbanl)             ! ratio of building height to street width
    real(r8), intent(in)  :: sdif(num_urbanl, numrad)           ! diffuse solar radiation incident on horizontal surface
    real(r8), intent(out) :: sdif_road(num_urbanl, numrad)      ! diffuse solar radiation incident on road
    real(r8), intent(out) :: sdif_sunwall(num_urbanl, numrad)   ! diffuse solar radiation (per unit wall area) incident on sunlit wall
    real(r8), intent(out) :: sdif_shadewall(num_urbanl, numrad) ! diffuse solar radiation (per unit wall area) incident on shaded wall
!
! local pointers to original implicit in arguments (clmtype)
!
    real(r8), pointer :: vf_sr(:)     ! view factor of sky for road
    real(r8), pointer :: vf_sw(:)     ! view factor of sky for one wall
!
! !CALLED FROM:
! subroutine UrbanAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: l, fl, ib       ! indices      
    real(r8) :: err(num_urbanl) ! energy conservation error (W/m**2)
    real(r8) :: swall_projected ! diffuse solar radiation (per unit ground area) incident on wall (W/m**2)
!-----------------------------------------------------------------------

    ! Assign landunit level pointer

    vf_sr              => clm3%g%l%lps%vf_sr
    vf_sw              => clm3%g%l%lps%vf_sw

    do ib = 1, numrad

       ! diffuse solar and conservation check. need to convert wall fluxes to ground area

       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          sdif_road(fl,ib)      = sdif(fl,ib) * vf_sr(l)
          sdif_sunwall(fl,ib)   = sdif(fl,ib) * vf_sw(l) 
          sdif_shadewall(fl,ib) = sdif(fl,ib) * vf_sw(l) 

          swall_projected = (sdif_shadewall(fl,ib) + sdif_sunwall(fl,ib)) * canyon_hwr(fl)
          err(fl) = sdif(fl,ib) - (sdif_road(fl,ib) + swall_projected)
       end do

       ! error check

       do l = 1, num_urbanl
          if (abs(err(l)) > 0.001_r8) then
             write (iulog,*) 'urban diffuse solar radiation balance error',err(l) 
             write (iulog,*) 'clm model is stopping'
             call endrun
          endif
       end do

    end do

  end subroutine incident_diffuse

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: net_solar
!
! !INTERFACE:
  subroutine net_solar (lbl, ubl, num_urbanl, filter_urbanl, coszen, canyon_hwr, wtroad_perv, sdir, sdif, &
                        alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir, &
                        alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif, &
                        sdir_road, sdir_sunwall, sdir_shadewall,  &
                        sdif_road, sdif_sunwall, sdif_shadewall,  &
                        sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir, sref_roof_dir, &
                        sref_improad_dif, sref_perroad_dif, sref_sunwall_dif, sref_shadewall_dif, sref_roof_dif)
!
! !DESCRIPTION: 
! Solar radiation absorbed by road and both walls in urban canyon allowing 
! for multiple reflection.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
!
!
! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: lbl, ubl                               ! landunit-index bounds
    integer , intent(in)  :: num_urbanl                             ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(ubl-lbl+1)               ! urban landunit filter
    real(r8), intent(in)  :: coszen(num_urbanl)                     ! cosine solar zenith angle
    real(r8), intent(in)  :: canyon_hwr(num_urbanl)                 ! ratio of building height to street width
    real(r8), intent(in)  :: wtroad_perv(num_urbanl)                ! weight of pervious road wrt total road
    real(r8), intent(in)  :: sdir(num_urbanl, numrad)               ! direct beam solar radiation incident on horizontal surface
    real(r8), intent(in)  :: sdif(num_urbanl, numrad)               ! diffuse solar radiation on horizontal surface
    real(r8), intent(in)  :: alb_improad_dir(num_urbanl, numrad)    ! direct impervious road albedo
    real(r8), intent(in)  :: alb_perroad_dir(num_urbanl, numrad)    ! direct pervious road albedo
    real(r8), intent(in)  :: alb_wall_dir(num_urbanl, numrad)       ! direct  wall albedo
    real(r8), intent(in)  :: alb_roof_dir(num_urbanl, numrad)       ! direct  roof albedo
    real(r8), intent(in)  :: alb_improad_dif(num_urbanl, numrad)    ! diffuse impervious road albedo
    real(r8), intent(in)  :: alb_perroad_dif(num_urbanl, numrad)    ! diffuse pervious road albedo
    real(r8), intent(in)  :: alb_wall_dif(num_urbanl, numrad)       ! diffuse wall albedo
    real(r8), intent(in)  :: alb_roof_dif(num_urbanl, numrad)       ! diffuse roof albedo
    real(r8), intent(in)  :: sdir_road(num_urbanl, numrad)          ! direct beam solar radiation incident on road per unit incident flux
    real(r8), intent(in)  :: sdir_sunwall(num_urbanl, numrad)       ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8), intent(in)  :: sdir_shadewall(num_urbanl, numrad)     ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8), intent(in)  :: sdif_road(num_urbanl, numrad)          ! diffuse solar radiation incident on road per unit incident flux
    real(r8), intent(in)  :: sdif_sunwall(num_urbanl, numrad)       ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8), intent(in)  :: sdif_shadewall(num_urbanl, numrad)     ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8), intent(inout) :: sref_improad_dir(num_urbanl, numrad)   ! direct  solar rad reflected by impervious road (per unit ground area) per unit incident flux
    real(r8), intent(inout) :: sref_perroad_dir(num_urbanl, numrad)   ! direct  solar rad reflected by pervious road (per unit ground area) per unit incident flux
    real(r8), intent(inout) :: sref_improad_dif(num_urbanl, numrad)   ! diffuse solar rad reflected by impervious road (per unit ground area) per unit incident flux
    real(r8), intent(inout) :: sref_perroad_dif(num_urbanl, numrad)   ! diffuse solar rad reflected by pervious road (per unit ground area) per unit incident flux
    real(r8), intent(inout) :: sref_sunwall_dir(num_urbanl, numrad)   ! direct solar  rad reflected by sunwall (per unit wall area) per unit incident flux
    real(r8), intent(inout) :: sref_sunwall_dif(num_urbanl, numrad)   ! diffuse solar rad reflected by sunwall (per unit wall area) per unit incident flux
    real(r8), intent(inout) :: sref_shadewall_dir(num_urbanl, numrad) ! direct solar  rad reflected by shadewall (per unit wall area) per unit incident flux
    real(r8), intent(inout) :: sref_shadewall_dif(num_urbanl, numrad) ! diffuse solar rad reflected by shadewall (per unit wall area) per unit incident flux
    real(r8), intent(inout) :: sref_roof_dir(num_urbanl, numrad)      ! direct  solar rad reflected by roof (per unit ground area) per unit incident flux
    real(r8), intent(inout) :: sref_roof_dif(num_urbanl, numrad)      ! diffuse solar rad reflected by roof (per unit ground area)  per unit incident flux
!
! local pointers to original implicit in arguments (clmtype)
!
    real(r8), pointer :: vf_sr(:)     ! view factor of sky for road
    real(r8), pointer :: vf_wr(:)     ! view factor of one wall for road
    real(r8), pointer :: vf_sw(:)     ! view factor of sky for one wall
    real(r8), pointer :: vf_rw(:)     ! view factor of road for one wall
    real(r8), pointer :: vf_ww(:)     ! view factor of opposing wall for one wall
    real(r8), pointer :: sabs_roof_dir(:,:)      ! direct  solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_roof_dif(:,:)      ! diffuse solar absorbed  by roof per unit ground area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dir(:,:)   ! direct  solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_sunwall_dif(:,:)   ! diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dir(:,:) ! direct  solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_shadewall_dif(:,:) ! diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
    real(r8), pointer :: sabs_improad_dir(:,:)   ! direct  solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_improad_dif(:,:)   ! diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dir(:,:)   ! direct  solar absorbed  by pervious road per unit ground area per unit incident flux
    real(r8), pointer :: sabs_perroad_dif(:,:)   ! diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
!
! !CALLED FROM:
! subroutine UrbanAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES
!EOP
!
    real(r8) :: wtroad_imperv(num_urbanl)           ! weight of impervious road wrt total road
    real(r8) :: sabs_canyon_dir(num_urbanl)         ! direct solar rad absorbed by canyon per unit incident flux
    real(r8) :: sabs_canyon_dif(num_urbanl)         ! diffuse solar rad absorbed by canyon per unit incident flux
    real(r8) :: sref_canyon_dir(num_urbanl)         ! direct solar reflected by canyon per unit incident flux 
    real(r8) :: sref_canyon_dif(num_urbanl)         ! diffuse solar reflected by canyon per unit incident flux

    real(r8) :: improad_a_dir(num_urbanl)           ! absorbed direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_a_dif(num_urbanl)           ! absorbed diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dir(num_urbanl)           ! reflected direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dif(num_urbanl)           ! reflected diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_sky_dir(num_urbanl)       ! improad_r_dir to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dir(num_urbanl)   ! improad_r_dir to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dir(num_urbanl) ! improad_r_dir to shaded wall per unit incident flux
    real(r8) :: improad_r_sky_dif(num_urbanl)       ! improad_r_dif to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dif(num_urbanl)   ! improad_r_dif to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dif(num_urbanl) ! improad_r_dif to shaded wall per unit incident flux

    real(r8) :: perroad_a_dir(num_urbanl)           ! absorbed direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_a_dif(num_urbanl)           ! absorbed diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dir(num_urbanl)           ! reflected direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dif(num_urbanl)           ! reflected diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_sky_dir(num_urbanl)       ! perroad_r_dir to sky per unit incident flux 
    real(r8) :: perroad_r_sunwall_dir(num_urbanl)   ! perroad_r_dir to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dir(num_urbanl) ! perroad_r_dir to shaded wall per unit incident flux
    real(r8) :: perroad_r_sky_dif(num_urbanl)       ! perroad_r_dif to sky per unit incident flux
    real(r8) :: perroad_r_sunwall_dif(num_urbanl)   ! perroad_r_dif to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dif(num_urbanl) ! perroad_r_dif to shaded wall per unit incident flux

    real(r8) :: road_a_dir(num_urbanl)              ! absorbed direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_a_dif(num_urbanl)              ! absorbed diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dir(num_urbanl)              ! reflected direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dif(num_urbanl)              ! reflected diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_sky_dir(num_urbanl)          ! road_r_dir to sky per unit incident flux
    real(r8) :: road_r_sunwall_dir(num_urbanl)      ! road_r_dir to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dir(num_urbanl)    ! road_r_dir to shaded wall per unit incident flux
    real(r8) :: road_r_sky_dif(num_urbanl)          ! road_r_dif to sky per unit incident flux
    real(r8) :: road_r_sunwall_dif(num_urbanl)      ! road_r_dif to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dif(num_urbanl)    ! road_r_dif to shaded wall per unit incident flux

    real(r8) :: sunwall_a_dir(num_urbanl)           ! absorbed direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_a_dif(num_urbanl)           ! absorbed diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dir(num_urbanl)           ! reflected direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dif(num_urbanl)           ! reflected diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_sky_dir(num_urbanl)       ! sunwall_r_dir to sky per unit incident flux
    real(r8) :: sunwall_r_road_dir(num_urbanl)      ! sunwall_r_dir to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dir(num_urbanl) ! sunwall_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: sunwall_r_sky_dif(num_urbanl)       ! sunwall_r_dif to sky per unit incident flux
    real(r8) :: sunwall_r_road_dif(num_urbanl)      ! sunwall_r_dif to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dif(num_urbanl) ! sunwall_r_dif to opposing (shaded) wall per unit incident flux

    real(r8) :: shadewall_a_dir(num_urbanl)         ! absorbed direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_a_dif(num_urbanl)         ! absorbed diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dir(num_urbanl)         ! reflected direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dif(num_urbanl)         ! reflected diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_sky_dir(num_urbanl)     ! shadewall_r_dir to sky per unit incident flux
    real(r8) :: shadewall_r_road_dir(num_urbanl)    ! shadewall_r_dir to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dir(num_urbanl) ! shadewall_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: shadewall_r_sky_dif(num_urbanl)     ! shadewall_r_dif to sky per unit incident flux
    real(r8) :: shadewall_r_road_dif(num_urbanl)    ! shadewall_r_dif to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dif(num_urbanl) ! shadewall_r_dif to opposing (sunlit) wall per unit incident flux

    real(r8) :: canyon_alb_dir(num_urbanl)          ! direct canyon albedo
    real(r8) :: canyon_alb_dif(num_urbanl)          ! diffuse canyon albedo

    real(r8) :: stot(num_urbanl)                    ! sum of radiative terms
    real(r8) :: stot_dir(num_urbanl)                ! sum of direct radiative terms
    real(r8) :: stot_dif(num_urbanl)                ! sum of diffuse radiative terms

    integer  :: l,fl,ib                             ! indices
    integer  :: iter_dir,iter_dif                   ! iteration counter
    real(r8) :: crit                                ! convergence criterion
    real(r8) :: err                                 ! energy conservation error
    integer  :: pass
    integer, parameter :: n = 50                    ! number of interations
    real(r8) :: sabs_road                           ! temporary for absorption over road
    real(r8) :: sref_road                           ! temporary for reflected over road
    real(r8), parameter :: errcrit  = .00001_r8     ! error criteria
!-----------------------------------------------------------------------

    ! Assign landunit level pointer

    vf_sr              => clm3%g%l%lps%vf_sr
    vf_wr              => clm3%g%l%lps%vf_wr
    vf_sw              => clm3%g%l%lps%vf_sw
    vf_rw              => clm3%g%l%lps%vf_rw
    vf_ww              => clm3%g%l%lps%vf_ww
    sabs_roof_dir      => clm3%g%l%lps%sabs_roof_dir
    sabs_roof_dif      => clm3%g%l%lps%sabs_roof_dif
    sabs_sunwall_dir   => clm3%g%l%lps%sabs_sunwall_dir
    sabs_sunwall_dif   => clm3%g%l%lps%sabs_sunwall_dif
    sabs_shadewall_dir => clm3%g%l%lps%sabs_shadewall_dir
    sabs_shadewall_dif => clm3%g%l%lps%sabs_shadewall_dif
    sabs_improad_dir   => clm3%g%l%lps%sabs_improad_dir
    sabs_improad_dif   => clm3%g%l%lps%sabs_improad_dif
    sabs_perroad_dir   => clm3%g%l%lps%sabs_perroad_dir
    sabs_perroad_dif   => clm3%g%l%lps%sabs_perroad_dif

    ! Calculate impervious road

    do l = 1,num_urbanl 
       wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
    end do

    do ib = 1,numrad
       do fl = 1,num_urbanl
        if (coszen(fl) .gt. 0._r8) then
          l = filter_urbanl(fl)

          ! initial absorption and reflection for road and both walls. 
          ! distribute reflected radiation to sky, road, and walls 
          ! according to appropriate view factor. radiation reflected to 
          ! road and walls will undergo multiple reflections within the canyon. 
          ! do separately for direct beam and diffuse radiation.
          
          ! direct beam

          road_a_dir(fl)              = 0.0_r8
          road_r_dir(fl)              = 0.0_r8
          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             improad_a_dir(fl)           = (1._r8-alb_improad_dir(fl,ib)) * sdir_road(fl,ib) 
             improad_r_dir(fl)           =     alb_improad_dir(fl,ib)  * sdir_road(fl,ib) 
             improad_r_sky_dir(fl)       = improad_r_dir(fl) * vf_sr(l)
             improad_r_sunwall_dir(fl)   = improad_r_dir(fl) * vf_wr(l)
             improad_r_shadewall_dir(fl) = improad_r_dir(fl) * vf_wr(l)
             road_a_dir(fl)              = road_a_dir(fl) + improad_a_dir(fl)*wtroad_imperv(fl)
             road_r_dir(fl)              = road_r_dir(fl) + improad_r_dir(fl)*wtroad_imperv(fl)
          end if

          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             perroad_a_dir(fl)           = (1._r8-alb_perroad_dir(fl,ib)) * sdir_road(fl,ib) 
             perroad_r_dir(fl)           =     alb_perroad_dir(fl,ib)  * sdir_road(fl,ib) 
             perroad_r_sky_dir(fl)       = perroad_r_dir(fl) * vf_sr(l)
             perroad_r_sunwall_dir(fl)   = perroad_r_dir(fl) * vf_wr(l)
             perroad_r_shadewall_dir(fl) = perroad_r_dir(fl) * vf_wr(l)
             road_a_dir(fl)              = road_a_dir(fl) + perroad_a_dir(fl)*wtroad_perv(fl)
             road_r_dir(fl)              = road_r_dir(fl) + perroad_r_dir(fl)*wtroad_perv(fl)
          end if

          road_r_sky_dir(fl)          = road_r_dir(fl) * vf_sr(l)
          road_r_sunwall_dir(fl)      = road_r_dir(fl) * vf_wr(l)
          road_r_shadewall_dir(fl)    = road_r_dir(fl) * vf_wr(l)

          sunwall_a_dir(fl)           = (1._r8-alb_wall_dir(fl,ib)) * sdir_sunwall(fl,ib)
          sunwall_r_dir(fl)           =     alb_wall_dir(fl,ib)  * sdir_sunwall(fl,ib)
          sunwall_r_sky_dir(fl)       = sunwall_r_dir(fl) * vf_sw(l)
          sunwall_r_road_dir(fl)      = sunwall_r_dir(fl) * vf_rw(l)
          sunwall_r_shadewall_dir(fl) = sunwall_r_dir(fl) * vf_ww(l)

          shadewall_a_dir(fl)         = (1._r8-alb_wall_dir(fl,ib)) * sdir_shadewall(fl,ib)
          shadewall_r_dir(fl)         =     alb_wall_dir(fl,ib)  * sdir_shadewall(fl,ib)
          shadewall_r_sky_dir(fl)     = shadewall_r_dir(fl) * vf_sw(l)
          shadewall_r_road_dir(fl)    = shadewall_r_dir(fl) * vf_rw(l)
          shadewall_r_sunwall_dir(fl) = shadewall_r_dir(fl) * vf_ww(l)

          ! diffuse

          road_a_dif(fl)              = 0.0_r8
          road_r_dif(fl)              = 0.0_r8
          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             improad_a_dif(fl)           = (1._r8-alb_improad_dif(fl,ib)) * sdif_road(fl,ib) 
             improad_r_dif(fl)           =     alb_improad_dif(fl,ib)  * sdif_road(fl,ib) 
             improad_r_sky_dif(fl)       = improad_r_dif(fl) * vf_sr(l)
             improad_r_sunwall_dif(fl)   = improad_r_dif(fl) * vf_wr(l)
             improad_r_shadewall_dif(fl) = improad_r_dif(fl) * vf_wr(l)
             road_a_dif(fl)              = road_a_dif(fl) + improad_a_dif(fl)*wtroad_imperv(fl)
             road_r_dif(fl)              = road_r_dif(fl) + improad_r_dif(fl)*wtroad_imperv(fl)
          end if

          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             perroad_a_dif(fl)           = (1._r8-alb_perroad_dif(fl,ib)) * sdif_road(fl,ib) 
             perroad_r_dif(fl)           =     alb_perroad_dif(fl,ib)  * sdif_road(fl,ib) 
             perroad_r_sky_dif(fl)       = perroad_r_dif(fl) * vf_sr(l)
             perroad_r_sunwall_dif(fl)   = perroad_r_dif(fl) * vf_wr(l)
             perroad_r_shadewall_dif(fl) = perroad_r_dif(fl) * vf_wr(l)
             road_a_dif(fl)              = road_a_dif(fl) + perroad_a_dif(fl)*wtroad_perv(fl)
             road_r_dif(fl)              = road_r_dif(fl) + perroad_r_dif(fl)*wtroad_perv(fl)
          end if

          road_r_sky_dif(fl)          = road_r_dif(fl) * vf_sr(l)
          road_r_sunwall_dif(fl)      = road_r_dif(fl) * vf_wr(l)
          road_r_shadewall_dif(fl)    = road_r_dif(fl) * vf_wr(l)

          sunwall_a_dif(fl)           = (1._r8-alb_wall_dif(fl,ib)) * sdif_sunwall(fl,ib)
          sunwall_r_dif(fl)           =     alb_wall_dif(fl,ib)  * sdif_sunwall(fl,ib)
          sunwall_r_sky_dif(fl)       = sunwall_r_dif(fl) * vf_sw(l)
          sunwall_r_road_dif(fl)      = sunwall_r_dif(fl) * vf_rw(l)
          sunwall_r_shadewall_dif(fl) = sunwall_r_dif(fl) * vf_ww(l)

          shadewall_a_dif(fl)         = (1._r8-alb_wall_dif(fl,ib)) * sdif_shadewall(fl,ib)
          shadewall_r_dif(fl)         =     alb_wall_dif(fl,ib)  * sdif_shadewall(fl,ib)
          shadewall_r_sky_dif(fl)     = shadewall_r_dif(fl) * vf_sw(l)
          shadewall_r_road_dif(fl)    = shadewall_r_dif(fl) * vf_rw(l) 
          shadewall_r_sunwall_dif(fl) = shadewall_r_dif(fl) * vf_ww(l) 

          ! initialize sum of direct and diffuse solar absorption and reflection for road and both walls

          if ( wtroad_imperv(fl) > 0.0_r8 ) sabs_improad_dir(l,ib)   = improad_a_dir(fl)
          if ( wtroad_perv(fl)   > 0.0_r8 ) sabs_perroad_dir(l,ib)   = perroad_a_dir(fl)
          sabs_sunwall_dir(l,ib)   = sunwall_a_dir(fl)
          sabs_shadewall_dir(l,ib) = shadewall_a_dir(fl)

          if ( wtroad_imperv(fl) > 0.0_r8 ) sabs_improad_dif(l,ib)   = improad_a_dif(fl)
          if ( wtroad_perv(fl)   > 0.0_r8 ) sabs_perroad_dif(l,ib)   = perroad_a_dif(fl)
          sabs_sunwall_dif(l,ib)   = sunwall_a_dif(fl)
          sabs_shadewall_dif(l,ib) = shadewall_a_dif(fl)

          if ( wtroad_imperv(fl) > 0.0_r8 ) sref_improad_dir(fl,ib)   = improad_r_sky_dir(fl) 
          if ( wtroad_perv(fl)   > 0.0_r8 ) sref_perroad_dir(fl,ib)   = perroad_r_sky_dir(fl) 
          sref_sunwall_dir(fl,ib)   = sunwall_r_sky_dir(fl) 
          sref_shadewall_dir(fl,ib) = shadewall_r_sky_dir(fl) 

          if ( wtroad_imperv(fl) > 0.0_r8 ) sref_improad_dif(fl,ib)   = improad_r_sky_dif(fl)
          if ( wtroad_perv(fl)   > 0.0_r8 ) sref_perroad_dif(fl,ib)   = perroad_r_sky_dif(fl)
          sref_sunwall_dif(fl,ib)   = sunwall_r_sky_dif(fl)
          sref_shadewall_dif(fl,ib) = shadewall_r_sky_dif(fl)
        endif

       end do

       ! absorption and reflection for walls and road with multiple reflections
       ! (i.e., absorb and reflect initial reflection in canyon and allow for 
       ! subsequent scattering)
       !
       ! (1) absorption and reflection of scattered solar radiation
       !     road: reflected fluxes from walls need to be projected to ground area
       !     wall: reflected flux from road needs to be projected to wall area
       !
       ! (2) add absorbed radiation for ith reflection to total absorbed
       !
       ! (3) distribute reflected radiation to sky, road, and walls according to view factors
       !
       ! (4) add solar reflection to sky for ith reflection to total reflection
       !
       ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
       !     small convergence criteria is required to ensure solar radiation is conserved
       !
       ! do separately for direct beam and diffuse
       
       do fl = 1,num_urbanl
        if (coszen(fl) .gt. 0._r8) then
          l = filter_urbanl(fl)

          ! reflected direct beam

          do iter_dir = 1, n
             ! step (1)

             stot(fl) = (sunwall_r_road_dir(fl) + shadewall_r_road_dir(fl))*canyon_hwr(fl)

             road_a_dir(fl) = 0.0_r8
             road_r_dir(fl) = 0.0_r8
             if ( wtroad_imperv(fl) > 0.0_r8 ) then
                improad_a_dir(fl) = (1._r8-alb_improad_dir(fl,ib)) * stot(fl) 
                improad_r_dir(fl) =     alb_improad_dir(fl,ib)  * stot(fl) 
                road_a_dir(fl)    = road_a_dir(fl) + improad_a_dir(fl)*wtroad_imperv(fl)
                road_r_dir(fl)    = road_r_dir(fl) + improad_r_dir(fl)*wtroad_imperv(fl)
             end if
             if ( wtroad_perv(fl)   > 0.0_r8 ) then
                perroad_a_dir(fl) = (1._r8-alb_perroad_dir(fl,ib)) * stot(fl) 
                perroad_r_dir(fl) =     alb_perroad_dir(fl,ib)  * stot(fl) 
                road_a_dir(fl)    = road_a_dir(fl) + perroad_a_dir(fl)*wtroad_perv(fl)
                road_r_dir(fl)    = road_r_dir(fl) + perroad_r_dir(fl)*wtroad_perv(fl)
             end if

             stot(fl) = road_r_sunwall_dir(fl)/canyon_hwr(fl) + shadewall_r_sunwall_dir(fl)
             sunwall_a_dir(fl) = (1._r8-alb_wall_dir(fl,ib)) * stot(fl)
             sunwall_r_dir(fl) =     alb_wall_dir(fl,ib)  * stot(fl)

             stot(fl) = road_r_shadewall_dir(fl)/canyon_hwr(fl) + sunwall_r_shadewall_dir(fl)
             shadewall_a_dir(fl) = (1._r8-alb_wall_dir(fl,ib)) * stot(fl)
             shadewall_r_dir(fl) =     alb_wall_dir(fl,ib)  * stot(fl)

             ! step (2)

             if ( wtroad_imperv(fl) > 0.0_r8 ) sabs_improad_dir(l,ib)   = sabs_improad_dir(l,ib)   + improad_a_dir(fl)
             if ( wtroad_perv(fl)   > 0.0_r8 ) sabs_perroad_dir(l,ib)   = sabs_perroad_dir(l,ib)   + perroad_a_dir(fl)
             sabs_sunwall_dir(l,ib)   = sabs_sunwall_dir(l,ib)   + sunwall_a_dir(fl)
             sabs_shadewall_dir(l,ib) = sabs_shadewall_dir(l,ib) + shadewall_a_dir(fl)

             ! step (3)

             if ( wtroad_imperv(fl) > 0.0_r8 ) then
                improad_r_sky_dir(fl)       = improad_r_dir(fl) * vf_sr(l)
                improad_r_sunwall_dir(fl)   = improad_r_dir(fl) * vf_wr(l)
                improad_r_shadewall_dir(fl) = improad_r_dir(fl) * vf_wr(l)
             end if

             if ( wtroad_perv(fl)   > 0.0_r8 ) then
                perroad_r_sky_dir(fl)       = perroad_r_dir(fl) * vf_sr(l)
                perroad_r_sunwall_dir(fl)   = perroad_r_dir(fl) * vf_wr(l)
                perroad_r_shadewall_dir(fl) = perroad_r_dir(fl) * vf_wr(l)
             end if

             road_r_sky_dir(fl)          = road_r_dir(fl) * vf_sr(l)
             road_r_sunwall_dir(fl)      = road_r_dir(fl) * vf_wr(l)
             road_r_shadewall_dir(fl)    = road_r_dir(fl) * vf_wr(l)

             sunwall_r_sky_dir(fl)       = sunwall_r_dir(fl) * vf_sw(l)
             sunwall_r_road_dir(fl)      = sunwall_r_dir(fl) * vf_rw(l)
             sunwall_r_shadewall_dir(fl) = sunwall_r_dir(fl) * vf_ww(l)

             shadewall_r_sky_dir(fl)     = shadewall_r_dir(fl) * vf_sw(l)
             shadewall_r_road_dir(fl)    = shadewall_r_dir(fl) * vf_rw(l)
             shadewall_r_sunwall_dir(fl) = shadewall_r_dir(fl) * vf_ww(l)

             ! step (4)

             if ( wtroad_imperv(fl) > 0.0_r8 ) sref_improad_dir(fl,ib)   = sref_improad_dir(fl,ib) + improad_r_sky_dir(fl)
             if ( wtroad_perv(fl)   > 0.0_r8 ) sref_perroad_dir(fl,ib)   = sref_perroad_dir(fl,ib) + perroad_r_sky_dir(fl)
             sref_sunwall_dir(fl,ib)   = sref_sunwall_dir(fl,ib) + sunwall_r_sky_dir(fl)
             sref_shadewall_dir(fl,ib) = sref_shadewall_dir(fl,ib) + shadewall_r_sky_dir(fl)

             ! step (5)

             crit = max(road_a_dir(fl), sunwall_a_dir(fl), shadewall_a_dir(fl))
             if (crit < errcrit) exit
          end do
          if (iter_dir >= n) then
             write (iulog,*) 'urban net solar radiation error: no convergence, direct beam'
             write (iulog,*) 'clm model is stopping'
             call endrun
          endif
 
          ! reflected diffuse
       
          do iter_dif = 1, n
             ! step (1)

             stot(fl) = (sunwall_r_road_dif(fl) + shadewall_r_road_dif(fl))*canyon_hwr(fl)
             road_a_dif(fl)    = 0.0_r8
             road_r_dif(fl)    = 0.0_r8
             if ( wtroad_imperv(fl) > 0.0_r8 ) then
                improad_a_dif(fl) = (1._r8-alb_improad_dif(fl,ib)) * stot(fl) 
                improad_r_dif(fl) =     alb_improad_dif(fl,ib)  * stot(fl) 
                road_a_dif(fl)    = road_a_dif(fl) + improad_a_dif(fl)*wtroad_imperv(fl)
                road_r_dif(fl)    = road_r_dif(fl) + improad_r_dif(fl)*wtroad_imperv(fl)
             end if
             if ( wtroad_perv(fl)   > 0.0_r8 ) then
                perroad_a_dif(fl) = (1._r8-alb_perroad_dif(fl,ib)) * stot(fl) 
                perroad_r_dif(fl) =     alb_perroad_dif(fl,ib)  * stot(fl) 
                road_a_dif(fl)    = road_a_dif(fl) + perroad_a_dif(fl)*wtroad_perv(fl)
                road_r_dif(fl)    = road_r_dif(fl) + perroad_r_dif(fl)*wtroad_perv(fl)
             end if

             stot(fl) = road_r_sunwall_dif(fl)/canyon_hwr(fl) + shadewall_r_sunwall_dif(fl)
             sunwall_a_dif(fl) = (1._r8-alb_wall_dif(fl,ib)) * stot(fl)
             sunwall_r_dif(fl) =     alb_wall_dif(fl,ib)  * stot(fl)

             stot(fl) = road_r_shadewall_dif(fl)/canyon_hwr(fl) + sunwall_r_shadewall_dif(fl)
             shadewall_a_dif(fl) = (1._r8-alb_wall_dif(fl,ib)) * stot(fl)
             shadewall_r_dif(fl) =     alb_wall_dif(fl,ib)  * stot(fl)

             ! step (2)

             if ( wtroad_imperv(fl) > 0.0_r8 ) sabs_improad_dif(l,ib)   = sabs_improad_dif(l,ib)   + improad_a_dif(fl)
             if ( wtroad_perv(fl)   > 0.0_r8 ) sabs_perroad_dif(l,ib)   = sabs_perroad_dif(l,ib)   + perroad_a_dif(fl)
             sabs_sunwall_dif(l,ib)   = sabs_sunwall_dif(l,ib)   + sunwall_a_dif(fl)
             sabs_shadewall_dif(l,ib) = sabs_shadewall_dif(l,ib) + shadewall_a_dif(fl)

             ! step (3)

             if ( wtroad_imperv(fl) > 0.0_r8 ) then
                improad_r_sky_dif(fl)       = improad_r_dif(fl) * vf_sr(l)
                improad_r_sunwall_dif(fl)   = improad_r_dif(fl) * vf_wr(l)
                improad_r_shadewall_dif(fl) = improad_r_dif(fl) * vf_wr(l)
             end if

             if ( wtroad_perv(fl)   > 0.0_r8 ) then
                perroad_r_sky_dif(fl)       = perroad_r_dif(fl) * vf_sr(l)
                perroad_r_sunwall_dif(fl)   = perroad_r_dif(fl) * vf_wr(l)
                perroad_r_shadewall_dif(fl) = perroad_r_dif(fl) * vf_wr(l)
             end if

             road_r_sky_dif(fl)          = road_r_dif(fl) * vf_sr(l)
             road_r_sunwall_dif(fl)      = road_r_dif(fl) * vf_wr(l)
             road_r_shadewall_dif(fl)    = road_r_dif(fl) * vf_wr(l)

             sunwall_r_sky_dif(fl)       = sunwall_r_dif(fl) * vf_sw(l)
             sunwall_r_road_dif(fl)      = sunwall_r_dif(fl) * vf_rw(l)
             sunwall_r_shadewall_dif(fl) = sunwall_r_dif(fl) * vf_ww(l)

             shadewall_r_sky_dif(fl)     = shadewall_r_dif(fl) * vf_sw(l)
             shadewall_r_road_dif(fl)    = shadewall_r_dif(fl) * vf_rw(l)
             shadewall_r_sunwall_dif(fl) = shadewall_r_dif(fl) * vf_ww(l)

             ! step (4)

             if ( wtroad_imperv(fl) > 0.0_r8 ) sref_improad_dif(fl,ib)   = sref_improad_dif(fl,ib)   + improad_r_sky_dif(fl)
             if ( wtroad_perv(fl)   > 0.0_r8 ) sref_perroad_dif(fl,ib)   = sref_perroad_dif(fl,ib)   + perroad_r_sky_dif(fl)
             sref_sunwall_dif(fl,ib)   = sref_sunwall_dif(fl,ib)   + sunwall_r_sky_dif(fl)
             sref_shadewall_dif(fl,ib) = sref_shadewall_dif(fl,ib) + shadewall_r_sky_dif(fl)

             ! step (5)

             crit = max(road_a_dif(fl), sunwall_a_dif(fl), shadewall_a_dif(fl))
             if (crit < errcrit) exit
          end do
          if (iter_dif >= n) then
             write (iulog,*) 'urban net solar radiation error: no convergence, diffuse'
             write (iulog,*) 'clm model is stopping'
             call endrun()
          endif

          ! total reflected by canyon - sum of solar reflection to sky from canyon.
          ! project wall fluxes to horizontal surface
          
          sref_canyon_dir(fl) = 0.0_r8
          sref_canyon_dif(fl) = 0.0_r8
          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             sref_canyon_dir(fl) = sref_canyon_dir(fl) + sref_improad_dir(fl,ib)*wtroad_imperv(fl)
             sref_canyon_dif(fl) = sref_canyon_dif(fl) + sref_improad_dif(fl,ib)*wtroad_imperv(fl)
          end if
          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             sref_canyon_dir(fl) = sref_canyon_dir(fl) + sref_perroad_dir(fl,ib)*wtroad_perv(fl)
             sref_canyon_dif(fl) = sref_canyon_dif(fl) + sref_perroad_dif(fl,ib)*wtroad_perv(fl)
          end if
          sref_canyon_dir(fl) = sref_canyon_dir(fl) + (sref_sunwall_dir(fl,ib) + sref_shadewall_dir(fl,ib))*canyon_hwr(fl)
          sref_canyon_dif(fl) = sref_canyon_dif(fl) + (sref_sunwall_dif(fl,ib) + sref_shadewall_dif(fl,ib))*canyon_hwr(fl)

          ! total absorbed by canyon. project wall fluxes to horizontal surface

          sabs_canyon_dir(fl) = 0.0_r8
          sabs_canyon_dif(fl) = 0.0_r8
          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             sabs_canyon_dir(fl) = sabs_canyon_dir(fl) + sabs_improad_dir(l,ib)*wtroad_imperv(fl)
             sabs_canyon_dif(fl) = sabs_canyon_dif(fl) + sabs_improad_dif(l,ib)*wtroad_imperv(fl)
          end if
          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             sabs_canyon_dir(fl) = sabs_canyon_dir(fl) + sabs_perroad_dir(l,ib)*wtroad_perv(fl)
             sabs_canyon_dif(fl) = sabs_canyon_dif(fl) + sabs_perroad_dif(l,ib)*wtroad_perv(fl)
          end if
          sabs_canyon_dir(fl) = sabs_canyon_dir(fl) + (sabs_sunwall_dir(l,ib) + sabs_shadewall_dir(l,ib))*canyon_hwr(fl)
          sabs_canyon_dif(fl) = sabs_canyon_dif(fl) + (sabs_sunwall_dif(l,ib) + sabs_shadewall_dif(l,ib))*canyon_hwr(fl)

          ! conservation check. note: previous conservation checks confirm partioning of total direct
          ! beam and diffuse radiation from atmosphere to road and walls is conserved as
          !    sdir (from atmosphere) = sdir_road + (sdir_sunwall + sdir_shadewall)*canyon_hwr
          !    sdif (from atmosphere) = sdif_road + (sdif_sunwall + sdif_shadewall)*canyon_hwr

          stot_dir(fl) = sdir_road(fl,ib) + (sdir_sunwall(fl,ib) + sdir_shadewall(fl,ib))*canyon_hwr(fl)
          stot_dif(fl) = sdif_road(fl,ib) + (sdif_sunwall(fl,ib) + sdif_shadewall(fl,ib))*canyon_hwr(fl)

          err = stot_dir(fl) + stot_dif(fl) &
                - (sabs_canyon_dir(fl) + sabs_canyon_dif(fl) + sref_canyon_dir(fl) + sref_canyon_dif(fl))
          if (abs(err) > 0.001_r8 ) then
             write(iulog,*)'urban net solar radiation balance error for ib=',ib,' err= ',err
             write(iulog,*)' l= ',l,' ib= ',ib 
             write(iulog,*)' stot_dir        = ',stot_dir(fl)
             write(iulog,*)' stot_dif        = ',stot_dif(fl)
             write(iulog,*)' sabs_canyon_dir = ',sabs_canyon_dir(fl)
             write(iulog,*)' sabs_canyon_dif = ',sabs_canyon_dif(fl)
             write(iulog,*)' sref_canyon_dir = ',sref_canyon_dir(fl)
             write(iulog,*)' sref_canyon_dif = ',sref_canyon_dir(fl)
             write(iulog,*) 'clm model is stopping'
             call endrun()
          endif

          ! canyon albedo

          canyon_alb_dif(fl) = sref_canyon_dif(fl) / max(stot_dif(fl), 1.e-06_r8)
          canyon_alb_dir(fl) = sref_canyon_dir(fl) / max(stot_dir(fl), 1.e-06_r8)
        end if

       end do   ! end of landunit loop

       ! Refected and absorbed solar radiation per unit incident radiation for roof

       do fl = 1,num_urbanl
        if (coszen(fl) .gt. 0._r8) then
          l = filter_urbanl(fl)
          sref_roof_dir(fl,ib) = alb_roof_dir(fl,ib) * sdir(fl,ib)
          sref_roof_dif(fl,ib) = alb_roof_dif(fl,ib) * sdif(fl,ib)
          sabs_roof_dir(l,ib) = sdir(fl,ib) - sref_roof_dir(fl,ib)
          sabs_roof_dif(l,ib) = sdif(fl,ib) - sref_roof_dif(fl,ib)
        end if
       end do

    end do   ! end of radiation band loop

  end subroutine net_solar

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: net_longwave
!
! !INTERFACE:
  subroutine net_longwave (lbl, ubl, num_urbanl, filter_urbanl, canyon_hwr, wtroad_perv, &
                           lwdown, em_roof, em_improad, em_perroad, em_wall, &
                           t_roof,  t_improad, t_perroad, t_sunwall, t_shadewall, &
                           lwnet_roof, lwnet_improad, lwnet_perroad, lwnet_sunwall, lwnet_shadewall, lwnet_canyon, &
                           lwup_roof, lwup_improad, lwup_perroad, lwup_sunwall, lwup_shadewall, lwup_canyon)
!
! !DESCRIPTION: 
! Net longwave radiation for road and both walls in urban canyon allowing for 
! multiple reflection. Also net longwave radiation for urban roof. 
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_varcon   , only : sb
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: num_urbanl                 ! number of urban landunits
    integer,  intent(in)  :: lbl, ubl                   ! landunit-index bounds
    integer , intent(in)  :: filter_urbanl(ubl-lbl+1)   ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr(num_urbanl)     ! ratio of building height to street width
    real(r8), intent(in)  :: wtroad_perv(num_urbanl)     ! weight of pervious road wrt total road

    real(r8), intent(in)  :: lwdown(num_urbanl)          ! atmospheric longwave radiation (W/m**2)
    real(r8), intent(in)  :: em_roof(num_urbanl)         ! roof emissivity
    real(r8), intent(in)  :: em_improad(num_urbanl)      ! impervious road emissivity
    real(r8), intent(in)  :: em_perroad(num_urbanl)      ! pervious road emissivity
    real(r8), intent(in)  :: em_wall(num_urbanl)         ! wall emissivity

    real(r8), intent(in)  :: t_roof(num_urbanl)          ! roof temperature (K)
    real(r8), intent(in)  :: t_improad(num_urbanl)       ! impervious road temperature (K)
    real(r8), intent(in)  :: t_perroad(num_urbanl)       ! ervious road temperature (K)
    real(r8), intent(in)  :: t_sunwall(num_urbanl)       ! sunlit wall temperature (K)
    real(r8), intent(in)  :: t_shadewall(num_urbanl)     ! shaded wall temperature (K)

    real(r8), intent(out) :: lwnet_roof(num_urbanl)      ! net (outgoing-incoming) longwave radiation, roof (W/m**2)
    real(r8), intent(out) :: lwnet_improad(num_urbanl)   ! net (outgoing-incoming) longwave radiation, impervious road (W/m**2)
    real(r8), intent(out) :: lwnet_perroad(num_urbanl)   ! net (outgoing-incoming) longwave radiation, pervious road (W/m**2)
    real(r8), intent(out) :: lwnet_sunwall(num_urbanl)   ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8), intent(out) :: lwnet_shadewall(num_urbanl) ! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8), intent(out) :: lwnet_canyon(num_urbanl)    ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)

    real(r8), intent(out) :: lwup_roof(num_urbanl)       ! upward longwave radiation, roof (W/m**2)
    real(r8), intent(out) :: lwup_improad(num_urbanl)    ! upward longwave radiation, impervious road (W/m**2)
    real(r8), intent(out) :: lwup_perroad(num_urbanl)    ! upward longwave radiation, pervious road (W/m**2)
    real(r8), intent(out) :: lwup_sunwall(num_urbanl)    ! upward longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8), intent(out) :: lwup_shadewall(num_urbanl)  ! upward longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8), intent(out) :: lwup_canyon(num_urbanl)     ! upward longwave radiation for canyon, per unit ground area (W/m**2)
!
! local pointers to original implicit in arguments (clmtype)
!
    real(r8), pointer :: vf_sr(:)     ! view factor of sky for road
    real(r8), pointer :: vf_wr(:)     ! view factor of one wall for road
    real(r8), pointer :: vf_sw(:)     ! view factor of sky for one wall
    real(r8), pointer :: vf_rw(:)     ! view factor of road for one wall
    real(r8), pointer :: vf_ww(:)     ! view factor of opposing wall for one wall
!
! !CALLED FROM:
! subroutine UrbanRadiation in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: lwdown_road(num_urbanl)         ! atmospheric longwave radiation for total road (W/m**2)
    real(r8) :: lwdown_sunwall(num_urbanl)      ! atmospheric longwave radiation (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: lwdown_shadewall(num_urbanl)    ! atmospheric longwave radiation (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: lwtot(num_urbanl)               ! incoming longwave radiation (W/m**2)

    real(r8) :: improad_a(num_urbanl)           ! absorbed longwave for improad (W/m**2)
    real(r8) :: improad_r(num_urbanl)           ! reflected longwave for improad (W/m**2)
    real(r8) :: improad_r_sky(num_urbanl)       ! improad_r to sky (W/m**2)
    real(r8) :: improad_r_sunwall(num_urbanl)   ! improad_r to sunlit wall (W/m**2)
    real(r8) :: improad_r_shadewall(num_urbanl) ! improad_r to shaded wall (W/m**2)
    real(r8) :: improad_e(num_urbanl)           ! emitted longwave for improad (W/m**2)
    real(r8) :: improad_e_sky(num_urbanl)       ! improad_e to sky (W/m**2)
    real(r8) :: improad_e_sunwall(num_urbanl)   ! improad_e to sunlit wall (W/m**2)
    real(r8) :: improad_e_shadewall(num_urbanl) ! improad_e to shaded wall (W/m**2)

    real(r8) :: perroad_a(num_urbanl)           ! absorbed longwave for perroad (W/m**2)
    real(r8) :: perroad_r(num_urbanl)           ! reflected longwave for perroad (W/m**2)
    real(r8) :: perroad_r_sky(num_urbanl)       ! perroad_r to sky (W/m**2)
    real(r8) :: perroad_r_sunwall(num_urbanl)   ! perroad_r to sunlit wall (W/m**2)
    real(r8) :: perroad_r_shadewall(num_urbanl) ! perroad_r to shaded wall (W/m**2)
    real(r8) :: perroad_e(num_urbanl)           ! emitted longwave for perroad (W/m**2)
    real(r8) :: perroad_e_sky(num_urbanl)       ! perroad_e to sky (W/m**2)
    real(r8) :: perroad_e_sunwall(num_urbanl)   ! perroad_e to sunlit wall (W/m**2)
    real(r8) :: perroad_e_shadewall(num_urbanl) ! perroad_e to shaded wall (W/m**2)

    real(r8) :: road_a(num_urbanl)              ! absorbed longwave for total road (W/m**2)
    real(r8) :: road_r(num_urbanl)              ! reflected longwave for total road (W/m**2)
    real(r8) :: road_r_sky(num_urbanl)          ! total road_r to sky (W/m**2)
    real(r8) :: road_r_sunwall(num_urbanl)      ! total road_r to sunlit wall (W/m**2)
    real(r8) :: road_r_shadewall(num_urbanl)    ! total road_r to shaded wall (W/m**2)
    real(r8) :: road_e(num_urbanl)              ! emitted longwave for total road (W/m**2)
    real(r8) :: road_e_sky(num_urbanl)          ! total road_e to sky (W/m**2)
    real(r8) :: road_e_sunwall(num_urbanl)      ! total road_e to sunlit wall (W/m**2)
    real(r8) :: road_e_shadewall(num_urbanl)    ! total road_e to shaded wall (W/m**2)

    real(r8) :: sunwall_a(num_urbanl)           ! absorbed longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r(num_urbanl)           ! reflected longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r_sky(num_urbanl)       ! sunwall_r to sky (W/m**2)
    real(r8) :: sunwall_r_road(num_urbanl)      ! sunwall_r to road (W/m**2)
    real(r8) :: sunwall_r_shadewall(num_urbanl) ! sunwall_r to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_e(num_urbanl)           ! emitted longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_e_sky(num_urbanl)       ! sunwall_e to sky (W/m**2)
    real(r8) :: sunwall_e_road(num_urbanl)      ! sunwall_e to road (W/m**2)
    real(r8) :: sunwall_e_shadewall(num_urbanl) ! sunwall_e to opposing (shaded) wall (W/m**2)

    real(r8) :: shadewall_a(num_urbanl)         ! absorbed longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r(num_urbanl)         ! reflected longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r_sky(num_urbanl)     ! shadewall_r to sky (W/m**2)
    real(r8) :: shadewall_r_road(num_urbanl)    ! shadewall_r to road (W/m**2)
    real(r8) :: shadewall_r_sunwall(num_urbanl) ! shadewall_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_e(num_urbanl)         ! emitted longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_e_sky(num_urbanl)     ! shadewall_e to sky (W/m**2)
    real(r8) :: shadewall_e_road(num_urbanl)    ! shadewall_e to road (W/m**2)
    real(r8) :: shadewall_e_sunwall(num_urbanl) ! shadewall_e to opposing (sunlit) wall (W/m**2)
    integer  :: l,fl,iter                       ! indices
    integer, parameter  :: n = 50               ! number of interations
    real(r8) :: crit                            ! convergence criterion (W/m**2)
    real(r8) :: err                             ! energy conservation error (W/m**2)
    real(r8) :: wtroad_imperv(num_urbanl)       ! weight of impervious road wrt total road
!-----------------------------------------------------------------------

    ! Assign landunit level pointer

    vf_sr              => clm3%g%l%lps%vf_sr
    vf_wr              => clm3%g%l%lps%vf_wr
    vf_sw              => clm3%g%l%lps%vf_sw
    vf_rw              => clm3%g%l%lps%vf_rw
    vf_ww              => clm3%g%l%lps%vf_ww

    ! Calculate impervious road

    do l = 1,num_urbanl 
       wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
    end do

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       ! atmospheric longwave radiation incident on walls and road in urban canyon.
       ! check for conservation (need to convert wall fluxes to ground area).
       ! lwdown (from atmosphere) = lwdown_road + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr
       
       lwdown_road(fl)      = lwdown(fl) * vf_sr(l) 
       lwdown_sunwall(fl)   = lwdown(fl) * vf_sw(l)
       lwdown_shadewall(fl) = lwdown(fl) * vf_sw(l) 

       err = lwdown(fl) - (lwdown_road(fl) + (lwdown_shadewall(fl) + lwdown_sunwall(fl))*canyon_hwr(fl))
       if (abs(err) > 0.10_r8 ) then
          write (iulog,*) 'urban incident atmospheric longwave radiation balance error',err
          write (iulog,*) 'clm model is stopping'
          call endrun
       endif
    end do

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       ! initial absorption, reflection, and emission for road and both walls. 
       ! distribute reflected and emitted radiation to sky, road, and walls according 
       ! to appropriate view factor. radiation reflected to road and walls will
       ! undergo multiple reflections within the canyon.

       road_a(fl)              = 0.0_r8
       road_r(fl)              = 0.0_r8
       road_e(fl)              = 0.0_r8
       if ( wtroad_imperv(fl) > 0.0_r8 ) then
          improad_a(fl)           =     em_improad(fl)  * lwdown_road(fl) 
          improad_r(fl)           = (1._r8-em_improad(fl)) * lwdown_road(fl) 
          improad_r_sky(fl)       = improad_r(fl) * vf_sr(l)
          improad_r_sunwall(fl)   = improad_r(fl) * vf_wr(l)
          improad_r_shadewall(fl) = improad_r(fl) * vf_wr(l)
          improad_e(fl)           = em_improad(fl) * sb * (t_improad(fl)**4) 
          improad_e_sky(fl)       = improad_e(fl) * vf_sr(l)
          improad_e_sunwall(fl)   = improad_e(fl) * vf_wr(l)
          improad_e_shadewall(fl) = improad_e(fl) * vf_wr(l)
          road_a(fl)              = road_a(fl) + improad_a(fl)*wtroad_imperv(fl)
          road_r(fl)              = road_r(fl) + improad_r(fl)*wtroad_imperv(fl)
          road_e(fl)              = road_e(fl) + improad_e(fl)*wtroad_imperv(fl)
       end if

       if ( wtroad_perv(fl)   > 0.0_r8 ) then
          perroad_a(fl)           =     em_perroad(fl)  * lwdown_road(fl)
          perroad_r(fl)           = (1._r8-em_perroad(fl)) * lwdown_road(fl)
          perroad_r_sky(fl)       = perroad_r(fl) * vf_sr(l)
          perroad_r_sunwall(fl)   = perroad_r(fl) * vf_wr(l)
          perroad_r_shadewall(fl) = perroad_r(fl) * vf_wr(l)
          perroad_e(fl)           = em_perroad(fl) * sb * (t_perroad(fl)**4) 
          perroad_e_sky(fl)       = perroad_e(fl) * vf_sr(l)
          perroad_e_sunwall(fl)   = perroad_e(fl) * vf_wr(l)
          perroad_e_shadewall(fl) = perroad_e(fl) * vf_wr(l)
          road_a(fl)              = road_a(fl) + perroad_a(fl)*wtroad_perv(fl)
          road_r(fl)              = road_r(fl) + perroad_r(fl)*wtroad_perv(fl)
          road_e(fl)              = road_e(fl) + perroad_e(fl)*wtroad_perv(fl)
       end if

       road_r_sky(fl)          = road_r(fl) * vf_sr(l)
       road_r_sunwall(fl)      = road_r(fl) * vf_wr(l)
       road_r_shadewall(fl)    = road_r(fl) * vf_wr(l)
       road_e_sky(fl)          = road_e(fl) * vf_sr(l)
       road_e_sunwall(fl)      = road_e(fl) * vf_wr(l)
       road_e_shadewall(fl)    = road_e(fl) * vf_wr(l)

       sunwall_a(fl)           = em_wall(fl) * lwdown_sunwall(fl)
       sunwall_r(fl)           = (1._r8-em_wall(fl)) * lwdown_sunwall(fl)
       sunwall_r_sky(fl)       = sunwall_r(fl) * vf_sw(l)
       sunwall_r_road(fl)      = sunwall_r(fl) * vf_rw(l)
       sunwall_r_shadewall(fl) = sunwall_r(fl) * vf_ww(l)
       sunwall_e(fl)           = em_wall(fl) * sb * (t_sunwall(fl)**4) 
       sunwall_e_sky(fl)       = sunwall_e(fl) * vf_sw(l)
       sunwall_e_road(fl)      = sunwall_e(fl) * vf_rw(l)
       sunwall_e_shadewall(fl) = sunwall_e(fl) * vf_ww(l)

       shadewall_a(fl)         = em_wall(fl) * lwdown_shadewall(fl)
       shadewall_r(fl)         = (1._r8-em_wall(fl)) * lwdown_shadewall(fl)
       shadewall_r_sky(fl)     = shadewall_r(fl) * vf_sw(l)
       shadewall_r_road(fl)    = shadewall_r(fl) * vf_rw(l)
       shadewall_r_sunwall(fl) = shadewall_r(fl) * vf_ww(l)
       shadewall_e(fl)         = em_wall(fl) * sb * (t_shadewall(fl)**4) 
       shadewall_e_sky(fl)     = shadewall_e(fl) * vf_sw(l)
       shadewall_e_road(fl)    = shadewall_e(fl) * vf_rw(l)
       shadewall_e_sunwall(fl) = shadewall_e(fl) * vf_ww(l)

       ! initialize sum of net and upward longwave radiation for road and both walls

       if ( wtroad_imperv(fl) > 0.0_r8 ) lwnet_improad(fl)   = improad_e(fl)   - improad_a(fl)
       if ( wtroad_perv(fl)   > 0.0_r8 ) lwnet_perroad(fl)   = perroad_e(fl)   - perroad_a(fl)
       lwnet_sunwall(fl)   = sunwall_e(fl)   - sunwall_a(fl)
       lwnet_shadewall(fl) = shadewall_e(fl) - shadewall_a(fl)

       if ( wtroad_imperv(fl) > 0.0_r8 ) lwup_improad(fl)   = improad_r_sky(fl)   + improad_e_sky(fl)
       if ( wtroad_perv(fl)   > 0.0_r8 ) lwup_perroad(fl)   = perroad_r_sky(fl)   + perroad_e_sky(fl)
       lwup_sunwall(fl)   = sunwall_r_sky(fl)   + sunwall_e_sky(fl)
       lwup_shadewall(fl) = shadewall_r_sky(fl) + shadewall_e_sky(fl)

    end do

    ! now account for absorption and reflection within canyon of fluxes from road and walls 
    ! allowing for multiple reflections
    !
    ! (1) absorption and reflection. note: emission from road and walls absorbed by walls and roads
    !     only occurs in first iteration. zero out for later iterations.
    !
    !     road: fluxes from walls need to be projected to ground area
    !     wall: fluxes from road need to be projected to wall area
    !
    ! (2) add net longwave for ith reflection to total net longwave
    !
    ! (3) distribute reflected radiation to sky, road, and walls according to view factors
    !
    ! (4) add upward longwave radiation to sky from road and walls for ith reflection to total
    !
    ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
    !     small convergence criteria is required to ensure radiation is conserved

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       do iter = 1, n
          ! step (1)

          lwtot(fl) =  (sunwall_r_road(fl) + sunwall_e_road(fl)  &
                       + shadewall_r_road(fl) + shadewall_e_road(fl))*canyon_hwr(fl)
          road_a(fl)    = 0.0_r8
          road_r(fl)    = 0.0_r8
          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             improad_r(fl) = (1._r8-em_improad(fl)) * lwtot(fl)
             improad_a(fl) =     em_improad(fl)  * lwtot(fl)
             road_a(fl)    = road_a(fl) + improad_a(fl)*wtroad_imperv(fl)
             road_r(fl)    = road_r(fl) + improad_r(fl)*wtroad_imperv(fl)
          end if
          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             perroad_r(fl) = (1._r8-em_perroad(fl)) * lwtot(fl) 
             perroad_a(fl) =     em_perroad(fl)  * lwtot(fl) 
             road_a(fl)    = road_a(fl) + perroad_a(fl)*wtroad_perv(fl)
             road_r(fl)    = road_r(fl) + perroad_r(fl)*wtroad_perv(fl)
          end if

          lwtot(fl) = (road_r_sunwall(fl) + road_e_sunwall(fl))/canyon_hwr(fl) &
                      + (shadewall_r_sunwall(fl) + shadewall_e_sunwall(fl))
          sunwall_a(fl) =     em_wall(fl)  * lwtot(fl)
          sunwall_r(fl) = (1._r8-em_wall(fl)) * lwtot(fl)

          lwtot(fl) = (road_r_shadewall(fl) + road_e_shadewall(fl))/canyon_hwr(fl) &
                      + (sunwall_r_shadewall(fl) + sunwall_e_shadewall(fl))
          shadewall_a(fl) =     em_wall(fl)  * lwtot(fl)
          shadewall_r(fl) = (1._r8-em_wall(fl)) * lwtot(fl)

          sunwall_e_road(fl)      = 0._r8
          shadewall_e_road(fl)    = 0._r8
          road_e_sunwall(fl)      = 0._r8
          shadewall_e_sunwall(fl) = 0._r8
          road_e_shadewall(fl)    = 0._r8
          sunwall_e_shadewall(fl) = 0._r8

          ! step (2)

          if ( wtroad_imperv(fl) > 0.0_r8 ) lwnet_improad(fl)   = lwnet_improad(fl)   - improad_a(fl)
          if ( wtroad_perv(fl)   > 0.0_r8 ) lwnet_perroad(fl)   = lwnet_perroad(fl)   - perroad_a(fl)
          lwnet_sunwall(fl)   = lwnet_sunwall(fl)   - sunwall_a(fl)
          lwnet_shadewall(fl) = lwnet_shadewall(fl) - shadewall_a(fl)

          ! step (3)

          if ( wtroad_imperv(fl) > 0.0_r8 ) then
             improad_r_sky(fl)       = improad_r(fl) * vf_sr(l)
             improad_r_sunwall(fl)   = improad_r(fl) * vf_wr(l)
             improad_r_shadewall(fl) = improad_r(fl) * vf_wr(l)
          end if

          if ( wtroad_perv(fl)   > 0.0_r8 ) then
             perroad_r_sky(fl)       = perroad_r(fl) * vf_sr(l)
             perroad_r_sunwall(fl)   = perroad_r(fl) * vf_wr(l)
             perroad_r_shadewall(fl) = perroad_r(fl) * vf_wr(l)
          end if

          road_r_sky(fl)          = road_r(fl) * vf_sr(l)
          road_r_sunwall(fl)      = road_r(fl) * vf_wr(l)
          road_r_shadewall(fl)    = road_r(fl) * vf_wr(l)

          sunwall_r_sky(fl)       = sunwall_r(fl) * vf_sw(l)
          sunwall_r_road(fl)      = sunwall_r(fl) * vf_rw(l)
          sunwall_r_shadewall(fl) = sunwall_r(fl) * vf_ww(l)

          shadewall_r_sky(fl)     = shadewall_r(fl) * vf_sw(l)
          shadewall_r_road(fl)    = shadewall_r(fl) * vf_rw(l)
          shadewall_r_sunwall(fl) = shadewall_r(fl) * vf_ww(l)

          ! step (4)

          if ( wtroad_imperv(fl) > 0.0_r8 ) lwup_improad(fl)   = lwup_improad(fl)   + improad_r_sky(fl)
          if ( wtroad_perv(fl)   > 0.0_r8 ) lwup_perroad(fl)   = lwup_perroad(fl)   + perroad_r_sky(fl)
          lwup_sunwall(fl)   = lwup_sunwall(fl)   + sunwall_r_sky(fl)
          lwup_shadewall(fl) = lwup_shadewall(fl) + shadewall_r_sky(fl)

          ! step (5)

          crit = max(road_a(fl), sunwall_a(fl), shadewall_a(fl))
          if (crit < .001_r8) exit
       end do
       if (iter >= n) then
          write (iulog,*) 'urban net longwave radiation error: no convergence'
          write (iulog,*) 'clm model is stopping'
          call endrun
       endif

       ! total net longwave radiation for canyon. project wall fluxes to horizontal surface

       lwnet_canyon(fl) = 0.0_r8
       if ( wtroad_imperv(fl) > 0.0_r8 ) lwnet_canyon(fl) = lwnet_canyon(fl) + lwnet_improad(fl)*wtroad_imperv(fl)
       if ( wtroad_perv(fl)   > 0.0_r8 ) lwnet_canyon(fl) = lwnet_canyon(fl) + lwnet_perroad(fl)*wtroad_perv(fl)
       lwnet_canyon(fl) = lwnet_canyon(fl) + (lwnet_sunwall(fl) + lwnet_shadewall(fl))*canyon_hwr(fl)

       ! total emitted longwave for canyon. project wall fluxes to horizontal

       lwup_canyon(fl) = 0.0_r8
       if( wtroad_imperv(fl) > 0.0_r8 ) lwup_canyon(fl) = lwup_canyon(fl) + lwup_improad(fl)*wtroad_imperv(fl)
       if( wtroad_perv(fl)   > 0.0_r8 ) lwup_canyon(fl) = lwup_canyon(fl) + lwup_perroad(fl)*wtroad_perv(fl)
       lwup_canyon(fl) = lwup_canyon(fl) + (lwup_sunwall(fl) + lwup_shadewall(fl))*canyon_hwr(fl)

       ! conservation check. note: previous conservation check confirms partioning of incident
       ! atmospheric longwave radiation to road and walls is conserved as
       ! lwdown (from atmosphere) = lwdown_improad + lwdown_perroad + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

       err = lwnet_canyon(fl) - (lwup_canyon(fl) - lwdown(fl))
       if (abs(err) > .10_r8 ) then
          write (iulog,*) 'urban net longwave radiation balance error',err
          write (iulog,*) 'clm model is stopping'
          call endrun()
       end if

    end do

    ! Net longwave radiation for roof
    
    do l = 1,num_urbanl
       lwup_roof(l) = em_roof(l)*sb*(t_roof(l)**4) + (1._r8-em_roof(l))*lwdown(l)
       lwnet_roof(l) = lwup_roof(l) - lwdown(l)
    end do

  end subroutine net_longwave

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanClumpInit
!
! !INTERFACE:
  subroutine UrbanClumpInit()
!
! !DESCRIPTION: 
! Initialize urban surface dataset variables
!
! !USES:
    use clmtype
    use clm_varcon   , only : spval, icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, udens_base
    use decompMod    , only : get_proc_clumps, ldecomp
    use filterMod    , only : filter
    use UrbanInputMod, only : urbinp
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein 04/2003
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: coli(:)         ! beginning column index for landunit 
    integer , pointer :: colf(:)         ! ending column index for landunit
    integer , pointer :: lgridcell(:)    ! gridcell of corresponding landunit
    integer , pointer :: ctype(:)        ! column type
    integer , pointer :: udenstype(:)    ! urban density type
!
!
! !OTHER LOCAL VARIABLES
!EOP
!
    integer :: nc,fl,ib,l,c,p,g          ! indices
    integer :: nclumps                   ! number of clumps on processor 
    integer :: num_urbanl                ! number of per-clump urban landunits
    integer :: ier                       ! error status
    integer :: dindx                     ! urban density type index
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit-level)

    coli         => clm3%g%l%coli
    colf         => clm3%g%l%colf
    lgridcell    => clm3%g%l%gridcell
    udenstype    => clm3%g%l%udenstype

    ! Assign local pointers to derived type members (column-level)

    ctype      => clm3%g%l%c%itype

    ! Allocate memory 

    nclumps = get_proc_clumps()
    allocate(urban_clump(nclumps), stat=ier)
    if (ier /= 0) then
       write (iulog,*) 'UrbanInit: allocation error for urban clumps'; call endrun()
    end if

    ! Loop over all clumps on this processor

    do nc = 1, nclumps

       ! Determine number of urban landunits in clump

       num_urbanl = filter(nc)%num_urbanl

       ! Consistency check for urban columns

       do fl = 1,num_urbanl
          l = filter(nc)%urbanl(fl)
          do c = coli(l),colf(l)
             if ( ctype(c) /= icol_roof .and.  &
	          ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall .and. &
	          ctype(c) /= icol_road_perv .and.  ctype(c) /= icol_road_imperv) then
                write(iulog,*)'error in urban column types for landunit = ',l
                write(iulog,*)'ctype= ',ctype(c)
                call endrun()
             endif
          end do
       end do

       ! Allocate memory for urban clump clumponents

       if (num_urbanl > 0) then
          allocate( urban_clump(nc)%canyon_hwr        (num_urbanl),        &
	            urban_clump(nc)%wtroad_perv       (num_urbanl),        &
                    urban_clump(nc)%ht_roof           (num_urbanl),        &
                    urban_clump(nc)%wtlunit_roof      (num_urbanl),        &
                    urban_clump(nc)%wind_hgt_canyon   (num_urbanl),        &
                    urban_clump(nc)%em_roof           (num_urbanl),        &
                    urban_clump(nc)%em_improad        (num_urbanl),        &
                    urban_clump(nc)%em_perroad        (num_urbanl),        &
                    urban_clump(nc)%em_wall           (num_urbanl),        &
                    urban_clump(nc)%alb_roof_dir      (num_urbanl,numrad), &
                    urban_clump(nc)%alb_roof_dif      (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_improad_dir   (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_perroad_dir   (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_improad_dif   (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_perroad_dif   (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_wall_dir      (num_urbanl,numrad), &        
                    urban_clump(nc)%alb_wall_dif      (num_urbanl,numrad), stat=ier )
          if (ier /= 0) then
             write(iulog,*)'UrbanClumpInit: allocation error for urban derived type'; call endrun()
          endif
       end if

       ! Set constants in derived type values for urban clump

       do fl = 1,num_urbanl
          l = filter(nc)%urbanl(fl)
          g = clm3%g%l%gridcell(l)
          dindx = udenstype(l) - udens_base
          urban_clump(nc)%canyon_hwr     (fl) = urbinp%canyon_hwr     (g,dindx)
          urban_clump(nc)%wtroad_perv    (fl) = urbinp%wtroad_perv    (g,dindx)
          urban_clump(nc)%ht_roof        (fl) = urbinp%ht_roof        (g,dindx)
          urban_clump(nc)%wtlunit_roof   (fl) = urbinp%wtlunit_roof   (g,dindx)
          urban_clump(nc)%wind_hgt_canyon(fl) = urbinp%wind_hgt_canyon(g,dindx)
          do ib = 1,numrad
             urban_clump(nc)%alb_roof_dir   (fl,ib) = urbinp%alb_roof_dir   (g,dindx,ib)
             urban_clump(nc)%alb_roof_dif   (fl,ib) = urbinp%alb_roof_dif   (g,dindx,ib)
             urban_clump(nc)%alb_improad_dir(fl,ib) = urbinp%alb_improad_dir(g,dindx,ib)
             urban_clump(nc)%alb_perroad_dir(fl,ib) = urbinp%alb_perroad_dir(g,dindx,ib)
             urban_clump(nc)%alb_improad_dif(fl,ib) = urbinp%alb_improad_dif(g,dindx,ib)
             urban_clump(nc)%alb_perroad_dif(fl,ib) = urbinp%alb_perroad_dif(g,dindx,ib)
             urban_clump(nc)%alb_wall_dir   (fl,ib) = urbinp%alb_wall_dir   (g,dindx,ib)
             urban_clump(nc)%alb_wall_dif   (fl,ib) = urbinp%alb_wall_dif   (g,dindx,ib)
          end do
          urban_clump(nc)%em_roof   (fl) = urbinp%em_roof   (g,dindx)
          urban_clump(nc)%em_improad(fl) = urbinp%em_improad(g,dindx)
          urban_clump(nc)%em_perroad(fl) = urbinp%em_perroad(g,dindx)
          urban_clump(nc)%em_wall   (fl) = urbinp%em_wall   (g,dindx)
       end do
    end do   ! end of loop over clumps

  end subroutine UrbanClumpInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanFluxes
!
! !INTERFACE:
  subroutine UrbanFluxes (nc, lbp, ubp, lbl, ubl, lbc, ubc, &
                          num_nourbanl, filter_nourbanl, &
                          num_urbanl, filter_urbanl, &
                          num_urbanc, filter_urbanc, &
                          num_urbanp, filter_urbanp)
!
! !DESCRIPTION: 
! Turbulent and momentum fluxes from urban canyon (consisting of roof, sunwall, 
! shadewall, pervious and impervious road).

! !USES:
    use clmtype
    use clm_varcon         , only : cpair, vkc, spval, icol_roof, icol_sunwall, &
                                    icol_shadewall, icol_road_perv, icol_road_imperv, &
                                    grav, pondmx_urban, rpi, rgas, &
                                    ht_wasteheat_factor, ac_wasteheat_factor, &
                                    wasteheat_limit
    use filterMod          , only : filter
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat
    use clm_varpar         , only : maxpatch_urb, nlevurb, nlevgrnd
    use clm_time_manager   , only : get_curr_date, get_step_size, get_nstep
    use clm_atmlnd         , only : clm_a2l

!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nc                         ! clump index
    integer, intent(in)  :: lbp, ubp                   ! pft-index bounds
    integer, intent(in)  :: lbl, ubl                   ! landunit-index bounds
    integer, intent(in)  :: lbc, ubc                   ! column-index bounds
    integer , intent(in) :: num_nourbanl               ! number of non-urban landunits in clump
    integer , intent(in) :: filter_nourbanl(ubl-lbl+1) ! non-urban landunit filter
    integer , intent(in) :: num_urbanl                 ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(ubl-lbl+1)   ! urban landunit filter
    integer , intent(in) :: num_urbanc                 ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(ubc-lbc+1)   ! urban column filter
    integer , intent(in) :: num_urbanp                 ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(ubp-lbp+1)   ! urban pft filter
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Author: Keith Oleson 10/2005
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments (urban clump)
!
    real(r8), pointer :: ht_roof(:)         ! height of urban roof (m)
    real(r8), pointer :: wtlunit_roof(:)    ! weight of roof with respect to landunit
    real(r8), pointer :: canyon_hwr(:)      ! ratio of building height to street width
    real(r8), pointer :: wtroad_perv(:)     ! weight of pervious road wrt total road
    real(r8), pointer :: wind_hgt_canyon(:) ! height above road at which wind in canyon is to be computed (m)
!
! local pointers to original implicit in arguments (clmtype)
!
    real(r8), pointer :: forc_u(:)     ! atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)     ! atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_rho(:)   ! density (kg/m**3)
    real(r8), pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft-level (m)
    real(r8), pointer :: forc_hgt_t_pft(:) ! observational height of temperature at pft-level (m)
    real(r8), pointer :: forc_q(:)     ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_t(:)     ! atmospheric temperature (K)
    real(r8), pointer :: forc_th(:)    ! atmospheric potential temperature (K)
    real(r8), pointer :: forc_pbot(:)  ! atmospheric pressure (Pa)

    real(r8), pointer :: z_0_town(:)   ! momentum roughness length of urban landunit (m)
    real(r8), pointer :: z_d_town(:)   ! displacement height of urban landunit (m)

    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
    integer , pointer :: plandunit(:)  ! pft's landunit index
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: coli(:)       ! beginning column index for landunit 
    integer , pointer :: colf(:)       ! ending column index for landunit
    integer , pointer :: pfti(:)       ! beginning pft index for landunit 
    integer , pointer :: pftf(:)       ! ending pft index for landunit

    real(r8), pointer :: taf(:)        ! urban canopy air temperature (K)
    real(r8), pointer :: qaf(:)        ! urban canopy air specific humidity (kg/kg)
    integer , pointer :: npfts(:)      ! landunit's number of pfts (columns)
    real(r8), pointer :: t_grnd(:)     ! ground surface temperature (K)
    real(r8), pointer :: qg(:)         ! specific humidity at ground surface (kg/kg)
    real(r8), pointer :: htvp(:)       ! latent heat of evaporation (/sublimation) (J/kg)
    real(r8), pointer :: dqgdT(:)      ! temperature derivative of "qg"
    real(r8), pointer :: eflx_traffic(:)        ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_traffic_factor(:) ! multiplicative urban traffic factor for sensible heat flux
    real(r8), pointer :: eflx_wasteheat(:)      ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac(:)   ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: t_soisno(:,:)          ! soil temperature (K)
    real(r8), pointer :: eflx_urban_ac(:)     ! urban air conditioning flux (W/m**2)
    real(r8), pointer :: eflx_urban_heat(:)   ! urban heating flux (W/m**2)
    real(r8), pointer :: londeg(:)            ! longitude (degrees)
    real(r8), pointer :: h2osoi_ice(:,:)      ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)      ! liquid water (kg/m2)
    real(r8), pointer :: frac_sno(:)          ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: snow_depth(:)            ! snow height (m)
    real(r8), pointer :: h2osno(:)            ! snow water (mm H2O)
    integer , pointer :: snl(:)               ! number of snow layers
    real(r8), pointer :: rootr_road_perv(:,:) ! effective fraction of roots in each soil layer for urban pervious road
    real(r8), pointer :: soilalpha_u(:)       ! Urban factor that reduces ground saturated specific humidity (-)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: dlrad(:)         ! downward longwave radiation below the canopy (W/m**2)
    real(r8), pointer :: ulrad(:)         ! upward longwave radiation above the canopy (W/m**2)
    real(r8), pointer :: cgrnds(:)        ! deriv, of soil sensible heat flux wrt soil temp (W/m**2/K)
    real(r8), pointer :: cgrndl(:)        ! deriv of soil latent heat flux wrt soil temp (W/m**2/K)
    real(r8), pointer :: cgrnd(:)         ! deriv. of soil energy flux wrt to soil temp (W/m**2/K)
    real(r8), pointer :: taux(:)          ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy(:)          ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_u(:) ! urban total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:) ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_veg(:) ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: t_ref2m(:)       ! 2 m height surface air temperature (K)
    real(r8), pointer :: q_ref2m(:)       ! 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: t_ref2m_u(:)     ! Urban 2 m height surface air temperature (K)
    real(r8), pointer :: t_veg(:)         ! vegetation temperature (K)
    real(r8), pointer :: ram1(:)          ! aerodynamical resistance (s/m)
    real(r8), pointer :: rootr(:,:)       ! effective fraction of roots in each soil layer
    real(r8), pointer :: psnsun(:)        ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha(:)        ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: t_building(:)    ! internal building temperature (K)
    real(r8), pointer :: rh_ref2m(:)      ! 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m_u(:)    ! Urban 2 m height surface relative humidity (%)
    real(r8), pointer :: eflx_sh_snow(:)          ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_soil(:)          ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_h2osfc(:)        ! sensible heat flux from soil (W/m**2) [+ to atm]
!
!
! !OTHER LOCAL VARIABLES
!EOP
!
    character(len=*), parameter :: sub="UrbanFluxes"
    integer  :: fp,fc,fl,f,p,c,l,g,j,pi,i     ! indices
                                 
    real(r8) :: canyontop_wind(num_urbanl)    ! wind at canyon top (m/s) 
    real(r8) :: canyon_u_wind(num_urbanl)     ! u-component of wind speed inside canyon (m/s)
    real(r8) :: canyon_wind(num_urbanl)       ! net wind speed inside canyon (m/s)
    real(r8) :: canyon_resistance(num_urbanl) ! resistance to heat and moisture transfer from canyon road/walls to canyon air (s/m)

    real(r8) :: ur(lbl:ubl)        ! wind speed at reference height (m/s)
    real(r8) :: ustar(lbl:ubl)     ! friction velocity (m/s)
    real(r8) :: ramu(lbl:ubl)      ! aerodynamic resistance (s/m)
    real(r8) :: rahu(lbl:ubl)      ! thermal resistance (s/m)
    real(r8) :: rawu(lbl:ubl)      ! moisture resistance (s/m)
    real(r8) :: temp1(lbl:ubl)     ! relation for potential temperature profile
    real(r8) :: temp12m(lbl:ubl)   ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(lbl:ubl)     ! relation for specific humidity profile
    real(r8) :: temp22m(lbl:ubl)   ! relation for specific humidity profile applied at 2-m
    real(r8) :: thm_g(lbl:ubl)     ! intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv_g(lbl:ubl)     ! virtual potential temperature (K)
    real(r8) :: dth(lbl:ubl)       ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(lbl:ubl)       ! diff of humidity between ref. height and surface
    real(r8) :: zldis(lbl:ubl)     ! reference height "minus" zero displacement height (m)
    real(r8) :: um(lbl:ubl)        ! wind speed including the stablity effect (m/s)
    real(r8) :: obu(lbl:ubl)       ! Monin-Obukhov length (m)
    real(r8) :: taf_numer(lbl:ubl) ! numerator of taf equation (K m/s)
    real(r8) :: taf_denom(lbl:ubl) ! denominator of taf equation (m/s)
    real(r8) :: qaf_numer(lbl:ubl) ! numerator of qaf equation (kg m/kg s)
    real(r8) :: qaf_denom(lbl:ubl) ! denominator of qaf equation (m/s)
    real(r8) :: wtas(lbl:ubl)      ! sensible heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wtaq(lbl:ubl)      ! latent heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wts_sum(lbl:ubl)   ! sum of wtas, wtus_roof, wtus_road_perv, wtus_road_imperv, wtus_sunwall, wtus_shadewall
    real(r8) :: wtq_sum(lbl:ubl)   ! sum of wtaq, wtuq_roof, wtuq_road_perv, wtuq_road_imperv, wtuq_sunwall, wtuq_shadewall
    real(r8) :: beta(lbl:ubl)      ! coefficient of convective velocity
    real(r8) :: zii(lbl:ubl)       ! convective boundary layer height (m)

    real(r8) :: fm(lbl:ubl)        ! needed for BGC only to diagnose 10m wind speed

    real(r8) :: wtus(lbc:ubc)      ! sensible heat conductance for urban columns (m/s)
    real(r8) :: wtuq(lbc:ubc)      ! latent heat conductance for urban columns (m/s)

    integer  :: iter               ! iteration index
    real(r8) :: dthv               ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: tstar              ! temperature scaling parameter
    real(r8) :: qstar              ! moisture scaling parameter
    real(r8) :: thvstar            ! virtual potential temperature scaling parameter
    real(r8) :: wtus_roof(lbl:ubl)        ! sensible heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtuq_roof(lbl:ubl)        ! latent heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtus_road_perv(lbl:ubl)   ! sensible heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_perv(lbl:ubl)   ! latent heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtus_road_imperv(lbl:ubl) ! sensible heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_imperv(lbl:ubl) ! latent heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtus_sunwall(lbl:ubl)     ! sensible heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtuq_sunwall(lbl:ubl)     ! latent heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtus_shadewall(lbl:ubl)   ! sensible heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: wtuq_shadewall(lbl:ubl)   ! latent heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: t_sunwall_innerl(lbl:ubl)   ! temperature of inner layer of sunwall (K)
    real(r8) :: t_shadewall_innerl(lbl:ubl) ! temperature of inner layer of shadewall (K)
    real(r8) :: t_roof_innerl(lbl:ubl)      ! temperature of inner layer of roof (K)
    real(r8) :: lngth_roof                  ! length of roof (m)
    real(r8) :: wc                     ! convective velocity (m/s)
    real(r8) :: zeta                   ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: eflx_sh_grnd_scale(lbp:ubp)  ! scaled sensible heat flux from ground (W/m**2) [+ to atm] 
    real(r8) :: qflx_evap_soi_scale(lbp:ubp) ! scaled soil evaporation (mm H2O/s) (+ = to atm) 
    real(r8) :: eflx_wasteheat_roof(lbl:ubl) ! sensible heat flux from urban heating/cooling sources of waste heat for roof (W/m**2)
    real(r8) :: eflx_wasteheat_sunwall(lbl:ubl) ! sensible heat flux from urban heating/cooling sources of waste heat for sunwall (W/m**2)
    real(r8) :: eflx_wasteheat_shadewall(lbl:ubl) ! sensible heat flux from urban heating/cooling sources of waste heat for shadewall (W/m**2)
    real(r8) :: eflx_heat_from_ac_roof(lbl:ubl) ! sensible heat flux put back into canyon due to heat removal by AC for roof (W/m**2)
    real(r8) :: eflx_heat_from_ac_sunwall(lbl:ubl) ! sensible heat flux put back into canyon due to heat removal by AC for sunwall (W/m**2)
    real(r8) :: eflx_heat_from_ac_shadewall(lbl:ubl) ! sensible heat flux put back into canyon due to heat removal by AC for shadewall (W/m**2)
    real(r8) :: eflx(lbl:ubl)                     ! total sensible heat flux for error check (W/m**2)
    real(r8) :: qflx(lbl:ubl)                     ! total water vapor flux for error check (kg/m**2/s)
    real(r8) :: eflx_scale(lbl:ubl)               ! sum of scaled sensible heat fluxes for urban columns for error check (W/m**2)
    real(r8) :: qflx_scale(lbl:ubl)               ! sum of scaled water vapor fluxes for urban columns for error check (kg/m**2/s)
    real(r8) :: eflx_err(lbl:ubl)                 ! sensible heat flux error (W/m**2)
    real(r8) :: qflx_err(lbl:ubl)                 ! water vapor flux error (kg/m**2/s)
    real(r8) :: fwet_roof                         ! fraction of roof surface that is wet (-)
    real(r8) :: fwet_road_imperv                  ! fraction of impervious road surface that is wet (-)

    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: local_secp1(lbl:ubl)   ! seconds into current date in local time (sec)
    real(r8) :: dtime                  ! land model time step (sec)
    integer  :: year,month,day,secs    ! calendar info for current time step
    logical  :: found                  ! flag in search loop
    integer  :: indexl                 ! index of first found in search loop
    integer  :: nstep                  ! time step number
    real(r8) :: z_d_town_loc(lbl:ubl)  ! temporary copy
    real(r8) :: z_0_town_loc(lbl:ubl)  ! temporary copy
    real(r8), parameter :: lapse_rate = 0.0098_r8     ! Dry adiabatic lapse rate (K/m)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m

!-----------------------------------------------------------------------

    ! Assign pointers into module urban clumps

    if ( num_urbanl > 0 )then
       ht_roof            => urban_clump(nc)%ht_roof
       wtlunit_roof       => urban_clump(nc)%wtlunit_roof
       canyon_hwr         => urban_clump(nc)%canyon_hwr
       wtroad_perv        => urban_clump(nc)%wtroad_perv
       wind_hgt_canyon    => urban_clump(nc)%wind_hgt_canyon
    end if

    ! Assign local pointers to multi-level derived type members (gridcell level)

    forc_t     => clm_a2l%forc_t
    forc_th    => clm_a2l%forc_th
    forc_u     => clm_a2l%forc_u
    forc_v     => clm_a2l%forc_v
    forc_rho   => clm_a2l%forc_rho
    forc_q     => clm_a2l%forc_q
    forc_pbot  => clm_a2l%forc_pbot
    londeg     => clm3%g%londeg

    ! Assign local pointers to derived type members (landunit level)

    pfti                => clm3%g%l%pfti
    pftf                => clm3%g%l%pftf
    coli                => clm3%g%l%coli
    colf                => clm3%g%l%colf
    lgridcell           => clm3%g%l%gridcell
    z_0_town            => clm3%g%l%z_0_town
    z_d_town            => clm3%g%l%z_d_town
    taf                 => clm3%g%l%lps%taf
    qaf                 => clm3%g%l%lps%qaf
    npfts               => clm3%g%l%npfts
    eflx_traffic        => clm3%g%l%lef%eflx_traffic
    eflx_traffic_factor => clm3%g%l%lef%eflx_traffic_factor
    eflx_wasteheat      => clm3%g%l%lef%eflx_wasteheat
    eflx_heat_from_ac   => clm3%g%l%lef%eflx_heat_from_ac
    t_building          => clm3%g%l%lps%t_building

    ! Assign local pointers to derived type members (column level)

    ctype              => clm3%g%l%c%itype
    t_grnd             => clm3%g%l%c%ces%t_grnd
    qg                 => clm3%g%l%c%cws%qg
    htvp               => clm3%g%l%c%cps%htvp
    dqgdT              => clm3%g%l%c%cws%dqgdT
    t_soisno           => clm3%g%l%c%ces%t_soisno
    eflx_urban_ac      => clm3%g%l%c%cef%eflx_urban_ac
    eflx_urban_heat    => clm3%g%l%c%cef%eflx_urban_heat
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    frac_sno           => clm3%g%l%c%cps%frac_sno
    snow_depth             => clm3%g%l%c%cps%snow_depth
    h2osno             => clm3%g%l%c%cws%h2osno
    snl                => clm3%g%l%c%cps%snl
    rootr_road_perv    => clm3%g%l%c%cps%rootr_road_perv
    soilalpha_u        => clm3%g%l%c%cws%soilalpha_u

    ! Assign local pointers to derived type members (pft level)

    pgridcell      => clm3%g%l%c%p%gridcell
    pcolumn        => clm3%g%l%c%p%column
    plandunit      => clm3%g%l%c%p%landunit
    ram1           => clm3%g%l%c%p%pps%ram1
    dlrad          => clm3%g%l%c%p%pef%dlrad
    ulrad          => clm3%g%l%c%p%pef%ulrad
    cgrnds         => clm3%g%l%c%p%pef%cgrnds
    cgrndl         => clm3%g%l%c%p%pef%cgrndl
    cgrnd          => clm3%g%l%c%p%pef%cgrnd
    taux           => clm3%g%l%c%p%pmf%taux
    tauy           => clm3%g%l%c%p%pmf%tauy
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_sh_tot_u  => clm3%g%l%c%p%pef%eflx_sh_tot_u
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
    qflx_evap_veg  => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_evap_tot  => clm3%g%l%c%p%pwf%qflx_evap_tot
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    q_ref2m        => clm3%g%l%c%p%pes%q_ref2m
    t_ref2m_u      => clm3%g%l%c%p%pes%t_ref2m_u
    t_veg          => clm3%g%l%c%p%pes%t_veg
    rootr          => clm3%g%l%c%p%pps%rootr
    psnsun         => clm3%g%l%c%p%pcf%psnsun
    psnsha         => clm3%g%l%c%p%pcf%psnsha
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
    forc_hgt_t_pft => clm3%g%l%c%p%pps%forc_hgt_t_pft
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
    forc_hgt_t_pft => clm3%g%l%c%p%pps%forc_hgt_t_pft
    rh_ref2m => clm3%g%l%c%p%pes%rh_ref2m
    rh_ref2m_u     => clm3%g%l%c%p%pes%rh_ref2m_u

    eflx_sh_snow   => clm3%g%l%c%p%pef%eflx_sh_snow
    eflx_sh_soil   => clm3%g%l%c%p%pef%eflx_sh_soil
    eflx_sh_h2osfc => clm3%g%l%c%p%pef%eflx_sh_h2osfc
    ! Define fields that appear on the restart file for non-urban landunits 

    do fl = 1,num_nourbanl
       l = filter_nourbanl(fl)
       taf(l) = spval
       qaf(l) = spval
    end do

    ! Get time step
    nstep = get_nstep()

    ! Set constants (same as in Biogeophysics1Mod)
    beta(:) = 1._r8             ! Should be set to the same values as in Biogeophysics1Mod
    zii(:)  = 1000._r8          ! Should be set to the same values as in Biogeophysics1Mod

    ! Get current date
    dtime = get_step_size()
    call get_curr_date (year, month, day, secs)
    
    ! Compute canyontop wind using Masson (2000)

    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)

       local_secp1(l)        = secs + nint((londeg(g)/degpsec)/dtime)*dtime
       local_secp1(l)        = mod(local_secp1(l),isecspday)

       ! Error checks

       if (ht_roof(fl) - z_d_town(l) <= z_0_town(l)) then
          write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
          write (iulog,*) 'h_r - z_d <= z_0'
          write (iulog,*) 'ht_roof, z_d_town, z_0_town: ', ht_roof(fl), z_d_town(l), &
                       z_0_town(l)
          write (iulog,*) 'clm model is stopping'
          call endrun()
       end if
       if (forc_hgt_u_pft(pfti(l)) - z_d_town(l) <= z_0_town(l)) then
          write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
          write (iulog,*) 'h_u - z_d <= z_0'
          write (iulog,*) 'forc_hgt_u_pft, z_d_town, z_0_town: ', forc_hgt_u_pft(pfti(l)), z_d_town(l), &
                       z_0_town(l)
          write (iulog,*) 'clm model is stopping'
          call endrun()
       end if

       ! Magnitude of atmospheric wind

       ur(l) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))

       ! Canyon top wind

       canyontop_wind(fl) = ur(l) * &
                           log( (ht_roof(fl)-z_d_town(l)) / z_0_town(l) ) / &
                           log( (forc_hgt_u_pft(pfti(l))-z_d_town(l)) / z_0_town(l) )

       ! U component of canyon wind 

       if (canyon_hwr(fl) < 0.5_r8) then  ! isolated roughness flow
         canyon_u_wind(fl) = canyontop_wind(fl) * exp( -0.5_r8*canyon_hwr(fl)* &
                            (1._r8-(wind_hgt_canyon(fl)/ht_roof(fl))) )
       else if (canyon_hwr(fl) < 1.0_r8) then ! wake interference flow
         canyon_u_wind(fl) = canyontop_wind(fl) * (1._r8+2._r8*(2._r8/rpi - 1._r8)* &
                            (ht_roof(fl)/(ht_roof(fl)/canyon_hwr(fl)) - 0.5_r8)) * &
                            exp(-0.5_r8*canyon_hwr(fl)*(1._r8-(wind_hgt_canyon(fl)/ht_roof(fl))))
       else  ! skimming flow
         canyon_u_wind(fl) = canyontop_wind(fl) * (2._r8/rpi) * &
                            exp(-0.5_r8*canyon_hwr(fl)*(1._r8-(wind_hgt_canyon(fl)/ht_roof(fl))))
       end if

    end do

! Compute fluxes - Follows CLM approach for bare soils (Oleson et al 2004)

    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)

       thm_g(l) = forc_t(g) + lapse_rate*forc_hgt_t_pft(pfti(l))
       thv_g(l) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))
       dth(l)   = thm_g(l)-taf(l)
       dqh(l)   = forc_q(g)-qaf(l)
       dthv     = dth(l)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(l)
       zldis(l) = forc_hgt_u_pft(pfti(l)) - z_d_town(l)

       ! Initialize Monin-Obukhov length and wind speed including convective velocity

       call MoninObukIni(ur(l), thv_g(l), dthv, zldis(l), z_0_town(l), um(l), obu(l))

    end do

    ! Initialize conductances
    wtus_roof(:)        = 0._r8
    wtus_road_perv(:)   = 0._r8
    wtus_road_imperv(:) = 0._r8
    wtus_sunwall(:)     = 0._r8
    wtus_shadewall(:)   = 0._r8
    wtuq_roof(:)        = 0._r8
    wtuq_road_perv(:)   = 0._r8
    wtuq_road_imperv(:) = 0._r8
    wtuq_sunwall(:)     = 0._r8
    wtuq_shadewall(:)   = 0._r8

   ! Make copies so that array sections are not passed in function calls to friction velocity

    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       z_d_town_loc(l) = z_d_town(l)
       z_0_town_loc(l) = z_0_town(l)
    end do

    ! Start stability iteration

    do iter = 1,niters

     ! Get friction velocity, relation for potential
     ! temperature and humidity profiles of surface boundary layer.

      if (num_urbanl .gt. 0) then
         call FrictionVelocity(lbl, ubl, num_urbanl, filter_urbanl, &
                               z_d_town_loc, z_0_town_loc, z_0_town_loc, z_0_town_loc, &
                               obu, iter, ur, um, ustar, &
                               temp1, temp2, temp12m, temp22m, fm, landunit_index=.true.)
      end if

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         ! Determine aerodynamic resistance to fluxes from urban canopy air to
         ! atmosphere

         ramu(l) = 1._r8/(ustar(l)*ustar(l)/um(l))
         rahu(l) = 1._r8/(temp1(l)*ustar(l))
         rawu(l) = 1._r8/(temp2(l)*ustar(l))

         ! Determine magnitude of canyon wind by using horizontal wind determined
         ! previously and vertical wind from friction velocity (Masson 2000)

         canyon_wind(fl) = sqrt(canyon_u_wind(fl)**2._r8 + ustar(l)**2._r8)

         ! Determine canyon_resistance (currently this single resistance determines the
         ! resistance from urban surfaces (roof, pervious and impervious road, sunlit and
         ! shaded walls) to urban canopy air, since it is only dependent on wind speed
         ! Also from Masson 2000.

         canyon_resistance(fl) = cpair * forc_rho(g) / (11.8_r8 + 4.2_r8*canyon_wind(fl))

      end do

      ! This is the first term in the equation solutions for urban canopy air temperature
      ! and specific humidity (numerator) and is a landunit quantity
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         taf_numer(l) = thm_g(l)/rahu(l)
         taf_denom(l) = 1._r8/rahu(l)
         qaf_numer(l) = forc_q(g)/rawu(l)
         qaf_denom(l) = 1._r8/rawu(l)

         ! First term needed for derivative of heat fluxes
         wtas(l) = 1._r8/rahu(l)
         wtaq(l) = 1._r8/rawu(l)

      end do


      ! Gather other terms for other urban columns for numerator and denominator of
      ! equations for urban canopy air temperature and specific humidity

      do pi = 1,maxpatch_urb
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if ( pi <= npfts(l) ) then
               c = coli(l) + pi - 1

               if (ctype(c) == icol_roof) then

                 ! scaled sensible heat conductance
                 wtus(c) = wtlunit_roof(fl)/canyon_resistance(fl)
                 ! unscaled sensible heat conductance
                 wtus_roof(l) = 1._r8/canyon_resistance(fl)

                 if (snow_depth(c) > 0._r8) then
                   fwet_roof = min(snow_depth(c)/0.05_r8, 1._r8)
                 else
                   fwet_roof = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
                   fwet_roof = min(fwet_roof,1._r8)
                 end if
                 if (qaf(l) > qg(c)) then 
                   fwet_roof = 1._r8
                 end if
                 ! scaled latent heat conductance
                 wtuq(c)      = fwet_roof*(wtlunit_roof(fl)/canyon_resistance(fl))
                 ! unscaled latent heat conductance
                 wtuq_roof(l) = fwet_roof*(1._r8/canyon_resistance(fl))

                 ! wasteheat from heating/cooling
                 if (trim(urban_hac) == urban_wasteheat_on) then
                    eflx_wasteheat_roof(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                                             ht_wasteheat_factor * eflx_urban_heat(c)
                 else
                    eflx_wasteheat_roof(l) = 0._r8
                 end if

                 ! If air conditioning on, always replace heat removed with heat into canyon
                 if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                    eflx_heat_from_ac_roof(l) = abs(eflx_urban_ac(c))
                 else
                    eflx_heat_from_ac_roof(l) = 0._r8
                 end if

               else if (ctype(c) == icol_road_perv) then

                 ! scaled sensible heat conductance
                 wtus(c) = wtroad_perv(fl)*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled sensible heat conductance
                 if (wtroad_perv(fl) > 0._r8) then
                   wtus_road_perv(l) = 1._r8/canyon_resistance(fl)
                 else
                   wtus_road_perv(l) = 0._r8
                 end if

                 ! scaled latent heat conductance
                 wtuq(c) = wtroad_perv(fl)*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled latent heat conductance
                 if (wtroad_perv(fl) > 0._r8) then
                   wtuq_road_perv(l) = 1._r8/canyon_resistance(fl)
                 else
                   wtuq_road_perv(l) = 0._r8
                 end if

               else if (ctype(c) == icol_road_imperv) then

                 ! scaled sensible heat conductance
                 wtus(c) = (1._r8-wtroad_perv(fl))*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled sensible heat conductance
                 if ((1._r8-wtroad_perv(fl)) > 0._r8) then
                   wtus_road_imperv(l) = 1._r8/canyon_resistance(fl)
                 else
                   wtus_road_imperv(l) = 0._r8
                 end if

                 if (snow_depth(c) > 0._r8) then
                   fwet_road_imperv = min(snow_depth(c)/0.05_r8, 1._r8)
                 else
                   fwet_road_imperv = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
                   fwet_road_imperv = min(fwet_road_imperv,1._r8)
                 end if
                 if (qaf(l) > qg(c)) then 
                   fwet_road_imperv = 1._r8
                 end if
                 ! scaled latent heat conductance
                 wtuq(c) = fwet_road_imperv*(1._r8-wtroad_perv(fl))*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled latent heat conductance
                 if ((1._r8-wtroad_perv(fl)) > 0._r8) then
                   wtuq_road_imperv(l) = fwet_road_imperv*(1._r8/canyon_resistance(fl))
                 else
                   wtuq_road_imperv(l) = 0._r8
                 end if

               else if (ctype(c) == icol_sunwall) then

                 ! scaled sensible heat conductance
                 wtus(c)         = canyon_hwr(fl)*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled sensible heat conductance
                 wtus_sunwall(l) = 1._r8/canyon_resistance(fl)

                 ! scaled latent heat conductance
                 wtuq(c)         = 0._r8
                 ! unscaled latent heat conductance
                 wtuq_sunwall(l) = 0._r8

                 ! wasteheat from heating/cooling
                 if (trim(urban_hac) == urban_wasteheat_on) then
                    eflx_wasteheat_sunwall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                                                ht_wasteheat_factor * eflx_urban_heat(c)
                 else
                    eflx_wasteheat_sunwall(l) = 0._r8
                 end if

                 ! If air conditioning on, always replace heat removed with heat into canyon
                 if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                    eflx_heat_from_ac_sunwall(l) = abs(eflx_urban_ac(c))
                 else
                    eflx_heat_from_ac_sunwall(l) = 0._r8
                 end if

               else if (ctype(c) == icol_shadewall) then

                 ! scaled sensible heat conductance
                 wtus(c) = canyon_hwr(fl)*(1._r8-wtlunit_roof(fl))/canyon_resistance(fl)
                 ! unscaled sensible heat conductance
                 wtus_shadewall(l) = 1._r8/canyon_resistance(fl)

                 ! scaled latent heat conductance
                 wtuq(c)           = 0._r8
                 ! unscaled latent heat conductance
                 wtuq_shadewall(l) = 0._r8

                 ! wasteheat from heating/cooling
                 if (trim(urban_hac) == urban_wasteheat_on) then
                    eflx_wasteheat_shadewall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                                                  ht_wasteheat_factor * eflx_urban_heat(c)
                 else
                    eflx_wasteheat_shadewall(l) = 0._r8
                 end if

                 ! If air conditioning on, always replace heat removed with heat into canyon
                 if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                    eflx_heat_from_ac_shadewall(l) = abs(eflx_urban_ac(c))
                 else
                    eflx_heat_from_ac_shadewall(l) = 0._r8
                 end if
               else
                 write(iulog,*) 'c, ctype, pi = ', c, ctype(c), pi
                 write(iulog,*) 'Column indices for: shadewall, sunwall, road_imperv, road_perv, roof: '
                 write(iulog,*) icol_shadewall, icol_sunwall, icol_road_imperv, icol_road_perv, icol_roof
                 call endrun( sub//':: ERROR, ctype out of range' )
               end if

               taf_numer(l) = taf_numer(l) + t_grnd(c)*wtus(c)
               taf_denom(l) = taf_denom(l) + wtus(c)
               qaf_numer(l) = qaf_numer(l) + qg(c)*wtuq(c)
               qaf_denom(l) = qaf_denom(l) + wtuq(c)

            end if
         end do
      end do

      ! Calculate new urban canopy air temperature and specific humidity

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         ! Total waste heat and heat from AC is sum of heat for walls and roofs
         ! accounting for different surface areas
         eflx_wasteheat(l) = wtlunit_roof(fl)*eflx_wasteheat_roof(l) + &
            (1._r8-wtlunit_roof(fl))*(canyon_hwr(fl)*(eflx_wasteheat_sunwall(l) + &
            eflx_wasteheat_shadewall(l)))

         ! Limit wasteheat to ensure that we don't get any unrealistically strong 
         ! positive feedbacks due to AC in a warmer climate
         eflx_wasteheat(l) = min(eflx_wasteheat(l),wasteheat_limit)

         eflx_heat_from_ac(l) = wtlunit_roof(fl)*eflx_heat_from_ac_roof(l) + &
            (1._r8-wtlunit_roof(fl))*(canyon_hwr(fl)*(eflx_heat_from_ac_sunwall(l) + &
            eflx_heat_from_ac_shadewall(l)))

         ! Calculate traffic heat flux
         ! Only comes from impervious road
         eflx_traffic(l) = (1._r8-wtlunit_roof(fl))*(1._r8-wtroad_perv(fl))* &
                           eflx_traffic_factor(l)

         taf(l) = taf_numer(l)/taf_denom(l)
         qaf(l) = qaf_numer(l)/qaf_denom(l)

         wts_sum(l) = wtas(l) + wtus_roof(l) + wtus_road_perv(l) + &
                      wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)

         wtq_sum(l) = wtaq(l) + wtuq_roof(l) + wtuq_road_perv(l) + &
                      wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)

      end do

      ! This section of code is not required if niters = 1
      ! Determine stability using new taf and qaf
      ! TODO: Some of these constants replicate what is in FrictionVelocity and BareGround fluxes should consildate. EBK
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         dth(l) = thm_g(l)-taf(l)
         dqh(l) = forc_q(g)-qaf(l)
         tstar = temp1(l)*dth(l)
         qstar = temp2(l)*dqh(l)
         thvstar = tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
         zeta = zldis(l)*vkc*grav*thvstar/(ustar(l)**2*thv_g(l))

         if (zeta >= 0._r8) then                   !stable
            zeta = min(2._r8,max(zeta,0.01_r8))
            um(l) = max(ur(l),0.1_r8)
         else                                      !unstable
            zeta = max(-100._r8,min(zeta,-0.01_r8))
            wc = beta(l)*(-grav*ustar(l)*thvstar*zii(l)/thv_g(l))**0.333_r8
            um(l) = sqrt(ur(l)*ur(l) + wc*wc)
         end if

         obu(l) = zldis(l)/zeta
      end do

    end do   ! end iteration

! Determine fluxes from canyon surfaces

    do f = 1, num_urbanp

       p = filter_urbanp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ram1(p) = ramu(l)  !pass value to global variable

       ! Upward and downward canopy longwave are zero

       ulrad(p)  = 0._r8
       dlrad(p)  = 0._r8

       ! Derivative of sensible and latent heat fluxes with respect to 
       ! ground temperature

       if (ctype(c) == icol_roof) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_road_perv(l) +  &
                     wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_roof(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_road_perv(l) +  &
                     wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_roof(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_road_perv) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_road_perv(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                     wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_road_perv(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_road_imperv) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_road_imperv(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                     wtuq_road_perv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_road_imperv(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_sunwall) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_road_imperv(l) + wtus_shadewall(l)) * &
                     (wtus_sunwall(l)/wts_sum(l))
         cgrndl(p) = 0._r8
       else if (ctype(c) == icol_shadewall) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_road_imperv(l) + wtus_sunwall(l)) * &
                     (wtus_shadewall(l)/wts_sum(l))
         cgrndl(p) = 0._r8
       end if
       cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

       ! Surface fluxes of momentum, sensible and latent heat

       taux(p)          = -forc_rho(g)*forc_u(g)/ramu(l)
       tauy(p)          = -forc_rho(g)*forc_v(g)/ramu(l)

       ! Use new canopy air temperature
       dth(l) = taf(l) - t_grnd(c)

       if (ctype(c) == icol_roof) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_roof(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_road_perv) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_perv(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_road_imperv) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_imperv(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_sunwall) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_sunwall(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_shadewall) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_shadewall(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       end if

       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       eflx_sh_tot_u(p) = eflx_sh_tot(p)

       dqh(l) = qaf(l) - qg(c)

       if (ctype(c) == icol_roof) then
         qflx_evap_soi(p) = -forc_rho(g)*wtuq_roof(l)*dqh(l)
       else if (ctype(c) == icol_road_perv) then
         ! Evaporation assigned to soil term if dew or snow
         ! or if no liquid water available in soil column
         if (dqh(l) > 0._r8 .or. frac_sno(c) > 0._r8 .or. soilalpha_u(c) .le. 0._r8) then
            qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_perv(l)*dqh(l)
            qflx_tran_veg(p) = 0._r8
         ! Otherwise, evaporation assigned to transpiration term
         else
            qflx_evap_soi(p) = 0._r8
            qflx_tran_veg(p) = -forc_rho(g)*wtuq_road_perv(l)*dqh(l)
         end if
         qflx_evap_veg(p) = qflx_tran_veg(p)
       else if (ctype(c) == icol_road_imperv) then
         qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_imperv(l)*dqh(l)
       else if (ctype(c) == icol_sunwall) then
         qflx_evap_soi(p) = 0._r8
       else if (ctype(c) == icol_shadewall) then
         qflx_evap_soi(p) = 0._r8
       end if

       ! SCALED sensible and latent heat flux for error check
       eflx_sh_grnd_scale(p)  = -forc_rho(g)*cpair*wtus(c)*dth(l)
       qflx_evap_soi_scale(p) = -forc_rho(g)*wtuq(c)*dqh(l)

    end do

    ! Check to see that total sensible and latent heat equal the sum of
    ! the scaled heat fluxes above
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)
       eflx(l)       = -(forc_rho(g)*cpair/rahu(l))*(thm_g(l) - taf(l))
       qflx(l)       = -(forc_rho(g)/rawu(l))*(forc_q(g) - qaf(l))
       eflx_scale(l) = sum(eflx_sh_grnd_scale(pfti(l):pftf(l)))
       qflx_scale(l) = sum(qflx_evap_soi_scale(pfti(l):pftf(l)))
       eflx_err(l)   = eflx_scale(l) - eflx(l)
       qflx_err(l)   = qflx_scale(l) - qflx(l)
    end do 

    found = .false.
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       if (abs(eflx_err(l)) > 0.01_r8) then
          found = .true.
          indexl = l
          exit
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  Total sensible heat does not equal sum of scaled heat fluxes for urban columns ',&
            ' nstep = ',nstep,' indexl= ',indexl,' eflx_err= ',eflx_err(indexl)
       if (abs(eflx_err(indexl)) > .01_r8) then
          write(iulog,*)'clm model is stopping - error is greater than .01 W/m**2'
          write(iulog,*)'eflx_scale    = ',eflx_scale(indexl)
          write(iulog,*)'eflx_sh_grnd_scale: ',eflx_sh_grnd_scale(pfti(indexl):pftf(indexl))
          write(iulog,*)'eflx          = ',eflx(indexl)
          call endrun
       end if
    end if

    found = .false.
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       ! 4.e-9 kg/m**2/s = 0.01 W/m**2
       if (abs(qflx_err(l)) > 4.e-9_r8) then
          found = .true.
          indexl = l
          exit
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  Total water vapor flux does not equal sum of scaled water vapor fluxes for urban columns ',&
            ' nstep = ',nstep,' indexl= ',indexl,' qflx_err= ',qflx_err(indexl)
       if (abs(qflx_err(indexl)) > 4.e-9_r8) then
          write(iulog,*)'clm model is stopping - error is greater than 4.e-9 kg/m**2/s'
          write(iulog,*)'qflx_scale    = ',qflx_scale(indexl)
          write(iulog,*)'qflx          = ',qflx(indexl)
          call endrun
       end if
    end if

    ! Gather terms required to determine internal building temperature

    do pi = 1,maxpatch_urb
       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          if ( pi <= npfts(l) ) then
             c = coli(l) + pi - 1

             if (ctype(c) == icol_roof) then
                t_roof_innerl(l) = t_soisno(c,nlevurb)
             else if (ctype(c) == icol_sunwall) then
                t_sunwall_innerl(l) = t_soisno(c,nlevurb)
             else if (ctype(c) == icol_shadewall) then
                t_shadewall_innerl(l) = t_soisno(c,nlevurb)
             end if

          end if
       end do
    end do

    ! Calculate internal building temperature
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
     
       lngth_roof = (ht_roof(fl)/canyon_hwr(fl))*wtlunit_roof(fl)/(1._r8-wtlunit_roof(fl))
       t_building(l) = (ht_roof(fl)*(t_shadewall_innerl(l) + t_sunwall_innerl(l)) &
                        +lngth_roof*t_roof_innerl(l))/(2._r8*ht_roof(fl)+lngth_roof)
    end do

    ! No roots for urban except for pervious road

    do j = 1, nlevgrnd
      do f = 1, num_urbanp
         p = filter_urbanp(f)
         c = pcolumn(p)
         if (ctype(c) == icol_road_perv) then
            rootr(p,j) = rootr_road_perv(c,j)
         else
            rootr(p,j) = 0._r8
         end if
       end do
    end do

    do f = 1, num_urbanp

       p = filter_urbanp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Use urban canopy air temperature and specific humidity to represent 
       ! 2-m temperature and humidity

       t_ref2m(p) = taf(l)
       q_ref2m(p) = qaf(l)
       t_ref2m_u(p) = taf(l)

       ! 2 m height relative humidity

       call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
       rh_ref2m_u(p) = rh_ref2m(p)

       ! Variables needed by history tape

       t_veg(p) = forc_t(g)

       ! Add the following to avoid NaN

       psnsun(p)                   = 0._r8
       psnsha(p)                   = 0._r8

    end do

  end subroutine UrbanFluxes

end module UrbanMod

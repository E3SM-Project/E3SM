module UrbanAlbedoMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use elm_varpar        , only : numrad
  use elm_varcon        , only : isecspday, degpsec, namel
  use elm_varctl        , only : iulog
  use UrbanParamsType   , only : urbanparams_type
  use SolarAbsorbedType , only : solarabs_type
  use SurfaceAlbedoType , only : surfalb_type
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_ws
  use VegetationType    , only : veg_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanAlbedo       ! Urban physics - albedos
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: SnowAlbedo       ! Snow albedos
  private :: incident_direct  ! Direct beam solar rad incident on walls and road in urban canyon
  private :: incident_diffuse ! Diffuse solar rad incident on walls and road in urban canyon
  private :: net_solar        ! Solar radiation absorbed by road and both walls in urban canyon
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanAlbedo (num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, solarabs_vars, surfalb_vars)
    !
    ! !DESCRIPTION:
    ! Determine urban landunit component albedos
    !
    ! Note that this is called with the "inactive_and_active" version of the filters, because
    ! the variables computed here are needed over inactive points that might later become
    ! active (due to landuse change). Thus, this routine cannot depend on variables that are
    ! only computed over active points.
    !
    ! !USES:
    !$acc routine seq
    use elm_varcon    , only : sb
    use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon , only : icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    !type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_urbanl       ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in)    :: num_urbanc       ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:) ! urban column filter
    integer                , intent(in)    :: num_urbanp       ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:) ! urban pft filter
    type(urbanparams_type) , intent(inout) :: urbanparams_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(surfalb_type)     , intent(inout) :: surfalb_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fl,fp,fc,g,l,p,c,ib                          ! indices
    integer  :: ic                                           ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: num_solar                                    ! counter
    real(r8) :: coszen              ! cosine solar zenith angle for next time step (landunit)
    real(r8) :: zen                 ! solar zenith angle (radians)
    real(r8) :: sdir                ! direct beam solar radiation on horizontal surface
    real(r8) :: sdif                ! diffuse solar radiation on horizontal surface
    real(r8) :: sdir_road           ! direct beam solar radiation incident on road
    real(r8) :: sdif_road           ! diffuse solar radiation incident on road
    real(r8) :: sdir_sunwall        ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdif_sunwall        ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdir_shadewall      ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_shadewall      ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: albsnd_roof         ! snow albedo for roof (direct)
    real(r8) :: albsni_roof         ! snow albedo for roof (diffuse)
    real(r8) :: albsnd_improad      ! snow albedo for impervious road (direct)
    real(r8) :: albsni_improad      ! snow albedo for impervious road (diffuse)
    real(r8) :: albsnd_perroad      ! snow albedo for pervious road (direct)
    real(r8) :: albsni_perroad      ! snow albedo for pervious road (diffuse)
    real(r8) :: alb_roof_dir_s      ! direct roof albedo with snow effects
    real(r8) :: alb_roof_dif_s      ! diffuse roof albedo with snow effects
    real(r8) :: alb_improad_dir_s   ! direct impervious road albedo with snow effects
    real(r8) :: alb_perroad_dir_s   ! direct pervious road albedo with snow effects
    real(r8) :: alb_improad_dif_s   ! diffuse impervious road albedo with snow effects
    real(r8) :: alb_perroad_dif_s   ! diffuse pervious road albedo with snow effects
    real(r8) :: sref_roof_dir      ! direct  solar reflected by roof per unit ground area per unit incident flux
    real(r8) :: sref_roof_dif      ! diffuse solar reflected by roof per unit ground area per unit incident flux
    real(r8) :: sref_sunwall_dir   ! direct  solar reflected by sunwall per unit wall area per unit incident flux
    real(r8) :: sref_sunwall_dif   ! diffuse solar reflected by sunwall per unit wall area per unit incident flux
    real(r8) :: sref_shadewall_dir ! direct  solar reflected by shadewall per unit wall area per unit incident flux
    real(r8) :: sref_shadewall_dif ! diffuse solar reflected by shadewall per unit wall area per unit incident flux
    real(r8) :: sref_improad_dir   ! direct  solar reflected by impervious road per unit ground area per unit incident flux
    real(r8) :: sref_improad_dif   ! diffuse solar reflected by impervious road per unit ground area per unit incident flux
    real(r8) :: sref_perroad_dir   ! direct  solar reflected by pervious road per unit ground area per unit incident flux
    real(r8) :: sref_perroad_dif   ! diffuse solar reflected by pervious road per unit ground area per unit incident flux
    !-----------------------------------------------------------------------

    associate(                                           &
         ctype              => col_pp%itype            , & ! Input:  [integer (:)    ]  column type
         coli               => lun_pp%coli             , & ! Input:  [integer (:)    ]  beginning column index for landunit
         canyon_hwr         => lun_pp%canyon_hwr       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width
         wtroad_perv        => lun_pp%wtroad_perv      , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road

         frac_sno           => col_ws%frac_sno         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)

         alb_roof_dir       => urbanparams_vars%alb_roof_dir        , & ! Output: [real(r8) (:,:) ]  direct roof albedo
         alb_roof_dif       => urbanparams_vars%alb_roof_dif        , & ! Output: [real(r8) (:,:) ]  diffuse roof albedo
         alb_improad_dir    => urbanparams_vars%alb_improad_dir     , & ! Output: [real(r8) (:,:) ]  direct impervious road albedo
         alb_improad_dif    => urbanparams_vars%alb_improad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse imprevious road albedo
         alb_perroad_dir    => urbanparams_vars%alb_perroad_dir     , & ! Output: [real(r8) (:,:) ]  direct pervious road albedo
         alb_perroad_dif    => urbanparams_vars%alb_perroad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse pervious road albedo
         alb_wall_dir       => urbanparams_vars%alb_wall_dir        , & ! Output: [real(r8) (:,:) ]  direct wall albedo
         alb_wall_dif       => urbanparams_vars%alb_wall_dif        , & ! Output: [real(r8) (:,:) ]  diffuse wall albedo

         sabs_roof_dir      => solarabs_vars%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      => solarabs_vars%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   => solarabs_vars%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   => solarabs_vars%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir => solarabs_vars%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif => solarabs_vars%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   => solarabs_vars%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   => solarabs_vars%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   => solarabs_vars%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   => solarabs_vars%sabs_perroad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux

         fabd               => surfalb_vars%fabd_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         fabd_sun           => surfalb_vars%fabd_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
         fabd_sha           => surfalb_vars%fabd_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
         fabi               => surfalb_vars%fabi_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
         fabi_sun           => surfalb_vars%fabi_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
         fabi_sha           => surfalb_vars%fabi_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
         ftdd               => surfalb_vars%ftdd_patch              , & ! Output:  [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
         ftid               => surfalb_vars%ftid_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
         ftii               => surfalb_vars%ftii_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flux
         albgrd             => surfalb_vars%albgrd_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (direct)
         albgri             => surfalb_vars%albgri_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (diffuse)
         albd               => surfalb_vars%albd_patch              , & ! Output  [real(r8) (:,:) ]  urban pft surface albedo (direct)
         albi               => surfalb_vars%albi_patch           & ! Output: [real(r8) (:,:) ]  urban pft surface albedo (diffuse)

         )

      ! ----------------------------------------------------------------------------
      ! Solar declination and cosine solar zenith angle and zenith angle for
      ! next time step
      ! ----------------------------------------------------------------------------

      ! do fl = 1,num_urbanl
      !    l = filter_urbanl(fl)
      !    g = lun_pp%gridcell(l)
      !    coszen(l) = surfalb_vars%coszen_col(coli(l)) ! Assumes coszen for each column are the same
      !    zen(l)    = acos(coszen(l))
      ! end do

       ! Initialize output because solar radiation only done if coszen > 0
       ! output is all global variables
      do ib = 1, numrad
         do fc = 1,num_urbanc
            c = filter_urbanc(fc)
            albgrd(c,ib) = 0._r8
            albgri(c,ib) = 0._r8
         end do

         do fp = 1,num_urbanp
            p = filter_urbanp(fp)
            l = veg_pp%landunit(p)
            c = veg_pp%column(p)
            coszen = surfalb_vars%coszen_col(coli(l))

            ! Initialize direct and diffuse albedo such that if the Sun is below
            ! the horizon, p2g scaling returns an albedo of 1.0.

            if (col_pp%itype(c) == icol_sunwall) then
               albd(p,ib) = 1._r8 / (3.0_r8 * canyon_hwr(l))
               albi(p,ib) = 1._r8 / (3.0_r8 * canyon_hwr(l))
            else if (col_pp%itype(c) == icol_shadewall) then
               albd(p,ib) = 1._r8 / (3.0_r8 * canyon_hwr(l))
               albi(p,ib) = 1._r8 / (3.0_r8 * canyon_hwr(l))
            else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
               albd(p,ib) = 1._r8 / (3.0_r8)
               albi(p,ib) = 1._r8 / (3.0_r8)
            else if (col_pp%itype(c) == icol_roof) then
               albd(p,ib) = 1._r8
               albi(p,ib) = 1._r8
            endif

            fabd(p,ib)     = 0._r8
            fabd_sun(p,ib) = 0._r8
            fabd_sha(p,ib) = 0._r8
            fabi(p,ib)     = 0._r8
            fabi_sun(p,ib) = 0._r8
            fabi_sha(p,ib) = 0._r8
            if (coszen > 0._r8) then
               ftdd(p,ib)  = 1._r8
            else
               ftdd(p,ib)  = 0._r8
            end if
            ftid(p,ib)     = 0._r8
            if (coszen > 0._r8) then
               ftii(p,ib)  = 1._r8
            else
               ftii(p,ib)  = 0._r8
            end if
         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Urban Code
      ! ----------------------------------------------------------------------------
      num_solar = 0
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         if (surfalb_vars%coszen_col(coli(l)) > 0._r8) num_solar = num_solar + 1
      end do

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


         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Only do the rest if all coszen are positive
      ! ----------------------------------------------------------------------------

      if (num_solar > 0)then

         ! Set constants - solar fluxes are per unit incoming flux

         do ib = 1,numrad
            do fl = 1,num_urbanl
               l = filter_urbanl(fl)
               coszen = surfalb_vars%coszen_col(coli(l)) ! Assumes coszen for each column are the same
               zen    = acos(coszen)
               sdir = 1._r8
               sdif = 1._r8

               ! Initialize direct and diffuse albedos such that if the Sun is below
               ! the horizon, p2g scaling returns an albedo of 1.0.

               sref_sunwall_dir   = 1._r8 / (3.0_r8 * canyon_hwr(l))
               sref_sunwall_dif   = 1._r8 / (3.0_r8 * canyon_hwr(l))
               sref_shadewall_dir = 1._r8 / (3.0_r8 * canyon_hwr(l))
               sref_shadewall_dif = 1._r8 / (3.0_r8 * canyon_hwr(l))
               sref_improad_dir   = 1._r8 / (3.0_r8)
               sref_improad_dif   = 1._r8 / (3.0_r8)
               sref_perroad_dir   = 1._r8 / (3.0_r8)
               sref_perroad_dif   = 1._r8 / (3.0_r8)
               sref_roof_dir      = 1._r8
               sref_roof_dif      = 1._r8

               ! Incident direct beam radiation for
               ! (a) roof and (b) road and both walls in urban canyon

               !if (coszen > 0._r8) then
               !if (num_urbanl > 0) then !!NOTE: Urban albedo is only called if num_urbanl >0??
                  ! Incident direct beam radiation for
                  ! (a) roof and (b) road and both walls in urban canyon

                  call incident_direct ( &
                       canyon_hwr(l), &
                       coszen    , &
                       zen       , &
                       sdir         , &
                       sdir_road    , &
                       sdir_sunwall , &
                       sdir_shadewall  )

                 ! Incident diffuse radiation for
                 ! (a) roof and (b) road and both walls in urban canyon.
                 call incident_diffuse (  l, &
                      canyon_hwr(l)  , &
                      sdif           , &
                      sdif_road      , &
                      sdif_sunwall   , &
                      sdif_shadewall , &
                      urbanparams_vars)

                ! Get snow albedos for roof and impervious and pervious road
                  ic = 0
                  !NOTE: Get rid of SnowAlbedo routine and merge
                  !! with below?
                  call SnowAlbedo(l, ib,&
                      coszen, &
                      ic, &
                      albsnd_roof    , &
                      albsnd_improad , &
                      albsnd_perroad  )

                  ic = 1
                  call SnowAlbedo( l,ib, &
                      coszen      , &
                      ic             , &
                      albsni_roof    , &
                      albsni_improad , &
                      albsni_perroad  )

               !end if

               ! Combine snow-free and snow albedos
               do c = lun_pp%coli(l), lun_pp%colf(l)
                 if (ctype(c) == icol_roof) then
                    alb_roof_dir_s = alb_roof_dir(l,ib)*(1._r8-frac_sno(c))  &
                         + albsnd_roof * frac_sno(c)
                    alb_roof_dif_s = alb_roof_dif(l,ib)*(1._r8-frac_sno(c))  &
                         + albsni_roof * frac_sno(c)

                 else if (ctype(c) == icol_road_imperv) then
                    alb_improad_dir_s = alb_improad_dir(l,ib)*(1._r8-frac_sno(c))  &
                         + albsnd_improad * frac_sno(c)
                    alb_improad_dif_s = alb_improad_dif(l,ib)*(1._r8-frac_sno(c))  &
                         + albsni_improad * frac_sno(c)

                 else if (ctype(c) == icol_road_perv) then
                    alb_perroad_dir_s = alb_perroad_dir(l,ib)*(1._r8-frac_sno(c))  &
                         + albsnd_perroad * frac_sno(c)
                    alb_perroad_dif_s = alb_perroad_dif(l,ib)*(1._r8-frac_sno(c))  &
                         + albsni_perroad * frac_sno(c)
                 end if
               end do !end of column loop


               ! Reflected and absorbed solar radiation per unit incident radiation
               ! for road and both walls in urban canyon allowing for multiple reflection
               ! Reflected and absorbed solar radiation per unit incident radiation for roof

               !if (num_urbanl > 0) then
                  call net_solar (l, ib, &
                       coszen             , &
                       canyon_hwr         (l), &
                       wtroad_perv        (l), &
                       sdir                  , &
                       sdif                  , &
                       alb_improad_dir_s  , &
                       alb_perroad_dir_s  , &
                       alb_wall_dir       (l, ib), &
                       alb_roof_dir_s     , &
                       alb_improad_dif_s  , &
                       alb_perroad_dif_s  , &
                       alb_wall_dif       (l, ib), &
                       alb_roof_dif_s     , &
                       sdir_road          , &
                       sdir_sunwall       , &
                       sdir_shadewall     , &
                       sdif_road          , &
                       sdif_sunwall       , &
                       sdif_shadewall     , &
                       sref_improad_dir   , &
                       sref_perroad_dir   , &
                       sref_sunwall_dir   , &
                       sref_shadewall_dir , &
                       sref_roof_dir      , &
                       sref_improad_dif   , &
                       sref_perroad_dif   , &
                       sref_sunwall_dif   , &
                       sref_shadewall_dif , &
                       sref_roof_dif      , &
                       urbanparams_vars, solarabs_vars)
               !end if

               do c = lun_pp%coli(l), lun_pp%colf(l)
                  if (ctype(c) == icol_roof) then
                     albgrd(c,ib) = sref_roof_dir
                     albgri(c,ib) = sref_roof_dif
                  else if (ctype(c) == icol_sunwall) then
                     albgrd(c,ib) = sref_sunwall_dir
                     albgri(c,ib) = sref_sunwall_dif
                  else if (ctype(c) == icol_shadewall) then
                     albgrd(c,ib) = sref_shadewall_dir
                     albgri(c,ib) = sref_shadewall_dif
                  else if (ctype(c) == icol_road_perv) then
                     albgrd(c,ib) = sref_perroad_dir
                     albgri(c,ib) = sref_perroad_dif
                  else if (ctype(c) == icol_road_imperv) then
                     albgrd(c,ib) = sref_improad_dir
                     albgri(c,ib) = sref_improad_dif
                  endif
               end do
            end do !loop over fl
         end do !loop over radiation band




         ! ----------------------------------------------------------------------------
         ! Map urban output to surfalb_vars components
         ! ----------------------------------------------------------------------------

         !  Set albgrd and albgri (ground albedos) and albd and albi (surface albedos)

         do ib = 1,numrad

            do fp = 1,num_urbanp
               p = filter_urbanp(fp)
               c = veg_pp%column(p)
               albd(p,ib) = albgrd(c,ib)
               albi(p,ib) = albgri(c,ib)
            end do
         end do
      end if

    end associate

  end subroutine UrbanAlbedo

  !-----------------------------------------------------------------------
  subroutine SnowAlbedo (l , ib, coszen, ind , &
       albsn_roof, albsn_improad, albsn_perroad )
    !
    ! !DESCRIPTION:
    ! Determine urban snow albedos
    ! ! !USES:
    !$acc routine seq
    use column_varcon, only : icol_roof, icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    integer , value, intent(in) :: l, ib     ! landunit index
    real(r8), value, intent(in) :: coszen    ! cosine solar zenith angle [landunit]
    integer , value, intent(in) :: ind       ! 0=direct beam, 1=diffuse radiation
    real(r8),  intent(out):: albsn_roof      ! roof snow albedo by waveband [landunit, numrad]
    real(r8),  intent(out):: albsn_improad   ! impervious road snow albedo by waveband [landunit, numrad]
    real(r8),  intent(out):: albsn_perroad   ! pervious road snow albedo by waveband [landunit, numrad]
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c ! indices
    integer  :: itype
    !
    ! These values are derived from Marshall (1989) assuming soot content of 1.5e-5
    ! (three times what LSM uses globally). Note that snow age effects are ignored here.
    real(r8), parameter :: snal0 = 0.66_r8 ! vis albedo of urban snow
    real(r8), parameter :: snal1 = 0.56_r8 ! nir albedo of urban snow
    real(r8) :: snal
    !-----------------------------------------------------------------------

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR
    if (ib == 1) then
      snal = snal0
    else
      snal = snal1
    end if
    associate(   &
         h2osno =>  col_ws%h2osno & ! Input:  [real(r8) (:) ]  snow water (mm H2O)
         )

      do c = lun_pp%coli(l), lun_pp%colf(l)

         if (coszen > 0._r8 .and. h2osno(c) > 0._r8) then
            if (col_pp%itype(c) == icol_roof) then
                 albsn_roof = snal
            else if (col_pp%itype(c) == icol_road_imperv) then
                albsn_improad = snal
            else if (col_pp%itype(c) == icol_road_perv) then
                albsn_perroad = snal
            end if
         else
            if (col_pp%itype(c) == icol_roof) then
               albsn_roof = 0._r8
            else if (col_pp%itype(c) == icol_road_imperv) then
               albsn_improad = 0._r8
            else if (col_pp%itype(c) == icol_road_perv) then
               albsn_perroad = 0._r8
            end if
         end if
      end do

    end associate

  end subroutine SnowAlbedo

  !-----------------------------------------------------------------------
  subroutine incident_direct ( canyon_hwr, coszen, zen , &
       sdir, sdir_road, sdir_sunwall, sdir_shadewall)
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
      !$acc routine seq
    use elm_varcon, only : rpi
    !
    ! !ARGUMENTS:
    !type(bounds_type), intent(in) :: bounds
    !integer , intent(in)  :: num_urbanl                          ! number of urban landunits
    !integer , intent(in)  :: filter_urbanl(:)                    ! urban landunit filter
    real(r8),value, intent(in)  :: canyon_hwr     ! ratio of building height to street width [landunit]
    real(r8),value, intent(in)  :: coszen         ! cosine solar zenith angle [landunit]
    real(r8),value, intent(in)  :: zen            ! solar zenith angle (radians) [landunit]
    real(r8),value, intent(in)  :: sdir           ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(out) :: sdir_road      ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(out) :: sdir_sunwall   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(out) :: sdir_shadewall ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    !
    ! !LOCAL VARIABLES:
    !integer  :: fl,l,i,ib                   ! indices
    logical, parameter  :: numchk = .false. ! true => perform numerical check of analytical solution
    real(r8) :: theta0                      ! critical canyon orientation for which road is no longer illuminated
    real(r8) :: tanzen                      ! tan(zenith angle)
    real(r8) :: swall_projected             ! direct beam solar radiation (per unit ground area) incident on wall
    real(r8) :: err1                        ! energy conservation error
    real(r8) :: err2                        ! energy conservation error
    real(r8) :: err3                        ! energy conservation error
    real(r8) :: sumr                        ! sum of sroad for each orientation (0 <= theta <= pi/2)
    real(r8) :: sumw                        ! sum of swall for each orientation (0 <= theta <= pi/2)
    real(r8) :: num                         ! number of orientations
    real(r8) :: theta                       ! canyon orientation relative to sun (0 <= theta <= pi/2)
    real(r8) :: zen0                        ! critical solar zenith angle for which sun begins to illuminate road
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    !
    ! do fl = 1,num_urbanl
    !    l = filter_urbanl(fl)
    !    if (coszen(l) > 0._r8) then
    !       theta0(l) = asin(min( (1._r8/(canyon_hwr(l)*tan(max(zen(l),0.000001_r8)))), 1._r8 ))
    !       tanzen(l) = tan(zen(l))
    !    end if
    ! end do

          !l = filter_urbanl(fl)
      if (coszen > 0._r8) then

         theta0 = asin(min( (1._r8/(canyon_hwr*tan(max(zen,0.000001_r8)))), 1._r8 ))
         tanzen = tan(zen)
         sdir_shadewall = 0._r8

         ! incident solar radiation on wall and road integrated over all canyon orientations (0 <= theta <= pi/2)

         sdir_road = sdir *                                    &
              (2._r8*theta0/rpi - 2./rpi*canyon_hwr*tanzen*(1._r8-cos(theta0)))
         sdir_sunwall = 2._r8 * sdir * ((1._r8/canyon_hwr)* &
              (0.5_r8-theta0/rpi) + (1._r8/rpi)*tanzen*(1._r8-cos(theta0)))

         ! conservation check for road and wall. need to use wall fluxes converted to ground area

         swall_projected = ( sdir_shadewall + sdir_sunwall ) * canyon_hwr
         err1 = sdir - (sdir_road + swall_projected)
      else
         sdir_road    = 0._r8
         sdir_sunwall = 0._r8
         sdir_shadewall  = 0._r8
      endif

       ! do fl = 1,num_urbanl
       !    l = filter_urbanl(fl)
       !    if (coszen(l) > 0._r8) then
       !       if (abs(err1(l)) > 0.001_r8) then
       !          !#py write (iulog,*) 'urban direct beam solar radiation balance error',err1(l)
       !          !#py write (iulog,*) 'clm model is stopping'
       !          !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       !       endif
       !    endif
       ! end do

       ! numerical check of analytical solution
       ! sum sroad and swall over all canyon orientations (0 <= theta <= pi/2)

       ! if (numchk) then
       !    do fl = 1,num_urbanl
       !       l = filter_urbanl(fl)
       !       if (coszen(l) > 0._r8) then
       !          sumr = 0._r8
       !          sumw = 0._r8
       !          num  = 0._r8
       !          do i = 1, 9000
       !             theta = i/100._r8 * rpi/180._r8
       !             zen0 = atan(1._r8/(canyon_hwr(l)*sin(theta)))
       !             if (zen(l) >= zen0) then
       !                sumr = sumr + 0._r8
       !                sumw = sumw + sdir(l,ib) / canyon_hwr(l)
       !             else
       !                sumr = sumr + sdir(l,ib) * (1._r8-canyon_hwr(l)*sin(theta)*tanzen(l))
       !                sumw = sumw + sdir(l,ib) * sin(theta)*tanzen(l)
       !             end if
       !             num = num + 1._r8
       !          end do
       !          err2 = sumr/num - sdir_road(l,ib)
       !          err3 = sumw/num - sdir_sunwall(l,ib)
       !       endif
       !    end do
          ! do fl = 1,num_urbanl
          !    l = filter_urbanl(fl)
          !    if (coszen(l) > 0._r8) then
          !       if (abs(err2(l)) > 0.0006_r8 ) then
          !          !#py write (iulog,*) 'urban road incident direct beam solar radiation error',err2(l)
          !          !#py write (iulog,*) 'clm model is stopping'
          !          !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
          !       endif
          !       if (abs(err3(l)) > 0.0006_r8 ) then
          !          !#py write (iulog,*) 'urban wall incident direct beam solar radiation error',err3(l)
          !          !#py write (iulog,*) 'clm model is stopping'
          !          !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
          !       end if
          !    end if
          ! end do
       ! end if

  end subroutine incident_direct

  !-----------------------------------------------------------------------
  subroutine incident_diffuse (l, canyon_hwr, &
       sdif, sdif_road, sdif_sunwall, sdif_shadewall, &
       urbanparams_vars)
    !
    ! !DESCRIPTION:
    ! Diffuse solar radiation incident on walls and road in urban canyon
    ! Conservation check: Total incoming diffuse
    ! (sdif) = sdif_road + (sdif_shadewall + sdif_sunwall)*canyon_hwr
    ! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
    !
    ! !ARGUMENTS:
    !$acc routine seq
    integer  , value      , intent(in)  :: l
    real(r8) , value      , intent(in)  :: canyon_hwr        ! ratio of building height to street width [landunit]
    real(r8) , value      , intent(in)  :: sdif              ! diffuse solar radiation incident on horizontal surface [landunit, numrad]
    real(r8)              , intent(out) :: sdif_road         ! diffuse solar radiation incident on road [landunit, numrad]
    real(r8)              , intent(out) :: sdif_sunwall      ! diffuse solar radiation (per unit wall area) incident on sunlit wall [landunit, numrad]
    real(r8)              , intent(out) :: sdif_shadewall    ! diffuse solar radiation (per unit wall area) incident on shaded wall [landunit, numrad]
    type(urbanparams_type), intent(in)  :: urbanparams_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fl, ib       ! indices
    real(r8) :: err             ! energy conservation error (W/m**2)
    real(r8) :: swall_projected ! diffuse solar radiation (per unit ground area) incident on wall (W/m**2)
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    associate(                            &
         vf_sr =>    urbanparams_vars%vf_sr , & ! Input:  [real(r8) (:) ]  view factor of sky for road
         vf_sw =>    urbanparams_vars%vf_sw   & ! Input:  [real(r8) (:) ]  view factor of sky for one wall
         )

         ! diffuse solar and conservation check. need to convert wall fluxes to ground area

         sdif_road       = sdif * vf_sr(l)
         sdif_sunwall    = sdif * vf_sw(l)
         sdif_shadewall  = sdif * vf_sw(l)

         swall_projected = (sdif_shadewall  + sdif_sunwall ) * canyon_hwr
         err = sdif - (sdif_road + swall_projected)

         ! error check

         ! do fl = 1, num_urbanl
         !    l = filter_urbanl(fl)
         !    if (abs(err(l)) > 0.001_r8) then
         !       !#py write (iulog,*) 'urban diffuse solar radiation balance error',err(l)
         !       !#py write (iulog,*) 'clm model is stopping'
         !       !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         !    endif
         ! end do

    end associate

  end subroutine incident_diffuse

  !-----------------------------------------------------------------------
  subroutine net_solar (l, ib, coszen, canyon_hwr, wtroad_perv, sdir, sdif                  , &
       alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir                            , &
       alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif                            , &
       sdir_road, sdir_sunwall, sdir_shadewall,                                                  &
       sdif_road, sdif_sunwall, sdif_shadewall,                                                  &
       sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir, sref_roof_dir , &
       sref_improad_dif, sref_perroad_dif, sref_sunwall_dif, sref_shadewall_dif, sref_roof_dif , &
       urbanparams_vars, solarabs_vars)
    !
    ! !DESCRIPTION:
    ! Solar radiation absorbed by road and both walls in urban canyon allowing
    ! for multiple reflection.
    !
    ! !ARGUMENTS:
      !$acc routine seq
    !type (bounds_type), intent(in) :: bounds
    !integer , intent(in)    :: num_urbanl                               ! number of urban landunits
    !integer , intent(in)    :: filter_urbanl(:)                         ! urban landunit filter
    integer , value, intent(in) :: l, ib
    real(r8),value, intent(in)  :: coszen              ! cosine solar zenith angle [landunit]
    real(r8),value, intent(in)  :: canyon_hwr          ! ratio of building height to street width [landunit]
    real(r8),value, intent(in)  :: wtroad_perv         ! weight of pervious road wrt total road [landunit]
    real(r8),value, intent(in)  :: sdir                ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8),value, intent(in)  :: sdif                ! diffuse solar radiation on horizontal surface [landunit, numrad]
    real(r8),value, intent(in)  :: alb_improad_dir     ! direct impervious road albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_perroad_dir     ! direct pervious road albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_wall_dir        ! direct  wall albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_roof_dir        ! direct  roof albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_improad_dif     ! diffuse impervious road albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_perroad_dif     ! diffuse pervious road albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_wall_dif        ! diffuse wall albedo [landunit, numrad]
    real(r8),value, intent(in)  :: alb_roof_dif        ! diffuse roof albedo [landunit, numrad]
    real(r8),value, intent(in)  :: sdir_road           ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8),value, intent(in)  :: sdir_sunwall        ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8),value, intent(in)  :: sdir_shadewall      ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8),value, intent(in)  :: sdif_road           ! diffuse solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8),value, intent(in)  :: sdif_sunwall        ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8),value, intent(in)  :: sdif_shadewall      ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dir    ! direct  solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dir    ! direct  solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dif    ! diffuse solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dif    ! diffuse solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dir    ! direct solar  rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dif    ! diffuse solar rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dir  ! direct solar  rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dif  ! diffuse solar rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dir       ! direct  solar rad reflected by roof (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dif       ! diffuse solar rad reflected by roof (per unit ground area)  per unit incident flux [landunit, numrad]
    type(urbanparams_type), intent(in)    :: urbanparams_vars
    type(solarabs_type)   , intent(inout) :: solarabs_vars
    !
    ! !LOCAL VARIABLES
    real(r8) :: wtroad_imperv         ! weight of impervious road wrt total road
    real(r8) :: sabs_canyon_dir         ! direct solar rad absorbed by canyon per unit incident flux
    real(r8) :: sabs_canyon_dif         ! diffuse solar rad absorbed by canyon per unit incident flux
    real(r8) :: sref_canyon_dir         ! direct solar reflected by canyon per unit incident flux
    real(r8) :: sref_canyon_dif         ! diffuse solar reflected by canyon per unit incident flux

    real(r8) :: improad_a_dir           ! absorbed direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_a_dif           ! absorbed diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dir           ! reflected direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dif           ! reflected diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_sky_dir       ! improad_r_dir to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dir   ! improad_r_dir to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dir ! improad_r_dir to shaded wall per unit incident flux
    real(r8) :: improad_r_sky_dif       ! improad_r_dif to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dif   ! improad_r_dif to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dif ! improad_r_dif to shaded wall per unit incident flux

    real(r8) :: perroad_a_dir           ! absorbed direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_a_dif           ! absorbed diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dir           ! reflected direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dif           ! reflected diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_sky_dir       ! perroad_r_dir to sky per unit incident flux
    real(r8) :: perroad_r_sunwall_dir   ! perroad_r_dir to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dir ! perroad_r_dir to shaded wall per unit incident flux
    real(r8) :: perroad_r_sky_dif       ! perroad_r_dif to sky per unit incident flux
    real(r8) :: perroad_r_sunwall_dif   ! perroad_r_dif to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dif ! perroad_r_dif to shaded wall per unit incident flux

    real(r8) :: road_a_dir              ! absorbed direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_a_dif              ! absorbed diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dir              ! reflected direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dif              ! reflected diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_sky_dir          ! road_r_dir to sky per unit incident flux
    real(r8) :: road_r_sunwall_dir      ! road_r_dir to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dir    ! road_r_dir to shaded wall per unit incident flux
    real(r8) :: road_r_sky_dif          ! road_r_dif to sky per unit incident flux
    real(r8) :: road_r_sunwall_dif      ! road_r_dif to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dif    ! road_r_dif to shaded wall per unit incident flux

    real(r8) :: sunwall_a_dir           ! absorbed direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_a_dif           ! absorbed diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dir           ! reflected direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dif           ! reflected diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_sky_dir       ! sunwall_r_dir to sky per unit incident flux
    real(r8) :: sunwall_r_road_dir      ! sunwall_r_dir to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dir ! sunwall_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: sunwall_r_sky_dif       ! sunwall_r_dif to sky per unit incident flux
    real(r8) :: sunwall_r_road_dif      ! sunwall_r_dif to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dif ! sunwall_r_dif to opposing (shaded) wall per unit incident flux

    real(r8) :: shadewall_a_dir         ! absorbed direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_a_dif         ! absorbed diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dir         ! reflected direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dif         ! reflected diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_sky_dir     ! shadewall_r_dir to sky per unit incident flux
    real(r8) :: shadewall_r_road_dir    ! shadewall_r_dir to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dir ! shadewall_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: shadewall_r_sky_dif     ! shadewall_r_dif to sky per unit incident flux
    real(r8) :: shadewall_r_road_dif    ! shadewall_r_dif to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dif ! shadewall_r_dif to opposing (sunlit) wall per unit incident flux

    real(r8) :: canyon_alb_dir          ! direct canyon albedo
    real(r8) :: canyon_alb_dif          ! diffuse canyon albedo

    real(r8) :: stot                ! sum of radiative terms
    real(r8) :: stot_dir            ! sum of direct radiative terms
    real(r8) :: stot_dif            ! sum of diffuse radiative terms
    integer  :: iter_dir,iter_dif                ! iteration counter
    real(r8) :: crit                             ! convergence criterion
    real(r8) :: err                              ! energy conservation error
    integer  :: pass
    integer, parameter :: n = 50                 ! number of interations
    real(r8) :: sabs_road                        ! temporary for absorption over road
    real(r8) :: sref_road                        ! temporary for reflected over road
    real(r8), parameter :: errcrit  = .00001_r8  ! error criteria
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes

    associate(                                                           &
         vf_sr              =>    urbanparams_vars%vf_sr               , & ! Input:  [real(r8) (:)   ]  view factor of sky for road
         vf_wr              =>    urbanparams_vars%vf_wr               , & ! Input:  [real(r8) (:)   ]  view factor of one wall for road
         vf_sw              =>    urbanparams_vars%vf_sw               , & ! Input:  [real(r8) (:)   ]  view factor of sky for one wall
         vf_rw              =>    urbanparams_vars%vf_rw               , & ! Input:  [real(r8) (:)   ]  view factor of road for one wall
         vf_ww              =>    urbanparams_vars%vf_ww               , & ! Input:  [real(r8) (:)   ]  view factor of opposing wall for one wall

         sabs_roof_dir      =>    solarabs_vars%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      =>    solarabs_vars%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_vars%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_vars%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_vars%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_vars%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_vars%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   =>    solarabs_vars%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_vars%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_vars%sabs_perroad_dif_lun     & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         )

      ! Calculate impervious road


          ! initial absorption and reflection for road and both walls.
          ! distribute reflected radiation to sky, road, and walls
          ! according to appropriate view factor. radiation reflected to
          ! road and walls will undergo multiple reflections within the canyon.
          ! do separately for direct beam and diffuse radiation.

          ! direct beam
          wtroad_imperv = 1._r8 - wtroad_perv
          road_a_dir              = 0.0_r8
          road_r_dir              = 0.0_r8
          improad_a_dir           = (1._r8-alb_improad_dir) * sdir_road
          improad_r_dir           =     alb_improad_dir  * sdir_road
          improad_r_sky_dir       = improad_r_dir * vf_sr(l)
          improad_r_sunwall_dir   = improad_r_dir * vf_wr(l)
          improad_r_shadewall_dir = improad_r_dir * vf_wr(l)
          road_a_dir              = road_a_dir + improad_a_dir * wtroad_imperv
          road_r_dir              = road_r_dir + improad_r_dir * wtroad_imperv

          perroad_a_dir           = (1._r8-alb_perroad_dir) * sdir_road
          perroad_r_dir           =     alb_perroad_dir * sdir_road
          perroad_r_sky_dir       = perroad_r_dir * vf_sr(l)
          perroad_r_sunwall_dir   = perroad_r_dir * vf_wr(l)
          perroad_r_shadewall_dir = perroad_r_dir * vf_wr(l)
          road_a_dir              = road_a_dir + perroad_a_dir * wtroad_perv
          road_r_dir              = road_r_dir + perroad_r_dir * wtroad_perv

          road_r_sky_dir          = road_r_dir * vf_sr(l)
          road_r_sunwall_dir      = road_r_dir * vf_wr(l)
          road_r_shadewall_dir    = road_r_dir * vf_wr(l)

          sunwall_a_dir           = (1._r8-alb_wall_dir ) * sdir_sunwall
          sunwall_r_dir           =     alb_wall_dir  * sdir_sunwall
          sunwall_r_sky_dir       = sunwall_r_dir * vf_sw(l)
          sunwall_r_road_dir      = sunwall_r_dir * vf_rw(l)
          sunwall_r_shadewall_dir = sunwall_r_dir * vf_ww(l)

          shadewall_a_dir         = (1._r8-alb_wall_dir) * sdir_shadewall
          shadewall_r_dir         =     alb_wall_dir  * sdir_shadewall
          shadewall_r_sky_dir     = shadewall_r_dir * vf_sw(l)
          shadewall_r_road_dir    = shadewall_r_dir * vf_rw(l)
          shadewall_r_sunwall_dir = shadewall_r_dir * vf_ww(l)

          ! diffuse

          road_a_dif              = 0.0_r8
          road_r_dif              = 0.0_r8
          improad_a_dif           = (1._r8-alb_improad_dif) * sdif_road
          improad_r_dif           =     alb_improad_dif  * sdif_road
          improad_r_sky_dif       = improad_r_dif * vf_sr(l)
          improad_r_sunwall_dif   = improad_r_dif * vf_wr(l)
          improad_r_shadewall_dif = improad_r_dif * vf_wr(l)
          road_a_dif              = road_a_dif + improad_a_dif * wtroad_imperv
          road_r_dif              = road_r_dif + improad_r_dif * wtroad_imperv

          perroad_a_dif           = (1._r8-alb_perroad_dif) * sdif_road
          perroad_r_dif           =     alb_perroad_dif  * sdif_road
          perroad_r_sky_dif       = perroad_r_dif * vf_sr(l)
          perroad_r_sunwall_dif   = perroad_r_dif * vf_wr(l)
          perroad_r_shadewall_dif = perroad_r_dif * vf_wr(l)
          road_a_dif              = road_a_dif + perroad_a_dif * wtroad_perv
          road_r_dif              = road_r_dif + perroad_r_dif * wtroad_perv

          road_r_sky_dif          = road_r_dif * vf_sr(l)
          road_r_sunwall_dif      = road_r_dif * vf_wr(l)
          road_r_shadewall_dif    = road_r_dif * vf_wr(l)

          sunwall_a_dif           = (1._r8-alb_wall_dif) * sdif_sunwall
          sunwall_r_dif           =     alb_wall_dif  * sdif_sunwall
          sunwall_r_sky_dif       = sunwall_r_dif * vf_sw(l)
          sunwall_r_road_dif      = sunwall_r_dif * vf_rw(l)
          sunwall_r_shadewall_dif = sunwall_r_dif * vf_ww(l)

          shadewall_a_dif         = (1._r8-alb_wall_dif) * sdif_shadewall
          shadewall_r_dif         =     alb_wall_dif  * sdif_shadewall
          shadewall_r_sky_dif     = shadewall_r_dif * vf_sw(l)
          shadewall_r_road_dif    = shadewall_r_dif * vf_rw(l)
          shadewall_r_sunwall_dif = shadewall_r_dif * vf_ww(l)

          ! initialize sum of direct and diffuse solar absorption and reflection for road and both walls

          sabs_improad_dir(l,ib)   = improad_a_dir
          sabs_perroad_dir(l,ib)   = perroad_a_dir
          sabs_sunwall_dir(l,ib)   = sunwall_a_dir
          sabs_shadewall_dir(l,ib) = shadewall_a_dir

          sabs_improad_dif(l,ib)   = improad_a_dif
          sabs_perroad_dif(l,ib)   = perroad_a_dif
          sabs_sunwall_dif(l,ib)   = sunwall_a_dif
          sabs_shadewall_dif(l,ib) = shadewall_a_dif

          sref_improad_dir   = improad_r_sky_dir
          sref_perroad_dir   = perroad_r_sky_dir
          sref_sunwall_dir   = sunwall_r_sky_dir
          sref_shadewall_dir = shadewall_r_sky_dir

          sref_improad_dif   = improad_r_sky_dif
          sref_perroad_dif   = perroad_r_sky_dif
          sref_sunwall_dif   = sunwall_r_sky_dif
          sref_shadewall_dif = shadewall_r_sky_dif


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

         ! reflected direct beam

         do iter_dir = 1, n
           ! step (1)
           stot = (sunwall_r_road_dir + shadewall_r_road_dir) * canyon_hwr
           !
           road_a_dir = 0.0_r8
           road_r_dir = 0.0_r8
           improad_a_dir = (1._r8-alb_improad_dir) * stot
           improad_r_dir =     alb_improad_dir  * stot
           road_a_dir    = road_a_dir + improad_a_dir * wtroad_imperv
           road_r_dir    = road_r_dir + improad_r_dir * wtroad_imperv
           perroad_a_dir = (1._r8-alb_perroad_dir) * stot
           perroad_r_dir =     alb_perroad_dir  * stot
           road_a_dir    = road_a_dir + perroad_a_dir * wtroad_perv
           road_r_dir    = road_r_dir + perroad_r_dir * wtroad_perv

           stot = road_r_sunwall_dir/canyon_hwr + shadewall_r_sunwall_dir
           sunwall_a_dir = (1._r8-alb_wall_dir) * stot
           sunwall_r_dir =     alb_wall_dir * stot

           stot = road_r_shadewall_dir/canyon_hwr + sunwall_r_shadewall_dir
           shadewall_a_dir = (1._r8-alb_wall_dir) * stot
           shadewall_r_dir =     alb_wall_dir  * stot

           ! step (2)

           sabs_improad_dir(l,ib)   = sabs_improad_dir(l,ib)   + improad_a_dir
           sabs_perroad_dir(l,ib)   = sabs_perroad_dir(l,ib)   + perroad_a_dir
           sabs_sunwall_dir(l,ib)   = sabs_sunwall_dir(l,ib)   + sunwall_a_dir
           sabs_shadewall_dir(l,ib) = sabs_shadewall_dir(l,ib) + shadewall_a_dir

           ! step (3)

           improad_r_sky_dir       = improad_r_dir * vf_sr(l)
           improad_r_sunwall_dir   = improad_r_dir * vf_wr(l)
           improad_r_shadewall_dir = improad_r_dir * vf_wr(l)

           perroad_r_sky_dir       = perroad_r_dir * vf_sr(l)
           perroad_r_sunwall_dir   = perroad_r_dir * vf_wr(l)
           perroad_r_shadewall_dir = perroad_r_dir * vf_wr(l)

           road_r_sky_dir          = road_r_dir * vf_sr(l)
           road_r_sunwall_dir      = road_r_dir * vf_wr(l)
           road_r_shadewall_dir    = road_r_dir * vf_wr(l)

           sunwall_r_sky_dir       = sunwall_r_dir * vf_sw(l)
           sunwall_r_road_dir      = sunwall_r_dir * vf_rw(l)
           sunwall_r_shadewall_dir = sunwall_r_dir * vf_ww(l)

           shadewall_r_sky_dir     = shadewall_r_dir * vf_sw(l)
           shadewall_r_road_dir    = shadewall_r_dir * vf_rw(l)
           shadewall_r_sunwall_dir = shadewall_r_dir * vf_ww(l)

           ! step (4)

           sref_improad_dir   = sref_improad_dir + improad_r_sky_dir
           sref_perroad_dir   = sref_perroad_dir + perroad_r_sky_dir
           sref_sunwall_dir   = sref_sunwall_dir + sunwall_r_sky_dir
           sref_shadewall_dir = sref_shadewall_dir + shadewall_r_sky_dir

           ! step (5)

           crit = max(road_a_dir, sunwall_a_dir, shadewall_a_dir)
           if (crit < errcrit) exit

         end do
         if (iter_dir >= n) then
           print *, 'urban net solar radiation error: no convergence, direct beam'
           stop
           !#py write (iulog,*) 'clm model is stopping'
           !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         endif

         ! reflected diffuse

         do iter_dif = 1, n
           ! step (1)

           stot = (sunwall_r_road_dif + shadewall_r_road_dif ) * canyon_hwr
           road_a_dif    = 0.0_r8
           road_r_dif    = 0.0_r8
           improad_a_dif = (1._r8-alb_improad_dif) * stot
           improad_r_dif =     alb_improad_dif  * stot
           road_a_dif    = road_a_dif + improad_a_dif * wtroad_imperv
           road_r_dif    = road_r_dif + improad_r_dif * wtroad_imperv
           perroad_a_dif = (1._r8-alb_perroad_dif) * stot
           perroad_r_dif =     alb_perroad_dif  * stot
           road_a_dif    = road_a_dif + perroad_a_dif * wtroad_perv
           road_r_dif    = road_r_dif + perroad_r_dif * wtroad_perv

           stot = road_r_sunwall_dif/canyon_hwr + shadewall_r_sunwall_dif
           sunwall_a_dif = (1._r8-alb_wall_dif) * stot
           sunwall_r_dif =     alb_wall_dif  * stot

           stot = road_r_shadewall_dif/canyon_hwr + sunwall_r_shadewall_dif
           shadewall_a_dif = (1._r8-alb_wall_dif) * stot
           shadewall_r_dif =     alb_wall_dif  * stot

           ! step (2)

           sabs_improad_dif(l,ib)   = sabs_improad_dif(l,ib)   + improad_a_dif
           sabs_perroad_dif(l,ib)   = sabs_perroad_dif(l,ib)   + perroad_a_dif
           sabs_sunwall_dif(l,ib)   = sabs_sunwall_dif(l,ib)   + sunwall_a_dif
           sabs_shadewall_dif(l,ib) = sabs_shadewall_dif(l,ib) + shadewall_a_dif

           ! step (3)

           improad_r_sky_dif        = improad_r_dif  * vf_sr(l)
           improad_r_sunwall_dif    = improad_r_dif  * vf_wr(l)
           improad_r_shadewall_dif  = improad_r_dif  * vf_wr(l)

           perroad_r_sky_dif        = perroad_r_dif  * vf_sr(l)
           perroad_r_sunwall_dif    = perroad_r_dif  * vf_wr(l)
           perroad_r_shadewall_dif  = perroad_r_dif  * vf_wr(l)

           road_r_sky_dif          = road_r_dif * vf_sr(l)
           road_r_sunwall_dif      = road_r_dif * vf_wr(l)
           road_r_shadewall_dif    = road_r_dif * vf_wr(l)

           sunwall_r_sky_dif       = sunwall_r_dif * vf_sw(l)
           sunwall_r_road_dif      = sunwall_r_dif * vf_rw(l)
           sunwall_r_shadewall_dif = sunwall_r_dif * vf_ww(l)

           shadewall_r_sky_dif     = shadewall_r_dif * vf_sw(l)
           shadewall_r_road_dif    = shadewall_r_dif * vf_rw(l)
           shadewall_r_sunwall_dif = shadewall_r_dif * vf_ww(l)

           ! step (4)

           sref_improad_dif   = sref_improad_dif  + improad_r_sky_dif
           sref_perroad_dif   = sref_perroad_dif  + perroad_r_sky_dif
           sref_sunwall_dif   = sref_sunwall_dif  + sunwall_r_sky_dif
           sref_shadewall_dif = sref_shadewall_dif + shadewall_r_sky_dif

           ! step (5)

           crit = max(road_a_dif, sunwall_a_dif, shadewall_a_dif)
           if (crit < errcrit) exit
         end do

         if (iter_dif >= n) then
           print *, 'urban net solar radiation error: no convergence, diffuse'
           stop
           !#py write (iulog,*) 'clm model is stopping'
           !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         endif

         ! total reflected by canyon - sum of solar reflection to sky from canyon.
         ! project wall fluxes to horizontal surface

         sref_canyon_dir = 0.0_r8
         sref_canyon_dif = 0.0_r8
         sref_canyon_dir = sref_canyon_dir + sref_improad_dir * wtroad_imperv
         sref_canyon_dif = sref_canyon_dif + sref_improad_dif * wtroad_imperv
         sref_canyon_dir = sref_canyon_dir + sref_perroad_dir * wtroad_perv
         sref_canyon_dif = sref_canyon_dif + sref_perroad_dif * wtroad_perv
         sref_canyon_dir = sref_canyon_dir + (sref_sunwall_dir + sref_shadewall_dir)*canyon_hwr
         sref_canyon_dif = sref_canyon_dif + (sref_sunwall_dif + sref_shadewall_dif)*canyon_hwr

         ! total absorbed by canyon. project wall fluxes to horizontal surface

         sabs_canyon_dir = 0.0_r8
         sabs_canyon_dif = 0.0_r8
         sabs_canyon_dir = sabs_canyon_dir + sabs_improad_dir(l,ib)*wtroad_imperv
         sabs_canyon_dif = sabs_canyon_dif + sabs_improad_dif(l,ib)*wtroad_imperv
         sabs_canyon_dir = sabs_canyon_dir + sabs_perroad_dir(l,ib)*wtroad_perv
         sabs_canyon_dif = sabs_canyon_dif + sabs_perroad_dif(l,ib)*wtroad_perv
         sabs_canyon_dir = sabs_canyon_dir + (sabs_sunwall_dir(l,ib) + sabs_shadewall_dir(l,ib))*canyon_hwr
         sabs_canyon_dif = sabs_canyon_dif + (sabs_sunwall_dif(l,ib) + sabs_shadewall_dif(l,ib))*canyon_hwr

         ! conservation check. note: previous conservation checks confirm partioning of total direct
         ! beam and diffuse radiation from atmosphere to road and walls is conserved as
         !    sdir (from atmosphere) = sdir_road + (sdir_sunwall + sdir_shadewall)*canyon_hwr
         !    sdif (from atmosphere) = sdif_road + (sdif_sunwall + sdif_shadewall)*canyon_hwr

         stot_dir = sdir_road + (sdir_sunwall + sdir_shadewall)*canyon_hwr
         stot_dif = sdif_road + (sdif_sunwall + sdif_shadewall)*canyon_hwr
         !
         err = stot_dir + stot_dif &
              - (sabs_canyon_dir + sabs_canyon_dif + sref_canyon_dir + sref_canyon_dif)
         if (abs(err) > 0.001_r8 ) then
                  print *, 'urban net solar radiation balance error for ib=',ib,' err= ',err
                  stop
                  !#py write(iulog,*)' l= ',l,' ib= ',ib
                  !#py write(iulog,*)' stot_dir        = ',stot_dir(l)
                  !#py write(iulog,*)' stot_dif        = ',stot_dif(l)
                  !#py write(iulog,*)' sabs_canyon_dir = ',sabs_canyon_dir(l)
                  !#py write(iulog,*)' sabs_canyon_dif = ',sabs_canyon_dif(l)
                  !#py write(iulog,*)' sref_canyon_dir = ',sref_canyon_dir(l)
                  !#py write(iulog,*)' sref_canyon_dif = ',sref_canyon_dir(l)
                  !#py write(iulog,*) 'clm model is stopping'
                  !#py !#py call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         endif

         ! canyon albedo  -- Not Used ??
         canyon_alb_dif = sref_canyon_dif / max(stot_dif, 1.e-06_r8)
         canyon_alb_dir = sref_canyon_dir / max(stot_dir, 1.e-06_r8)


         ! Refected and absorbed solar radiation per unit incident radiation for roof
         sref_roof_dir = alb_roof_dir * sdir
         sref_roof_dif = alb_roof_dif * sdir
         sabs_roof_dir(l,ib) = sdir - sref_roof_dir
         sabs_roof_dif(l,ib) = sdif - sref_roof_dif

    end associate

  end subroutine net_solar

end module UrbanAlbedoMod

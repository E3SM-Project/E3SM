module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  ! ! USES:
  use clm_varctl , only : paramfile, iulog, use_ed, use_cn
  use spmdMod    , only : masterproc
  use fileutils  , only : getfil
  use ncdio_pio  , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio  , only : file_desc_t , ncd_inqdid, ncd_inqdlen

  implicit none
  private
  !
  public :: readParameters
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParameters (nutrient_competition_method)
    !
    ! ! USES:
    use EDSharedParamsMod                 , only : EDParamsReadShared
    use EDParamsMod                       , only : EDParamsRead 
    use SFParamsMod                       , only : SFParamsRead
    use CNSharedParamsMod                 , only : CNParamsReadShared
    use CNGapMortalityMod                 , only : readCNGapMortParams                    => readParams
    use CNMRespMod                        , only : readCNMRespParams                      => readParams
    use CNPhenologyMod                    , only : readCNPhenolParams                     => readParams
    use SoilBiogeochemCompetitionMod      , only : readSoilBiogeochemCompetitionParams    => readParams
    use SoilBiogeochemNLeachingMod        , only : readSoilBiogeochemNLeachingParams      => readParams
    use SoilBiogeochemNitrifDenitrifMod   , only : readSoilBiogeochemNitrifDenitrifParams => readParams
    use SoilBiogeochemLittVertTranspMod   , only : readSoilBiogeochemLittVertTranspParams => readParams
    use SoilBiogeochemPotentialMod        , only : readSoilBiogeochemPotentialParams      => readParams
    use SoilBiogeochemDecompMod           , only : readSoilBiogeochemDecompParams         => readParams
    use SoilBiogeochemDecompCascadeBGCMod , only : readSoilBiogeochemDecompBgcParams      => readParams
    use SoilBiogeochemDecompCascadeCNMod  , only : readSoilBiogeochemDecompCnParams       => readParams
    use ch4Mod                            , only : readCH4Params                          => readParams
    use NutrientCompetitionMethodMod      , only : nutrient_competition_method_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    character(len=32)  :: subname = 'readParameters'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: reading ED '//' parameters '
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    if (use_ed) then
       call EDParamsReadShared(ncid)
       call EDParamsRead(ncid)
       call SFParamsRead(ncid)
    end if

    if (use_cn) then
       call CNParamsReadShared(ncid)
       call nutrient_competition_method%readParams(ncid)
       call readCNGapMortParams(ncid)
       call readCNMRespParams(ncid)
       call readCNPhenolParams(ncid)
    end if

    call readSoilBiogeochemCompetitionParams(ncid)
    call readSoilBiogeochemDecompBgcParams(ncid)
    call readSoilBiogeochemDecompCnParams(ncid)
    call readSoilBiogeochemDecompParams(ncid)
    call readSoilBiogeochemLittVertTranspParams(ncid)
    call readSoilBiogeochemNitrifDenitrifParams(ncid)
    call readSoilBiogeochemNLeachingParams(ncid)
    call readSoilBiogeochemPotentialParams(ncid)

    call readCH4Params (ncid)

    call ncd_pio_closefile(ncid)

  end subroutine readParameters

end module readParamsMod

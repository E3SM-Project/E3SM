module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  !
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varctl   , only: use_cn, use_century_decomp, use_nitrif_denitrif, &
                           use_lch4
  implicit none
  save
  private
  !
  public :: readParameters
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParameters ()
    !
    implicit none
    !-----------------------------------------------------------------------

    call CNParamsReadFile()

  end subroutine readParameters
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine CNParamsReadFile ()
    !
    ! read CN and BGC shared parameters
    !
    use CNAllocationMod         , only : readCNAllocParams
    use CNDecompMod             , only : readCNDecompParams
    use CNDecompCascadeBGCMod   , only : readCNDecompBgcParams
    use CNDecompCascadeCNMod    , only : readCNDecompCnParams
    use CNPhenologyMod          , only : readCNPhenolParams
    use CNMRespMod              , only : readCNMRespParams
    use CNNDynamicsMod          , only : readCNNDynamicsParams
    use CNGapMortalityMod       , only : readCNGapMortParams 
    use CNNitrifDenitrifMod     , only : readCNNitrifDenitrifParams
    use CNSoilLittVertTranspMod , only : readCNSoilLittVertTranspParams
    use CNSharedParamsMod       , only : CNParamsReadShared,CNParamsShareInst
    use ch4Mod                  , only : readCH4Params
    use clm_varctl              , only : paramfile, iulog
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use abortutils              , only : endrun
    use ncdio_pio               , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNParamsReadShared'
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'readParamsMod.F90::'//trim(subname)//' :: reading CN '//&
          'and BGC parameter file'
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    !
    ! some parameters (eg. organic_max) are used in non-CN, non-BGC cases
    !
    call CNParamsReadShared(ncid)

    if (use_cn) then
       !
       ! populate each module with private parameters
       !
       call readCNAllocParams(ncid)
       call readCNDecompParams(ncid)
       if (use_century_decomp) then
          call readCNDecompBgcParams(ncid)
       else
          call readCNDecompCnParams(ncid)
       end if
       call readCNPhenolParams(ncid)
       call readCNMRespParams (ncid)
       call readCNNDynamicsParams (ncid)
       call readCNGapMortParams (ncid)
       if (use_nitrif_denitrif) then
          call readCNNitrifDenitrifParams(ncid)
       end if

       call readCNSoilLittVertTranspParams(ncid)

       if (use_lch4) then
          call readCH4Params (ncid)
       end if

       !
       ! close CN params file
       !
       call ncd_pio_closefile(ncid)

    end if

  end subroutine CNParamsReadFile

end module readParamsMod

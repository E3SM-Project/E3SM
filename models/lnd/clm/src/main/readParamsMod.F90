module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  use clm_varctl   , only: use_cn, use_century_decomp, use_nitrif_denitrif, &
                           use_lch4
  implicit none
  save
  private
  !
  public :: readSharedParameters
  public :: readPrivateParameters
  
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readSharedParameters ()
    !
    implicit none
    !-----------------------------------------------------------------------

    call CNParamsSharedReadFile()
    ! calls for ED parameters
    call EDParamsReadFile()

  end subroutine readSharedParameters
  !-----------------------------------------------------------------------
  subroutine EDParamsReadFile ()
    !
    ! read ED shared parameters
    !
    use EDParamsMod             , only : EDParamsRead 
    use SFParamsMod             , only : SFParamsRead
    use clm_varctl              , only : paramfile, iulog, use_ed
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'EDParamsReadFile'
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: reading ED '//&
          ' parameters '
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    if ( use_ed ) then
       call EDParamsRead(ncid)
       call SFParamsRead(ncid)
    endif
    !
    ! close ED params file
    !
    call ncd_pio_closefile(ncid)

  end subroutine EDParamsReadFile

  !-----------------------------------------------------------------------
  subroutine CNParamsSharedReadFile ()
    !
    ! read CN and BGC shared parameters
    !

    use CNSharedParamsMod       , only : CNParamsReadShared

    use clm_varctl              , only : paramfile, iulog
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen                                       
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNParamsSharedReadFile'
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


  end subroutine CNParamsSharedReadFile

  !-----------------------------------------------------------------------
  subroutine readPrivateParameters
    ! read CN and BGC shared parameters
    !
    use CNAllocationBetrMod     , only : readCNAllocBetrParams
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
    use ch4Mod                  , only : readCH4Params
    use clm_varctl              , only : paramfile, iulog, use_betr
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
    use tracer_varcon           , only : is_active_betr_bgc                                         
    use betr_initializeMod      , only : bgc_reaction, betrtracer_vars
    
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'readPrivateParameters'
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
    
    if(use_betr)then
      call bgc_reaction%readParams(ncid, betrtracer_vars)   
    endif
    
    if (use_cn) then
       !
       ! populate each module with private parameters
       !       
       if (is_active_betr_bgc)then

          call readCNAllocBetrParams(ncid) 

       else
          call readCNAllocParams(ncid)

          call readCNDecompParams(ncid)
          if (use_century_decomp) then
            call readCNDecompBgcParams(ncid)
          else
            call readCNDecompCnParams(ncid)
          end if
          if (use_nitrif_denitrif) then
            call readCNNitrifDenitrifParams(ncid)
          end if

          call readCNSoilLittVertTranspParams(ncid)

          if (use_lch4) then
            call readCH4Params (ncid)
          end if
       endif
       
       call readCNPhenolParams(ncid)
       call readCNMRespParams (ncid)
       call readCNNDynamicsParams (ncid)
       call readCNGapMortParams (ncid)

    end if

    !
    ! close CN params file
    !
    call ncd_pio_closefile(ncid)

  end subroutine readPrivateParameters
end module readParamsMod

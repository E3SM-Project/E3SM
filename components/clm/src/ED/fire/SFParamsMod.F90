module SFParamsMod
   !
   ! module that deals with reading the SF parameter file
   !
   use shr_kind_mod      , only: r8 => shr_kind_r8
   use EDtypesMod        , only: NLSC,NFSC,NCWD

   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   real(r8),protected :: SF_val_fdi_a
   real(r8),protected :: SF_val_fdi_b
   real(r8),protected :: SF_val_fdi_alpha
   real(r8),protected :: SF_val_miner_total
   real(r8),protected :: SF_val_fuel_energy
   real(r8),protected :: SF_val_part_dens
   real(r8),protected :: SF_val_miner_damp
   real(r8),protected :: SF_val_max_durat
   real(r8),protected :: SF_val_durat_slope
   real(r8),protected :: SF_val_alpha_SH
   real(r8),protected :: SF_val_alpha_FMC(NLSC)
   real(r8),protected :: SF_val_CWD_frac(NCWD)
   real(r8),protected :: SF_val_max_decomp(NLSC)
   real(r8),protected :: SF_val_SAV(NFSC)
   real(r8),protected :: SF_val_FBD(NFSC)
   real(r8),protected :: SF_val_min_moisture(NFSC)
   real(r8),protected :: SF_val_mid_moisture(NFSC)
   real(r8),protected :: SF_val_low_moisture_C(NFSC)
   real(r8),protected :: SF_val_low_moisture_S(NFSC)
   real(r8),protected :: SF_val_mid_moisture_C(NFSC)
   real(r8),protected :: SF_val_mid_moisture_S(NFSC)
  
   character(len=20),parameter :: SF_name_fdi_a = "fdi_a"
   character(len=20),parameter :: SF_name_fdi_b = "fdi_b"
   character(len=20),parameter :: SF_name_fdi_alpha = "fdi_alpha"
   character(len=20),parameter :: SF_name_miner_total = "miner_total"
   character(len=20),parameter :: SF_name_fuel_energy = "fuel_energy"
   character(len=20),parameter :: SF_name_part_dens = "part_dens"
   character(len=20),parameter :: SF_name_miner_damp = "miner_damp"
   character(len=20),parameter :: SF_name_max_durat = "max_durat"
   character(len=20),parameter :: SF_name_durat_slope = "durat_slope"
   character(len=20),parameter :: SF_name_alpha_SH = "alpha_SH"
   character(len=20),parameter :: SF_name_alpha_FMC = "alpha_FMC"
   character(len=20),parameter :: SF_name_CWD_frac = "CWD_frac"
   character(len=20),parameter :: SF_name_max_decomp = "max_decomp"
   character(len=20),parameter :: SF_name_SAV = "SAV"
   character(len=20),parameter :: SF_name_FBD = "FBD"
   character(len=20),parameter :: SF_name_min_moisture = "min_moisture"
   character(len=20),parameter :: SF_name_mid_moisture = "mid_moisture"
   character(len=20),parameter :: SF_name_low_moisture_C = "low_moisture_C"
   character(len=20),parameter :: SF_name_low_moisture_S = "low_moisture_S"
   character(len=20),parameter :: SF_name_mid_moisture_C = "mid_moisture_C"
   character(len=20),parameter :: SF_name_mid_moisture_S = "mid_moisture_S"

   public :: SFParamsRead
  
contains
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine SFParamsRead(ncid)
     !
     ! calls to initialize parameter instance and do ncdio read
     !
     use ncdio_pio    , only : file_desc_t
     
     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     call SFParamsReadLocal(ncid)

  end subroutine SFParamsRead
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine SFParamsReadLocal(ncid)
     !
     ! read the netcdf file and populate internalInstScalar
     !
     use ncdio_pio         , only : file_desc_t
     use paramUtilMod      , only : readNcdio

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     ! local vars
     character(len=32)  :: subname = 'SFParamsReadLocal::'

     !
     ! call read function
     !

      call readNcdio(ncid = ncid, &
         varName=SF_name_fdi_a, &
         callingName=subname, &
         retVal=SF_val_fdi_a)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_fdi_b, &
         callingName=subname, &
         retVal=SF_val_fdi_b)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_fdi_alpha, &
         callingName=subname, &
         retVal=SF_val_fdi_alpha)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_miner_total, &
         callingName=subname, &
         retVal=SF_val_miner_total)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_fuel_energy, &
         callingName=subname, &
         retVal=SF_val_fuel_energy)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_part_dens, &
         callingName=subname, &
         retVal=SF_val_part_dens)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_miner_damp, &
         callingName=subname, &
         retVal=SF_val_miner_damp)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_max_durat, &
         callingName=subname, &
         retVal=SF_val_max_durat)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_durat_slope, &
         callingName=subname, &
         retVal=SF_val_durat_slope)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_alpha_SH, &
         callingName=subname, &
         retVal=SF_val_alpha_SH)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_alpha_FMC, &
         callingName=subname, &
         retVal=SF_val_alpha_FMC)

      call readNcdio(ncid = ncid, &
         varName=SF_name_CWD_frac, &
         callingName=subname, &
         retVal=SF_val_CWD_frac)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_max_decomp, &
         callingName=subname, &
         retVal=SF_val_max_decomp)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_SAV, &
         callingName=subname, &
         retVal=SF_val_SAV)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_FBD, &
         callingName=subname, &
         retVal=SF_val_FBD)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_min_moisture, &
         callingName=subname, &
         retVal=SF_val_min_moisture)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture, &
         callingName=subname, &
         retVal=SF_val_mid_moisture)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_low_moisture_C, &
         callingName=subname, &
         retVal=SF_val_low_moisture_C)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_low_moisture_S, &
         callingName=subname, &
         retVal=SF_val_low_moisture_S)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture_C, &
         callingName=subname, &
         retVal=SF_val_mid_moisture_C)
  
      call readNcdio(ncid = ncid, &
         varName=SF_name_mid_moisture_S, &
         callingName=subname, &
         retVal=SF_val_mid_moisture_S)
  
  end subroutine SFParamsReadLocal
  !-----------------------------------------------------------------------

end module SFParamsMod

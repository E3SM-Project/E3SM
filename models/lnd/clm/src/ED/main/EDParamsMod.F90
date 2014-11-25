module EDParamsMod
   !
   ! module that deals with reading the ED parameter file
   !
   use shr_kind_mod      , only: r8 => shr_kind_r8

   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be 
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   real(r8),protected :: ED_val_grass_spread
   real(r8),protected :: ED_val_comp_excln
   real(r8),protected :: ED_val_stress_mort
   real(r8),protected :: ED_val_dispersal
   real(r8),protected :: ED_val_grperc
   real(r8),protected :: ED_val_maxspread
   real(r8),protected :: ED_val_minspread
   real(r8),protected :: ED_val_init_litter
   real(r8),protected :: ED_val_nfires
   real(r8),protected :: ED_val_understorey_death
   real(r8),protected :: ED_val_profile_tol
   real(r8),protected :: ED_val_ag_biomass
  
   character(len=20),parameter :: ED_name_grass_spread = "grass_spread"
   character(len=20),parameter :: ED_name_comp_excln = "comp_excln"
   character(len=20),parameter :: ED_name_stress_mort = "stress_mort"
   character(len=20),parameter :: ED_name_dispersal = "dispersal"
   character(len=20),parameter :: ED_name_grperc = "grperc"
   character(len=20),parameter :: ED_name_maxspread = "maxspread"
   character(len=20),parameter :: ED_name_minspread = "minspread"
   character(len=20),parameter :: ED_name_init_litter = "init_litter"
   character(len=20),parameter :: ED_name_nfires = "nfires"
   character(len=20),parameter :: ED_name_understorey_death = "understorey_death"
   character(len=20),parameter :: ED_name_profile_tol = "profile_tol"
   character(len=20),parameter :: ED_name_ag_biomass= "ag_biomass"   
   
   public :: EDParamsRead
  
contains

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine EDParamsRead(ncid)
     !
     ! calls to initialize parameter instance and do ncdio read
     !
     use ncdio_pio    , only : file_desc_t
     
     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     call EDParamsReadLocal(ncid)

  end subroutine EDParamsRead
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine EDParamsReadLocal(ncid)
     !
     ! read the netcdf file and populate internalInstScalar
     !
     use ncdio_pio         , only : file_desc_t
     use paramUtilMod      , only : readNcdio

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id

     ! local vars
     character(len=32)  :: subname = 'EDParamsReadLocal::'

     !
     ! call read function
     !

      call readNcdio(ncid = ncid, &
         varName=ED_name_grass_spread, &
         callingName=subname, &
         retVal=ED_val_grass_spread)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_comp_excln, &
         callingName=subname, &
         retVal=ED_val_comp_excln)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_stress_mort, &
         callingName=subname, &
         retVal=ED_val_stress_mort)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_dispersal, &
         callingName=subname, &
         retVal=ED_val_dispersal)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_grperc, &
         callingName=subname, &
         retVal=ED_val_grperc)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_maxspread, &
         callingName=subname, &
         retVal=ED_val_maxspread)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_minspread, &
         callingName=subname, &
         retVal=ED_val_minspread)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_init_litter, &
         callingName=subname, &
         retVal=ED_val_init_litter)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_nfires, &
         callingName=subname, &
         retVal=ED_val_nfires)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_understorey_death, &
         callingName=subname, &
         retVal=ED_val_understorey_death)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_profile_tol, &
         callingName=subname, &
         retVal=ED_val_profile_tol)
  
      call readNcdio(ncid = ncid, &
         varName=ED_name_ag_biomass, &
         callingName=subname, &
         retVal=ED_val_ag_biomass)
  
  end subroutine EDParamsReadLocal
  !-----------------------------------------------------------------------

end module EDParamsMod

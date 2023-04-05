module FireMethodType

   !---------------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Abstract base class for functions to implement fire model and data  for 
   ! FATES.  This could be adopted by elm fire module in the future if desired.
   ! Adapted by Greg Lemieux (LBL) based on work by Erik Kluzek (NCAR).  
   !
   ! !USES:
   implicit none
   private

   ! !PUBLIC TYPES:
   public :: fire_method_type
 
   type, abstract :: fire_method_type
    contains
 
      ! Initialize the fire datasets
      procedure(FireInit_interface)   , public, deferred :: FireInit
 
      ! Interpolate the fire datasets
      procedure(FireInterp_interface) , public, deferred :: FireInterp

   end type fire_method_type
 
   abstract interface
 
      ! Note: The following code is adapted based on what Bill Sacks has done for soil water retention curve
      ! polymorphism. Therefore, I also keep some suggestions he gave there.
      !
      ! - Make the interfaces contain all possible inputs that are needed by any
      !   implementation; each implementation will then ignore the inputs it doesn't need.
      !
      ! - For inputs that are needed only by particular implementations - and particularly
      !   for inputs that are constant in time 
      !   pass these into the constructor, and save pointers to these inputs as components
      !   of the child type that needs them. Then they aren't needed as inputs to the
      !   individual routines, allowing the interfaces for these routines to remain more
      !   consistent between different implementations.
      !
      !---------------------------------------------------------------------------
      ! subroutine FireInit_interface(this, bounds, cnstate_vars )
      subroutine FireInit_interface(this, bounds, NLFilename )
        !
        ! !DESCRIPTION:
        ! Initialize Fire datasets
        !
        ! USES
        use decompMod              , only : bounds_type
        ! use CNStateType            , only : cnstate_type
        import :: fire_method_type
        ! !ARGUMENTS:
        class(fire_method_type)     :: this
        type(bounds_type), intent(in) :: bounds
        character(len=*), intent(in) :: NLFilename ! Namelist filename
        ! type(cnstate_type), intent(in) :: cnstate_vars

      end subroutine FireInit_interface

      !---------------------------------------------------------------------------      
      subroutine FireInterp_interface(this, bounds)
        !
        ! !DESCRIPTION:
        ! Interpolate Fire datasets
        !
        ! USES
        use decompMod              , only : bounds_type
        import :: fire_method_type
        ! !ARGUMENTS:
        class(fire_method_type)     :: this
        type(bounds_type), intent(in) :: bounds
    
      end subroutine FireInterp_interface
    
      !---------------------------------------------------------------------------
 
   end interface
 
 end module FireMethodType
 

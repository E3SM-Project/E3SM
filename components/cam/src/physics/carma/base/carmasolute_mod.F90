!! The CARMASOLUTE module contains configuration information about a solute used by CARMA.
!!
!!  @version May-2009 
!!  @author  Chuck Bardeen 
module carmasolute_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMASOLUTE_Create
  public CARMASOLUTE_Destroy
  public CARMASOLUTE_Get
  public CARMASOLUTE_Print

contains

  !! Defines a solute used by CARMA for nucleation and growth of cloud and 
  !! aerosol particles.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_AddGas
  !! @see CARMASOLUTE_Destroy
 subroutine CARMASOLUTE_Create(carma, isolute, name, ions, wtmol, rho, rc, shortname)
    type(carma_type), intent(inout)       :: carma               !! the carma object
    integer, intent(in)                   :: isolute             !! the solute index
    character(*), intent(in)              :: name                !! the solute name, maximum of 255 characters
		integer, intent(in)                   :: ions                !! Number of ions solute dissociates into
    real(kind=f), intent(in)              :: wtmol               !! the solute molecular weight [g/mol]
		real(kind=f), intent(in)              :: rho                 !! Mass density of solute
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    character(*), optional, intent(in)    :: shortname           !! the solute shortname, maximum of 6 characters

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough solutes allocated.
    if (isolute > carma%f_NSOLUTE) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMASOLUTE_Create:: ERROR - The specifed solute (", &
        isolute, ") is larger than the number of solutes (", carma%f_NSOLUTE, ")."
      rc = RC_ERROR
      return
    end if

    ! Save off the settings.
    carma%f_solute(isolute)%f_name        = name
    carma%f_solute(isolute)%f_ions        = ions
    carma%f_solute(isolute)%f_wtmol       = wtmol
    carma%f_solute(isolute)%f_rho         = rho
    
    
    ! Defaults for optional parameters
    carma%f_solute(isolute)%f_shortname   = ""
    
    ! Set optional parameters.
    if (present(shortname))  carma%f_solute(isolute)%f_shortname    = shortname

    return
  end subroutine CARMASOLUTE_Create
    

  !! Deallocates the memory associated with a CARMASOLUTE object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMASOLUTE_Create
  subroutine CARMASOLUTE_Destroy(carma, isolute, rc)
    type(carma_type), intent(inout)             :: carma         !! the carma object
    integer, intent(in)                         :: isolute       !! the solute index
    integer, intent(out)                        :: rc            !! return code, negative indicates failure

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough solutes allocated.
    if (isolute > carma%f_NSOLUTE) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMASOLUTE_Destroy:: ERROR - The specifed solute (", &
        isolute, ") is larger than the number of solutes (", carma%f_NSOLUTE, ")."
      rc = RC_ERROR
      return
    end if

    return
  end subroutine CARMASOLUTE_Destroy


  !! Gets information about a solute.
  !!
  !! The group name and other properties are available after a call to
  !! CARMASOLUTE_Create().
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMASOLUTE_Create
  !! @see CARMA_GetGas
  subroutine CARMASOLUTE_Get(carma, isolute, rc, name, shortname, ions, wtmol, rho)
    type(carma_type), intent(in)                :: carma        !! the carma object
    integer, intent(in)                         :: isolute      !! the solute index
    integer, intent(out)                        :: rc           !! return code, negative indicates failure
    character(len=*), optional, intent(out)     :: name         !! the solute name
    character(len=*), optional, intent(out)     :: shortname    !! the solute short name
		integer, optional, intent(out)              :: ions         !! Number of ions solute dissociates into
    real(kind=f), optional, intent(out)         :: wtmol        !! the solute molecular weight [g/mol]
		real(kind=f), optional, intent(out)         :: rho          !! Mass density of solute
    
    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough solutes allocated.
    if (isolute > carma%f_NSOLUTE) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMASOLUTE_Get:: ERROR - The specifed solute (", &
        isolute, ") is larger than the number of solutes (", carma%f_NSOLUTE, ")."
      rc = RC_ERROR
      return
    end if

    ! Return any requested properties of the group.
    if (present(name))         name         = carma%f_solute(isolute)%f_name
    if (present(shortname))    shortname    = carma%f_solute(isolute)%f_shortname
    if (present(ions))         ions         = carma%f_solute(isolute)%f_ions
    if (present(wtmol))        wtmol        = carma%f_solute(isolute)%f_wtmol
    if (present(rho))          rho          = carma%f_solute(isolute)%f_rho
        
    return
  end subroutine CARMASOLUTE_Get
  
  
  !! Prints information about a solute.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMASOLUTE_Get
  subroutine CARMASOLUTE_Print(carma, isolute, rc)
    type(carma_type), intent(in)              :: carma         !! the carma object
    integer, intent(in)                       :: isolute       !! the solute index
    integer, intent(out)                      :: rc            !! return code, negative indicates failure
    
    ! Local variables
    character(len=CARMA_NAME_LEN)             :: name          !! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname     !! shortname
		integer                                   :: ions          !! Number of ions solute dissociates into
    real(kind=f)                              :: wtmol         !! the solute molecular weight [g/mol]
		real(kind=f)                              :: rho           !! Mass density of solute

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%f_do_print) then
      call CARMASOLUTE_Get(carma, isolute, rc, name=name, shortname=shortname, ions=ions, wtmol=wtmol, rho=rho) 
      if (rc < 0) return

    
      write(carma%f_LUNOPRT,*) "    name          : ", trim(name)
      write(carma%f_LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%f_LUNOPRT,*) "    ions          : ", ions
      write(carma%f_LUNOPRT,*) "    wtmol         : ", wtmol, " (g/mol)"
      write(carma%f_LUNOPRT,*) "    rho           : ", rho, " (g/cm3)"
    end if
    
    return
  end subroutine CARMASOLUTE_Print
end module

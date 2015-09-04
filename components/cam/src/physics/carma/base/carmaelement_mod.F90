!! The CARMAELEMENT module contains configuration information about a particle
!! element used by CARMA.
!!
!!  @version March-2010
!!  @author  Chuck Bardeen 
module CARMAELEMENT_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAELEMENT_Create
  public CARMAELEMENT_Destroy
  public CARMAELEMENT_Get
  public CARMAELEMENT_Print

contains

  !! Defines a gas used by CARMA for nucleation and growth of cloud and 
  !! aerosol particles.
  !!
  !! NOTE: The element density can be specifeid per bin using rhobin; however,
  !! if only the bulk density is provided (rho) then the same value will be used
  !! for all bins. The bulk density allows for backward compatability and ease of
  !! configuration. If rhobin is provided, then rho is ignored.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMA_AddGas
  !! @see CARMAELEMENT_Destroy
 subroutine CARMAELEMENT_Create(carma, ielement, igroup, name, rho, itype, icomposition, rc, &
              shortname, isolute, rhobin, arat)
    type(carma_type), intent(inout)       :: carma               !! the carma object
    integer, intent(in)                   :: ielement            !! the element index
    integer, intent(in)                   :: igroup              !! Group to which the element belongs
    character(*), intent(in)              :: name                !! the element name, maximum of 255 characters
    real(kind=f), intent(in)              :: rho                 !! bulk mass density of particle element [g/cm^3]
    integer, intent(in)                   :: itype               !! Particle type specification
    integer, intent(in)                   :: icomposition        !! Particle compound specification
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    character(*), optional, intent(in)    :: shortname           !! the element shortname, maximum of 6 characters
    integer, optional, intent(in)         :: isolute             !! Index of solute for the particle element
    real(kind=f), optional, intent(in)    :: rhobin(carma%f_NBIN)  !! mass density per bin of particle element [g/cm^3]
    real(kind=f), optional, intent(in)    :: arat(carma%f_NBIN)    !! projected area ratio

    ! Local variables
    integer                               :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough elements allocated.
    if (ielement > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    allocate( &
      carma%f_element(ielement)%f_rho(carma%f_NBIN), &
      stat=ier) 
    if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Add: ERROR allocating, status=", ier
      rc = RC_ERROR
      return
    end if

    ! Save off the settings.
    carma%f_element(ielement)%f_igroup       = igroup
    carma%f_element(ielement)%f_name         = name
    carma%f_element(ielement)%f_rho(:)       = rho
    carma%f_element(ielement)%f_itype        = itype
    carma%f_element(ielement)%f_icomposition = icomposition
    
    
    ! Defaults for optional parameters
    carma%f_element(ielement)%f_shortname   = ""
    carma%f_element(ielement)%f_isolute     = 0
    
    ! Set optional parameters.
    if (present(shortname))  carma%f_element(ielement)%f_shortname    = shortname
    if (present(isolute)) then
    
      ! Make sure there are enough solutes allocated.
      if (isolute > carma%f_NSOLUTE) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Create:: ERROR - The specifed solute (", &
          isolute, ") is larger than the number of solutes (", carma%f_NSOLUTE, ")."
        rc = RC_ERROR
        return
      end if

      carma%f_element(ielement)%f_isolute      = isolute
    end if
    if (present(rhobin)) carma%f_element(ielement)%f_rho(:) = rhobin(:)
    
    ! If the area ratio is specfied (usually along with rhobin), then set this
    ! for the group.
    if (present(arat)) carma%f_group(igroup)%f_arat(:) = arat(:)
    
    ! Keep track of the fact that another element has been added to the group.
    carma%f_group(igroup)%f_nelem = carma%f_group(igroup)%f_nelem + 1
    
    return
  end subroutine CARMAELEMENT_Create
    

  !! Deallocates the memory associated with a CARMAELEMENT object.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Create
  subroutine CARMAELEMENT_Destroy(carma, ielement, rc)
    type(carma_type), intent(inout)        :: carma         !! the carma object
    integer, intent(in)                    :: ielement      !! the element index
    integer, intent(out)                   :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                               :: ier

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough elements allocated.
    if (ielement > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Destroy:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (allocated(carma%f_element(ielement)%f_rho)) then
      deallocate( &
        carma%f_element(ielement)%f_rho, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif

    return
  end subroutine CARMAELEMENT_Destroy


  !! Gets information about a particle element.
  !!
  !! The group name and other properties are available after a call to
  !! CARMAELEMENT_Create().
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Create
  !! @see CARMA_GetElement
  subroutine CARMAELEMENT_Get(carma, ielement, rc, igroup, name, shortname, rho, itype, icomposition, isolute)
    type(carma_type), intent(in)                :: carma           !! the carma object
    integer, intent(in)                         :: ielement        !! the element index
    integer, intent(out)                        :: rc              !! return code, negative indicates failure
    integer, optional, intent(out)              :: igroup          !! Group to which the element belongs
    character(len=*), optional, intent(out)     :: name            !! the element name
    character(len=*), optional, intent(out)     :: shortname       !! the element short name
    real(kind=f), optional, intent(out)         :: rho(carma%f_NBIN) !! Mass density of particle element [g/cm^3]
    integer, optional, intent(out)              :: itype           !! Particle type specification
    integer, optional, intent(out)              :: icomposition    !! Particle compound specification
    integer, optional, intent(out)              :: isolute         !! Index of solute for the particle element
    
    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough elements allocated.
    if (ielement > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAELEMENT_Get:: ERROR - The specifed element (", &
        ielement, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if

    ! Return any requested properties of the group.
    if (present(igroup))       igroup       = carma%f_element(ielement)%f_igroup
    if (present(name))         name         = carma%f_element(ielement)%f_name
    if (present(shortname))    shortname    = carma%f_element(ielement)%f_shortname
    if (present(rho))          rho(:)       = carma%f_element(ielement)%f_rho(:)
    if (present(itype))        itype        = carma%f_element(ielement)%f_itype
    if (present(icomposition)) icomposition = carma%f_element(ielement)%f_icomposition
    if (present(isolute))      isolute      = carma%f_element(ielement)%f_isolute
        
    return
  end subroutine CARMAELEMENT_Get
  
  
  !! Prints information about an element.
  !!
  !! @author  Chuck Bardeen
  !! @version March-2010
  !!
  !! @see CARMAELEMENT_Get
  subroutine CARMAELEMENT_Print(carma, ielement, rc)
    type(carma_type), intent(in)              :: carma         !! the carma object
    integer, intent(in)                       :: ielement      !! the element index
    integer, intent(out)                      :: rc            !! return code, negative indicates failure
    
    ! Local variables
    character(len=CARMA_NAME_LEN)             :: name             ! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname        ! shortname
    real(kind=f)                              :: rho(carma%f_NBIN)  ! density (g/cm3)
    integer                                   :: igroup           ! Group to which the element belongs
    integer                                   :: itype            ! Particle type specification
    integer                                   :: icomposition     ! Particle compound specification
    integer                                   :: isolute          ! Index of solute for the particle element

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%f_do_print) then
      call CARMAELEMENT_Get(carma, ielement, rc, name=name, shortname=shortname, igroup=igroup, &
                            itype=itype, icomposition=icomposition, rho=rho, isolute=isolute) 
      if (rc < 0) return

    
      write(carma%f_LUNOPRT,*) "    name          : ", trim(name)
      write(carma%f_LUNOPRT,*) "    igroup        : ", igroup
      write(carma%f_LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%f_LUNOPRT,*) "    rho           : ", rho, " (g/cm3)"

      select case(itype)
        case (I_INVOLATILE)
          write(carma%f_LUNOPRT,*) "    itype         :    involatile"
        case (I_VOLATILE)
          write(carma%f_LUNOPRT,*) "    itype         :    volatile"
        case (I_COREMASS)
          write(carma%f_LUNOPRT,*) "    itype         :    core mass"
        case (I_VOLCORE)
          write(carma%f_LUNOPRT,*) "    itype         :    volatile core"
        case (I_CORE2MOM)
          write(carma%f_LUNOPRT,*) "    itype         :    core mass - second moment"
        case default
          write(carma%f_LUNOPRT,*) "    itype         :    unknown, ", itype
      end select

      write(carma%f_LUNOPRT,*) "    icomposition  : ", icomposition
      write(carma%f_LUNOPRT,*) "    isolute       : ", isolute
    end if
    
    return
  end subroutine CARMAELEMENT_Print
end module

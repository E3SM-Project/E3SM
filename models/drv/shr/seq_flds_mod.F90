module seq_flds_mod

  !====================================================================
  !   New standardized naming convention
  !====================================================================
  !  
  !  ---------
  !  definitions:
  !  ---------
  !  state-prefix
  !    first 3 characters: Sx_, Sa_, Si_, Sl_, So_ 
  !    one letter indices: x,a,l,i,o,g,r 
  !    x => coupler (mapping, merging, atm/ocn flux calc done on coupler procs)
  !    a => atm
  !    l => lnd
  !    i => ice
  !    o => ocn
  !    g => glc
  !    r => rof
  !    w => wav
  !
  !  state-name 
  !    what follows state prefix
  !
  !  flux-prefix
  !    first 5 characters: Flmn__ 
  !    lm => between components l and m
  !    n  => computed by component n
  !    example: Fioi => ice/ocn flux computed by ice
  !    example: Fall => atm/lnd flux computed by lnd
  !    If flux prefix has first letter of P (so first five characters are PFlmn_)
  !    then flux is passed straight through without scaling by the corresponding fraction)
  !    
  !  flux-name
  !    what follows flux-prefix
  !
  !  ---------
  !  rules:
  !  ---------
  !  1) states: 
  !     a) atm attributes fields that HAVE a state-prefix of Sx_ in seq_flds_x2a_states
  !        rule: will merge all identical values of the state-names from
  !           seq_flds_i2x_states 
  !           seq_flds_l2x_states 
  !           seq_flds_o2x_states 
  !           seq_flds_xao_states
  !         to obtain output state-name in seq_flds_x2a_states
  !  
  !        rule: to merge input states that originate in the 
  !           lnd (l2x_a) will be scaled by the lndfrac
  !           ice (i2x_a) will be scaled by the icefrac
  !           cpl (xao_a) will be scaled by the ocnfrac
  !           ocn (o2x_a) will be scaled by the ocnfrac
  !  
  !        example: 
  !           seq_flds_l2x_states = "Sl_t"
  !           seq_flds_i2x_states = "Si_t"
  !           seq_flds_o2x_states = "So_t"
  !           seq_flds_x2a_states = "Sx_t" 
  !           attribute fields Sl_t, Si_t, So_t, in 
  !           attribute vectors l2x_a, i2x_a, o2x_a will be
  !           merged to obtain attribute Sx_t in attribute vector x2a_a
  !  
  !     b) atm attribute fields that DO NOT HAVE a state-prefix of Sx_ in seq_flds_x2a_states
  !        rule: copy directly all variables that identical state-prefix 
  !               AND state-name in
  !           seq_flds_i2x_states and seq_flds_x2a_states
  !           seq_flds_l2x_states and seq_flds_x2a_states
  !           seq_flds_o2x_states and seq_flds_x2a_states
  !           seq_flds_xao_states and seq_flds_x2a_states
  !  
  !        example 
  !           seq_flds_i2x_states = ":Si_snowh"
  !           seq_flds_x2a_states = ":Si_snowh"
  !           attribute field of Si_snowh in i2x_a will be copied to 
  !           attribute field Si_snowh in x2a_a
  !  
  !  2) fluxes: 
  !     rule: will merge all identical values of the flux-names from
  !         seq_flds_i2x_states 
  !         seq_flds_l2x_states 
  !         seq_flds_o2x_states 
  !         seq_flds_xao_states
  !       to obtain output state-name in seq_flds_x2a_states
  !  
  !     rule: input flux fields that originate in the 
  !         lnd (l2x_a) will be scaled by the lndfrac
  !         ice (i2x_a) will be scaled by the icefrac
  !            - ignore all fluxes that are ice/ocn fluxes (e.g. Fioi_)
  !         cpl (xao_a) will be scaled by the ocnfrac
  !         ocn (o2x_a) will be scaled by the ocnfrac+icefrac
  !
  !====================================================================
  !
  !   New user specified fields
  ! 
  !====================================================================
  ! New fields that are user specidied can be added as namelist variables
  ! by the user in the cpl namelist seq_flds_user using the namelist variable
  ! array cplflds_customs. The user specified new fields must follow the
  ! above naming convention.
  ! As an example, say you want to add a new state 'foo' that is passed 
  ! from the land to the atm - you would do this as follows 
  !    &seq_flds_user
  !       cplflds_custom = 'Sa_foo->a2x', 'Sa_foo->x2a'
  !    /
  ! This would add the field 'Sa_foo' to the character strings defining the 
  ! attribute vectors a2x and x2a. It is assumed that code would need to be 
  ! introduced in the atm and land components to deal with this new attribute 
  ! vector field.
  ! Currently, the only way to add this is to edit $CASEROOT/user_nl_cpl
  !====================================================================
  !
  !   Coupler fields use cases
  !
  !====================================================================
  ! Previously, new fields that were needed to be passed between components
  ! for certain compsets were specified by cpp-variables. This has been 
  ! modified to now be use cases. The use cases are specified in the 
  ! namelist cpl_flds_inparm and are currently triggered by the xml 
  ! variables CCSM_VOC, CCSM_BGC and GLC_NEC.  
  !====================================================================
 
   use shr_kind_mod,   only : CX => shr_kind_CX, CXX => shr_kind_CXX
   use shr_sys_mod,    only : shr_sys_abort
   use seq_drydep_mod, only : seq_drydep_init, seq_drydep_read, lnd_drydep
   use seq_comm_mct,   only : seq_comm_iamroot, seq_comm_setptrs, logunit
   use shr_megan_mod,  only : shr_megan_readnl, shr_megan_mechcomps_n
   use shr_carma_mod,  only : shr_carma_readnl

   implicit none
   public
   save

   interface seq_flds_lookup; module procedure &
     seq_flds_esmf_metadata_get
   end interface

   integer, parameter, private :: CSS = 256  ! use longer short character
   integer, parameter, private :: CLL = 1024
   character(len=CXX) :: seq_drydep_fields   ! List of dry-deposition fields
   character(len=CXX) :: megan_voc_fields    ! List of MEGAN VOC emission fields
   character(len=CX)  :: carma_fields        ! List of CARMA fields from lnd->atm
   integer            :: seq_flds_glc_nec    ! number of glc elevation classes

   !----------------------------------------------------------------------------
   ! metadata
   !----------------------------------------------------------------------------

   character(len=*),parameter :: undef     = 'undefined'
   integer         ,parameter :: nmax      = 1000        ! maximum number of entries in lookup_entry
   integer                    :: n_entries = 0           ! actual number of entries in lookup_entry
   character(len=80), dimension(nmax, 4) :: lookup_entry = undef

   !----------------------------------------------------------------------------
   ! for the domain
   !----------------------------------------------------------------------------

   character(CXX) :: seq_flds_dom_coord 
   character(CXX) :: seq_flds_dom_other

   !----------------------------------------------------------------------------
   ! state + flux fields
   !----------------------------------------------------------------------------

   character(CXX) :: seq_flds_a2x_states 
   character(CXX) :: seq_flds_a2x_fluxes 
   character(CXX) :: seq_flds_x2a_states 
   character(CXX) :: seq_flds_x2a_fluxes

   character(CXX) :: seq_flds_i2x_states 
   character(CXX) :: seq_flds_i2x_fluxes 
   character(CXX) :: seq_flds_x2i_states 
   character(CXX) :: seq_flds_x2i_fluxes

   character(CXX) :: seq_flds_l2x_states 
   character(CXX) :: seq_flds_l2x_fluxes 
   character(CXX) :: seq_flds_x2l_states 
   character(CXX) :: seq_flds_x2l_fluxes

   character(CXX) :: seq_flds_o2x_states 
   character(CXX) :: seq_flds_o2x_fluxes 
   character(CXX) :: seq_flds_x2o_states 
   character(CXX) :: seq_flds_x2o_fluxes

   character(CXX) :: seq_flds_g2x_states 
   character(CXX) :: seq_flds_g2x_fluxes 
   character(CXX) :: seq_flds_x2g_states 
   character(CXX) :: seq_flds_x2g_fluxes

   character(CXX) :: seq_flds_w2x_states 
   character(CXX) :: seq_flds_w2x_fluxes 
   character(CXX) :: seq_flds_x2w_states 
   character(CXX) :: seq_flds_x2w_fluxes

   character(CXX) :: seq_flds_xao_albedo
   character(CXX) :: seq_flds_xao_states 
   character(CXX) :: seq_flds_xao_fluxes

   character(CXX) :: seq_flds_r2x_states 
   character(CXX) :: seq_flds_r2x_fluxes
   character(CXX) :: seq_flds_x2r_states 
   character(CXX) :: seq_flds_x2r_fluxes

   !----------------------------------------------------------------------------
   ! combined state/flux fields
   !----------------------------------------------------------------------------

   character(CXX) :: seq_flds_dom_fields 
   character(CXX) :: seq_flds_a2x_fields 
   character(CXX) :: seq_flds_x2a_fields 
   character(CXX) :: seq_flds_i2x_fields 
   character(CXX) :: seq_flds_x2i_fields 
   character(CXX) :: seq_flds_l2x_fields 
   character(CXX) :: seq_flds_x2l_fields 
   character(CXX) :: seq_flds_o2x_fields 
   character(CXX) :: seq_flds_x2o_fields 
   character(CXX) :: seq_flds_xao_fields 
   character(CXX) :: seq_flds_r2x_fields
   character(CXX) :: seq_flds_x2r_fields
   character(CXX) :: seq_flds_g2x_fields 
   character(CXX) :: seq_flds_x2g_fields 
   character(CXX) :: seq_flds_w2x_fields 
   character(CXX) :: seq_flds_x2w_fields 

   !----------------------------------------------------------------------------
   ! component names
   !----------------------------------------------------------------------------

   character(32) :: atmname='atm'
   character(32) :: ocnname='ocn'
   character(32) :: icename='ice'
   character(32) :: lndname='lnd'
   character(32) :: glcname='glc'
   character(32) :: wavname='wav'
   character(32) :: rofname='rof'

!----------------------------------------------------------------------------
 contains
!----------------------------------------------------------------------------

   subroutine seq_flds_set(nmlfile, ID)

! !USES:
     use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit
     use shr_string_mod, only : shr_string_listIntersect
     use shr_mpi_mod,    only : shr_mpi_bcast

! !INPUT/OUTPUT PARAMETERS:
     character(len=*), intent(in) :: nmlfile   ! Name-list filename
     integer         , intent(in) :: ID        ! seq_comm ID

     !----- local -----
     integer :: mpicom             ! MPI communicator
     integer :: ierr               ! I/O error code
     integer :: unitn              ! Namelist unit number to read

     character(len=CSS) :: attname
     character(len=CSS) :: units
     character(len=CSS) :: longname
     character(len=CSS) :: stdname
     integer            :: num
     character(len=  2) :: cnum
     character(len=CSS) :: name

     character(CXX) :: dom_coord  = ''
     character(CXX) :: dom_other  = ''

     character(CXX) :: a2x_states = ''
     character(CXX) :: a2x_fluxes = ''
     character(CXX) :: x2a_states = ''
     character(CXX) :: x2a_fluxes = ''
     character(CXX) :: i2x_states = ''
     character(CXX) :: i2x_fluxes = ''
     character(CXX) :: x2i_states = ''
     character(CXX) :: x2i_fluxes = ''
     character(CXX) :: l2x_states = ''
     character(CXX) :: l2x_fluxes = ''
     character(CXX) :: x2l_states = ''
     character(CXX) :: x2l_fluxes = ''
     character(CXX) :: o2x_states = ''
     character(CXX) :: o2x_fluxes = ''
     character(CXX) :: x2o_states = ''
     character(CXX) :: x2o_fluxes = ''
     character(CXX) :: g2x_states = ''
     character(CXX) :: g2x_fluxes = ''
     character(CXX) :: x2g_states = ''
     character(CXX) :: x2g_fluxes = ''
     character(CXX) :: xao_albedo = ''
     character(CXX) :: xao_states = ''
     character(CXX) :: xao_fluxes = ''
     character(CXX) :: r2x_states = ''
     character(CXX) :: r2x_fluxes = ''
     character(CXX) :: x2r_states = ''
     character(CXX) :: x2r_fluxes = ''
     character(CXX) :: w2x_states = ''
     character(CXX) :: w2x_fluxes = ''
     character(CXX) :: x2w_states = ''
     character(CXX) :: x2w_fluxes = ''

     character(CXX) :: stringtmp  = ''

     !------ namelist -----
     character(len=CSS)  :: fldname, fldflow
     logical :: is_state, is_flux
     integer :: i,n

     ! use cases namelists
     logical :: flds_co2a 
     logical :: flds_co2b 
     logical :: flds_co2c 
     logical :: flds_co2_dmsa 
     integer :: glc_nec

     namelist /seq_cplflds_inparm/  &
          flds_co2a, flds_co2b, flds_co2c, flds_co2_dmsa, glc_nec

     ! user specified new fields
     integer,  parameter :: nfldmax = 200
     character(len=CLL)  :: cplflds_custom(nfldmax) = ''

     namelist /seq_cplflds_userspec/ &          
          cplflds_custom

     character(len=*),parameter :: subname = '(seq_flds_set) '

!-------------------------------------------------------------------------------

     call seq_comm_setptrs(ID,mpicom=mpicom)

     !---------------------------------------------------------------------------
     ! Read in namelist for use cases
     !---------------------------------------------------------------------------
     ! TODO: permit duplicates to occur - then check for this in seq_flds_add
     ! TODO: add entries for lookup entry table for custom fields 
     !---------------------------------------------------------------------------
     
     if (seq_comm_iamroot(ID)) then
        flds_co2a = .false.
        flds_co2b = .false.
        flds_co2c = .false.
        flds_co2_dmsa = .false.
        glc_nec   = 0

        unitn = shr_file_getUnit()
        write(logunit,"(A)") subname//': read seq_cplflds_inparm namelist from: '&
             //trim(nmlfile)
        open( unitn, file=trim(nmlfile), status='old' )
        ierr = 1
        do while( ierr /= 0 )
           read(unitn,nml=seq_cplflds_inparm,iostat=ierr)
           if (ierr < 0) then
              call shr_sys_abort( &
                   subname//"ERROR: namelist read returns an EOF or EOR condition" )
           end if
        end do
        close(unitn)
        call shr_file_freeUnit( unitn )
     end if
     call shr_mpi_bcast(flds_co2a    , mpicom)
     call shr_mpi_bcast(flds_co2b    , mpicom)
     call shr_mpi_bcast(flds_co2c    , mpicom)
     call shr_mpi_bcast(flds_co2_dmsa, mpicom)
     call shr_mpi_bcast(glc_nec      , mpicom)
     seq_flds_glc_nec = glc_nec

     !---------------------------------------------------------------------------
     ! Read in namelists for user specified new fields
     !---------------------------------------------------------------------------
     ! TODO: permit duplicates to occur - then check for this in seq_flds_add
     ! TODO: add entries for lookup entry table for custom fields 
     !---------------------------------------------------------------------------

     if (seq_comm_iamroot(ID)) then
        cplflds_custom(:) = ' '

        unitn = shr_file_getUnit()
        write(logunit,"(A)") subname//': read seq_cplflds_userspec namelist from: '&
             //trim(nmlfile)
        open( unitn, file=trim(nmlfile), status='old' )
        ierr = 1
        do while( ierr /= 0 )
           read(unitn,nml=seq_cplflds_userspec,iostat=ierr)
           if (ierr < 0) then
              call shr_sys_abort( &
                   subname//"ERROR: namelist read returns an EOF or EOR condition" )
           end if
        end do
        close(unitn)
        call shr_file_freeUnit( unitn )
     end if
     do n = 1, nfldmax
        call shr_mpi_bcast(cplflds_custom(n), mpicom)
     end do
        
     ! add customized fields through coupler

     do n = 1,nfldmax
        if (cplflds_custom(n) /= ' ') then
           i = scan(cplflds_custom(n),'->')
           fldname = trim(adjustl(cplflds_custom(n)(:i-1)))
           fldflow = trim(adjustl(cplflds_custom(n)(i+2:)))

           if (fldname(1:1) == 'S') then
              is_state = .true.
              is_flux  = .false.
           else if (fldname (1:1) == 'F')  then
              is_state = .false.
              is_flux  = .true.
           else if (fldname (1:2) == 'PF') then
              is_state = .false.
              is_flux  = .true.
           else
              write(logunit,*) subname//'ERROR: fldname must start with S,F,P, not ',trim(fldname)
              call shr_sys_abort(subname//"ERROR: fldname must start with S, F, or P")
           end if

           select case (trim(fldflow))
           case('a2x')
              if (is_state) call seq_flds_add(a2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(a2x_fluxes,trim(fldname))
           case('x2a')
              if (is_state) call seq_flds_add(x2a_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2a_fluxes,trim(fldname))
           case('l2x')
              if (is_state) call seq_flds_add(l2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(l2x_fluxes,trim(fldname))
           case('x2l')
              if (is_state) call seq_flds_add(x2l_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2l_fluxes,trim(fldname))
           case('r2x')
              if (is_state) call seq_flds_add(r2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(r2x_fluxes,trim(fldname))
           case('x2r')
              if (is_state) call seq_flds_add(x2r_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2r_fluxes,trim(fldname))
           case('i2x')
              if (is_state) call seq_flds_add(i2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(i2x_fluxes,trim(fldname))
           case('x2i')
              if (is_state) call seq_flds_add(x2i_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2i_fluxes,trim(fldname))
           case('o2x')
              if (is_state) call seq_flds_add(o2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(o2x_fluxes,trim(fldname))
           case('x2o')
              if (is_state) call seq_flds_add(x2o_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2o_fluxes,trim(fldname))
           case('g2x')
              if (is_state) call seq_flds_add(g2x_states,trim(fldname))
              if (is_flux ) call seq_flds_add(g2x_fluxes,trim(fldname))
           case('x2g')
              if (is_state) call seq_flds_add(x2g_states,trim(fldname))
              if (is_flux ) call seq_flds_add(x2g_fluxes,trim(fldname))
           case default
              write(logunit,*) subname//'ERROR: ',trim(cplflds_custom(n)),&
                   ' not a recognized value'
              call shr_sys_abort()
           end select
        else
           exit
        end if
     end do

     !----------------------------------------------------------
     ! domain coordinates
     !----------------------------------------------------------

     call seq_flds_add(dom_coord,'lat')
     longname = ''
     stdname  = 'latitude'
     units    = 'degrees north'
     attname  = 'lat' 
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_coord,'lon')
     longname = ''
     stdname  = 'longitude'
     units    = 'degrees east'
     attname  = 'lon' 
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'area')
     longname = ''
     stdname  = 'cell area'
     units    = 'radian^2'
     attname  = 'area' 
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'aream')
     longname = ''
     stdname  = 'cell area from mapping file'
     units    = 'radian^2'
     attname  = 'aream'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'mask')
     longname = ''
     stdname  = 'mask'
     units    = 'unitless'
     attname  = 'mask'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'frac')
     longname = 'area_fraction'
     stdname  = 'area fraction'
     units    = 'unitless'
     attname  = 'frac' 
     call metadata_set(attname, longname, stdname, units)

     !----------------------------------------------------------
     ! states/fluxes from atm 
     !----------------------------------------------------------
     
     ! height at the lowest model level (m)
     call seq_flds_add(a2x_states,"Sa_z")
     call seq_flds_add(x2l_states,"Sa_z")    
     call seq_flds_add(x2i_states,"Sa_z")
     longname = 'Height at the lowest model level'
     stdname  = 'height'
     units    = 'm'
     attname  = 'Sa_z'
     call metadata_set(attname, longname, stdname, units)

     ! zonal wind at the lowest model level (m/s)
     call seq_flds_add(a2x_states,"Sa_u")
     call seq_flds_add(x2l_states,"Sa_u")
     call seq_flds_add(x2i_states,"Sa_u")
     call seq_flds_add(x2w_states,"Sa_u")
     longname = 'Zonal wind at the lowest model level'
     stdname  = 'eastward_wind'
     units    = 'm s-1'
     attname  = 'Sa_u'
     call metadata_set(attname, longname, stdname, units)

     ! meridional wind at the lowest model level (m/s)
     call seq_flds_add(a2x_states,"Sa_v")
     call seq_flds_add(x2l_states,"Sa_v")
     call seq_flds_add(x2i_states,"Sa_v")
     call seq_flds_add(x2w_states,"Sa_v")
     longname = 'Meridional wind at the lowest model level'
     stdname  = 'northward_wind'
     units    = 'm s-1'
     attname  = 'Sa_v'
     call metadata_set(attname, longname, stdname, units)

     ! temperature at the lowest model level (K)
     call seq_flds_add(a2x_states,"Sa_tbot")
     call seq_flds_add(x2l_states,"Sa_tbot")
     call seq_flds_add(x2i_states,"Sa_tbot")
     call seq_flds_add(x2w_states,"Sa_tbot")
     longname = 'Temperature at the lowest model level'
     stdname  = 'air_temperature'
     units    = 'K'
     attname  = 'Sa_tbot'
     call metadata_set(attname, longname, stdname, units)

     ! potential temperature at the lowest model level (K)
     call seq_flds_add(a2x_states,"Sa_ptem")
     call seq_flds_add(x2l_states,"Sa_ptem")
     call seq_flds_add(x2i_states,"Sa_ptem")
     longname = 'Potential temperature at the lowest model level'
     stdname  = 'air_potential_temperature'
     units    = 'K'
     attname  = 'Sa_ptem'
     call metadata_set(attname, longname, stdname, units)

     ! ppecific humidity at the lowest model level (kg/kg)
     call seq_flds_add(a2x_states,"Sa_shum")
     call seq_flds_add(x2l_states,"Sa_shum")
     call seq_flds_add(x2i_states,"Sa_shum")
     longname = 'Specific humidity at the lowest model level'
     stdname  = 'specific_humidity'
     units    = 'kg kg-1'
     attname  = 'Sa_shum'
     call metadata_set(attname, longname, stdname, units)

     ! pressure at the lowest model level (Pa)
     call seq_flds_add(a2x_states,"Sa_pbot")
     call seq_flds_add(x2l_states,"Sa_pbot")
     call seq_flds_add(x2i_states,"Sa_pbot")
     longname = 'Pressure at the lowest model level'
     stdname  = 'air_pressure'
     units    = 'Pa'
     attname  = 'Sa_pbot'
     call metadata_set(attname, longname, stdname, units)

     ! air density at the lowest model level (kg/m**3)
     call seq_flds_add(a2x_states,"Sa_dens")
     call seq_flds_add(x2i_states,"Sa_dens")    
     longname = 'Density at the lowest model level'
     stdname  = 'air_density'
     units    = 'kg m-3'
     attname  = 'Sa_dens'
     call metadata_set(attname, longname, stdname, units)

     ! convective precipitation rate
     ! large-scale (stable) snow rate (water equivalent)
     call seq_flds_add(a2x_fluxes,"Faxa_rainc")
     call seq_flds_add(a2x_fluxes,"Faxa_rainl")
     call seq_flds_add(x2l_fluxes,"Faxa_rainc")
     call seq_flds_add(x2l_fluxes,"Faxa_rainl")
     call seq_flds_add(x2i_fluxes,"Faxa_rain" )
     call seq_flds_add(x2o_fluxes,"Faxa_rain" )
     units    = 'kg m-2 s-1'
     longname = 'Convective precipitation rate'
     stdname  = 'convective_precipitation_flux'
     attname  = 'Faxa_rainc'
     call metadata_set(attname, longname, stdname, units)
     longname = 'Large-scale (stable) precipitation rate'
     stdname  = 'large_scale_precipitation_flux'
     attname  = 'Faxa_rainl' 
     call metadata_set(attname, longname, stdname, units)
     longname = 'Water flux due to rain'
     stdname  = 'rainfall_flux'
     attname  = 'Faxa_rain'  
     call metadata_set(attname, longname, stdname, units)

     ! convective snow rate (water equivalent)
     ! large-scale (stable) snow rate (water equivalent)
     call seq_flds_add(a2x_fluxes,"Faxa_snowc")
     call seq_flds_add(a2x_fluxes,"Faxa_snowl")
     call seq_flds_add(x2l_fluxes,"Faxa_snowc")
     call seq_flds_add(x2l_fluxes,"Faxa_snowl")
     call seq_flds_add(x2i_fluxes,"Faxa_snow" )
     call seq_flds_add(x2o_fluxes,"Faxa_snow" )
     units    = 'kg m-2 s-1'
     longname = 'Convective snow rate (water equivalent)'
     stdname  = 'convective_snowfall_flux'
     attname  = 'Faxa_snowc'
     call metadata_set(attname, longname, stdname, units)
     longname = 'Large-scale (stable) snow rate (water equivalent)'
     stdname  = 'large_scale_snowfall_flux'
     attname  = 'Faxa_snowl' 
     call metadata_set(attname, longname, stdname, units)
     longname = 'Water flux due to snow'
     stdname  = 'surface_snow_melt_flux'
     attname  = 'Faxa_snow'  
     call metadata_set(attname, longname, stdname, units)
     
     ! total precipitation to ocean
     call seq_flds_add(x2o_fluxes,"Faxa_prec")  ! derived rain+snow
     longname = 'Water flux (rain+snow)'
     stdname  = 'precipitation_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_prec'
     call metadata_set(attname, longname, stdname, units)

     ! downward longwave heat flux (W/m**2)
     call seq_flds_add(a2x_fluxes,"Faxa_lwdn")
     call seq_flds_add(x2l_fluxes,"Faxa_lwdn")
     call seq_flds_add(x2i_fluxes,"Faxa_lwdn")
     call seq_flds_add(x2o_fluxes,"Faxa_lwdn")
     longname = 'Downward longwave heat flux'
     stdname  = 'downwelling_longwave_flux'
     units    = 'W m-2'
     attname  = 'Faxa_lwdn'
     call metadata_set(attname, longname, stdname, units)

     ! direct near-infrared incident solar radiation
     call seq_flds_add(a2x_fluxes,"Faxa_swndr")
     call seq_flds_add(x2i_fluxes,"Faxa_swndr")
     call seq_flds_add(x2l_fluxes,"Faxa_swndr")
     longname = 'Direct near-infrared incident solar radiation'
     stdname  = 'surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation'
     units    = 'W m-2'
     attname  = 'Faxa_swndr'
     call metadata_set(attname, longname, stdname, units)

     ! direct visible incident solar radiation
     call seq_flds_add(a2x_fluxes,"Faxa_swvdr")
     call seq_flds_add(x2i_fluxes,"Faxa_swvdr")
     call seq_flds_add(x2l_fluxes,"Faxa_swvdr")
     longname = 'Direct visible incident solar radiation'
     stdname  = 'surface_downward_direct_shortwave_flux_due_to_visible_radiation'
     units    = 'W m-2'
     attname  = 'Faxa_swvdr'
     call metadata_set(attname, longname, stdname, units)

     ! diffuse near-infrared incident solar radiation
     call seq_flds_add(a2x_fluxes,"Faxa_swndf")
     call seq_flds_add(x2i_fluxes,"Faxa_swndf")
     call seq_flds_add(x2l_fluxes,"Faxa_swndf")
     longname = 'Diffuse near-infrared incident solar radiation'
     stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation'
     units    = 'W m-2'
     attname  = 'Faxa_swndf'
     call metadata_set(attname, longname, stdname, units)

     ! diffuse visible incident solar radiation
     call seq_flds_add(a2x_fluxes,"Faxa_swvdf")
     call seq_flds_add(x2i_fluxes,"Faxa_swvdf")
     call seq_flds_add(x2l_fluxes,"Faxa_swvdf")
     longname = 'Diffuse visible incident solar radiation'
     stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_visible_radiation'
     units    = 'W m-2'
     attname  = 'Faxa_swvdf'
     call metadata_set(attname, longname, stdname, units)

     ! Net shortwave radiation
     call seq_flds_add(a2x_fluxes,"Faxa_swnet") ! diagnostic 
     call seq_flds_add(l2x_fluxes,"Fall_swnet") ! diagnostic 
     call seq_flds_add(i2x_fluxes,"Faii_swnet") ! diagnostic

     call seq_flds_add(i2x_fluxes,"Fioi_swpen") ! used for Foxx_swnet below
     call seq_flds_add(x2o_fluxes,"Foxx_swnet") ! derived using albedos, Faxa_swxxx and swpen
     units    = 'W m-2'
     longname = 'Net shortwave radiation'
     stdname  = 'surface_net_shortwave_flux'
     attname  = 'Faxa_swnet' 
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Fall_swnet' 
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_swnet' 
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_swnet' 
     call metadata_set(attname, longname, stdname, units)
     longname = 'Net shortwave radiation penetrating into ice and ocean'
     stdname  = 'net_downward_shortwave_flux_in_sea_ice_due_to_penetration'
     attname  = 'Fioi_swpen' 
     call metadata_set(attname, longname, stdname, units)

     ! Black Carbon hydrophilic dry deposition  
     call seq_flds_add(a2x_fluxes,"Faxa_bcphidry" )
     call seq_flds_add(x2i_fluxes,"Faxa_bcphidry" )
     call seq_flds_add(x2l_fluxes,"Faxa_bcphidry" )
     call seq_flds_add(x2o_fluxes,"Faxa_bcphidry"   )
     longname = 'Hydrophylic black carbon dry deposition flux'
     stdname  = 'dry_deposition_flux_of_hydrophylic_black_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_bcphidry'
     call metadata_set(attname, longname, stdname, units)

     ! Black Carbon hydrophobic dry deposition   
     call seq_flds_add(a2x_fluxes,"Faxa_bcphodry" )
     call seq_flds_add(x2i_fluxes,"Faxa_bcphodry" )
     call seq_flds_add(x2l_fluxes,"Faxa_bcphodry" )
     call seq_flds_add(x2o_fluxes,"Faxa_bcphodry")
     longname = 'Hydrophobic black carbon dry deposition flux'
     stdname  = 'dry_deposition_flux_of_hydrophobic_black_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_bcphodry'
     call metadata_set(attname, longname, stdname, units)

     ! Black Carbon hydrophilic wet deposition   
     call seq_flds_add(a2x_fluxes,"Faxa_bcphiwet" )
     call seq_flds_add(x2i_fluxes,"Faxa_bcphiwet" )
     call seq_flds_add(x2l_fluxes,"Faxa_bcphiwet" )
     call seq_flds_add(x2o_fluxes,"Faxa_bcphiwet" )
     longname = 'Hydrophylic black carbon wet deposition flux'
     stdname  = 'wet_deposition_flux_of_hydrophylic_black_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_bcphiwet'
     call metadata_set(attname, longname, stdname, units)

     ! Organic Carbon hydrophilic dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_ocphidry" )
     call seq_flds_add(x2i_fluxes,"Faxa_ocphidry" )
     call seq_flds_add(x2l_fluxes,"Faxa_ocphidry" )
     call seq_flds_add(x2o_fluxes,"Faxa_ocphidry" )
     longname = 'Hydrophylic organic carbon dry deposition flux'
     stdname  = 'dry_deposition_flux_of_hydrophylic_organic_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_ocphidry'
     call metadata_set(attname, longname, stdname, units)

     ! Organic Carbon hydrophobic dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_ocphodry" )
     call seq_flds_add(x2i_fluxes,"Faxa_ocphodry" )
     call seq_flds_add(x2l_fluxes,"Faxa_ocphodry" )
     call seq_flds_add(x2o_fluxes,"Faxa_ocphodry" )
     longname = 'Hydrophobic organic carbon dry deposition flux'
     stdname  = 'dry_deposition_flux_of_hydrophobic_organic_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_ocphodry'
     call metadata_set(attname, longname, stdname, units)

     ! Organic Carbon hydrophilic wet deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_ocphiwet" )
     call seq_flds_add(x2i_fluxes,"Faxa_ocphiwet" )
     call seq_flds_add(x2l_fluxes,"Faxa_ocphiwet" )
     call seq_flds_add(x2o_fluxes,"Faxa_ocphiwet" )
     longname = 'Hydrophylic organic carbon wet deposition flux'
     stdname  = 'wet_deposition_flux_of_hydrophylic_organic_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_ocphiwet'
     call metadata_set(attname, longname, stdname, units)

     ! Size 1 dust -- wet deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstwet1"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstwet1"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstwet1"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstwet1"  )
     longname = 'Dust wet deposition flux (size 1)'
     stdname  = 'wet_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstwet1'
     call metadata_set(attname, longname, stdname, units)

     ! Size 2 dust -- wet deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstwet2"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstwet2"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstwet2"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstwet2"  )
     longname = 'Dust wet deposition flux (size 2)'
     stdname  = 'wet_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstwet2'
     call metadata_set(attname, longname, stdname, units)

     ! Size 3 dust -- wet deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstwet3"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstwet3"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstwet3"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstwet3"  )
     longname = 'Dust wet deposition flux (size 3)'
     stdname  = 'wet_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstwet3'
     call metadata_set(attname, longname, stdname, units)

     ! Size 4 dust -- wet deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstwet4"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstwet4"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstwet4"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstwet4"  )
     longname = 'Dust wet deposition flux (size 4)'
     stdname  = 'wet_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstwet4'
     call metadata_set(attname, longname, stdname, units)

     ! Size 1 dust -- dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstdry1"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstdry1"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstdry1"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstdry1"  )
     longname = 'Dust dry deposition flux (size 1)'
     stdname  = 'dry_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstdry1'
     call metadata_set(attname, longname, stdname, units)

     ! Size 2 dust -- dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstdry2"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstdry2"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstdry2"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstdry2"  )
     longname = 'Dust dry deposition flux (size 2)'
     stdname  = 'dry_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstdry2'
     call metadata_set(attname, longname, stdname, units)

     ! Size 3 dust -- dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstdry3"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstdry3"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstdry3"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstdry3"  )
     longname = 'Dust dry deposition flux (size 3)'
     stdname  = 'dry_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstdry3'
     call metadata_set(attname, longname, stdname, units)

     ! Size 4 dust -- dry deposition 
     call seq_flds_add(a2x_fluxes,"Faxa_dstdry4"  )
     call seq_flds_add(x2i_fluxes,"Faxa_dstdry4"  )
     call seq_flds_add(x2l_fluxes,"Faxa_dstdry4"  )
     call seq_flds_add(x2o_fluxes,"Faxa_dstdry4"  )
     longname = 'Dust dry deposition flux (size 4)'
     stdname  = 'dry_deposition_flux_of_dust'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_dstdry4'
     call metadata_set(attname, longname, stdname, units)

     !----------------------------------------------------------
     ! states/fluxes to atm (and ocean)
     !----------------------------------------------------------

     ! land/sea-ice/ocean fractions
     call seq_flds_add(x2a_states,'Sf_lfrac')
     call seq_flds_add(x2a_states,'Sf_ifrac')
     call seq_flds_add(x2a_states,'Sf_ofrac')
     longname = 'Surface land fraction'
     stdname  = 'land_area_fraction'
     units    = 'unitless'
     attname  = 'Sf_lfrac'
     call metadata_set(attname, longname, stdname, units)
     longname = 'Surface ice fraction'
     stdname  = 'sea_ice_area_fraction'
     attname  = 'Sf_ifrac'
     call metadata_set(attname, longname, stdname, units)
     longname = 'Surface ocean fraction'
     stdname  = 'sea_area_fraction'
     attname  = 'Sf_ofrac'
     call metadata_set(attname, longname, stdname, units)

     ! Direct albedo (visible radiation)
     call seq_flds_add(i2x_states,"Si_avsdr")
     call seq_flds_add(l2x_states,"Sl_avsdr")
     call seq_flds_add(xao_albedo,"So_avsdr") 
     call seq_flds_add(x2a_states,"Sx_avsdr")
     longname = 'Direct albedo (visible radiation)'
     stdname  = 'surface_direct_albedo_due_to_visible_radiation'
     units    = 'unitless'
     attname  = 'Si_avsdr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_avsdr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_avsdr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_avsdr'
     call metadata_set(attname, longname, stdname, units)

     ! Direct albedo (near-infrared radiation)
     call seq_flds_add(i2x_states,"Si_anidr")
     call seq_flds_add(l2x_states,"Sl_anidr")
     call seq_flds_add(xao_albedo,"So_anidr") 
     call seq_flds_add(x2a_states,"Sx_anidr")
     longname = 'Direct albedo (near-infrared radiation)'
     stdname  = 'surface_direct_albedo_due_to_near_infrared_radiation'
     units    = 'unitless'
     attname  = 'Si_anidr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_anidr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_anidr'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_anidr'
     call metadata_set(attname, longname, stdname, units)

     ! Diffuse albedo (visible radiation)
     call seq_flds_add(i2x_states,"Si_avsdf")
     call seq_flds_add(l2x_states,"Sl_avsdf")
     call seq_flds_add(xao_albedo,"So_avsdf") 
     call seq_flds_add(x2a_states,"Sx_avsdf")
     longname = 'Diffuse albedo (visible radiation)'
     stdname  = 'surface_diffuse_albedo_due_to_visible_radiation'
     units    = 'unitless'
     attname  = 'Si_avsdf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_avsdf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_avsdf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_avsdf'
     call metadata_set(attname, longname, stdname, units)

     ! Diffuse albedo (near-infrared radiation)
     call seq_flds_add(i2x_states,"Si_anidf")
     call seq_flds_add(l2x_states,"Sl_anidf")
     call seq_flds_add(xao_albedo,"So_anidf") 
     call seq_flds_add(x2a_states,"Sx_anidf")
     longname = 'Diffuse albedo (near-infrared radiation)'
     stdname  = 'surface_diffuse_albedo_due_to_near_infrared_radiation'
     units    = 'unitless'
     attname  = 'Si_anidf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_anidf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_anidf'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_anidf'
     call metadata_set(attname, longname, stdname, units)

     ! Reference temperature at 2 meters
     call seq_flds_add(l2x_states,"Sl_tref")
     call seq_flds_add(i2x_states,"Si_tref")
     call seq_flds_add(xao_states,"So_tref")
     call seq_flds_add(x2a_states,"Sx_tref")
     longname = 'Reference temperature at 2 meters'
     stdname  = 'air_temperature'
     units    = 'K'
     attname  = 'Si_tref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_tref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_tref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_tref'
     call metadata_set(attname, longname, stdname, units)

     ! Reference specific humidity at 2 meters
     call seq_flds_add(l2x_states,"Sl_qref")
     call seq_flds_add(i2x_states,"Si_qref")
     call seq_flds_add(xao_states,"So_qref")
     call seq_flds_add(x2a_states,"Sx_qref")
     longname = 'Reference specific humidity at 2 meters'
     stdname  = 'specific_humidity'
     units    = 'kg kg-1'
     attname  = 'Si_qref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_qref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_qref'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_qref'
     call metadata_set(attname, longname, stdname, units)

     ! Surface temperature
     call seq_flds_add(l2x_states,"Sl_t")
     call seq_flds_add(i2x_states,"Si_t")
     call seq_flds_add(x2a_states,"So_t")
     call seq_flds_add(x2a_states,"Sx_t")
     longname = 'Surface temperature'
     stdname  = 'surface_temperature'
     units    = 'K'
     attname  = 'Si_t'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sl_t'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'So_t'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Sx_t'
     call metadata_set(attname, longname, stdname, units)

     ! Surface friction velocity in land (land/atm only)
     call seq_flds_add(l2x_states,"Sl_fv")
     call seq_flds_add(x2a_states,"Sl_fv")
     longname = 'Surface fraction velocity in land'
     stdname  = 'fraction_velocity'
     units    = 'm s-1'
     attname  = 'Sl_fv'
     call metadata_set(attname, longname, stdname, units)

     ! Aerodynamical resistance (land/atm only)
     call seq_flds_add(l2x_states,"Sl_ram1")
     call seq_flds_add(x2a_states,"Sl_ram1")

     ! Surface snow water equivalent (land/atm only) 
     call seq_flds_add(l2x_states,"Sl_snowh")
     call seq_flds_add(x2a_states,"Sl_snowh")
     longname = 'Surface snow water equivalent'
     stdname  = 'surface_snow_water_equivalent'
     units    = 'm'
     attname  = 'Sl_snowh'
     call metadata_set(attname, longname, stdname, units)

     ! Surface snow depth (ice/atm only)
     call seq_flds_add(i2x_states,"Si_snowh")
     call seq_flds_add(x2a_states,"Si_snowh")
     longname = 'Surface snow depth'
     stdname  = 'surface_snow_thickness'
     units    = 'm'
     attname  = 'Si_snowh'
     call metadata_set(attname, longname, stdname, units)

     ! Surface saturation specific humidity in ocean (ocn/atm only)
     call seq_flds_add(xao_states,"So_ssq")
     call seq_flds_add(x2a_states,"So_ssq")
     longname = 'Surface saturation specific humidity in ocean'
     stdname  = 'specific_humidity_at_saturation'
     units    = 'kg kg-1'
     attname  = 'So_ssq'
     call metadata_set(attname, longname, stdname, units)

     ! Square of exch. coeff (tracers) (ocn/atm only)
     call seq_flds_add(xao_states,"So_re")
     call seq_flds_add(x2a_states,"So_re")
     longname = 'Square of exch. coeff (tracers)'
     stdname  = ''
     units    = ''
     attname  = 'So_re'
     call metadata_set(attname, longname, stdname, units)

     ! 10 meter wind
     call seq_flds_add(i2x_states,"Si_u10")
     call seq_flds_add(xao_states,"So_u10")
     call seq_flds_add(l2x_states,"Sl_u10")
     call seq_flds_add(x2a_states,"Sx_u10")
     longname = '10m wind'
     stdname  = '10m_wind'
     units    = 'm'
     attname  = 'u10'
     call metadata_set(attname, longname, stdname, units)

     ! Zonal surface stress"
     call seq_flds_add(l2x_fluxes,"Fall_taux")
     call seq_flds_add(xao_fluxes,"Faox_taux")
     call seq_flds_add(i2x_fluxes,"Faii_taux")
     call seq_flds_add(x2a_fluxes,"Faxx_taux")
     call seq_flds_add(i2x_fluxes,"Fioi_taux")
     call seq_flds_add(x2o_fluxes,"Foxx_taux")
     longname = 'Zonal surface stress'
     stdname  = 'surface_downward_eastward_stress'
     units    = 'N m-2'
     attname  = 'Fall_taux'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_taux'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_taux'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Fioi_taux'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_taux'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_taux'
     call metadata_set(attname, longname, stdname, units)

     ! Meridional surface stress
     call seq_flds_add(l2x_fluxes,"Fall_tauy")
     call seq_flds_add(xao_fluxes,"Faox_tauy")
     call seq_flds_add(i2x_fluxes,"Faii_tauy")
     call seq_flds_add(x2a_fluxes,"Faxx_tauy")
     call seq_flds_add(i2x_fluxes,"Fioi_tauy")
     call seq_flds_add(x2o_fluxes,"Foxx_tauy")
     longname = 'Meridional surface stress'
     stdname  = 'surface_downward_northward_stress'
     units    = 'N m-2'
     attname  = 'Fall_tauy'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_tauy'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_tauy'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Fioi_tauy'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_tauy'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_tauy'
     call metadata_set(attname, longname, stdname, units)

     ! Surface latent heat flux
     call seq_flds_add(l2x_fluxes,"Fall_lat")
     call seq_flds_add(xao_fluxes,"Faox_lat")
     call seq_flds_add(i2x_fluxes,"Faii_lat")
     call seq_flds_add(x2a_fluxes,"Faxx_lat")
     call seq_flds_add(x2o_fluxes,"Foxx_lat")
     longname = 'Surface latent heat flux'
     stdname  = 'surface_upward_latent_heat_flux'
     units    = 'W m-2'
     attname  = 'Fall_lat'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_lat'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_lat'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_lat'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_lat'
     call metadata_set(attname, longname, stdname, units)

     ! Surface sensible heat flux
     call seq_flds_add(l2x_fluxes,"Fall_sen")
     call seq_flds_add(xao_fluxes,"Faox_sen")
     call seq_flds_add(i2x_fluxes,"Faii_sen")
     call seq_flds_add(x2a_fluxes,"Faxx_sen")
     call seq_flds_add(x2o_fluxes,"Foxx_sen")
     longname = 'Sensible heat flux'
     stdname  = 'surface_upward_sensible_heat_flux'
     units    = 'W m-2'
     attname  = 'Fall_sen'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_sen'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_sen'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_sen'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_sen'
     call metadata_set(attname, longname, stdname, units)

     ! Surface upward longwave heat flux
     call seq_flds_add(l2x_fluxes,"Fall_lwup")
     call seq_flds_add(xao_fluxes,"Faox_lwup")
     call seq_flds_add(i2x_fluxes,"Faii_lwup")
     call seq_flds_add(x2a_fluxes,"Faxx_lwup")
     call seq_flds_add(x2o_fluxes,"Foxx_lwup")
     longname = 'Surface upward longwave heat flux'
     stdname  = 'surface_net_upward_longwave_flux'
     units    = 'W m-2'
     attname  = 'Fall_lwup'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_lwup'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_lwup'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_lwup'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_lwup'
     call metadata_set(attname, longname, stdname, units)

     ! Evaporation water flux
     call seq_flds_add(l2x_fluxes,"Fall_evap")
     call seq_flds_add(xao_fluxes,"Faox_evap")
     call seq_flds_add(i2x_fluxes,"Faii_evap")
     call seq_flds_add(x2a_fluxes,"Faxx_evap")
     call seq_flds_add(x2o_fluxes,"Foxx_evap")
     longname = 'Evaporation water flux'
     stdname  = 'water_evaporation_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fall_evap'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faox_evap'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faii_evap'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Faxx_evap'
     call metadata_set(attname, longname, stdname, units)

     ! Dust flux (particle bin number 1)
     call seq_flds_add(l2x_fluxes,"Fall_flxdst1")
     call seq_flds_add(x2a_fluxes,"Fall_flxdst1")
     longname = 'Dust flux (particle bin number 1)'
     stdname  = 'dust_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fall_flxdst1'
     call metadata_set(attname, longname, stdname, units)

     ! Dust flux (particle bin number 2)
     call seq_flds_add(l2x_fluxes,"Fall_flxdst2")
     call seq_flds_add(x2a_fluxes,"Fall_flxdst2")
     longname = 'Dust flux (particle bin number 2)'
     stdname  = 'dust_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fall_flxdst2'
     call metadata_set(attname, longname, stdname, units)

     ! Dust flux (particle bin number 3)
     call seq_flds_add(l2x_fluxes,"Fall_flxdst3")
     call seq_flds_add(x2a_fluxes,"Fall_flxdst3")
     longname = 'Dust flux (particle bin number 3)'
     stdname  = 'dust_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fall_flxdst3'
     call metadata_set(attname, longname, stdname, units)

     ! Dust flux (particle bin number 4)
     call seq_flds_add(l2x_fluxes,"Fall_flxdst4")
     call seq_flds_add(x2a_fluxes,"Fall_flxdst4")
     longname = 'Dust flux (particle bin number 4)'
     stdname  = 'dust_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fall_flxdst4'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! atm<->ocn only exchange
     !-----------------------------

     ! Sea level pressure (Pa)
     call seq_flds_add(a2x_states,"Sa_pslv")
     call seq_flds_add(x2o_states,"Sa_pslv")
     longname = 'Sea level pressure'
     stdname  = 'air_pressure_at_sea_level'
     units    = 'Pa'
     attname  = 'Sa_pslv'
     call metadata_set(attname, longname, stdname, units)

     ! Wind speed squared at 10 meters
     call seq_flds_add(xao_states,"So_duu10n")
     call seq_flds_add(x2o_states,"So_duu10n")
     longname = 'Wind speed squared at 10 meters'
     stdname  = 'square_of_wind_speed'
     units    = 'm2 s-2'
     attname  = 'So_duu10n'
     call metadata_set(attname, longname, stdname, units)

     ! Surface friction velocity in ocean
     call seq_flds_add(xao_states,"So_ustar")
     call seq_flds_add(x2a_states,"So_ustar")
     longname = 'Surface fraction velocity in ocean'
     stdname  = 'fraction_velocity'
     units    = 'm s-1'
     attname  = 'So_ustar'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! ice<->ocn only exchange
     !-----------------------------

     ! Fractional ice coverage wrt ocean
     call seq_flds_add(i2x_states,"Si_ifrac")
     call seq_flds_add(x2o_states,"Si_ifrac")
     call seq_flds_add(x2w_states,"Si_ifrac")
     longname = 'Fractional ice coverage wrt ocean'
     stdname  = 'sea_ice_area_fraction'
     units    = 'unitless'
     attname  = 'Si_ifrac'
     call metadata_set(attname, longname, stdname, units)

     ! Ocean freeze (q>0) or melt (q<0) potential
     call seq_flds_add(o2x_fluxes,"Fioo_q")
     call seq_flds_add(x2i_fluxes,"Fioo_q")
     longname = 'Ocean freeze (q>0) or melt (q<0) potential'
     stdname  = 'surface_snow_and_ice_melt_heat_flux'
     units    = 'W m-2'
     attname  = 'Fioo_q'
     call metadata_set(attname, longname, stdname, units)

     ! Heat flux from melting
     call seq_flds_add(i2x_fluxes,"Fioi_melth")
     call seq_flds_add(x2o_fluxes,"Fioi_melth")
     longname = 'Heat flux from melting'
     stdname  = 'surface_snow_melt_heat_flux'
     units    = 'W m-2'
     attname  = 'Faxa_melth'
     call metadata_set(attname, longname, stdname, units)

     ! Water flux from melting
     call seq_flds_add(i2x_fluxes,"Fioi_meltw")
     call seq_flds_add(x2o_fluxes,"Fioi_meltw")
     longname = 'Water flux due to melting'
     stdname  = 'surface_snow_melt_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_meltw'
     call metadata_set(attname, longname, stdname, units)

     ! Salt flux
     call seq_flds_add(i2x_fluxes,"Fioi_salt")
     call seq_flds_add(x2o_fluxes,"Fioi_salt")
     longname = 'Salt flux'
     stdname  = 'virtual_salt_flux_into_sea_water'
     units    = 'kg m-2 s-1'
     attname  = 'Faxa_salt'
     call metadata_set(attname, longname, stdname, units)

     ! Sea surface temperature
     call seq_flds_add(o2x_states,"So_t")
     call seq_flds_add(x2i_states,"So_t")    
     call seq_flds_add(x2w_states,"So_t")    

     ! Sea surface  salinity
     call seq_flds_add(o2x_states,"So_s")
     call seq_flds_add(x2i_states,"So_s")    
     longname = 'Sea surface salinity'
     stdname  = 'sea_surface_salinity'
     units    = 'g kg-1'
     attname  = 'So_s'
     call metadata_set(attname, longname, stdname, units)

     ! Zonal sea water velocity
     call seq_flds_add(o2x_states,"So_u")
     call seq_flds_add(x2i_states,"So_u")    
     call seq_flds_add(x2w_states,"So_u")    
     longname = 'Zonal sea water velocity'
     stdname  = 'eastward_sea_water_velocity'
     units    = 'm s-1'
     attname  = 'So_u'
     call metadata_set(attname, longname, stdname, units)

     ! Meridional sea water velocity
     call seq_flds_add(o2x_states,"So_v")
     call seq_flds_add(x2i_states,"So_v")    
     call seq_flds_add(x2w_states,"So_v")    
     longname = 'Meridional sea water velocity'
     stdname  = 'northward_sea_water_velocity'
     units    = 'm s-1'
     attname  = 'So_v'

     ! Zonal sea surface slope
     call seq_flds_add(o2x_states,"So_dhdx")
     call seq_flds_add(x2i_states,"So_dhdx") 
     longname = 'Zonal sea surface slope'
     stdname  = 'sea_surface_eastward_slope'
     units    = 'm m-1'
     attname  = 'So_dhdx'
     call metadata_set(attname, longname, stdname, units)

     ! Meridional sea surface slope
     call seq_flds_add(o2x_states,"So_dhdy")
     call seq_flds_add(x2i_states,"So_dhdy") 
     longname = 'Meridional sea surface slope'
     stdname  = 'sea_surface_northward_slope'
     units    = 'm m-1'
     attname  = 'So_dhdy'
     call metadata_set(attname, longname, stdname, units)

     ! Boundary Layer Depth
     call seq_flds_add(o2x_states,"So_bldepth")
     call seq_flds_add(x2w_states,"So_bldepth") 
     longname = 'Ocean Boundary Layer Depth'
     stdname  = 'ocean_boundary_layer_depth'
     units    = 'm'
     attname  = 'So_bldepth'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! lnd->rof exchange
     ! TODO: put in attributes below
     !-----------------------------

     call seq_flds_add(l2x_fluxes,'Flrl_rofl')
     call seq_flds_add(x2r_fluxes,'Flrl_rofl')
     longname = 'Water flux from land (liquid)'
     stdname  = 'water_flux_into_runoff'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofl'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(l2x_fluxes,'Flrl_rofi')
     call seq_flds_add(x2r_fluxes,'Flrl_rofi')
     longname = 'Water flux from land (frozen)'
     stdname  = 'frozen_water_flux_into_runoff'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofi'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! rof->ocn (runoff) and rof->lnd (flooding)
     !-----------------------------

     call seq_flds_add(r2x_fluxes,'Forr_rofl')
     call seq_flds_add(x2o_fluxes,'Foxx_rofl')
     longname = 'Water flux due to runoff (liquid)'
     stdname  = 'water_flux_into_sea_water'
     units    = 'kg m-2 s-1'
     attname  = 'Forr_rofl'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_rofl'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(r2x_fluxes,'Forr_rofi')
     call seq_flds_add(x2o_fluxes,'Foxx_rofi')
     longname = 'Water flux due to runoff (frozen)'
     stdname  = 'frozen_water_flux_into_sea_water'
     units    = 'kg m-2 s-1'
     attname  = 'Forr_rofi'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Foxx_rofi'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(r2x_fluxes,'Firr_rofi')
     call seq_flds_add(x2i_fluxes,'Fixx_rofi')
     longname = 'Water flux due to runoff (frozen)'
     stdname  = 'frozen_water_flux_into_sea_ice'
     units    = 'kg m-2 s-1'
     attname  = 'Firr_rofi'
     call metadata_set(attname, longname, stdname, units)
     attname  = 'Fixx_rofi'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(r2x_fluxes,'Flrr_flood')
     call seq_flds_add(x2l_fluxes,'Flrr_flood')
     longname = 'Waterrflux due to flooding'
     stdname  = 'flooding_water_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Flrr_flood'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(r2x_fluxes,'Flrr_volr')
     call seq_flds_add(x2l_fluxes,'Flrr_volr')
     longname = 'River channel water volume'
     stdname  = 'rtm_volr'
     units    = 'm'
     attname  = 'Flrr_volr'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! wav->ocn and ocn->wav
     !-----------------------------

     call seq_flds_add(w2x_states,'Sw_lamult')
     call seq_flds_add(x2o_states,'Sw_lamult')
     longname = 'Langmuir multiplier'
     stdname  = 'wave_model_langmuir_multiplier'
     units    = ''
     attname  = 'Sw_lamult'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(w2x_states,'Sw_ustokes')
     call seq_flds_add(x2o_states,'Sw_ustokes')
     longname = 'Stokes drift u component'
     stdname  = 'wave_model_stokes_drift_eastward_velocity'
     units    = 'm/s'
     attname  = 'Sw_ustokes'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(w2x_states,'Sw_vstokes')
     call seq_flds_add(x2o_states,'Sw_vstokes')
     longname = 'Stokes drift v component'
     stdname  = 'wave_model_stokes_drift_northward_velocity'
     units    = 'm/s'
     attname  = 'Sw_vstokes'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(w2x_states,'Sw_hstokes')
     call seq_flds_add(x2o_states,'Sw_hstokes')
     longname = 'Stokes drift depth'
     stdname  = 'wave_model_stokes_drift_depth'
     units    = 'm'
     attname  = 'Sw_hstokes'
     call metadata_set(attname, longname, stdname, units)

     !-----------------------------
     ! glc fields
     !-----------------------------

     name = 'Fogg_rofl'
     call seq_flds_add(g2x_fluxes,trim(name))
     longname = 'glc liquid runoff flux to ocean'
     stdname  = 'glacier_liquid_runoff_flux_to_ocean'
     units    = 'kg m-2 s-1'
     attname  = 'Fogg_rofl'
     call metadata_set(attname, longname, stdname, units)

     name = 'Fogg_rofi'
     call seq_flds_add(g2x_fluxes,trim(name))
     longname = 'glc frozen runoff flux to ocean'
     stdname  = 'glacier_frozen_runoff_flux_to_ocean'
     units    = 'kg m-2 s-1'
     attname  = 'Fogg_rofi'
     call metadata_set(attname, longname, stdname, units)

     name = 'Figg_rofi'
     call seq_flds_add(g2x_fluxes,trim(name))
     longname = 'glc frozen runoff_iceberg flux to ice'
     stdname  = 'glacier_frozen_runoff_flux_to_seaice'
     units    = 'kg m-2 s-1'
     attname  = 'Figg_rofi'
     call metadata_set(attname, longname, stdname, units)

     name = 'Sg_icemask'
     call seq_flds_add(g2x_states,trim(name))     
     call seq_flds_add(x2l_states,trim(name))
     longname = 'Ice sheet grid coverage on global grid'
     stdname  = 'ice_sheet_grid_mask'
     units    = 'unitless'
     attname  = 'Sg_icemask'
     call metadata_set(attname, longname, stdname, units)     

     name = 'Sg_icemask_coupled_fluxes'
     call seq_flds_add(g2x_states,trim(name))     
     call seq_flds_add(x2l_states,trim(name))
     longname = 'Ice sheet mask where we are potentially sending non-zero fluxes'
     stdname  = 'icemask_coupled_fluxes'
     units    = 'unitless'
     attname  = 'Sg_icemask_coupled_fluxes'
     call metadata_set(attname, longname, stdname, units)     

     ! If glc_nec > 0, then create coupling fields for all glc elevation classes
     ! (1:glc_nec) plus bare land (index 0). Note that, if glc_nec = 0, then we don't
     ! even need the bare land (0) index.
     if (seq_flds_glc_nec > 0) then
        do num = 0,seq_flds_glc_nec
           write(cnum,'(i2.2)') num

           ! glc fields: lnd->glc 

           name = 'Sl_tsrf' // cnum
           call seq_flds_add(l2x_states,trim(name))
           call seq_flds_add(x2g_states,trim(name))
           longname = 'Surface temperature  of glacier elevation class ' // cnum 
           stdname  = 'surface_temperature'
           units    = 'deg C'
           attname  = 'Sl_tsrf' // cnum
           call metadata_set(attname, longname, stdname, units)

           name = 'Sl_topo' // cnum
           call seq_flds_add(l2x_states,trim(name))
           call seq_flds_add(x2g_states,trim(name))
           longname = 'Surface height of glacier elevation class ' // cnum 
           stdname  = 'height'
           units    = 'm'
           attname  = 'Sl_topo' // cnum
           call metadata_set(attname, longname, stdname, units)

           name = 'Flgl_qice' // cnum
           call seq_flds_add(l2x_fluxes,trim(name))
           call seq_flds_add(x2g_fluxes,trim(name))
           longname = 'New glacier ice flux of elevation class ' // cnum
           stdname  = 'ice_flux_out_of_glacier'
           units    = 'kg m-2 s-1'
           attname  = 'Fgll_qice' // cnum
           call metadata_set(attname, longname, stdname, units)

           ! glc fields: glc->lnd 

           name = 'Sg_frac' // cnum
           call seq_flds_add(g2x_states,trim(name))
           call seq_flds_add(x2l_states,trim(name))
           longname = 'Fraction of glacier area of elevation class ' // cnum
           stdname  = 'glacier_area_fraction'
           units    = 'unitless'    
           attname  = 'Sg_frac' // cnum
           call metadata_set(attname, longname, stdname, units)

           name = 'Sg_topo' // cnum
           call seq_flds_add(g2x_states,trim(name))
           call seq_flds_add(x2l_states,trim(name))
           longname = 'Surface height of glacier of elevation class ' // cnum
           stdname  = 'height'
           units    = 'm'
           attname  = 'Sg_topo' // cnum
           call metadata_set(attname, longname, stdname, units)

           name = 'Flgg_hflx' // cnum
           call seq_flds_add(g2x_fluxes,trim(name))
           call seq_flds_add(x2l_fluxes,trim(name))
           longname = 'Downward heat flux from glacier interior of elevation class ' // cnum
           stdname  = 'downward_heat_flux_in_glacier'
           units    = 'W m-2'    
           attname  = 'Flgg_hflx' // cnum
           call metadata_set(attname, longname, stdname, units)
        end do
     end if

     if (flds_co2a) then

        call seq_flds_add(a2x_states, "Sa_co2prog")
        call seq_flds_add(x2l_states, "Sa_co2prog")
        longname = 'Prognostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2prog'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(a2x_states, "Sa_co2diag")
        call seq_flds_add(x2l_states, "Sa_co2diag")
        longname = 'Diagnostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2diag'
        call metadata_set(attname, longname, stdname, units)

     else if (flds_co2b) then

        call seq_flds_add(a2x_states,  "Sa_co2prog")
        call seq_flds_add(x2l_states,  "Sa_co2prog")
        longname = 'Prognostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2prog'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(a2x_states,  "Sa_co2diag")
        call seq_flds_add(x2l_states,  "Sa_co2diag")
        longname = 'Diagnostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2diag'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(l2x_fluxes,  "Fall_fco2_lnd")
        call seq_flds_add(x2a_fluxes,  "Fall_fco2_lnd")
        longname = 'Surface flux of CO2 from land'
        stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
        units    = 'moles m-2 s-1'
        attname  = 'Fall_fco2_lnd'
        call metadata_set(attname, longname, stdname, units)

     else if (flds_co2c) then

        call seq_flds_add(a2x_states, "Sa_co2prog")
        call seq_flds_add(x2l_states, "Sa_co2prog")
        call seq_flds_add(x2o_states, "Sa_co2prog")
        longname = 'Prognostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2prog'
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(a2x_states, "Sa_co2diag")
        call seq_flds_add(x2l_states, "Sa_co2diag")
        call seq_flds_add(x2o_states, "Sa_co2diag")
        longname = 'Diagnostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2diag'
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(l2x_fluxes, "Fall_fco2_lnd")
        call seq_flds_add(x2a_fluxes, "Fall_fco2_lnd")
        longname = 'Surface flux of CO2 from land'
        stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
        units    = 'moles m-2 s-1'
        attname  = 'Fall_foc2_lnd' 
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(o2x_fluxes, "Faoo_fco2_ocn")
        call seq_flds_add(x2a_fluxes, "Faoo_fco2_ocn")
        longname = 'Surface flux of CO2 from ocean'
        stdname  = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
        units    = 'moles m-2 s-1'
        attname  = 'Faoo_fco2_ocn' 
        call metadata_set(attname, longname, stdname, units)
        
     else if (flds_co2_dmsa) then

        call seq_flds_add(a2x_states, "Sa_co2prog")
        call seq_flds_add(x2l_states, "Sa_co2prog")
        longname = 'Prognostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2prog'
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(a2x_states, "Sa_co2diag")
        call seq_flds_add(x2l_states, "Sa_co2diag")
        longname = 'Diagnostic CO2 at the lowest model level'
        stdname  = ''
        units    = '1e-6 mol/mol'
        attname  = 'Sa_co2diag'
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(o2x_fluxes, "Faoo_fdms_ocn")
        call seq_flds_add(x2a_fluxes, "Faoo_fdms_ocn")
        longname = 'Surface flux of DMS'
        stdname  = 'surface_upward_flux_of_dimethyl_sulfide'
        units    = 'moles m-2 s-1'
        attname  = 'Faoo_fdms'     
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(l2x_fluxes, "Fall_fco2_lnd")
        call seq_flds_add(x2a_fluxes, "Fall_fco2_lnd")
        longname = 'Surface flux of CO2 from land'
        stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
        units    = 'moles m-2 s-1'
        attname  = 'Fall_foc2_lnd' 
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(o2x_fluxes, "Faoo_fco2_ocn")
        call seq_flds_add(x2a_fluxes, "Faoo_fco2_ocn")
        longname = 'Surface flux of CO2 from ocean'
        stdname  = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
        units    = 'moles m-2 s-1'
        attname  = 'Faoo_fco2_ocn' 
        call metadata_set(attname, longname, stdname, units)
        
     endif

     !-----------------------------------------------------------------------------
     ! Read namelist for CARMA 
     ! if carma_flds are specified then setup fields for CLM to CAM communication
     !-----------------------------------------------------------------------------

     call shr_carma_readnl(nlfilename='drv_flds_in', carma_fields=carma_fields)
     if (carma_fields /= ' ') then
        call seq_flds_add(l2x_fluxes, trim(carma_fields))
        call seq_flds_add(x2a_fluxes, trim(carma_fields))
     endif

     !-----------------------------------------------------------------------------
     ! Read namelist for MEGAN
     ! if MEGAN emission are specified then setup fields for CLM to CAM communication 
     ! (emissions fluxes)
     !-----------------------------------------------------------------------------

     call shr_megan_readnl(nlfilename='drv_flds_in', megan_fields=megan_voc_fields)
     if (shr_megan_mechcomps_n>0) then
        call seq_flds_add(l2x_fluxes, trim(megan_voc_fields))
        call seq_flds_add(x2a_fluxes, trim(megan_voc_fields))
     endif

     !-----------------------------------------------------------------------------
     ! Dry Deposition fields
     ! First read namelist and figure out the drydep field list to pass
     ! Then check if file exists and if not, n_drydep will be zero
     ! Then add dry deposition fields to land export and atmosphere import states
     ! Then initialize dry deposition fields
     ! Note: CAM and CLM will then call seq_drydep_setHCoeff
     !-----------------------------------------------------------------------------

     call seq_drydep_read(nlfilename="drv_flds_in", seq_drydep_fields=seq_drydep_fields)
     if ( lnd_drydep ) then
        call seq_flds_add(l2x_states, trim(seq_drydep_fields))
        call seq_flds_add(x2a_states, trim(seq_drydep_fields))
     endif
     call seq_drydep_init( )

     !----------------------------------------------------------------------------
     ! state + flux fields
     !----------------------------------------------------------------------------

     seq_flds_dom_coord  = trim(dom_coord )
     seq_flds_a2x_states = trim(a2x_states)
     seq_flds_x2a_states = trim(x2a_states)
     seq_flds_i2x_states = trim(i2x_states)
     seq_flds_x2i_states = trim(x2i_states)
     seq_flds_l2x_states = trim(l2x_states)
     seq_flds_x2l_states = trim(x2l_states)
     seq_flds_o2x_states = trim(o2x_states)
     seq_flds_x2o_states = trim(x2o_states)
     seq_flds_g2x_states = trim(g2x_states)
     seq_flds_x2g_states = trim(x2g_states)
     seq_flds_xao_states = trim(xao_states)
     seq_flds_xao_albedo = trim(xao_albedo)
     seq_flds_r2x_states = trim(r2x_states)
     seq_flds_x2r_states = trim(x2r_states)
     seq_flds_w2x_states = trim(w2x_states)
     seq_flds_x2w_states = trim(x2w_states)

     seq_flds_dom_other  = trim(dom_other )
     seq_flds_a2x_fluxes = trim(a2x_fluxes)
     seq_flds_x2a_fluxes = trim(x2a_fluxes)
     seq_flds_i2x_fluxes = trim(i2x_fluxes)
     seq_flds_x2i_fluxes = trim(x2i_fluxes)
     seq_flds_l2x_fluxes = trim(l2x_fluxes)
     seq_flds_x2l_fluxes = trim(x2l_fluxes)
     seq_flds_o2x_fluxes = trim(o2x_fluxes)
     seq_flds_x2o_fluxes = trim(x2o_fluxes)
     seq_flds_g2x_fluxes = trim(g2x_fluxes)
     seq_flds_x2g_fluxes = trim(x2g_fluxes)
     seq_flds_xao_fluxes = trim(xao_fluxes)
     seq_flds_r2x_fluxes = trim(r2x_fluxes)
     seq_flds_x2r_fluxes = trim(x2r_fluxes)
     seq_flds_w2x_fluxes = trim(w2x_fluxes)
     seq_flds_x2w_fluxes = trim(x2w_fluxes)

     if (seq_comm_iamroot(ID)) then
        write(logunit,"(A)") subname//': seq_flds_a2x_states= ',trim(seq_flds_a2x_states)
        write(logunit,"(A)") subname//': seq_flds_a2x_fluxes= ',trim(seq_flds_a2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2a_states= ',trim(seq_flds_x2a_states)
        write(logunit,"(A)") subname//': seq_flds_x2a_fluxes= ',trim(seq_flds_x2a_fluxes)
        write(logunit,"(A)") subname//': seq_flds_l2x_states= ',trim(seq_flds_l2x_states)
        write(logunit,"(A)") subname//': seq_flds_l2x_fluxes= ',trim(seq_flds_l2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2l_states= ',trim(seq_flds_x2l_states)
        write(logunit,"(A)") subname//': seq_flds_x2l_fluxes= ',trim(seq_flds_x2l_fluxes)
        write(logunit,"(A)") subname//': seq_flds_i2x_states= ',trim(seq_flds_i2x_states)
        write(logunit,"(A)") subname//': seq_flds_i2x_fluxes= ',trim(seq_flds_i2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2i_states= ',trim(seq_flds_x2i_states)
        write(logunit,"(A)") subname//': seq_flds_x2i_fluxes= ',trim(seq_flds_x2i_fluxes)
        write(logunit,"(A)") subname//': seq_flds_o2x_states= ',trim(seq_flds_o2x_states)
        write(logunit,"(A)") subname//': seq_flds_o2x_fluxes= ',trim(seq_flds_o2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2o_states= ',trim(seq_flds_x2o_states)
        write(logunit,"(A)") subname//': seq_flds_x2o_fluxes= ',trim(seq_flds_x2o_fluxes)
        write(logunit,"(A)") subname//': seq_flds_g2x_states= ',trim(seq_flds_g2x_states)
        write(logunit,"(A)") subname//': seq_flds_g2x_fluxes= ',trim(seq_flds_g2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2g_states= ',trim(seq_flds_x2g_states)
        write(logunit,"(A)") subname//': seq_flds_x2g_fluxes= ',trim(seq_flds_x2g_fluxes)
        write(logunit,"(A)") subname//': seq_flds_xao_states= ',trim(seq_flds_xao_states)
        write(logunit,"(A)") subname//': seq_flds_xao_fluxes= ',trim(seq_flds_xao_fluxes)
        write(logunit,"(A)") subname//': seq_flds_xao_albedo= ',trim(seq_flds_xao_albedo)
        write(logunit,"(A)") subname//': seq_flds_r2x_states= ',trim(seq_flds_r2x_states)
        write(logunit,"(A)") subname//': seq_flds_r2x_fluxes= ',trim(seq_flds_r2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2r_states= ',trim(seq_flds_x2r_states)
        write(logunit,"(A)") subname//': seq_flds_x2r_fluxes= ',trim(seq_flds_x2r_fluxes)
        write(logunit,"(A)") subname//': seq_flds_w2x_states= ',trim(seq_flds_w2x_states)
        write(logunit,"(A)") subname//': seq_flds_w2x_fluxes= ',trim(seq_flds_w2x_fluxes)
        write(logunit,"(A)") subname//': seq_flds_x2w_states= ',trim(seq_flds_x2w_states)
        write(logunit,"(A)") subname//': seq_flds_x2w_fluxes= ',trim(seq_flds_x2w_fluxes)
     end if

     call catFields(seq_flds_dom_fields, seq_flds_dom_coord , seq_flds_dom_other )
     call catFields(seq_flds_a2x_fields, seq_flds_a2x_states, seq_flds_a2x_fluxes)
     call catFields(seq_flds_x2a_fields, seq_flds_x2a_states, seq_flds_x2a_fluxes)
     call catFields(seq_flds_i2x_fields, seq_flds_i2x_states, seq_flds_i2x_fluxes)
     call catFields(seq_flds_x2i_fields, seq_flds_x2i_states, seq_flds_x2i_fluxes)
     call catFields(seq_flds_l2x_fields, seq_flds_l2x_states, seq_flds_l2x_fluxes)
     call catFields(seq_flds_x2l_fields, seq_flds_x2l_states, seq_flds_x2l_fluxes)
     call catFields(seq_flds_o2x_fields, seq_flds_o2x_states, seq_flds_o2x_fluxes)
     call catFields(seq_flds_x2o_fields, seq_flds_x2o_states, seq_flds_x2o_fluxes)
     call catFields(seq_flds_g2x_fields, seq_flds_g2x_states, seq_flds_g2x_fluxes)
     call catFields(seq_flds_x2g_fields, seq_flds_x2g_states, seq_flds_x2g_fluxes)
     call catFields(stringtmp          , seq_flds_xao_albedo, seq_flds_xao_states)
     call catFields(seq_flds_xao_fields, stringtmp          , seq_flds_xao_fluxes)
     call catFields(seq_flds_r2x_fields, seq_flds_r2x_states, seq_flds_r2x_fluxes)
     call catFields(seq_flds_x2r_fields, seq_flds_x2r_states, seq_flds_x2r_fluxes)
     call catFields(seq_flds_w2x_fields, seq_flds_w2x_states, seq_flds_w2x_fluxes)
     call catFields(seq_flds_x2w_fields, seq_flds_x2w_states, seq_flds_x2w_fluxes)

   end subroutine seq_flds_set

   !===============================================================================
   !BOP ===========================================================================
   !
   ! !IROUTINE: seq_flds_add
   !
   ! !DESCRIPTION:
   !  Returns new concatentated field list 
   !  in the output character string {\tt outfld}.
   !
   ! !REVISION HISTORY:
   !  2011-Nov-27  - M. Vertenstein - first version
   !
   ! !INTERFACE: ------------------------------------------------------------------

   subroutine seq_flds_add(outfld, str)

     ! !USES:

     ! !INPUT/OUTPUT PARAMETERS:

     character(len=*),intent(in)    :: str      ! string 
     character(len=*),intent(inout) :: outfld   ! output field name

     !EOP

     character(len=*),parameter :: subname = '(seq_flds_add) '
     !-------------------------------------------------------------------------------
     !
     !-------------------------------------------------------------------------------

     if (trim(outfld) == '') then
        outfld = trim(str)
     else
        outfld = trim(outfld)//':'//trim(str)
     end if
     if (len_trim(outfld) >= CXX) then
        write(logunit,*)'fields are = ',trim(outfld)
        write(logunit,*)'fields length = ',len_trim(outfld)  
        call shr_sys_abort(subname//'ERROR: maximum length of xxx_states or xxx_fluxes has been exceeded')
     end if

   end subroutine seq_flds_add

   !===============================================================================
   !BOP ===========================================================================
   !
   ! !IROUTINE: catFields
   !
   ! !DESCRIPTION:
   !  Returns {\tt nfld} concatentated field lists
   !  in the output character string {\tt outfield}.
   !
   ! !REVISION HISTORY:
   !  2003-Jan-24  - T. Craig - first version
   !
   ! !INTERFACE: ------------------------------------------------------------------

   subroutine catFields(outfield, str1, str2)

     ! !USES:

     ! !INPUT/OUTPUT PARAMETERS:

     character(len=*),intent(inout) :: outfield   ! output field name
     character(len=*),intent(in)    :: str1       ! string1 
     character(len=*),intent(in )   :: str2       ! string2

     !EOP

     character(len=*),parameter :: subname = '(seq_flds_catFields) '
     !-------------------------------------------------------------------------------
     !
     !-------------------------------------------------------------------------------

     outfield = ''
     if (len_trim(str1) > 0 .and. len_trim(str2) > 0) then
        if (len_trim(str1) + len_trim(str2) + 1 > len(outfield)) then
           call shr_sys_abort(subname//' ERROR: maximum length of string has been exceeded sum')
        endif
        outfield = trim(str1)//':'//trim(str2)
     else
        if (len_trim(str1) > 0) then
           if (len_trim(str1) > len(outfield)) then
              call shr_sys_abort(subname//' ERROR: maximum length of string has been exceeded str1')
           endif
           outfield = trim(str1)
        endif
        if (len_trim(str2) > 0) then
           if (len_trim(str2) > len(outfield)) then
              call shr_sys_abort(subname//' ERROR: maximum length of string has been exceeded str2')
           endif
           outfield = trim(str2)
        endif
     endif

   end subroutine catFields

   !===============================================================================
   !BOP ===========================================================================
   !
   ! !IROUTINE: seq_flds_getField
   !
   ! !DESCRIPTION:
   !  Returns {\tt nfld} element of the colon-delimited string {\tt cstring}
   !  in the output character string {\tt outfield}.
   !
   ! !REVISION HISTORY:
   !  2003-Jan-24  - T. Craig - first version
   !
   ! !INTERFACE: ------------------------------------------------------------------

   subroutine seq_flds_getField(outfield, nfld, cstring)

     ! !USES:
     use mct_mod

     ! !INPUT/OUTPUT PARAMETERS:

     character(len=*),intent(out) :: outfield   ! output field name
     integer         ,intent(in ) :: nfld       ! field number
     character(len=*),intent(in ) :: cstring    ! colon delimited field string

     !EOP

     type(mct_list)   :: mctIstr  ! mct list from input cstring
     type(mct_string) :: mctOStr  ! mct string for output outfield
     character(len=*),parameter :: subname = '(seq_flds_getField) '

     !-------------------------------------------------------------------------------
     !
     !-------------------------------------------------------------------------------

     outfield = ''

     call mct_list_init(mctIstr,cstring)
     call mct_list_get(mctOStr,nfld,mctIstr)
     outfield = mct_string_toChar(mctOStr)
     call mct_list_clean(mctIstr)
     call mct_string_clean(mctOStr)

   end subroutine seq_flds_getField

   !===============================================================================

   subroutine metadata_set(attname , longname, stdname , units   )

     ! !USES:
     implicit none

     ! !INPUT/OUTPUT PARAMETERS:
     character(len=*), intent(in) :: attname  
     character(len=*), intent(in) :: longname
     character(len=*), intent(in) :: stdname  
     character(len=*), intent(in) :: units    

     !EOP
     character(len=*),parameter :: subname = '(seq_flds_metadata_set) '

     n_entries = n_entries + 1
     if (n_entries > nmax) then
        write(logunit,*)'n_entries= ',n_entries,' nmax = ',nmax,' attname= ',trim(attname)
        call shr_sys_abort(subname//'ERROR: nmax fields in lookup_entry table exceeded') 
     end if

     lookup_entry(n_entries,1) = trim(attname )
     lookup_entry(n_entries,2) = trim(longname)
     lookup_entry(n_entries,3) = trim(stdname )
     lookup_entry(n_entries,4) = trim(units   )


   end subroutine metadata_set

   !===============================================================================

   subroutine seq_flds_esmf_metadata_get(shortname, longname, stdname, units)

     ! !USES:
     use shr_string_mod, only : shr_string_lastindex
     implicit none

     ! !INPUT/OUTPUT PARAMETERS:
     character(len=*), intent(in)  :: shortname 
     character(len=*),optional, intent(out) :: longname
     character(len=*),optional, intent(out) :: stdname  
     character(len=*),optional, intent(out) :: units    

     !EOP

     !--- local ---
     integer :: i,n
     character(len=80) :: llongname, lstdname, lunits, lshortname  ! local copies
     character(len=*),parameter :: undef = 'undefined'
     character(len=*),parameter :: unknown = 'unknown'
     logical :: found
     character(len=*),parameter :: subname = '(seq_flds_esmf_metadata_get) '

     !--- define field metadata (name, long_name, standard_name, units) ---

     llongname = trim(unknown)
     lstdname  = trim(unknown)
     lunits    = trim(unknown)

     found = .false.

     if (.not.found) then
        i = 1
        do while (i <= n_entries .and. .not.found)
           lshortname = trim(shortname)
           if (trim(lshortname) == trim(lookup_entry(i,1))) then
              llongname = trim(lookup_entry(i,2)) 
              lstdname  = trim(lookup_entry(i,3))
              lunits    = trim(lookup_entry(i,4))    
              found     =.true.
           end if
           i = i + 1
        end do
     endif

     if (.not.found) then
        i = 1
        do while (i <= n_entries .and. .not.found)
           n = shr_string_lastIndex(shortname,"_")
           lshortname = ""
           if (n < len_trim(shortname)) lshortname = shortname(n+1:len_trim(shortname))
           if (trim(lshortname) == trim(lookup_entry(i,1))) then
              llongname = trim(lookup_entry(i,2)) 
              lstdname  = trim(lookup_entry(i,3))
              lunits    = trim(lookup_entry(i,4))    
              found     = .true.
           end if
           i = i + 1
        end do
     endif

     if (present(longname)) then
        longname = trim(llongname)
     endif
     if (present(stdname))  then
        stdname = trim(lstdname)
     endif
     if (present(units)) then
        units = trim(lunits)
     endif

   end subroutine seq_flds_esmf_metadata_get

 end module seq_flds_mod


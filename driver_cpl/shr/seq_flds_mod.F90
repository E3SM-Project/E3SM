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
   use seq_drydep_mod, only : seq_drydep_init, seq_drydep_readnl, lnd_drydep
   use seq_comm_mct,   only : seq_comm_iamroot, seq_comm_setptrs, logunit
   use shr_megan_mod,  only : shr_megan_readnl, shr_megan_mechcomps_n
   use shr_fire_emis_mod,  only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
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
   character(len=CXX) :: fire_emis_fields    ! List of fire emission fields
   character(len=CX)  :: carma_fields        ! List of CARMA fields from lnd->atm
   integer            :: ice_ncat            ! number of sea ice thickness categories
   logical            :: seq_flds_i2o_per_cat! .true. if select per ice thickness category fields are passed from ice to ocean

   !----------------------------------------------------------------------------
   ! metadata
   !----------------------------------------------------------------------------

   character(len=*),parameter :: undef     = 'undefined'
   integer         ,parameter :: nmax      = 1000        ! maximum number of entries in lookup_entry
   integer                    :: n_entries = 0           ! actual number of entries in lookup_entry
   character(len=CSS), dimension(nmax, 4) :: lookup_entry = undef

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
   character(CXX) :: seq_flds_l2x_states_to_glc
   character(CXX) :: seq_flds_l2x_fluxes
   character(CXX) :: seq_flds_l2x_fluxes_to_glc
   character(CXX) :: seq_flds_x2l_states
   character(CXX) :: seq_flds_x2l_states_from_glc
   character(CXX) :: seq_flds_x2l_fluxes
   character(CXX) :: seq_flds_x2l_fluxes_from_glc

   character(CXX) :: seq_flds_o2x_states
   character(CXX) :: seq_flds_o2x_fluxes
   character(CXX) :: seq_flds_x2o_states
   character(CXX) :: seq_flds_x2o_fluxes

   character(CXX) :: seq_flds_g2x_states
   character(CXX) :: seq_flds_g2x_states_to_lnd
   character(CXX) :: seq_flds_g2x_fluxes
   character(CXX) :: seq_flds_g2x_fluxes_to_lnd
   character(CXX) :: seq_flds_x2g_states
   character(CXX) :: seq_flds_x2g_fluxes

   character(CXX) :: seq_flds_w2x_states
   character(CXX) :: seq_flds_w2x_fluxes
   character(CXX) :: seq_flds_x2w_states
   character(CXX) :: seq_flds_x2w_fluxes

   character(CXX) :: seq_flds_xao_albedo
   character(CXX) :: seq_flds_xao_states
   character(CXX) :: seq_flds_xao_fluxes
   character(CXX) :: seq_flds_xao_diurnl  ! for diurnal cycle

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
   character(CXX) :: seq_flds_l2x_fields_to_glc
   character(CXX) :: seq_flds_x2l_fields
   character(CXX) :: seq_flds_x2l_fields_from_glc
   character(CXX) :: seq_flds_o2x_fields
   character(CXX) :: seq_flds_x2o_fields
   character(CXX) :: seq_flds_xao_fields
   character(CXX) :: seq_flds_r2x_fields
   character(CXX) :: seq_flds_x2r_fields
   character(CXX) :: seq_flds_g2x_fields
   character(CXX) :: seq_flds_g2x_fields_to_lnd
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

   subroutine seq_flds_set(nmlfile, ID, infodata)

! !USES:
     use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit
     use shr_string_mod, only : shr_string_listIntersect
     use shr_mpi_mod,    only : shr_mpi_bcast
     use glc_elevclass_mod, only : glc_elevclass_init
     use seq_infodata_mod, only : seq_infodata_type, seq_infodata_getdata

! !INPUT/OUTPUT PARAMETERS:
     character(len=*), intent(in) :: nmlfile   ! Name-list filename
     integer         , intent(in) :: ID        ! seq_comm ID
     type(seq_infodata_type), intent(in) :: infodata

! !LOCAL VARIABLES:
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
     character(len=CSS) :: cime_model

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
     character(CXX) :: l2x_states_to_glc = ''
     character(CXX) :: l2x_fluxes = ''
     character(CXX) :: l2x_fluxes_to_glc = ''
     character(CXX) :: x2l_states = ''
     character(CXX) :: x2l_states_from_glc = ''
     character(CXX) :: x2l_fluxes = ''
     character(CXX) :: x2l_fluxes_from_glc = ''
     character(CXX) :: o2x_states = ''
     character(CXX) :: o2x_fluxes = ''
     character(CXX) :: x2o_states = ''
     character(CXX) :: x2o_fluxes = ''
     character(CXX) :: g2x_states = ''
     character(CXX) :: g2x_states_to_lnd = ''
     character(CXX) :: g2x_fluxes = ''
     character(CXX) :: g2x_fluxes_to_lnd = ''
     character(CXX) :: x2g_states = ''
     character(CXX) :: x2g_fluxes = ''
     character(CXX) :: xao_albedo = ''
     character(CXX) :: xao_states = ''
     character(CXX) :: xao_fluxes = ''
     character(CXX) :: xao_diurnl = ''
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
     logical :: flds_bgc
     logical :: flds_wiso
     integer :: glc_nec

     namelist /seq_cplflds_inparm/  &
          flds_co2a, flds_co2b, flds_co2c, flds_co2_dmsa, flds_wiso, glc_nec, &
          ice_ncat, seq_flds_i2o_per_cat, flds_bgc

     ! user specified new fields
     integer,  parameter :: nfldmax = 200
     character(len=CLL)  :: cplflds_custom(nfldmax) = ''

     namelist /seq_cplflds_userspec/ &
          cplflds_custom

     character(len=*),parameter :: subname = '(seq_flds_set) '

!-------------------------------------------------------------------------------

     call seq_comm_setptrs(ID,mpicom=mpicom)

     call seq_infodata_GetData(infodata, cime_model=cime_model)

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
        flds_bgc  = .false.
        flds_wiso = .false.
        glc_nec   = 0
        ice_ncat  = 1
        seq_flds_i2o_per_cat = .false.

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
     call shr_mpi_bcast(flds_bgc     , mpicom)
     call shr_mpi_bcast(flds_wiso    , mpicom)
     call shr_mpi_bcast(glc_nec      , mpicom)
     call shr_mpi_bcast(ice_ncat     , mpicom)
     call shr_mpi_bcast(seq_flds_i2o_per_cat, mpicom)

     call glc_elevclass_init(glc_nec)

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
     longname = 'latitude'
     stdname  = 'latitude'
     units    = 'degrees north'
     attname  = 'lat'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_coord,'lon')
     longname = 'longitude'
     stdname  = 'longitude'
     units    = 'degrees east'
     attname  = 'lon'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_coord,'hgt')
     longname = 'height'
     stdname  = 'height, depth, or levels'
     units    = 'unitless'
     attname  = 'hgt'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'area')
     longname = 'cell_area_model'
     stdname  = 'cell area from model'
     units    = 'radian^2'
     attname  = 'area'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'aream')
     longname = 'cell_area_mapping'
     stdname  = 'cell area from mapping file'
     units    = 'radian^2'
     attname  = 'aream'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'mask')
     longname = 'mask'
     stdname  = 'mask'
     units    = '1'
     attname  = 'mask'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(dom_other,'frac')
     longname = 'area_fraction'
     stdname  = 'area fraction'
     units    = '1'
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

     ! topographic height (m)
     call seq_flds_add(a2x_states,"Sa_topo")
     call seq_flds_add(x2l_states,"Sa_topo")
     longname = 'Surface height'
     stdname  = 'height'
     units    = 'm'
     attname  = 'Sa_topo'
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

     ! specific humidity at the lowest model level (kg/kg)
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
     if (cime_model == 'acme') then
        call seq_flds_add(x2o_states,"Sa_pbot")
     end if
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
     units    = '1'
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
     units    = '1'
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
     units    = '1'
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
     units    = '1'
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
     units    = '1'
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
     longname = 'aerodynamic resistance'
     stdname = 'aerodynamic_resistance'
     attname = 'SI_ram1'
     units = 's/m'
     call metadata_set(attname, longname, stdname, units)


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
     units    = '1'
     attname  = 'Si_ifrac'
     call metadata_set(attname, longname, stdname, units)

     if (trim(cime_model) == 'acme') then
        ! Sea ice basal pressure 
        call seq_flds_add(i2x_states,"Si_bpress")
        call seq_flds_add(x2o_states,"Si_bpress")
        longname = 'Sea ice basal pressure'
        stdname  = 'cice_basal_pressure'
        units    = 'Pa'
        attname  = 'Si_bpress'
        call metadata_set(attname, longname, stdname, units)
     end if

     ! Ocean melt and freeze potential 
     call seq_flds_add(o2x_fluxes,"Fioo_q")
     call seq_flds_add(x2i_fluxes,"Fioo_q")
     longname = 'Ocean melt and freeze potential'
     stdname  = 'surface_snow_and_ice_melt_heat_flux'
     units    = 'W m-2'
     attname  = 'Fioo_q'
     call metadata_set(attname, longname, stdname, units)

     if (trim(cime_model) == 'acme') then
        ! Ocean melt (q<0) potential
        call seq_flds_add(o2x_fluxes,"Fioo_meltp")
        call seq_flds_add(x2i_fluxes,"Fioo_meltp")
        longname = 'Ocean melt (q<0) potential'
        stdname  = 'surface_snow_and_ice_melt_heat_flux'
        units    = 'W m-2'
        attname  = 'Fioo_meltp'
        call metadata_set(attname, longname, stdname, units)
     end if

     if (trim(cime_model) == 'acme') then
        ! Ocean frazil production
        call seq_flds_add(o2x_fluxes,"Fioo_frazil")
        call seq_flds_add(x2i_fluxes,"Fioo_frazil")
        longname = 'Ocean frazil production'
        stdname  = 'ocean_frazil_ice_production'
        units    = 'kg m-2 s-1'
        attname  = 'Fioo_frazil'
        call metadata_set(attname, longname, stdname, units)
     end if

     ! Heat flux from melting
     call seq_flds_add(i2x_fluxes,"Fioi_melth")
     call seq_flds_add(x2o_fluxes,"Fioi_melth")
     longname = 'Heat flux from melting'
     stdname  = 'surface_snow_melt_heat_flux'
     units    = 'W m-2'
     attname  = 'Fioi_melth'
     call metadata_set(attname, longname, stdname, units)

     ! Water flux from melting
     call seq_flds_add(i2x_fluxes,"Fioi_meltw")
     call seq_flds_add(x2o_fluxes,"Fioi_meltw")
     longname = 'Water flux due to melting'
     stdname  = 'surface_snow_melt_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fioi_meltw'
     call metadata_set(attname, longname, stdname, units)

     ! Salt flux
     call seq_flds_add(i2x_fluxes,"Fioi_salt")
     call seq_flds_add(x2o_fluxes,"Fioi_salt")
     longname = 'Salt flux'
     stdname  = 'virtual_salt_flux_into_sea_water'
     units    = 'kg m-2 s-1'
     attname  = 'Fioi_salt'
     call metadata_set(attname, longname, stdname, units)

     ! Black Carbon hydrophilic deposition  
     call seq_flds_add(i2x_fluxes,"Fioi_bcphi" )
     call seq_flds_add(x2o_fluxes,"Fioi_bcphi"   )
     longname = 'Hydrophylic black carbon deposition flux'
     stdname  = 'deposition_flux_of_hydrophylic_black_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Fioi_bcphi'
     call metadata_set(attname, longname, stdname, units)

     ! Black Carbon hydrophobic deposition  
     call seq_flds_add(i2x_fluxes,"Fioi_bcpho" )
     call seq_flds_add(x2o_fluxes,"Fioi_bcpho"   )
     longname = 'Hydrophobic black carbon deposition flux'
     stdname  = 'deposition_flux_of_hydrophobic_black_carbon'
     units    = 'kg m-2 s-1'
     attname  = 'Fioi_bcpho'
     call metadata_set(attname, longname, stdname, units)

     ! Dust flux
     call seq_flds_add(i2x_fluxes,"Fioi_flxdst")
     call seq_flds_add(x2o_fluxes,"Fioi_flxdst")
     longname = 'Dust flux'
     stdname  = 'dust_flux'
     units    = 'kg m-2 s-1'
     attname  = 'Fioi_flxdst'
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

     call seq_flds_add(xao_states,"So_fswpen")
     call seq_flds_add(o2x_states,"So_fswpen")
     longname = 'Fraction of sw penetrating surface layer for diurnal cycle'
     stdname  = 'Fraction of sw penetrating surface layer for diurnal cycle'
     units    = '1'
     attname  = 'So_fswpen'
     call metadata_set(attname, longname, stdname, units)

     !------------------------------
     ! ice<->ocn only exchange - BGC
     !------------------------------
     if (trim(cime_model) == 'acme' .and. flds_bgc) then

        ! Ocean algae concentration 1 - diatoms?
        call seq_flds_add(o2x_states,"So_algae1")
        call seq_flds_add(x2i_states,"So_algae1")
        longname = 'Ocean algae concentration 1 - diatoms'
        stdname  = 'ocean_algae_conc_1'
        units    = 'mmol C m-3'
        attname  = 'So_algae1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean algae concentration 2 - flagellates?
        call seq_flds_add(o2x_states,"So_algae2")
        call seq_flds_add(x2i_states,"So_algae2")
        longname = 'Ocean algae concentration 2 - flagellates'
        stdname  = 'ocean_algae_conc_2'
        units    = 'mmol C m-3'
        attname  = 'So_algae2'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean algae concentration 3 - phaeocyctis?
        call seq_flds_add(o2x_states,"So_algae3")
        call seq_flds_add(x2i_states,"So_algae3")
        longname = 'Ocean algae concentration 3 - phaeocyctis'
        stdname  = 'ocean_algae_conc_3'
        units    = 'mmol C m-3'
        attname  = 'So_algae3'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved organic carbon concentration 1 - saccharides?
        call seq_flds_add(o2x_states,"So_doc1")
        call seq_flds_add(x2i_states,"So_doc1")
        longname = 'Ocean dissolved organic carbon concentration 1 - saccharides'
        stdname  = 'ocean_dissolved_organic_carbon_conc_1'
        units    = 'mmol C m-3'
        attname  = 'So_doc1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved organic carbon concentration 2 - lipids?
        call seq_flds_add(o2x_states,"So_doc2")
        call seq_flds_add(x2i_states,"So_doc2")
        longname = 'Ocean dissolved organic carbon concentration 2 - lipids'
        stdname  = 'ocean_dissolved_organic_carbon_conc_2'
        units    = 'mmol C m-3'
        attname  = 'So_doc2'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved organic carbon concentration 3 - tbd?
        call seq_flds_add(o2x_states,"So_doc3")
        call seq_flds_add(x2i_states,"So_doc3")
        longname = 'Ocean dissolved organic carbon concentration 3 - tbd'
        stdname  = 'ocean_dissolved_organic_carbon_conc_3'
        units    = 'mmol C m-3'
        attname  = 'So_doc3'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved inorganic carbon concentration 1
        call seq_flds_add(o2x_states,"So_dic1")
        call seq_flds_add(x2i_states,"So_dic1")
        longname = 'Ocean dissolved inorganic carbon concentration 1'
        stdname  = 'ocean_dissolved_inorganic_carbon_conc_1'
        units    = 'mmol C m-3'
        attname  = 'So_dic1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved organic nitrogen concentration 1
        call seq_flds_add(o2x_states,"So_don1")
        call seq_flds_add(x2i_states,"So_don1")
        longname = 'Ocean dissolved organic nitrogen concentration 1'
        stdname  = 'ocean_dissolved_organic_nitrogen_conc_1'
        units    = 'mmol N m-3'
        attname  = 'So_don1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean nitrate concentration
        call seq_flds_add(o2x_states,"So_no3")
        call seq_flds_add(x2i_states,"So_no3")
        longname = 'Ocean nitrate concentration'
        stdname  = 'ocean_nitrate_conc'
        units    = 'mmol N m-3'
        attname  = 'So_no3'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean silicate concentration
        call seq_flds_add(o2x_states,"So_sio3")
        call seq_flds_add(x2i_states,"So_sio3")
        longname = 'Ocean silicate concentration'
        stdname  = 'ocean_silicate_conc'
        units    = 'mmol Si m-3'
        attname  = 'So_sio3'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean ammonium concentration
        call seq_flds_add(o2x_states,"So_nh4")
        call seq_flds_add(x2i_states,"So_nh4")
        longname = 'Ocean ammonium concentration'
        stdname  = 'ocean_ammonium_conc'
        units    = 'mmol N m-3'
        attname  = 'So_nh4'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dimethyl sulfide (DMS) concentration
        call seq_flds_add(o2x_states,"So_dms")
        call seq_flds_add(x2i_states,"So_dms")
        longname = 'Ocean dimethyl sulfide concentration'
        stdname  = 'ocean_dimethyl_sulfide_conc'
        units    = 'mmol S m-3'
        attname  = 'So_dms'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dimethylsulphonio-propionate (DMSP) concentration
        call seq_flds_add(o2x_states,"So_dmsp")
        call seq_flds_add(x2i_states,"So_dmsp")
        longname = 'Ocean dimethylsulphonio-propionate concentration'
        stdname  = 'ocean_dimethylsulphoniopropionate_conc'
        units    = 'mmol S m-3'
        attname  = 'So_dmsp'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean DOCr concentration
        call seq_flds_add(o2x_states,"So_docr")
        call seq_flds_add(x2i_states,"So_docr")
        longname = 'Ocean DOCr concentration'
        stdname  = 'ocean_DOCr_conc'
        units    = 'mmol C m-3'
        attname  = 'So_docr'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean particulate iron concentration 1
        call seq_flds_add(o2x_states,"So_fep1")
        call seq_flds_add(x2i_states,"So_fep1")
        longname = 'Ocean particulate iron concentration 1'
        stdname  = 'ocean_particulate_iron_conc_1'
        units    = 'umol Fe m-3'
        attname  = 'So_fep1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean particulate iron concentration 2
        call seq_flds_add(o2x_states,"So_fep2")
        call seq_flds_add(x2i_states,"So_fep2")
        longname = 'Ocean particulate iron concentration 2'
        stdname  = 'ocean_particulate_iron_conc_2'
        units    = 'umol Fe m-3'
        attname  = 'So_fep2'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved iron concentration 1
        call seq_flds_add(o2x_states,"So_fed1")
        call seq_flds_add(x2i_states,"So_fed1")
        longname = 'Ocean dissolved iron concentration 1'
        stdname  = 'ocean_dissolved_iron_conc_1'
        units    = 'umol Fe m-3'
        attname  = 'So_fed1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean dissolved iron concentration 2
        call seq_flds_add(o2x_states,"So_fed2")
        call seq_flds_add(x2i_states,"So_fed2")
        longname = 'Ocean dissolved iron concentration 2'
        stdname  = 'ocean_dissolved_iron_conc_2'
        units    = 'umol Fe m-3'
        attname  = 'So_fed2'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 1
        call seq_flds_add(o2x_states,"So_zaer1")
        call seq_flds_add(x2i_states,"So_zaer1")
        longname = 'Ocean z-aerosol concentration 1'
        stdname  = 'ocean_z_aerosol_conc_1'
        units    = 'kg m-3'
        attname  = 'So_zaer1'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 2
        call seq_flds_add(o2x_states,"So_zaer2")
        call seq_flds_add(x2i_states,"So_zaer2")
        longname = 'Ocean z-aerosol concentration 2'
        stdname  = 'ocean_z_aerosol_conc_2'
        units    = 'kg m-3'
        attname  = 'So_zaer2'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 3
        call seq_flds_add(o2x_states,"So_zaer3")
        call seq_flds_add(x2i_states,"So_zaer3")
        longname = 'Ocean z-aerosol concentration 3'
        stdname  = 'ocean_z_aerosol_conc_3'
        units    = 'kg m-3'
        attname  = 'So_zaer3'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 4
        call seq_flds_add(o2x_states,"So_zaer4")
        call seq_flds_add(x2i_states,"So_zaer4")
        longname = 'Ocean z-aerosol concentration 4'
        stdname  = 'ocean_z_aerosol_conc_4'
        units    = 'kg m-3'
        attname  = 'So_zaer4'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 5
        call seq_flds_add(o2x_states,"So_zaer5")
        call seq_flds_add(x2i_states,"So_zaer5")
        longname = 'Ocean z-aerosol concentration 5'
        stdname  = 'ocean_z_aerosol_conc_5'
        units    = 'kg m-3'
        attname  = 'So_zaer5'
        call metadata_set(attname, longname, stdname, units)

        ! Ocean z-aerosol concentration 6
        call seq_flds_add(o2x_states,"So_zaer6")
        call seq_flds_add(x2i_states,"So_zaer6")
        longname = 'Ocean z-aerosol concentration 6'
        stdname  = 'ocean_z_aerosol_conc_6'
        units    = 'kg m-3'
        attname  = 'So_zaer6'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice algae flux 1 - diatoms?
        call seq_flds_add(i2x_fluxes,"Fioi_algae1")
        call seq_flds_add(x2o_fluxes,"Fioi_algae1")
        longname = 'Sea ice algae flux 1 - diatoms'
        stdname  = 'seaice_algae_flux_1'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_algae1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice algae flux 2 - flagellates?
        call seq_flds_add(i2x_fluxes,"Fioi_algae2")
        call seq_flds_add(x2o_fluxes,"Fioi_algae2")
        longname = 'Sea ice algae flux 2 - flagellates'
        stdname  = 'seaice_algae_flux_2'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_algae2'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice algae flux 3 - phaeocyctis?
        call seq_flds_add(i2x_fluxes,"Fioi_algae3")
        call seq_flds_add(x2o_fluxes,"Fioi_algae3")
        longname = 'Sea ice algae flux 3 - phaeocyctis'
        stdname  = '_algae_flux_3'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_algae3'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved organic carbon flux 1 - saccharides?
        call seq_flds_add(i2x_fluxes,"Fioi_doc1")
        call seq_flds_add(x2o_fluxes,"Fioi_doc1")
        longname = 'Sea ice dissolved organic carbon flux 1 - saccharides'
        stdname  = 'seaice_dissolved_organic_carbon_flux_1'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_doc1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved organic carbon flux 2 - lipids?
        call seq_flds_add(i2x_fluxes,"Fioi_doc2")
        call seq_flds_add(x2o_fluxes,"Fioi_doc2")
        longname = 'Sea ice dissolved organic carbon flux 2 - lipids'
        stdname  = 'seaice_dissolved_organic_carbon_flux_2'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_doc2'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved organic carbon flux 3 - tbd?
        call seq_flds_add(i2x_fluxes,"Fioi_doc3")
        call seq_flds_add(x2o_fluxes,"Fioi_doc3")
        longname = 'Sea ice dissolved organic carbon flux 3 - tbd'
        stdname  = 'seaice_dissolved_organic_carbon_flux_3'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_doc3'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved inorganic carbon flux 1
        call seq_flds_add(i2x_fluxes,"Fioi_dic1")
        call seq_flds_add(x2o_fluxes,"Fioi_dic1")
        longname = 'Sea ice dissolved inorganic carbon flux 1'
        stdname  = 'seaice_dissolved_inorganic_carbon_flux_1'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_dic1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved organic nitrogen flux 1
        call seq_flds_add(i2x_fluxes,"Fioi_don1")
        call seq_flds_add(x2o_fluxes,"Fioi_don1")
        longname = 'Sea ice dissolved organic nitrogen flux 1'
        stdname  = 'seaice_dissolved_organic_nitrogen_flux_1'
        units    = 'mmol N m-2 s-1'
        attname  = 'Fioi_don1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice nitrate flux
        call seq_flds_add(i2x_fluxes,"Fioi_no3")
        call seq_flds_add(x2o_fluxes,"Fioi_no3")
        longname = 'Sea ice nitrate flux'
        stdname  = 'seaice_nitrate_flux'
        units    = 'mmol N m-2 s-1'
        attname  = 'Fioi_no3'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice silicate flux
        call seq_flds_add(i2x_fluxes,"Fioi_sio3")
        call seq_flds_add(x2o_fluxes,"Fioi_sio3")
        longname = 'Sea ice silicate flux'
        stdname  = 'seaice_silicate_flux'
        units    = 'mmol Si m-2 s-1'
        attname  = 'Fioi_sio3'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice ammonium flux
        call seq_flds_add(i2x_fluxes,"Fioi_nh4")
        call seq_flds_add(x2o_fluxes,"Fioi_nh4")
        longname = 'Sea ice ammonium flux'
        stdname  = 'seaice_ammonium_flux'
        units    = 'mmol N m-2 s-1'
        attname  = 'Fioi_nh4'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dimethyl sulfide (DMS) flux
        call seq_flds_add(i2x_fluxes,"Fioi_dms")
        call seq_flds_add(x2o_fluxes,"Fioi_dms")
        longname = 'Sea ice dimethyl sulfide flux'
        stdname  = 'seaice_dimethyl_sulfide_flux'
        units    = 'mmol S m-2 s-1'
        attname  = 'Fioi_dms'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice DMSPp flux
        call seq_flds_add(i2x_fluxes,"Fioi_dmspp")
        call seq_flds_add(x2o_fluxes,"Fioi_dmspp")
        longname = 'Sea ice DSMPp flux'
        stdname  = 'seaice_DSMPp_flux'
        units    = 'mmol S m-2 s-1'
        attname  = 'Fioi_dmspp'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice DMSPd flux
        call seq_flds_add(i2x_fluxes,"Fioi_dmspd")
        call seq_flds_add(x2o_fluxes,"Fioi_dmspd")
        longname = 'Sea ice DSMPd flux'
        stdname  = 'seaice_DSMPd_flux'
        units    = 'mmol S m-2 s-1'
        attname  = 'Fioi_dmspd'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice DOCr flux
        call seq_flds_add(i2x_fluxes,"Fioi_docr")
        call seq_flds_add(x2o_fluxes,"Fioi_docr")
        longname = 'Sea ice DOCr flux'
        stdname  = 'seaice_DOCr_flux'
        units    = 'mmol C m-2 s-1'
        attname  = 'Fioi_docr'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice particulate iron flux 1
        call seq_flds_add(i2x_fluxes,"Fioi_fep1")
        call seq_flds_add(x2o_fluxes,"Fioi_fep1")
        longname = 'Sea ice particulate iron flux 1'
        stdname  = 'seaice_particulate_iron_flux_1'
        units    = 'umol Fe m-2 s-1'
        attname  = 'Fioi_fep1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice particulate iron flux 2
        call seq_flds_add(i2x_fluxes,"Fioi_fep2")
        call seq_flds_add(x2o_fluxes,"Fioi_fep2")
        longname = 'Sea ice particulate iron flux 2'
        stdname  = 'seaice_particulate_iron_flux_2'
        units    = 'umol Fe m-2 s-1'
        attname  = 'Fioi_fep2'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved iron flux 1
        call seq_flds_add(i2x_fluxes,"Fioi_fed1")
        call seq_flds_add(x2o_fluxes,"Fioi_fed1")
        longname = 'Sea ice dissolved iron flux 1'
        stdname  = 'seaice_dissolved_iron_flux_1'
        units    = 'umol Fe m-2 s-1'
        attname  = 'Fioi_fed1'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice dissolved iron flux 2
        call seq_flds_add(i2x_fluxes,"Fioi_fed2")
        call seq_flds_add(x2o_fluxes,"Fioi_fed2")
        longname = 'Sea ice dissolved iron flux 2'
        stdname  = 'seaice_dissolved_iron_flux_2'
        units    = 'umol Fe m-2 s-1'
        attname  = 'Fioi_fed2'
        call metadata_set(attname, longname, stdname, units)

        ! Sea ice iron dust
        call seq_flds_add(i2x_fluxes,"Fioi_dust1")
        call seq_flds_add(x2o_fluxes,"Fioi_dust1")
        longname = 'Sea ice iron dust 1'
        stdname  = 'seaice_iron_dust_1'
        units    = 'kg m-2 s-1'
        attname  = 'Fioi_dust1'
        call metadata_set(attname, longname, stdname, units)

     endif


     !-----------------------------
     ! lnd->rof exchange
     ! TODO: put in attributes below
     !-----------------------------

     call seq_flds_add(l2x_fluxes,'Flrl_rofsur')
     call seq_flds_add(x2r_fluxes,'Flrl_rofsur')
     longname = 'Water flux from land (liquid surface)'
     stdname  = 'water_flux_into_runoff_surface'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofsur'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(l2x_fluxes,'Flrl_rofgwl')
     call seq_flds_add(x2r_fluxes,'Flrl_rofgwl')
     longname = 'Water flux from land (liquid glacier, wetland, and lake)'
     stdname  = 'water_flux_into_runoff_from_gwl'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofgwl'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(l2x_fluxes,'Flrl_rofsub')
     call seq_flds_add(x2r_fluxes,'Flrl_rofsub')
     longname = 'Water flux from land (liquid subsurface)'
     stdname  = 'water_flux_into_runoff_subsurface'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofsub'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(l2x_fluxes,'Flrl_rofdto')
     call seq_flds_add(x2r_fluxes,'Flrl_rofdto')
     longname = 'Water flux from land direct to ocean'
     stdname  = 'water_flux_direct_to_ocean'
     units    = 'kg m-2 s-1'
     attname  = 'Flrl_rofdto'
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
     longname = 'River channel total water volume'
     stdname  = 'rtm_volr'
     units    = 'm'
     attname  = 'Flrr_volr'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(r2x_fluxes,'Flrr_volrmch')
     call seq_flds_add(x2l_fluxes,'Flrr_volrmch')
     longname = 'River channel main channel water volume'
     stdname  = 'rtm_volrmch'
     units    = 'm'
     attname  = 'Flrr_volrmch'
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
     ! New xao_states diagnostic
     ! fields for history output only
     !-----------------------------

     call seq_flds_add(xao_fluxes,"Faox_swdn")
     longname = 'Downward solar radiation'
     stdname  = 'surface_downward_shortwave_flux'
     units    = 'W m-2'
     attname  = 'Faox_swdn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_fluxes,"Faox_swup")
     longname = 'Upward solar radiation'
     stdname  = 'surface_upward_shortwave_flux'
     units    = 'W m-2'
     attname  = 'Faox_swup'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_tbulk_diurn")
     longname = 'atm/ocn flux temperature bulk'
     stdname  = 'aoflux_tbulk'
     units    = 'K'
     attname  = 'So_tbulk_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_tskin_diurn")
     longname = 'atm/ocn flux temperature skin'
     stdname  = 'aoflux_tskin'
     units    = 'K'
     attname  = 'So_tskin_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_tskin_night_diurn")
     longname = 'atm/ocn flux temperature skin at night'
     stdname  = 'aoflux_tskin_night'
     units    = 'K'
     attname  = 'So_tskin_night_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_tskin_day_diurn")
     longname = 'atm/ocn flux temperature skin at day'
     stdname  = 'aoflux_tskin_day'
     units    = 'K'
     attname  = 'So_tskin_day_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_cskin_diurn")
     longname = 'atm/ocn flux cool skin'
     stdname  = 'aoflux_cskin'
     units    = 'K'
     attname  = 'So_cskin_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_cskin_night_diurn")
     longname = 'atm/ocn flux cool skin at night'
     stdname  = 'aoflux_cskin_night'
     units    = 'K'
     attname  = 'So_cskin_night_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_warm_diurn")
     longname = 'atm/ocn flux warming'
     stdname  = 'aoflux_warm'
     units    = 'unitless'
     attname  = 'So_warm_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_salt_diurn")
     longname = 'atm/ocn flux salting'
     stdname  = 'aoflux_salt'
     units    = 'unitless'
     attname  = 'So_salt_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_speed_diurn")
     longname = 'atm/ocn flux speed'
     stdname  = 'aoflux_speed'
     units    = 'unitless'
     attname  = 'So_speed_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_regime_diurn")
     longname = 'atm/ocn flux regime'
     stdname  = 'aoflux_regime'
     units    = 'unitless'
     attname  = 'So_regime_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_warmmax_diurn")
     longname = 'atm/ocn flux warming dialy max'
     stdname  = 'aoflux_warmmax'
     units    = 'unitless'
     attname  = 'So_warmmax_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_windmax_diurn")
     longname = 'atm/ocn flux wind daily max'
     stdname  = 'aoflux_windmax'
     units    = 'unitless'
     attname  = 'So_windmax_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_qsolavg_diurn")
     longname = 'atm/ocn flux q-solar daily avg'
     stdname  = 'aoflux_qsolavg'
     units    = 'unitless'
     attname  = 'So_qsolavg_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_windavg_diurn")
     longname = 'atm/ocn flux wind daily avg'
     stdname  = 'aoflux_windavg'
     units    = 'unitless'
     attname  = 'So_windavg_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_warmmaxinc_diurn")
     longname = 'atm/ocn flux daily max increment'
     stdname  = 'aoflux_warmmaxinc'
     units    = 'unitless'
     attname  = 'So_warmmaxinc_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_windmaxinc_diurn")
     longname = 'atm/ocn flux wind daily max increment'
     stdname  = 'aoflux_windmaxinc'
     units    = 'unitless'
     attname  = 'So_windmaxinc_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_qsolinc_diurn")
     longname = 'atm/ocn flux q-solar increment'
     stdname  = 'aoflux_qsolinc'
     units    = 'unitless'
     attname  = 'So_qsolinc_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_windinc_diurn")
     longname = 'atm/ocn flux wind increment'
     stdname  = 'aoflux_windinc'
     units    = 'unitless'
     attname  = 'So_windinc_diurn'
     call metadata_set(attname, longname, stdname, units)

     call seq_flds_add(xao_diurnl,"So_ninc_diurn")
     longname = 'atm/ocn flux increment counter'
     stdname  = 'aoflux_ninc'
     units    = 'unitless'
     attname  = 'So_ninc_diurn'
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
     call seq_flds_add(g2x_states_to_lnd,trim(name))
     call seq_flds_add(x2l_states,trim(name))
     call seq_flds_add(x2l_states_from_glc,trim(name))
     longname = 'Ice sheet grid coverage on global grid'
     stdname  = 'ice_sheet_grid_mask'
     units    = '1'
     attname  = 'Sg_icemask'
     call metadata_set(attname, longname, stdname, units)

     name = 'Sg_icemask_coupled_fluxes'
     call seq_flds_add(g2x_states,trim(name))
     call seq_flds_add(g2x_states_to_lnd,trim(name))
     call seq_flds_add(x2l_states,trim(name))
     call seq_flds_add(x2l_states_from_glc,trim(name))
     longname = 'Ice sheet mask where we are potentially sending non-zero fluxes'
     stdname  = 'icemask_coupled_fluxes'
     units    = '1'
     attname  = 'Sg_icemask_coupled_fluxes'
     call metadata_set(attname, longname, stdname, units)

     ! glc fields with multiple elevation classes: lnd->glc
     !
     ! Note that these fields are sent in multiple elevation classes from lnd->cpl, but
     ! the fields sent from cpl->glc do NOT have elevation classes
     !
     ! Also note that we need to keep track of the l2x fields destined for glc in the
     ! additional variables, l2x_fluxes_to_glc and l2x_states_to_glc. This is needed so that
     ! we can set up an additional attribute vector holding accumulated quantities of just
     ! these fields. (We can't determine these field lists with a call to
     ! mct_aVect_initSharedFields, because the field names differ between l2x and x2g.)

     name = 'Flgl_qice'
     longname = 'New glacier ice flux'
     stdname  = 'ice_flux_out_of_glacier'
     units    = 'kg m-2 s-1'
     attname  = 'Fgll_qice'
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_fluxes)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_fluxes_to_glc, &
          additional_list = .true.)
     call seq_flds_add(x2g_fluxes,trim(name))
     call metadata_set(attname, longname, stdname, units)

     name = 'Sl_tsrf'
     longname = 'Surface temperature of glacier'
     stdname  = 'surface_temperature'
     units    = 'deg C'
     attname  = 'Sl_tsrf'
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_states)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_states_to_glc, &
          additional_list = .true.)
     call seq_flds_add(x2g_states,trim(name))
     call metadata_set(attname, longname, stdname, units)

     ! Sl_topo is sent from lnd -> cpl, but is NOT sent to glc (it is only used for the
     ! remapping in the coupler)
     name = 'Sl_topo'
     longname = 'Surface height'
     stdname  = 'height'
     units    = 'm'
     attname  = 'Sl_topo'
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_states)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, l2x_states_to_glc, &
          additional_list = .true.)

     ! glc fields with multiple elevation classes: glc->lnd
     !
     ! Note that the fields sent from glc->cpl do NOT have elevation classes, but the
     ! fields from cpl->lnd are broken into multiple elevation classes

     name = 'Sg_ice_covered'
     longname = 'Fraction of glacier area'
     stdname  = 'glacier_area_fraction'
     units    = '1'
     attname  = 'Sg_ice_covered'
     call seq_flds_add(g2x_states,trim(name))
     call seq_flds_add(g2x_states_to_lnd,trim(name))
     call metadata_set(attname, longname, stdname, units)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_states)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_states_from_glc, &
          additional_list = .true.)

     name = 'Sg_topo'
     longname = 'Surface height of glacier'
     stdname  = 'height'
     units    = 'm'
     attname  = 'Sg_topo'
     call seq_flds_add(g2x_states,trim(name))
     call seq_flds_add(g2x_states_to_lnd,trim(name))
     call metadata_set(attname, longname, stdname, units)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_states)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_states_from_glc, &
          additional_list = .true.)

     name = 'Flgg_hflx'
     longname = 'Downward heat flux from glacier interior'
     stdname  = 'downward_heat_flux_in_glacier'
     units    = 'W m-2'
     attname  = 'Flgg_hflx'
     call seq_flds_add(g2x_fluxes,trim(name))
     call seq_flds_add(g2x_fluxes_to_lnd,trim(name))
     call metadata_set(attname, longname, stdname, units)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_fluxes)
     call set_glc_elevclass_field(name, attname, longname, stdname, units, x2l_fluxes_from_glc, &
          additional_list = .true.)

     ! Done glc fields

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

     if (flds_wiso) then
        call seq_flds_add(o2x_states, "So_roce_16O")
        call seq_flds_add(x2i_states, "So_roce_16O")
        longname = 'Ratio of ocean surface level abund. H2_16O/H2O/Rstd'
        stdname  = ''
        units    = ' '
        attname  = 'So_roce_16O'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(o2x_states, "So_roce_18O")
        call seq_flds_add(x2i_states, "So_roce_18O")
        longname = 'Ratio of ocean surface level abund. H2_18O/H2O/Rstd'
        attname  = 'So_roce_18O'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(o2x_states, "So_roce_HDO")
        call seq_flds_add(x2i_states, "So_roce_HDO")
        longname = 'Ratio of ocean surface level abund. HDO/H2O/Rstd'
        attname  = 'So_roce_HDO'
        call metadata_set(attname, longname, stdname, units)
 
        !--------------------------------------------
        !Atmospheric specific humidty at lowest level: 
        !--------------------------------------------       

      ! specific humidity of H216O at the lowest model level (kg/kg)
        call seq_flds_add(a2x_states,"Sa_shum_16O")
        call seq_flds_add(x2l_states,"Sa_shum_16O")
        call seq_flds_add(x2i_states,"Sa_shum_16O")
        longname = 'Specific humidty of H216O at the lowest model level'
        stdname  = 'H216OV'
        units    = 'kg kg-1'
        attname  = 'Sa_shum_16O'
        call metadata_set(attname, longname, stdname, units)

       ! specific humidity of HD16O at the lowest model level (kg/kg)
        call seq_flds_add(a2x_states,"Sa_shum_HDO")
        call seq_flds_add(x2l_states,"Sa_shum_HDO")
        call seq_flds_add(x2i_states,"Sa_shum_HDO")
        longname = 'Specific humidty of HD16O at the lowest model level'
        stdname  = 'HD16OV'
        attname  = 'Sa_shum_HDO'
        call metadata_set(attname, longname, stdname, units)

       ! specific humidity of H218O at the lowest model level (kg/kg)
        call seq_flds_add(a2x_states,"Sa_shum_18O")
        call seq_flds_add(x2l_states,"Sa_shum_18O")
        call seq_flds_add(x2i_states,"Sa_shum_18O")
        longname = 'Specific humidty of H218O at the lowest model level'
        stdname  = 'H218OV'
        attname  = 'Sa_shum_18O'
        call metadata_set(attname, longname, stdname, units)

       ! Surface snow water equivalent (land/atm only) 
        call seq_flds_add(l2x_states,"Sl_snowh_16O")
        call seq_flds_add(l2x_states,"Sl_snowh_18O")
        call seq_flds_add(l2x_states,"Sl_snowh_HDO")
        call seq_flds_add(x2a_states,"Sl_snowh_16O")
        call seq_flds_add(x2a_states,"Sl_snowh_18O")
        call seq_flds_add(x2a_states,"Sl_snowh_HDO")
        longname = 'Isotopic surface snow water equivalent'
        stdname  = 'surface_snow_water_equivalent'
        units    = 'm'
        attname  = 'Sl_snowh_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_snowh_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_snowh_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_snowh_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_snowh_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_snowh_HDO'
        call metadata_set(attname, longname, stdname, units)

        !--------------
        !Isotopic Rain:
        !-------------- 

       !Isotopic Precipitation Fluxes:
        units    = 'kg m-2 s-1'
        call seq_flds_add(a2x_fluxes,"Faxa_rainc_16O")
        call seq_flds_add(a2x_fluxes,"Faxa_rainl_16O")
        call seq_flds_add(x2o_fluxes, "Faxa_rain_16O")
        call seq_flds_add(x2l_fluxes,"Faxa_rainc_16O")
        call seq_flds_add(x2l_fluxes,"Faxa_rainl_16O")
        call seq_flds_add(x2i_fluxes, "Faxa_rain_16O")
        longname = 'Water flux due to H216O rain' !equiv. to bulk
        stdname  = 'H2_16O_rainfall_flux'
        attname  = 'Faxa_rain_16O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H216O Convective precipitation rate'
        stdname  = 'H2_16O_convective_precipitation_flux'
        attname  = 'Faxa_rainc_16O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H216O Large-scale (stable) precipitation rate'
        stdname  = 'H2_16O_large_scale_precipitation_flux'
        attname  = 'Faxa_rainl_16O' 
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(a2x_fluxes,"Faxa_rainc_18O")
        call seq_flds_add(a2x_fluxes,"Faxa_rainl_18O")
        call seq_flds_add(x2o_fluxes, "Faxa_rain_18O")
        call seq_flds_add(x2l_fluxes,"Faxa_rainc_18O")
        call seq_flds_add(x2l_fluxes,"Faxa_rainl_18O")
        call seq_flds_add(x2i_fluxes, "Faxa_rain_18O")
        longname = 'Water flux due to H218O rain'
        stdname  = 'h2_18o_rainfall_flux'
        attname  = 'Faxa_rain_18O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H218O Convective precipitation rate'
        stdname  = 'H2_18O_convective_precipitation_flux'
        attname  = 'Faxa_rainc_18O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H218O Large-scale (stable) precipitation rate'
        stdname  = 'H2_18O_large_scale_precipitation_flux'
        attname  = 'Faxa_rainl_18O' 
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(a2x_fluxes,"Faxa_rainc_HDO")
        call seq_flds_add(a2x_fluxes,"Faxa_rainl_HDO")
        call seq_flds_add(x2o_fluxes, "Faxa_rain_HDO")
        call seq_flds_add(x2l_fluxes,"Faxa_rainc_HDO")
        call seq_flds_add(x2l_fluxes,"Faxa_rainl_HDO")
        call seq_flds_add(x2i_fluxes, "Faxa_rain_HDO")
        longname = 'Water flux due to HDO rain'
        stdname  = 'hdo_rainfall_flux'
        attname  = 'Faxa_rain_HDO' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'HDO Convective precipitation rate'
        stdname  = 'HDO_convective_precipitation_flux'
        attname  = 'Faxa_rainc_HDO' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'HDO Large-scale (stable) precipitation rate'
        stdname  = 'HDO_large_scale_precipitation_flux'
        attname  = 'Faxa_rainl_HDO' 
        call metadata_set(attname, longname, stdname, units)

        !-------------
        !Isotopic snow:
        !------------- 

        call seq_flds_add(a2x_fluxes,"Faxa_snowc_16O")
        call seq_flds_add(a2x_fluxes,"Faxa_snowl_16O")
        call seq_flds_add(x2o_fluxes, "Faxa_snow_16O")
        call seq_flds_add(x2l_fluxes,"Faxa_snowc_16O")
        call seq_flds_add(x2l_fluxes,"Faxa_snowl_16O")
        call seq_flds_add(x2i_fluxes, "Faxa_snow_16O")
        longname = 'Water equiv. H216O snow flux'
        stdname  = 'h2_16o_snowfall_flux'
        attname  = 'Faxa_snow_16O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H2_16O Convective snow rate (water equivalent)'
        stdname  = 'H2_16O_convective_snowfall_flux'
        attname  = 'Faxa_snowc_16O'
        call metadata_set(attname, longname, stdname, units)
        longname = 'H2_16O Large-scale (stable) snow rate (water equivalent)'
        stdname  = 'H2_16O_large_scale_snowfall_flux'
        attname  = 'Faxa_snowl_16O' 
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(a2x_fluxes,"Faxa_snowc_18O")
        call seq_flds_add(a2x_fluxes,"Faxa_snowl_18O")
        call seq_flds_add(x2o_fluxes, "Faxa_snow_18O")
        call seq_flds_add(x2l_fluxes,"Faxa_snowc_18O")
        call seq_flds_add(x2l_fluxes,"Faxa_snowl_18O")
        call seq_flds_add(x2i_fluxes, "Faxa_snow_18O")
        longname = 'Isotopic water equiv. snow flux of H218O'
        stdname  = 'h2_18o_snowfall_flux'
        attname  = 'Faxa_snow_18O' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'H2_18O Convective snow rate (water equivalent)'
        stdname  = 'H2_18O_convective_snowfall_flux'
        attname  = 'Faxa_snowc_18O'
        call metadata_set(attname, longname, stdname, units)
        longname = 'H2_18O Large-scale (stable) snow rate (water equivalent)'
        stdname  = 'H2_18O_large_scale_snowfall_flux'
        attname  = 'Faxa_snowl_18O' 
        call metadata_set(attname, longname, stdname, units)
        
        call seq_flds_add(a2x_fluxes,"Faxa_snowc_HDO")
        call seq_flds_add(a2x_fluxes,"Faxa_snowl_HDO")
        call seq_flds_add(x2o_fluxes, "Faxa_snow_HDO")
        call seq_flds_add(x2l_fluxes,"Faxa_snowc_HDO")
        call seq_flds_add(x2l_fluxes,"Faxa_snowl_HDO")
        call seq_flds_add(x2i_fluxes, "Faxa_snow_HDO")
        longname = 'Isotopic water equiv. snow flux of HDO'
        stdname  = 'hdo_snowfall_flux'
        attname  = 'Faxa_snow_HDO' 
        call metadata_set(attname, longname, stdname, units)
        longname = 'HDO Convective snow rate (water equivalent)'
        stdname  = 'HDO_convective_snowfall_flux'
        attname  = 'Faxa_snowc_HDO'
        call metadata_set(attname, longname, stdname, units)
        longname = 'HDO Large-scale (stable) snow rate (water equivalent)'
        stdname  = 'HDO_large_scale_snowfall_flux'
        attname  = 'Faxa_snowl_HDO' 
        call metadata_set(attname, longname, stdname, units)

        !----------------------------------
        !Isotopic precipitation (rain+snow):
        !----------------------------------  

        call seq_flds_add(x2o_fluxes,"Faxa_prec_16O")  ! derived rain+snow
        longname = 'Isotopic Water flux (rain+snow) for H2_16O'
        stdname  = 'h2_18o_precipitation_flux'
        attname  = 'Faxa_prec_16O'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(x2o_fluxes,"Faxa_prec_18O")  ! derived rain+snow
        longname = 'Isotopic Water flux (rain+snow) for H2_18O'
        stdname  = 'h2_18o_precipitation_flux'
        units    = 'kg m-2 s-1'
        attname  = 'Faxa_prec_18O'
        call metadata_set(attname, longname, stdname, units)

        call seq_flds_add(x2o_fluxes,"Faxa_prec_HDO")  ! derived rain+snow
        longname = 'Isotopic Water flux (rain+snow) for HD_O'
        stdname  = 'hdo_precipitation_flux'
        units    = 'kg m-2 s-1'
        attname  = 'Faxa_prec_HDO'
        call metadata_set(attname, longname, stdname, units)

        !-------------------------------------
        !Isotopic two meter reference humidity:
        !-------------------------------------

        ! H216O Reference specific humidity at 2 meters
        call seq_flds_add(l2x_states,"Sl_qref_16O")
        call seq_flds_add(i2x_states,"Si_qref_16O")
        call seq_flds_add(xao_states,"So_qref_16O")
        call seq_flds_add(x2a_states,"Sx_qref_16O")
        longname = 'Reference H216O specific humidity at 2 meters'
        stdname  = 'H216O_specific_humidity'
        units    = 'kg kg-1'
        attname  = 'Si_qref_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_qref_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'So_qref_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sx_qref_16O'
        call metadata_set(attname, longname, stdname, units)

        ! HD16O Reference specific humidity at 2 meters
        call seq_flds_add(l2x_states,"Sl_qref_HDO")
        call seq_flds_add(i2x_states,"Si_qref_HDO")
        call seq_flds_add(xao_states,"So_qref_HDO")
        call seq_flds_add(x2a_states,"Sx_qref_HDO")
        longname = 'Reference HD16O specific humidity at 2 meters'
        stdname  = 'HD16O_specific_humidity'
        units    = 'kg kg-1'
        attname  = 'Si_qref_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_qref_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'So_qref_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sx_qref_HDO'
        call metadata_set(attname, longname, stdname, units)

        ! H218O Reference specific humidity at 2 meters
        call seq_flds_add(l2x_states,"Sl_qref_18O")
        call seq_flds_add(i2x_states,"Si_qref_18O")
        call seq_flds_add(xao_states,"So_qref_18O")
        call seq_flds_add(x2a_states,"Sx_qref_18O")
        longname = 'Reference H218O specific humidity at 2 meters'
        stdname  = 'H218O_specific_humidity'
        units    = 'kg kg-1'
        attname  = 'Si_qref_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sl_qref_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'So_qref_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Sx_qref_18O'
        call metadata_set(attname, longname, stdname, units)

        !-------------------------
        !Isotopic Evaporation flux:
        !-------------------------

       ! H216O Evaporation water flux
        call seq_flds_add(l2x_fluxes,"Fall_evap_16O")
        call seq_flds_add(i2x_fluxes,"Faii_evap_16O") 
        call seq_flds_add(xao_fluxes,"Faox_evap_16O")
        call seq_flds_add(x2a_fluxes,"Faxx_evap_16O")
        call seq_flds_add(x2o_fluxes,"Foxx_evap_16O")
        longname = 'Evaporation H216O flux'
        stdname  = 'H216O_evaporation_flux'
        units    = 'kg m-2 s-1'
        attname  = 'Fall_evap_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faii_evap_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faox_evap_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faxx_evap_16O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Foxx_evap_16O'
        call metadata_set(attname, longname, stdname, units)

      ! HD16O Evaporation water flux
        call seq_flds_add(l2x_fluxes,"Fall_evap_HDO")
        call seq_flds_add(i2x_fluxes,"Faii_evap_HDO") 
        call seq_flds_add(xao_fluxes,"Faox_evap_HDO")
        call seq_flds_add(x2a_fluxes,"Faxx_evap_HDO")
        call seq_flds_add(x2o_fluxes,"Foxx_evap_HDO")
        longname = 'Evaporation HD16O flux'
        stdname  = 'HD16O_evaporation_flux'
        units    = 'kg m-2 s-1'
        attname  = 'Fall_evap_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faii_evap_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faox_evap_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Foxx_evap_HDO'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faxx_evap_HDO'
        call metadata_set(attname, longname, stdname, units)  

      ! H218O Evaporation water flux
        call seq_flds_add(l2x_fluxes,"Fall_evap_18O")
        call seq_flds_add(i2x_fluxes,"Faii_evap_18O")
        call seq_flds_add(xao_fluxes,"Faox_evap_18O")
        call seq_flds_add(x2a_fluxes,"Faxx_evap_18O")
        call seq_flds_add(x2o_fluxes,"Foxx_evap_18O")
        longname = 'Evaporation H218O flux'
        stdname  = 'H218O_evaporation_flux'
        units    = 'kg m-2 s-1'
        attname  = 'Fall_evap_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faii_evap_18O' 
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faox_evap_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Faxx_evap_18O'
        call metadata_set(attname, longname, stdname, units)
        attname  = 'Foxx_evap_18O'
        call metadata_set(attname, longname, stdname, units)

       !-----------------------------
       !Isotopic sea ice melting flux: 
       !-----------------------------

       ! H216O Water flux from melting
       units    = 'kg m-2 s-1'
       call seq_flds_add(i2x_fluxes,"Fioi_meltw_16O")
       call seq_flds_add(x2o_fluxes,"Fioi_meltw_16O")
       longname = 'H2_16O flux due to melting'
       stdname  = 'h2_16o_surface_snow_melt_flux'
       attname  = 'Fioi_meltw_16O'
       call metadata_set(attname, longname, stdname, units)

       ! H218O Water flux from melting
       call seq_flds_add(i2x_fluxes,"Fioi_meltw_18O")
       call seq_flds_add(x2o_fluxes,"Fioi_meltw_18O")
       longname = 'H2_18O flux due to melting'
       stdname  = 'h2_18o_surface_snow_melt_flux'
       attname  = 'Fioi_meltw_18O'
       call metadata_set(attname, longname, stdname, units)

       ! HDO Water flux from melting
       units    = 'kg m-2 s-1'
       call seq_flds_add(i2x_fluxes,"Fioi_meltw_HDO")
       call seq_flds_add(x2o_fluxes,"Fioi_meltw_HDO")
       longname = 'HDO flux due to melting'
       stdname  = 'hdo_surface_snow_melt_flux'
       attname  = 'Fioi_meltw_HDO'
       call metadata_set(attname, longname, stdname, units)

       !Iso-Runoff
       ! l2x, x2r
       units    = 'kg m-2 s-1'
       call seq_flds_add(l2x_fluxes,'Flrl_rofi_16O')
       call seq_flds_add(x2r_fluxes,'Flrl_rofi_16O')
       longname = 'H2_16O Water flux from land (frozen)'
       stdname  = 'H2_16O_frozen_water_flux_into_runoff'
       attname  = 'Flrl_rofi_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(l2x_fluxes,'Flrl_rofi_18O')
       call seq_flds_add(x2r_fluxes,'Flrl_rofi_18O')
       longname = 'H2_18O Water flux from land (frozen)'
       stdname  = 'H2_18O_frozen_water_flux_into_runoff'
       attname  = 'Flrl_rofi_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(l2x_fluxes,'Flrl_rofi_HDO')
       call seq_flds_add(x2r_fluxes,'Flrl_rofi_HDO')
       longname = 'HDO Water flux from land (frozen)'
       stdname  = 'HDO_frozen_water_flux_into_runoff'
       attname  = 'Flrl_rofi_HDO'
       call metadata_set(attname, longname, stdname, units)

       call seq_flds_add(l2x_fluxes,'Flrl_rofl_16O')
       call seq_flds_add(x2r_fluxes,'Flrl_rofl_16O')
       longname = 'H2_16O Water flux from land (liquid)'
       stdname  = 'H2_16O_liquid_water_flux_into_runoff'
       attname  = 'Flrl_rofl_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(l2x_fluxes,'Flrl_rofl_18O')
       call seq_flds_add(x2r_fluxes,'Flrl_rofl_18O')
       longname = 'H2_18O Water flux from land (liquid)'
       stdname  = 'H2_18O_liquid_water_flux_into_runoff'
       attname  = 'Flrl_rofl_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(l2x_fluxes,'Flrl_rofl_HDO')
       call seq_flds_add(x2r_fluxes,'Flrl_rofl_HDO')
       longname = 'HDO Water flux from land (liquid)'
       stdname  = 'HDO_liquid_water_flux_into_runoff'
       attname  = 'Flrl_rofl_HDO'
       call metadata_set(attname, longname, stdname, units)

       ! r2x, x2o
       call seq_flds_add(r2x_fluxes,'Forr_rofl_16O')
       call seq_flds_add(x2o_fluxes,'Foxx_rofl_16O')
       longname = 'H2_16O Water flux due to liq runoff '
       stdname  = 'H2_16O_water_flux_into_sea_water'
       attname  = 'Forr_rofl_16O'
       call metadata_set(attname, longname, stdname, units)
       attname  = 'Foxx_rofl_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Forr_rofl_18O')
       call seq_flds_add(x2o_fluxes,'Foxx_rofl_18O')
       longname = 'H2_18O Water flux due to liq runoff '
       stdname  = 'H2_18O_water_flux_into_sea_water'
       attname  = 'Forr_rofl_18O'
       call metadata_set(attname, longname, stdname, units)
       attname  = 'Foxx_rofl_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Forr_rofl_HDO')
       call seq_flds_add(x2o_fluxes,'Foxx_rofl_HDO')
       longname = 'HDO Water flux due to liq runoff '
       stdname  = 'HDO_water_flux_into_sea_water'
       attname  = 'Forr_rofl_HDO'
       call metadata_set(attname, longname, stdname, units)
       attname  = 'Foxx_rofl_HDO'
       call metadata_set(attname, longname, stdname, units)

       call seq_flds_add(r2x_fluxes,'Forr_rofi_16O')
       call seq_flds_add(x2o_fluxes,'Foxx_rofi_16O')
       longname = 'H2_16O Water flux due to ice runoff '
       stdname  = 'H2_16O_water_flux_into_sea_water'
       attname  = 'Forr_rofi_16O'
       call metadata_set(attname, longname, stdname, units)
       attname  = 'Foxx_rofi_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Forr_rofi_18O')
       call seq_flds_add(x2o_fluxes,'Foxx_rofi_18O')
       longname = 'H2_18O Water flux due to ice runoff '
       stdname  = 'H2_18O_water_flux_into_sea_water'
       attname  = 'Forr_rofi_18O'
       call metadata_set(attname, longname, stdname, units)
       attname  = 'Foxx_rofi_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Forr_rofi_HDO')
       call seq_flds_add(x2o_fluxes,'Foxx_rofi_HDO')
       longname = 'HDO Water flux due to ice runoff '
       stdname  = 'HDO_water_flux_into_sea_water'
       attname  = 'Forr_rofi_HDO'
       call metadata_set(attname, longname, stdname, units)

       ! r2x, x2l
       call seq_flds_add(r2x_fluxes,'Flrr_flood_16O')
       call seq_flds_add(x2l_fluxes,'Flrr_flood_16O')
       longname = 'H2_16O waterrflux due to flooding'
       stdname  = 'H2_16O_flodding_water_flux_back_to_land'
       attname  = 'Flrr_flood_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Flrr_flood_18O')
       call seq_flds_add(x2l_fluxes,'Flrr_flood_18O')
       longname = 'H2_18O waterrflux due to flooding'
       stdname  = 'H2_18O_flodding_water_flux_back_to_land'
       attname  = 'Flrr_flood_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_fluxes,'Flrr_flood_HDO')
       call seq_flds_add(x2l_fluxes,'Flrr_flood_HDO')
       longname = 'HDO Waterrflux due to flooding'
       stdname  = 'HDO_flodding_water_flux_back_to_land'
       attname  = 'Flrr_flood_HDO'
       call metadata_set(attname, longname, stdname, units)

       units    = 'm3'
       call seq_flds_add(r2x_states,'Flrr_volr_16O')
       call seq_flds_add(x2l_states,'Flrr_volr_16O')
       longname = 'H2_16O river channel water volume '
       stdname  = 'H2_16O_rtm_volr'
       attname  = 'Flrr_volr_16O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_states,'Flrr_volr_18O')
       call seq_flds_add(x2l_states,'Flrr_volr_18O')
       longname = 'H2_18O river channel water volume '
       stdname  = 'H2_18O_rtm_volr'
       attname  = 'Flrr_volr_18O'
       call metadata_set(attname, longname, stdname, units)
       call seq_flds_add(r2x_states,'Flrr_volr_HDO')
       call seq_flds_add(x2l_states,'Flrr_volr_HDO')
       longname = 'HDO river channel water volume '
       stdname  = 'HDO_rtm_volr'
       attname  = 'Flrr_volr_HDO'
       call metadata_set(attname, longname, stdname, units)

       ! call seq_flds_add(r2x_fluxes,'Flrr_flood_HDO')
       ! call seq_flds_add(x2l_fluxes,'Flrr_flood_HDO')
       ! longname = 'H2_18O Waterrflux due to flooding'
       ! stdname  = 'H2_18O_flodding_water_flux_back_to_land'
       ! attname  = 'Flrr_flood_18O'
       ! call metadata_set(attname, longname, stdname, units)
 
       !-----------------------------

     endif !Water isotopes

     !-----------------------------------------------------------------------------
     ! optional per thickness category fields
     !-----------------------------------------------------------------------------

     if (seq_flds_i2o_per_cat) then
        do num = 1, ice_ncat
           write(cnum,'(i2.2)') num

           ! Fractional ice coverage wrt ocean

           name = 'Si_ifrac_' // cnum
           call seq_flds_add(i2x_states,name)
           call seq_flds_add(x2o_states,name)
           longname = 'fractional ice coverage wrt ocean for thickness category ' // cnum
           stdname  = 'sea_ice_area_fraction'
           units    = '1'
           attname  = name
           call metadata_set(attname, longname, stdname, units)

           ! Net shortwave radiation

           name = 'PFioi_swpen_ifrac_' // cnum
           call seq_flds_add(i2x_fluxes,name)
           call seq_flds_add(x2o_fluxes,name)
           longname = 'net shortwave radiation penetrating into ice and ocean times ice fraction for thickness category ' // cnum
           stdname  = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_sea_ice_area_fraction'
           units    = 'W m-2'
           attname  = name
           call metadata_set(attname, longname, stdname, units)

        end do

        ! Fractional atmosphere coverage wrt ocean

        name = 'Sf_afrac'
        call seq_flds_add(x2o_states,name)
        longname = 'fractional atmosphere coverage wrt ocean'
        stdname  = 'atmosphere_area_fraction'
        units    = '1'
        attname  = name
        call metadata_set(attname, longname, stdname, units)

        name = 'Sf_afracr'
        call seq_flds_add(x2o_states,name)
        longname = 'fractional atmosphere coverage used in radiation computations wrt ocean'
        stdname  = 'atmosphere_area_fraction'
        units    = '1'
        attname  = name
        call metadata_set(attname, longname, stdname, units)

        ! Net shortwave radiation

        name = 'Foxx_swnet_afracr'
        call seq_flds_add(x2o_fluxes,name)
        longname = 'net shortwave radiation times atmosphere fraction'
        stdname = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_atmosphere_area_fraction'
        units = 'W m-2'
        attname = name
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
        longname = 'Volumetric soil water'
        stdname  = 'soil_water'
        units    = 'm3/m3'
        call metadata_set(carma_fields, longname, stdname, units)
     endif

     !-----------------------------------------------------------------------------
     ! Read namelist for MEGAN
     ! if MEGAN emission are specified then setup fields for CLM to CAM communication
     ! (emissions fluxes)
     !-----------------------------------------------------------------------------

     call shr_megan_readnl(nlfilename='drv_flds_in', ID=ID, megan_fields=megan_voc_fields)
     if (shr_megan_mechcomps_n>0) then
        call seq_flds_add(l2x_fluxes, trim(megan_voc_fields))
        call seq_flds_add(x2a_fluxes, trim(megan_voc_fields))
        longname = 'MEGAN emission fluxes'
        stdname  = 'megan_fluxes'
        units    = 'molecules/m2/sec'
        call metadata_set(megan_voc_fields, longname, stdname, units)
     endif

     !-----------------------------------------------------------------------------
     ! Read namelist for Fire Emissions
     ! if fire emission are specified then setup fields for CLM to CAM communication
     ! (emissions fluxes)
     !-----------------------------------------------------------------------------

     call shr_fire_emis_readnl(nlfilename='drv_flds_in', ID=ID, emis_fields=fire_emis_fields)
     if (shr_fire_emis_mechcomps_n>0) then
        call seq_flds_add(l2x_fluxes, trim(fire_emis_fields))
        call seq_flds_add(x2a_fluxes, trim(fire_emis_fields))
        longname = 'wild fire emission fluxes'
        stdname  = 'fire_emis'
        units    = 'kg/m2/sec'
        call metadata_set(fire_emis_fields, longname, stdname, units)

        call seq_flds_add(l2x_states, trim(shr_fire_emis_ztop_token))
        call seq_flds_add(x2a_states, trim(shr_fire_emis_ztop_token))
        longname = 'wild fire plume height'
        stdname  = 'fire_plume_top'
        units    = 'm'

        call metadata_set(shr_fire_emis_ztop_token, longname, stdname, units)

     endif

     !-----------------------------------------------------------------------------
     ! Dry Deposition fields
     ! First read namelist and figure out the drydep field list to pass
     ! Then check if file exists and if not, n_drydep will be zero
     ! Then add dry deposition fields to land export and atmosphere import states
     ! Then initialize dry deposition fields
     ! Note: CAM and CLM will then call seq_drydep_setHCoeff
     !-----------------------------------------------------------------------------

     call seq_drydep_readnl(nlfilename="drv_flds_in", ID=ID, seq_drydep_fields=seq_drydep_fields)
     if ( lnd_drydep ) then
        call seq_flds_add(l2x_states, seq_drydep_fields)
        call seq_flds_add(x2a_states, seq_drydep_fields)

        longname = 'dry deposition velocity'
        stdname  = 'drydep_vel'
        units    = 'cm/sec'
        call metadata_set(seq_drydep_fields, longname, stdname, units)

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
     seq_flds_l2x_states_to_glc = trim(l2x_states_to_glc)
     seq_flds_x2l_states = trim(x2l_states)
     seq_flds_x2l_states_from_glc = trim(x2l_states_from_glc)
     seq_flds_o2x_states = trim(o2x_states)
     seq_flds_x2o_states = trim(x2o_states)
     seq_flds_g2x_states = trim(g2x_states)
     seq_flds_g2x_states_to_lnd = trim(g2x_states_to_lnd)
     seq_flds_x2g_states = trim(x2g_states)
     seq_flds_xao_states = trim(xao_states)
     seq_flds_xao_albedo = trim(xao_albedo)
     seq_flds_xao_diurnl = trim(xao_diurnl)
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
     seq_flds_l2x_fluxes_to_glc = trim(l2x_fluxes_to_glc)
     seq_flds_x2l_fluxes = trim(x2l_fluxes)
     seq_flds_x2l_fluxes_from_glc = trim(x2l_fluxes_from_glc)
     seq_flds_o2x_fluxes = trim(o2x_fluxes)
     seq_flds_x2o_fluxes = trim(x2o_fluxes)
     seq_flds_g2x_fluxes = trim(g2x_fluxes)
     seq_flds_g2x_fluxes_to_lnd = trim(g2x_fluxes_to_lnd)
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
        write(logunit,"(A)") subname//': seq_flds_xao_diurnl= ',trim(seq_flds_xao_diurnl)
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
     call catFields(seq_flds_l2x_fields_to_glc, seq_flds_l2x_states_to_glc, seq_flds_l2x_fluxes_to_glc)
     call catFields(seq_flds_x2l_fields, seq_flds_x2l_states, seq_flds_x2l_fluxes)
     call catFields(seq_flds_x2l_fields_from_glc, seq_flds_x2l_states_from_glc, seq_flds_x2l_fluxes_from_glc)
     call catFields(seq_flds_o2x_fields, seq_flds_o2x_states, seq_flds_o2x_fluxes)
     call catFields(seq_flds_x2o_fields, seq_flds_x2o_states, seq_flds_x2o_fluxes)
     call catFields(seq_flds_g2x_fields, seq_flds_g2x_states, seq_flds_g2x_fluxes)
     call catFields(seq_flds_g2x_fields_to_lnd, seq_flds_g2x_states_to_lnd, seq_flds_g2x_fluxes_to_lnd)
     call catFields(seq_flds_x2g_fields, seq_flds_x2g_states, seq_flds_x2g_fluxes)
     call catFields(seq_flds_xao_fields, seq_flds_xao_albedo, seq_flds_xao_states)
     call catFields(stringtmp          , seq_flds_xao_fields, seq_flds_xao_fluxes)
     call catFields(seq_flds_xao_fields, stringtmp          , seq_flds_xao_diurnl)
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
! If the attname passed in contains colons it is assumed to be a list of fields
! all of which have the same names and units
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
     integer :: i, j

     i = index(attname,':')
     j=1

     do while(i>j .and. i<=len_trim(attname))
        n_entries = n_entries + 1
        lookup_entry(n_entries,1) = attname(j:i-1)
        lookup_entry(n_entries,2) = trim(longname)
        lookup_entry(n_entries,3) = trim(stdname )
        lookup_entry(n_entries,4) = trim(units   )
        j=i+1
        i =  index(attname(j:),':') + j - 1
     enddo
     n_entries = n_entries + 1
     i = len_trim(attname)
     lookup_entry(n_entries,1) = attname(j:i)
     lookup_entry(n_entries,2) = trim(longname)
     lookup_entry(n_entries,3) = trim(stdname )
     lookup_entry(n_entries,4) = trim(units   )




     if (n_entries .ge. nmax) then
        write(logunit,*)'n_entries= ',n_entries,' nmax = ',nmax,' attname= ',trim(attname)
        call shr_sys_abort(subname//'ERROR: nmax fields in lookup_entry table exceeded')
     end if

   end subroutine metadata_set

   !===============================================================================

   subroutine set_glc_elevclass_field(name, attname, longname, stdname, units, fieldlist, &
                                      additional_list)

     ! Sets a coupling field for all glc elevation classes (1:glc_nec) plus bare land
     ! (index 0).
     !
     ! Note that, if glc_nec = 0, then we don't create any coupling fields (not even the
     ! bare land (0) index)
     !
     ! Puts the coupling fields in the given fieldlist, and also does the appropriate
     ! metadata_set calls.
     !
     ! additional_list should be .false. (or absent) the first time this is called for a
     ! given set of coupling fields. However, if this same set of coupling fields is being
     ! added to multiple field lists, then additional_list should be set to true for the
     ! second and subsequent calls; in this case, the metadata_set calls are not done
     ! (because they have already been done).
     !
     ! name, attname and longname give the base name of the field; the elevation class
     ! index will be appended as a suffix

     ! !USES:
     use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_elevclass_as_string

     ! !INPUT/OUTPUT PARAMETERS:
     character(len=*), intent(in) :: name     ! base field name to add to fieldlist
     character(len=*), intent(in) :: attname  ! base field name for metadata
     character(len=*), intent(in) :: longname ! base long name for metadata
     character(len=*), intent(in) :: stdname  ! standard name for metadata
     character(len=*), intent(in) :: units    ! units for metadata
     character(len=*), intent(inout) :: fieldlist  ! field list into which the fields should be added

     logical, intent(in), optional :: additional_list  ! whether this is an additional list for the same set of coupling fields (see above for details; defaults to false)

     !EOP
     integer            :: num
     character(len= 16) :: cnum
     logical :: l_additional_list  ! local version of the optional additional_list argument

     l_additional_list = .false.
     if (present(additional_list)) then
        l_additional_list = additional_list
     end if

     if (glc_get_num_elevation_classes() > 0) then
        do num = 0, glc_get_num_elevation_classes()
           cnum = glc_elevclass_as_string(num)

           call seq_flds_add(fieldlist, trim(name) // trim(cnum))

           if (.not. l_additional_list) then
              call metadata_set(attname  = trim(attname) // trim(cnum), &
                   longname = trim(longname) // ' of elevation class ' // trim(cnum), &
                   stdname  = stdname, &
                   units    = units)
           end if
        end do
     end if
   end subroutine set_glc_elevclass_field

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
     character(len=CSS) :: llongname, lstdname, lunits, lshortname  ! local copies
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


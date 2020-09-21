module glc2lndMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from glc to clm.
  ! For now glc datais send and received on the lnd decomposition and grid.
  !
  ! The fields sent from the lnd component to the glc component via
  !  the coupler are labeled 's2x', or sno to coupler.
  ! The fields received by the lnd component from the glc component
  !  via the coupler are labeled 'x2s', or coupler to sno.
  ! 'Sno' is a misnomer in that the exchanged data are related to
  !  the ice beneath the snow, not the snow itself.  But by CESM convention,
  ! 'ice' refers to sea ice, not land ice.
  !
  ! !USES:
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar     , only : maxpatch_glcmec
  use elm_varctl     , only : iulog, glc_smb
  use abortutils     , only : endrun
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp 
  !
  ! !REVISION HISTORY:
  ! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! glc -> land variables structure
  type, public :: glc2lnd_type

     real(r8), pointer :: frac_grc    (:,:) => null()
     real(r8), pointer :: topo_grc    (:,:) => null()
     real(r8), pointer :: hflx_grc    (:,:) => null()

     ! Total ice sheet grid coverage mask (0-1)
     ! (true ice sheet mask, received from glc, in contrast to glcmask, which is just a
     ! guess available at initialization)
     real(r8), pointer :: icemask_grc (:)   => null()

     ! icemask_coupled_fluxes_grc is like icemask_grc, but the mask only contains icesheet
     ! points that potentially send non-zero fluxes to the coupler. i.e., it does not
     ! contain icesheets that are diagnostic only, because for those diagnostic ice sheets
     ! (which do not send calving fluxes to the coupler), we need to use the non-dynamic
     ! form of runoff routing in CLM in order to conserve water properly.
     !
     ! (However, note that this measure of "diagnostic-only" does not necessarily
     ! correspond to whether CLM is updating its glacier areas there - for example, we
     ! could theoretically have an icesheet whose areas are evolving, and CLM is updating
     ! its glacier areas to match, but where we're zeroing out the fluxes sent to the
     ! coupler, and so we're using the non-dynamic form of runoff routing in CLM.)
     real(r8), pointer :: icemask_coupled_fluxes_grc (:)  => null()

     ! Where we should do runof routing that is appropriate for having a dynamic icesheet underneath.
     logical , pointer :: glc_dyn_runoff_routing_grc (:) => null()

   contains

     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: update_glc2lnd

     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: check_glc2lnd_icemask  ! sanity-check icemask from GLC
     procedure, private :: check_glc2lnd_icemask_coupled_fluxes ! sanity-check icemask_coupled_fluxes from GLC
     procedure, private :: update_glc2lnd_dyn_runoff_routing ! update glc_dyn_runoff_routing field based on input from GLC
     procedure, private :: update_glc2lnd_fracs   ! update subgrid fractions based on input from GLC
     procedure, private :: update_glc2lnd_topo    ! update column-level topographic heights based on input from GLC

  end type glc2lnd_type

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize glc variables required by the land
    !
    ! !ARGUMENTS:
    class (glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg
    !------------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%frac_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%frac_grc    (:,:) = nan
    allocate(this%topo_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%topo_grc    (:,:) = nan
    allocate(this%hflx_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%hflx_grc    (:,:) = nan
    allocate(this%icemask_grc (begg:endg))                   ;   this%icemask_grc (:)   = nan
    allocate(this%icemask_coupled_fluxes_grc (begg:endg))    ;   this%icemask_coupled_fluxes_grc (:)   = nan
    allocate(this%glc_dyn_runoff_routing_grc (begg:endg))    ;   this%glc_dyn_runoff_routing_grc (:)   = .false.

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    use elm_varcon , only : spval
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg
    
    if (maxpatch_glcmec > 0) then
       this%icemask_grc(begg:endg) = spval
       call hist_addfld1d (fname='ICE_MASK',  units='unitless',  &
            avgflag='I', long_name='Ice sheet mask coverage', &
            ptr_gcell=this%icemask_grc)
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use domainMod      , only : ldomain
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg
    
    this%frac_grc(begg:endg, :) = 0.0_r8
    this%topo_grc(begg:endg, :) = 0.0_r8
    this%hflx_grc(begg:endg, :) = 0.0_r8

    ! glcmask (from a file) provides a rough guess of the icemask (from CISM); thus, in
    ! initialization, set icemask equal to glcmask; icemask will later get updated at
    ! the start of the run loop, as soon as we have data from CISM
    this%icemask_grc(begg:endg) = ldomain%glcmask(begg:endg)

    ! initialize icemask_coupled_fluxes to 0; this seems safest in case we aren't coupled
    ! to CISM (to ensure that we use the uncoupled form of runoff routing)
    this%icemask_coupled_fluxes_grc(begg:endg) = 0.0_r8
    
    call this%update_glc2lnd_dyn_runoff_routing(bounds)

  end subroutine InitCold


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write glc2lnd information to/from restart file.
    !
    ! !USES:
    use ncdio_pio , only : ncd_double, file_desc_t
    use decompMod , only : bounds_type
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) , intent(inout) :: this
    type(bounds_type)   , intent(in)    :: bounds 
    type(file_desc_t)   , intent(inout) :: ncid ! netcdf id
    character(len=*)    , intent(in)    :: flag ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    
    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='icemask', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='total ice-sheet grid coverage mask', units='fraction', &
         interpinic_flag='skip', readvar=readvar, data=this%icemask_grc)

  end subroutine Restart


  !-----------------------------------------------------------------------
  subroutine update_glc2lnd(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update values to derived-type CLM variables based on input from GLC (via the coupler)
    !
    ! icemask, icemask_coupled_fluxes, glc_dyn_runoff_routing, and topo are always updated
    ! (although note that this routine should only be called when
    ! create_glacier_mec_landunit is true, or some similar condition; this should be
    ! controlled in a conditional around the call to this routine); fracs are updated if
    ! glc_do_dynglacier is true
    !
    ! !USES:
    use elm_varctl , only : glc_do_dynglacier
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in)    :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'update_glc2lnd'
    !-----------------------------------------------------------------------

    ! Note that nothing is needed to update icemask or icemask_coupled_fluxes here,
    ! because these values have already been set in lnd_import_export. However, we do
    ! some sanity-checking of those fields here.
    call this%check_glc2lnd_icemask(bounds)
    call this%check_glc2lnd_icemask_coupled_fluxes(bounds)

    call this%update_glc2lnd_dyn_runoff_routing(bounds)

    if (glc_do_dynglacier) then
       call this%update_glc2lnd_fracs(bounds)
    end if

    call this%update_glc2lnd_topo(bounds)

  end subroutine update_glc2lnd

  !-----------------------------------------------------------------------
  subroutine check_glc2lnd_icemask(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do a sanity check on the icemask received from CISM via coupler.
    !
    ! !USES:
    use domainMod , only : ldomain
    use elm_varcon, only : nameg
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'check_glc2lnd_icemask'
    !-----------------------------------------------------------------------
    
    do g = bounds%begg, bounds%endg

       ! Ensure that icemask is a subset of glcmask. This is needed because we allocated
       ! memory based on glcmask, so it is a problem if the ice sheet tries to expand
       ! beyond the area defined by glcmask.

       if (this%icemask_grc(g) > 0._r8 .and. ldomain%glcmask(g) == 0._r8) then
          write(iulog,*) subname//' ERROR: icemask must be a subset of glcmask.'
          write(iulog,*) 'You can fix this problem by adding more grid cells'
          write(iulog,*) 'to the mask defined by the fglcmask file.'
          write(iulog,*) '(Change grid cells to 1 everywhere that CISM can operate.)'
          call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine check_glc2lnd_icemask

  !-----------------------------------------------------------------------
  subroutine check_glc2lnd_icemask_coupled_fluxes(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do a sanity check on the icemask_coupled_fluxes field received from CISM via coupler.
    !
    ! !USES:
    use elm_varcon, only : nameg
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'check_glc2lnd_icemask_coupled_fluxes'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg

       ! Ensure that icemask_coupled_fluxes is a subset of icemask. Although there
       ! currently is no code in CLM that depends on this relationship, it seems helpful
       ! to ensure that this intuitive relationship holds, so that code developed in the
       ! future can rely on it.

       if (this%icemask_coupled_fluxes_grc(g) > 0._r8 .and. this%icemask_grc(g) == 0._r8) then
          write(iulog,*) subname//' ERROR: icemask_coupled_fluxes must be a subset of icemask.'
          call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine check_glc2lnd_icemask_coupled_fluxes

  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_dyn_runoff_routing(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update glc_dyn_runoff_routing field based on updated icemask_coupled_fluxes field
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'update_glc2lnd_dyn_runoff_routing'
    !-----------------------------------------------------------------------

    ! Wherever we have an icesheet that is computing and sending fluxes to the coupler -
    ! which particularly means it is computing a calving flux - we will use the
    ! "glc_dyn_runoff_routing" scheme. In other places - including places where CISM is
    ! not running at all, as well as places where CISM is running in diagnostic-only mode
    ! and therefore is not sending a calving flux - we use the alternative
    ! glc_dyn_runoff_routing=false scheme. This is needed to conserve water correctly in
    ! the absence of a calving flux.

    do g = bounds%begg, bounds%endg
       if (this%icemask_coupled_fluxes_grc(g) > 0._r8) then
          this%glc_dyn_runoff_routing_grc(g) = .true.
       else
          this%glc_dyn_runoff_routing_grc(g) = .false.
       end if
    end do

  end subroutine update_glc2lnd_dyn_runoff_routing



  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_fracs(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update subgrid fractions based on input from GLC (via the coupler)
    !
    ! The weights updated here are some col_pp%wtlunit and lun_pp%wtgcell values
    !
    ! !USES:
    use elm_varcon        , only : ispval
    use landunit_varcon   , only : istice_mec
    use column_varcon     , only : col_itype_to_icemec_class
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,t,c                              ! indices
    real(r8):: area_ice_mec                     ! area of the ice_mec landunit
    integer :: l_ice_mec                        ! index of the ice_mec landunit
    integer :: icemec_class                     ! current icemec class (1..maxpatch_glcmec)
    logical :: frac_assigned(1:maxpatch_glcmec) ! whether this%frac has been assigned for each elevation class
    logical :: error                            ! if an error was found
    
    character(len=*), parameter :: subname = 'update_glc2lnd_fracs'
    !-----------------------------------------------------------------------
    
    do g = bounds%begg, bounds%endg
       ! Values from GLC are only valid within the icemask, so we only update CLM's areas there
       if (this%icemask_grc(g) > 0._r8) then
          do t = grc_pp%topi(g), grc_pp%topf(g)

             ! For now, all topounits on the gridcell will have identical settings for 
             ! the fractional area of the glacier landunit. This makes sense for max_topounits = 1,
             ! but should not break for max_topounits > 1. Revisit this code during the update for glacier_mec.
             ! Set total icemec landunit area
             area_ice_mec = sum(this%frac_grc(g, 1:maxpatch_glcmec))
             call set_landunit_weight(t, istice_mec, area_ice_mec)

             ! If new landunit area is greater than 0, then update column areas
             ! (If new landunit area is 0, col_pp%wtlunit is arbitrary, so we might as well keep the existing values)
             if (area_ice_mec > 0) then
                ! Determine index of the glc_mec landunit
                l_ice_mec = top_pp%landunit_indices(istice_mec, t)
                if (l_ice_mec == ispval) then
                   write(iulog,*) subname//' ERROR: no ice_mec landunit found within the icemask, for g = ', g
                   call endrun()
                end if
             
                frac_assigned(1:maxpatch_glcmec) = .false.
                do c = lun_pp%coli(l_ice_mec), lun_pp%colf(l_ice_mec)
                   icemec_class = col_itype_to_icemec_class(col_pp%itype(c))
                   col_pp%wtlunit(c) = this%frac_grc(g, icemec_class) / lun_pp%wttopounit(l_ice_mec)
                   frac_assigned(icemec_class) = .true.
                end do

                ! Confirm that all elevation classes that have non-zero area according to
                ! this%frac have been assigned to a column in CLM's data structures
                error = .false.
                do icemec_class = 1, maxpatch_glcmec
                   if (this%frac_grc(g, icemec_class) > 0._r8 .and. &
                        .not. frac_assigned(icemec_class)) then
                      error = .true.
                   end if
                end do
                if (error) then
                   write(iulog,*) subname//' ERROR: at least one glc_mec column has non-zero area from the coupler,'
                   write(iulog,*) 'but there was no slot in memory for this column; g = ', g
                   write(iulog,*) 'this%frac_grc(g, 1:maxpatch_glcmec) = ', &
                        this%frac_grc(g, 1:maxpatch_glcmec)
                   write(iulog,*) 'frac_assigned(1:maxpatch_glcmec) = ', &
                        frac_assigned(1:maxpatch_glcmec)
                   call endrun()
                end if  ! error
             end if  ! area_ice_mec > 0
          end do ! loop over topounits on this gridcell
       end if  ! this%icemask_grc(g) > 0
    end do  ! g

  end subroutine update_glc2lnd_fracs

  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_topo(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update column-level topographic heights based on input from GLC (via the coupler)
    !
    ! !USES:
    use landunit_varcon , only : istice_mec
    use column_varcon   , only : col_itype_to_icemec_class
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g      ! indices
    integer :: icemec_class ! current icemec class (1..maxpatch_glcmec)

    character(len=*), parameter :: subname = 'update_glc2lnd_topo'
    !-----------------------------------------------------------------------
    
    ! It is tempting to use the do_smb_c filter here, since we only need glc_topo inside
    ! this filter. But the problem with using the filter is that this routine is called
    ! before the filters are updated to reflect the updated weights. As long as
    ! glacier_mec, natural veg and any other landunit within the smb filter are always
    ! active, regardless of their weights, this isn't a problem. But we don't want to
    ! build in assumptions that those rules will be in place regarding active flags.
    ! Other ways around this problem would be:
    ! (1) Use the inactive_and_active filter  - but we're trying to avoid use of that
    !     filter if possible, because it can be confusing
    ! (2) Call this topo update routine later in the driver loop, after filters have been
    !     updated  - but that leads to greater complexity in the driver loop.
    ! So it seems simplest just to take the minor performance hit of setting glc_topo
    ! over all columns, even those outside the do_smb_c filter.
    
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       g = col_pp%gridcell(c)

       ! Values from GLC are only valid within the icemask, so we only update CLM's topo values there
       if (this%icemask_grc(g) > 0._r8) then
          if (lun_pp%itype(l) == istice_mec) then
             icemec_class = col_itype_to_icemec_class(col_pp%itype(c))
          else
             ! If not on a glaciated column, assign topography to the bare-land value determined by GLC.
             icemec_class = 0
          end if

          col_pp%glc_topo(c) = this%topo_grc(g, icemec_class)
       end if
    end do

  end subroutine update_glc2lnd_topo

end module glc2lndMod


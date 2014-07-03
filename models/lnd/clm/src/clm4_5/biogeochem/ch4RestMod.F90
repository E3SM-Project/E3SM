module ch4RestMod
#ifdef LCH4

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ch4RestMod
!
! !DESCRIPTION:
! Reads from or writes restart data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4Rest
!
! !REVISION HISTORY:
! 2009, August: Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4Rest
!
! !INTERFACE:
  subroutine ch4Rest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio_pio
    use decompMod     , only : get_proc_bounds
    use clm_varctl    , only : nsrest
    use clm_time_manager , only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
!
! !CALLED FROM:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! column ch4 state variable - conc_ch4_sat

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CONC_CH4_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methane soil concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CONC_CH4_SAT', data=cptr%cch4%conc_ch4_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - conc_ch4_unsat
    
    if (flag == 'define') then 
       call ncd_defvar(ncid=ncid, varname='CONC_CH4_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methane soil concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CONC_CH4_UNSAT', data=cptr%cch4%conc_ch4_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - conc_o2_sat
    
    if (flag == 'define') then 
       call ncd_defvar(ncid=ncid, varname='CONC_O2_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='oxygen soil concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CONC_O2_SAT', data=cptr%cch4%conc_o2_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - conc_o2_unsat

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CONC_O2_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='oxygen soil concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CONC_O2_UNSAT', data=cptr%cch4%conc_o2_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 flux variable - o2stress_sat (used in CNDecompCascade)

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='O2STRESS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='oxygen stress fraction', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='O2STRESS_SAT', data=cptr%cch4%o2stress_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 flux variable - o2stress_unsat (used in CNDecompCascade)

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='O2STRESS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='oxygen stress fraction', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='O2STRESS_UNSAT', data=cptr%cch4%o2stress_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - layer_sat_lag

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LAYER_SAT_LAG', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='lagged saturation status of layer in unsat. zone', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='LAYER_SAT_LAG', data=cptr%cch4%layer_sat_lag, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - qflx_surf_lag
 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='QFLX_SURF_LAG', xtype=ncd_double, &
            dim1name='column', long_name='time-lagged surface runoff', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='QFLX_SURF_LAG', data=cptr%cch4%qflx_surf_lag, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then 
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 state variable - finundated_lag

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FINUNDATED_LAG', xtype=ncd_double, &
            dim1name='column', long_name='time-lagged inundated fraction', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='FINUNDATED_LAG', data=cptr%cch4%finundated_lag, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then 
          if (is_restart()) call endrun()
       end if
    end if


    ! column ch4 state variable - fsat_bef
    ! fsat_bef = finundated except inside methane code
    ! Only necessary for bit-for-bit restarts

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FINUNDATED', xtype=ncd_double, &
            dim1name='column', long_name='inundated fraction', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='FINUNDATED', data=cptr%cch4%fsat_bef, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then 
          if (is_restart()) call endrun()
       end if
    end if

#ifdef CN

    ! column ch4 state variable - annavg_somhr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_somhr', xtype=ncd_double,  &
            dim1name='column',long_name='Annual Average SOMHR',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_somhr', data=cptr%cch4%annavg_somhr, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! column ch4 state variable - annavg_finrw
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_finrw', xtype=ncd_double,  &
            dim1name='column',long_name='Annual Average Respiration-Weighted FINUNDATED',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_finrw', data=cptr%cch4%annavg_finrw, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! column ch4 state variable - annsum_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_counter_ch4', xtype=ncd_double,  &
            dim1name='column',long_name='CH4 Ann. Sum Time Counter',units='s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_counter_ch4', data=cptr%cch4%annsum_counter, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! column ch4 state variable - tempavg_somhr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_somhr', xtype=ncd_double,  &
            dim1name='column',long_name='Temp. Average SOMHR',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_somhr', data=cptr%cch4%tempavg_somhr, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! column ch4 state variable - tempavg_finrwi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_finrw', xtype=ncd_double,  &
            dim1name='column',long_name='Temp. Average Respiration-Weighted FINUNDATED',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_finrw', data=cptr%cch4%tempavg_finrw, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - tempavg_agnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Temp. Average AGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_agnpp', data=pptr%pcf%tempavg_agnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - tempavg_bgnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Temp. Average BGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_bgnpp', data=pptr%pcf%tempavg_bgnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - annavg_agnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Ann. Average AGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_agnpp', data=pptr%pcf%annavg_agnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - annavg_bgnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Ann. Average BGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_bgnpp', data=pptr%pcf%annavg_bgnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

   ! column ch4 flux variable - o2_decomp_depth_sat (used in CNNitrifDenitrif)

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='O2_DECOMP_DEPTH_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='O2 consumption during decomposition', units='mol/m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='O2_DECOMP_DEPTH_SAT', data=cptr%cch4%o2_decomp_depth_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column ch4 flux variable - o2_decomp_depth_unsat (used in CNNitrifDenitrif)

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='O2_DECOMP_DEPTH_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='O2 consumption during decomposition', units='mol/m3/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='O2_DECOMP_DEPTH_UNSAT', data=cptr%cch4%o2_decomp_depth_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

#endif

    ! column ch4 state variable - lake_soilc

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LAKE_SOILC', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true.,&
            long_name='lake soil carbon concentration', units='g/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='LAKE_SOILC', data=cptr%cch4%lake_soilc, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then 
          if (is_restart()) call endrun()
       end if
    end if

  end subroutine ch4Rest

#endif

end module ch4RestMod

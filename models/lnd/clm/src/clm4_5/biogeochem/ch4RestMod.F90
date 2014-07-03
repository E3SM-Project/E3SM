module ch4RestMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes restart data
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  ! save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4Rest
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ch4Rest( bounds, ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/Write biogeophysics information to/from restart file.
    !
    ! !USES:
    use clmtype
    use clm_varctl, only : use_cn
    use ncdio_pio,  only : ncd_double 
    use pio,        only : file_desc_t
    use decompMod,  only : bounds_type
    use restUtilMod
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*),  intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    ! column ch4 state variable - conc_ch4_sat
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='methane soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=cch4%conc_ch4_sat)

    ! column ch4 state variable - conc_ch4_unsat
    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='methane soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=cch4%conc_ch4_unsat)

    ! column ch4 state variable - conc_o2_sat
    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=cch4%conc_o2_sat)

    ! column ch4 state variable - conc_o2_unsat
    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=cch4%conc_o2_unsat)

    ! column ch4 flux variable - o2stress_sat (used in CNDecompCascade)
    call restartvar(ncid=ncid, flag=flag, varname='O2STRESS_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen stress fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=cch4%o2stress_sat)

    ! column ch4 flux variable - o2stress_unsat (used in CNDecompCascade)
    call restartvar(ncid=ncid, flag=flag, varname='O2STRESS_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen stress fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=cch4%o2stress_unsat)

    ! column ch4 state variable - layer_sat_lag
    call restartvar(ncid=ncid, flag=flag, varname='LAYER_SAT_LAG', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='lagged saturation status of layer in unsat. zone', units='', &
         readvar=readvar, interpinic_flag='interp', data=cch4%layer_sat_lag)

    ! column ch4 state variable - qflx_surf_lag
    call restartvar(ncid=ncid, flag=flag, varname='QFLX_SURF_LAG', xtype=ncd_double, &
         dim1name='column', &
         long_name='time-lagged surface runoff', units='mm/s', &
         readvar=readvar, interpinic_flag='interp', data=cch4%qflx_surf_lag)

    ! column ch4 state variable - finundated_lag
    call restartvar(ncid=ncid, flag=flag, varname='FINUNDATED_LAG', xtype=ncd_double, &
         dim1name='column', &
         long_name='time-lagged inundated fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=cch4%finundated_lag)

    ! column ch4 state variable - fsat_bef  (finundated except inside methane code)
    ! only necessary for bit-for-bit restarts
    call restartvar(ncid=ncid, flag=flag, varname='FINUNDATED', xtype=ncd_double, &
            dim1name='column', &
            long_name='inundated fraction', units='', &
            readvar=readvar, interpinic_flag='interp', data=cch4%fsat_bef)
    ! Set finundated to fsat_bef when reading in restart or finidat
    cws%finundated(bounds%begc:bounds%endc) = cch4%fsat_bef(bounds%begc:bounds%endc)

    if (use_cn) then

       ! column ch4 state variable - annavg_somhr
       call restartvar(ncid=ncid, flag=flag, varname='annavg_somhr', xtype=ncd_double,  &
            dim1name='column',&
            long_name='Annual Average SOMHR',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=cch4%annavg_somhr)

       ! column ch4 state variable - annavg_finrw
       call restartvar(ncid=ncid, flag=flag, varname='annavg_finrw', xtype=ncd_double,  &
            dim1name='column',&
            long_name='Annual Average Respiration-Weighted FINUNDATED',units='', &
            readvar=readvar, interpinic_flag='interp', data=cch4%annavg_finrw)

       ! column ch4 state variable - annsum_counter
       call restartvar(ncid=ncid, flag=flag, varname='annsum_counter_ch4', xtype=ncd_double,  &
            dim1name='column',&
            long_name='CH4 Ann. Sum Time Counter',units='s', &
            readvar=readvar, interpinic_flag='interp', data=cch4%annsum_counter)

       ! column ch4 state variable - tempavg_somhr
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_somhr', xtype=ncd_double,  &
            dim1name='column',&
            long_name='Temp. Average SOMHR',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=cch4%tempavg_somhr)

       ! column ch4 state variable - tempavg_finrwi
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_finrw', xtype=ncd_double,  &
            dim1name='column',&
            long_name='Temp. Average Respiration-Weighted FINUNDATED',units='', &
            readvar=readvar, interpinic_flag='interp', data=cch4%tempavg_finrw)

       ! pft ch4 state variable - tempavg_agnpp
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=pcf%tempavg_agnpp)

       ! pft ch4 state variable - tempavg_bgnpp
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=pcf%tempavg_bgnpp)

       ! pft ch4 state variable - annavg_agnpp
       call restartvar(ncid=ncid, flag=flag, varname='annavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=pcf%annavg_agnpp)

       ! pft ch4 state variable - annavg_bgnpp
       call restartvar(ncid=ncid, flag=flag, varname='annavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=pcf%annavg_bgnpp)

       ! column ch4 flux variable - o2_decomp_depth_sat (used in CNNitrifDenitrif)
       call restartvar(ncid=ncid, flag=flag, varname='O2_DECOMP_DEPTH_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='O2 consumption during decomposition', units='mol/m3/s', &
            readvar=readvar, interpinic_flag='interp', data=cch4%o2_decomp_depth_sat)

       ! column ch4 flux variable - o2_decomp_depth_unsat (used in CNNitrifDenitrif)
       call restartvar(ncid=ncid, flag=flag, varname='O2_DECOMP_DEPTH_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='O2 consumption during decomposition', units='mol/m3/s', &
            readvar=readvar, interpinic_flag='interp', data=cch4%o2_decomp_depth_unsat)

    end if  ! end of use_cn if-block

    ! column ch4 state variable - lake_soilc
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_SOILC', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true.,&
         long_name='lake soil carbon concentration', units='g/m^3', &
         readvar=readvar, interpinic_flag='interp', data=cch4%lake_soilc)

  end subroutine ch4Rest

end module ch4RestMod

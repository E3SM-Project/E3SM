!
! output geos from HOMME's baroclinic instability
! this is used to generate test data for the other tool utilities
!
&ctl_nl
ne = 4
mesh_file = 'none'
test_case = 'dcmip2012_test4'
/

&vert_nl
/

&analysis_nl
tool = 'topo_gll_to_smoothed'
infilenames = ''
output_dir = "./"
output_timeunits=1
output_frequency=1
output_varnames1='PHIS'
output_type='netcdf'
io_stride = 16
/

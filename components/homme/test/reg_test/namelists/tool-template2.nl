&ctl_nl
ne = 0
mesh_file = 'mountain_10_x2.g'    
/

&vert_nl
/

&analysis_nl
tool = 'grid_template_tool'

!output_dir = "./" 
output_prefix     = "template2-"
output_timeunits=1
output_frequency=1
output_varnames1='area','corners','cv_lat','cv_lon'
output_type='netcdf'
!output_type='netcdf4p'  ! needed for ne1024
io_stride = 8
/

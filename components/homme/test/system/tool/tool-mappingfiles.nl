&ctl_nl
ne = 4
mesh_file='none'
/
&vert_nl
/
&analysis_nl
tool = 'gll_mapping_file'
infilenames = '9x16_scrip.nc'

output_dir = "movies/" 
output_timeunits=1                                                                                                
output_frequency=1                                                                                                
output_varnames1='area'
output_type='netcdf'                                                                                              
io_stride = 16                                                                                                    
/

&sw_model
   config_test_case = 5
   config_time_integration = 'RK4'
   config_dt = 172.8
   config_start_time = '0000-01-01_00:00:00'
   config_run_duration = '15_00:00:00'
   config_stats_interval = 0
   config_h_ScaleWithMesh = .false.
   config_h_mom_eddy_visc2  = 0.0
   config_h_mom_eddy_visc4  = 0.0
   config_h_tracer_eddy_diff2  = 0.0
   config_h_tracer_eddy_diff4  = 0.0
   config_thickness_adv_order = 2
   config_tracer_adv_order = 2
   config_positive_definite = .false.
   config_monotonic = .false.
   config_wind_stress = .false.
   config_bottom_drag = .false.
/
   config_stop_time  = '0000-01-16_00:00:00'

&io
   config_input_name = 'grid.nc'
   config_output_name = 'output.nc'
   config_restart_name = 'restart.nc'
   config_output_interval = '1_00:00:00'
   config_frames_per_outfile = 0
   config_pio_num_iotasks = 0
   config_pio_stride      = 1
/
&decomposition
   config_number_of_blocks = 0
   config_block_decomp_file_prefix = 'graph.info.part.'
   config_explicit_proc_decomp = .false.
   config_proc_decomp_file_prefix = 'graph.info.part.'
/
&restart
   config_restart_interval = '15_00:00:00'
   config_do_restart = .false.
/

&ctl_nl
NThreads      = 8
partmethod    = 4
topology      = "cube"
test_case = 'jw_baroclinic'
ne            = 30
qsize = 2 
ndays          = 30
statefreq     = 120
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = 90
energy_fixer  = 0
integration   = "explicit"
smooth        = 0.005               ! default = 0.005
nu            = 9.6e14
nu_s          = 9.6e14
nu_p          = 0
limiter_option = 4
hypervis_order = 2
hypervis_subcycle = 1 
 hypervis_subcycle_q = 1
u_perturb = 1 
rotate_grid = 0
/
&solver_nl
precon_method = "identity"
maxits        = 500
tol           = 1.e-9
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/camm-26.fbin.littleendian"
vfile_int     = "vcoord/cami-26.fbin.littleendian"
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
 output_timeunits=1
 output_frequency=1
 output_varnames1='ps'
 output_type='netcdf'
 io_stride=8
/
&analysis_nl
 output_timeunits=1,1,1
 output_frequency=0,0,0
 output_varnames1='ps'
 output_varnames2='zeta'
 output_varnames3='T'
/
&analysis_nlX2
 nlat=512
 nlon=1024
 output_timeunits=1,1,1,1
 output_frequency=1,1,1,5
 output_start_time=0,0,0,0
 output_end_time=30,30,30,30
 output_varnames1='ps'
 output_varnames2='zeta'
 output_varnames3='T'
 output_varnames4='u','v'
 num_io_procs=8
/
&analysis_nlXX2
 nlat=512
 nlon=1024
 output_timeunits=1,1,1,1
 output_frequency=1,1,1,5
 output_start_time=0,0,0,0
 output_end_time=30,30,30,30
 output_varnames1='ps'
 output_varnames2='zeta'
 output_varnames3='T','Q','Q2'
 output_varnames4='u','v'
 num_io_procs=8
/

&analysis_nlXX3
 gridtype=1
 nlat=91
 nlon=180
 output_timeunits=1,1,1,1
 output_frequency=1,1,1,5
 output_start_time=0,0,0,0
 output_end_time=30,30,30,30
 output_varnames1='ps'
 output_varnames2='zeta'
 output_varnames3='T','Q','Q2'
 output_varnames4='u','v'
 num_io_procs=8
/


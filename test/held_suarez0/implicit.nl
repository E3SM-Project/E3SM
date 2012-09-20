&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = 15
ndays          = 10
statefreq     = 480 
tracer_advection_formulation = 1
tstep_type    = 1
qsplit        = 1
qsize         = 1
ftype         = 0 
accumfreq     = -1
accumstart    = 300
accumstop     = 600
restartfreq   =  30
restartfile   = "restart/R0001"
restartdir    = "./restart/"
runtype       = 0
tstep         = 180
smooth         = 0.05
integration   = "full_imp"
nu            = 1e16
nu_s          = 1e16
nu_p          = 0
hypervis_order = 2
hypervis_subcycle = 1
/
&solver_nl
precon_method = "block_jacobi"
maxits        = 100
tol           = 1.e-12
debug_level   = 0
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/camm-26.fbin.littleendian"
vfile_int     = "../vcoord/cami-26.fbin.littleendian"
/
&analysis_nl
output_timeunits=1,1
output_frequency=6,6
output_start_time=200,200
output_end_time=30000,30000
output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
output_varnames2='DIFFT','DIFFU','DIFFV','CONVU','CONVV','FU','FV'
!io_stride = 8
!nlat=1024
!nlon=2048
/

&prof_inparm
profile_outpe_num = 100
/









&ctl_nl
NThreads      = 8
partmethod    = 4
topology      = "cube"
test_case     = "asp_rossby"
ne            = NE
qsize         = 0
tracer_advection_formulation = 1
ndays          = 30
statefreq     = SFREQ
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = TSTEP
integration   = "explicit"
smooth        = 0.005               ! default = 0.005
nu            = NU1
nu_s          = NU2
nu_p          = 0
hypervis_order = 2
hypervis_subcycle= -1
u_perturb      = 1
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
&analysis_nl
 gridtype=1
 output_timeunits=1
 output_frequency=1
 output_start_time=0
 output_end_time=30
 output_varnames1='ps','geos','T','u','v','geo','omega'
 io_stride=4
/

&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 1200
statefreq     = SFREQ
qsize         = 0
ftype         = 0 
restartfreq   =  50
restartfile   = "restart/R0001"
runtype       = RUNTYPE
tstep         = TSTEP
smooth        = .005
integration   = "explicit"
nu            = NU1
nu_s          = NU2
nu_p          = NUP
hypervis_order = ORDER
hypervis_subcycle = SUBCYCLE
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
filter_freq   = -1
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "/home/mataylo/homme/test/vcoord/sabm-20.ascii"
vfile_int     = "/home/mataylo/homme/test/vcoord/sabi-20.ascii"
/
&analysis_nl
output_timeunits=1,
output_frequency=1,
output_start_time=200
output_end_time=30000
output_varnames1='u','v','T','zeta','ps'
io_stride = 8
/











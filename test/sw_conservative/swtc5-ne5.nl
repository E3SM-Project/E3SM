&ctl_nl
NThreads      = 1
partmethod    = 0
topology      = "cube"
test_case     = "swtc5"
ne            = 5
ndays         = 15
statefreq     = 360
accumfreq     = 90
accumstart    = 300
accumstop     = 600
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = 120.
smooth        = 0.00
integration   = "explicit"
nu = 1e15
nu_s = 1e15
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
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
kcut_fm       = 2
wght_fm       = 0.10D0
/
&analysis_nl
output_timeunits=1,
output_frequency=5,
/












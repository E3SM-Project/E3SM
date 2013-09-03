&ctl_nl
NThreads      = 1
partmethod    = 0
topology      = "cube"
test_case     = "swtc5"
ne            = 5
ndays         = 15
accumfreq     = 90
accumstart    = 300
accumstop     = 600
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0

!statefreq     = 360
!tstep         = 180.     ! 200 also seems stable, 225 bad
!integration   = "explicit"

statefreq     = 60
tstep         = 450      ! 600 bad
integration   = "semi_imp"

smooth        = 0.05
nu = 0e5
nu_s = 0e5
hypervis_order = 0
hypervis_subcycle = 1
/
&solver_nl
!precon_method = "block_jacobi"
precon_method = "identity"
maxits        = 10000
tol           = 1.e-12
debug_level   = 0
/
&filter_nl
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 1
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
kcut_fm       = 2
wght_fm       = 0.10D0
/
&analysis_nl
output_timeunits=1,
output_frequency=5,
/












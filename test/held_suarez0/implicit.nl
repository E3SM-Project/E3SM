&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = 9
nmax          = 4
ndays         = 10
statefreq     = 1000
accumfreq     = -1
accumstart    = 300
accumstop     = 600
restartfreq   =  1
restartfile   = "./R0001"
runtype       = 0
tstep         = 100.0
integration   = "full_imp"
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
filter_freq   = 10
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
output_frequency=1,1
output_varnames1='u','v','T','zeta','ps','Q','DIFFT'
!output_varnames2='DIFFT','DIFFU','DIFFV','CONVU','CONVV','FU','FV'
io_stride = 8
/


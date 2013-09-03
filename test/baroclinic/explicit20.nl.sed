&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "baroclinic"
ne            = NE
qsize          = 0
ndays          = 12
statefreq     = SFREQ
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = TSTEP
integration   = "explicit"
smooth        = 0.05
nu            = NU1
nu_s          = NU2
hypervis_order = ORDER
hypervis_subcycle = SUBCYCLE
/
&solver_nl
precon_method = "identity"
maxits        = 50
tol           = 1.e-7
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/sabm-20.fbin.littleendian"
vfile_int     = "vcoord/sabi-20.fbin.littleendian"
/
&analysis_nl
 output_timeunits=1,
 output_frequency=4,
 output_varnames1='T','zeta'
 interp_type = 0          ! native high order
! interp_type = 1         ! bilinear
 output_type='netcdf'
 io_stride = 8
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

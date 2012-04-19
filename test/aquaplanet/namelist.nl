&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "aquaplanet"
ne            =  9
ndays         = 10
statefreq     = 300
accumfreq     = 200
accumstart    = 200 
accumstop     = 1200 
restartfreq   = 25 
! restartdir is the directory location for restart file output. 
! restartdir = "./restart/"
! restartfile is the full path name of the restart input file.
restartfile   = "./restart.one/R000000050"
runtype       = 0
tstep         = 30.0
integration   = "explicit"
smooth        = .05
nu            = 1.0e3
nu_s          = 1.0e3
moisture      = "vapor"
columnpackage = "emanuel"
/
&solver_nl
precon_method = "identity"
maxits        = 200 
tol           = 1.e-7
debug_level   = 0
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 1
filter_mu     = 0.1D0
filter_freq_advection   = 1
filter_mu_advection     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/sabm-18.fbin"
vfile_int     = "../vcoord/sabi-18.fbin"
/
&aquaplanet_nl
cool_ampl     = -1.0D0
cool_min      = 12.0D3
cool_max      = 15.0D3
qv_flag       = 0
qv_pert_flag  = 1
qv_pert_ampl  = 0.1D0
qv_pert_zmin  = 2.0D3
qv_pert_zmax  = 18.0D3
isrf_forc     = 1
h_dis         = 1000.0D0
Cdrag         = 0.001D0
wstar         = 1.0D0
tsurf         = 303.16D0
u0            = 0.D0
zabsampl      = 0.0D0
zabsmid       = 22.0D3
zabsmin       = 12.0D3
noisef = 2
/
&analysis_nl
!  currently up to 5 streams are allowed
!  output_stream_count=1
!
!  timunits: 0= steps, 1=days, 2=hours  
!  output_timeunits=1,2
!  output_start_time=0,1176,
!  output_end_time=-1,-1
!  output_frequency=1,1
!  output_dir ="./movies/"
!  
!  allowed variables are
!  'ps   ','zeta ','T    ','Th   ','u    ','v    ','ke   ','Q    ','prec ','accum',&
!  'udrag','vdrag','tsflx','qsflx','usfrc','vsfrc','tsfrc','qsfrc','geo  ','omega'
! 
!  output_varnames1-5
!
/

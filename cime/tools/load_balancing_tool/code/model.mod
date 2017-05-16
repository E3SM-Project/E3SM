# math model to find optimal allocation of cores for CESM components

param D >= 1, integer;       # number of data points
param CPUS >= 1, integer;    # number of nodes
param CPN >=1, integer;      # number of cores per node
param Tsync >= 0.0;          # time to sync ice and lnd in seconds
param Nlat >= 1, integer;    # max lat
param Nlev >= 1, integer;    # max lev
param Etarget >= 0.0;        # target Efficiency
param MinNodes >=1, integer; # minimum number of nodes
param MaxNodes >=1, integer; # maximum number of nodes

set M := {'lnd','ice','atm','ocn'};   # set of components
set DD := 1..D;         # set of data points

param rawx{M, DD};     # given by user
param rawy{M, DD};     # given by user

param x{DD};            # extracted automatically
param y{DD};            # extracted automatically

param A{M};            # best fit value
param B{M};            # best fit value
param C{M};            # best fit value
param K{M};            # best fit value

# special ordered set variables
var maxtasks >= 0;
var ntasks >= 0;
var nz >= 0;
var remainder >= 0;
var ny >= 0;
var taskcounter >= 0;
var ntasksrestrict >= 0;

# ... fitting parameters (bounds and initial values)
var a >= 0;
var b >= 0;
var c >= 0;
var k >= 0;
var eta_m1{M} >= 0;
var eta_m2{M} >= 0;
var eta_m3{M} >= 0;
var etaT >= 0;
var etaTi_m1 >= 0;

var n{M} >= 1, integer;

# efficiency var
var xcounter >= 1, integer;
var maxx >= 1, integer;
var fmod_xcounter >= 0;
var fmod_MinNodes >= 0;
var fmod_eff >= 0;

### special ordered sets
### Ocn
set OcnSet := 1..25;
param OcnPart{OcnSet};
var z2{OcnSet} binary;
subject to SOS2:           1 = sum{i in OcnSet} z2[i];
subject to DefNocn: n['ocn'] = sum{i in OcnSet} z2[i]*OcnPart[i];
### Atm
set AtmSet := 1..149;
param AtmPart{AtmSet};
var z1_2{AtmSet} binary;

# ... objective function is the least-squares error
minimize L2Error: sum{i in DD} (y[i] - a/x[i] - b*x[i]^c - k)^2;

minimize MaxTime: etaT;


### max(max(ice,lnd)+atm,ocn) model
# added suffix m1, meaning model 1, to all constraints
# time constraints
subject to
DefEta_m1{i in M}: eta_m1[i] = A[i]/n[i]  + B[i]*(n[i]^C[i]) + K[i];
DefEtaTi1_m1: etaTi_m1 >= eta_m1['ice'];
DefEtaTi2_m1: etaTi_m1 >= eta_m1['lnd'];
DefEta1_m1:   etaT  >= etaTi_m1 + eta_m1['atm'];
DefEtaT_m1:   etaT  >= eta_m1['ocn'];

# add constraint to force eta[ice] = eta[land]
EqualT1_m1: eta_m1['lnd'] >= eta_m1['ice'] - Tsync;
EqualT2_m1: eta_m1['lnd'] <= eta_m1['ice'] + Tsync;

# constrain number of nodes
TotalNumber_m1: n['atm'] + n['ocn'] <= CPUS;
IceLndNumber_m1: n['ice'] + n['lnd'] <= n['atm'];


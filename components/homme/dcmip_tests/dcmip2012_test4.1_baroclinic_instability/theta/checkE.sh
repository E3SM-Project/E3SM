#!/bin/tcsh 
#
#SBATCH -p ec
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 4
#SBATCH --time=0:10:00
#XXSBATCH -N 54
#XXSBATCH --time=12:00:00
#PBS -l walltime=10:00
#PBS -l nodes=4
#PBS -q acme
#
# hydrostatic x1:  4 nodes, 3.3min
#

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?PBS_NNODES} ) then   # redsky
    if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR     
    set NCPU = $PBS_NNODES
endif
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# hydrostatic
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30

#set case  = h-x1
#set case  = nh-x100
set case  = nh-x1000

\rm -f input.nl
mkdir restart


# run 15 days, eulerian,  (no restart)
#set namelist = ${case}-erun3.nl            
#\cp -f $namelist input.nl
#mpirun -np $NCPU $EXEC < input.nl

# run 15 days, write a restart file
set namelist = ${case}-erun1.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

# restart, run 2 timesteps with no viscosity:
set namelist = ${case}-erun2.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

exit 0

########################################
# h-x1 output:
########################################
#erun1:
KE,d/dt,diss:  0.16959893106811E+07  0.85016771568658E-01 -0.14983004751348E+00
IE,d/dt,diss:  0.26214194601295E+10 -0.11560603971835E+00  0.11924077936411E+00
 E,d/dt,diss:  0.26231154494402E+10 -0.30589268295853E-01

#erun2:
KE,d/dt,diss:  0.16960343913057E+07  0.22538985072169E+00 -0.18446424507157E-04
IE,d/dt,diss:  0.26214194150489E+10 -0.22538986682892E+00  0.18430317581225E-04
 E,d/dt,diss:  0.26231154494402E+10 -0.14305114746094E-07


########################################
# h-x1000 output:
########################################
#erun1:
KE,d/dt,diss:  0.17287572991040E+07  0.13452358263739E+03 -0.18908306774746E+03
IE,d/dt,diss:  0.18813163207871E+10 -0.28296764850616E+03  0.30175106945915E+02
PE,d/dt,diss:  0.73815345257935E+09  0.82550016442935E+02  0.92436756759263E+02
 E,d/dt,diss:  0.26211985306655E+10 -0.65894049008687E+02 -0.66471203625492E+02

#erun2: dt=.02
KE,d/dt,diss:  0.17287700957918E+07  0.32135602410417E+03  0.14425181334169E+01
IE,d/dt,diss:  0.18813163101971E+10 -0.26400336027145E+03  0.74164061239662E+00
PE,d/dt,diss:  0.73815345245899E+09 -0.28730630874634E+01  0.13618637544408E+00
 E,d/dt,diss:  0.26211985327519E+10  0.54479598999023E+02  0.23203433746311E+01




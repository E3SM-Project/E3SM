# 
# CMake initial cache file for Sandia's redsky
# configurd for redsky's netcdf-intel/4.1 module 
#
# With these modules one does not need binding options for mpiexec. 
#Currently Loaded Modules:
#  1) sems-env           4) sems-intel/18.0.0     7) sems-cmake/3.5.2       10) sems-netcdf/4.4.1/exo
#  2) intel/18.0.0.128   5) openmpi-intel/1.10    8) sems-git/2.10.1        11) totalview/2018.2.6
#  3) mkl/18.0.0.128     6) sems-openmpi/1.10.5   9) sems-hdf5/1.8.12/base
#
# However we may want to bring each machine file into sync with cime options for this machine. 
# Then mpiexec options (binding) may require changes.
# For example, right now in cime/config/e3sm/machines/config_machines.xml for sandiatoss3 
# these options are listed:
#  <mpirun mpilib="default">
#    <executable>mpiexec</executable>
#    <arguments>
#      <arg name="num_tasks"> --n {{ total_tasks }}</arg>
#      <arg name="tasks_per_node"> --map-by ppr:{{ tasks_per_numa }}:socket:PE=$ENV{OMP_NUM_THREADS} --bind-to core</arg>
#    </arguments>
#  </mpirun>
#  <mpirun mpilib="mpi-serial">
#    <executable></executable>
#  </mpirun>
#...
#    <modules>
#      <command name="purge"/>
#      <command name="load">sems-env</command>
#      <command name="load">sems-git</command>
#      <command name="load">sems-python/2.7.9</command>
#      <command name="load">sems-cmake</command>
#      <command name="load">gnu/4.9.2</command>
#      <command name="load">sems-intel/17.0.0</command>
#    </modules>
#    <modules mpilib="!mpi-serial">
#      <command name="load">sems-openmpi/1.10.5</command>
#      <command name="load">sems-netcdf/4.4.1/exo_parallel</command>
#    </modules>
#...


SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# Openmpi 1.8 only
SET (USE_MPI_OPTIONS " --map-by node:SPAN " CACHE FILEPATH "")

# Openmpi 1.6
#SET (USE_MPI_OPTIONS "-loadbalance" CACHE FILEPATH "")

# this is ignored if we use FORCE_Fortran_FLAGS
SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

# redsky upgrade 8/2017, need to load sems-netcdf module:
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{SEMS_NETCDF_ROOT} CACHE FILEPATH "")


SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")

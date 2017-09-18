#!/bin/bash
# script to make final modifications needed for Thwaites initial condition.
# run with one argument for the file to be modified.


# put it large negative SMB beyond current calving front
ncap2 -s "where(thickness(0,:)==0.0) sfcMassBal(0,:)=-10.0" $1 temp1.nc

# put small-ish beta value beyond current GL in case GL advances
ncap2 -s "where(thickness(0,:)*910.0/1028.0+bedTopography(0,:)<0.0) beta(0,:)=200.0" temp1.nc temp2.nc

cp temp2.nc $1
rm temp1.nc temp2.nc

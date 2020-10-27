#!/bin/sh 
# Test for missing r8
#


# Check physics
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/physics/cam
rc=$?
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/physics/waccm
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/physics/waccmx
rc=`expr $? + $rc`

#Check Chemistry
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/chemistry
rc=`expr $? + $rc`

#Check Dynamics
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/dynamics/se -s share
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/dynamics/fv
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/dynamics/eul
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/dynamics/sld
rc=`expr $? + $rc`

#Check other
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/advection
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/control
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/utils
rc=`expr $? + $rc`

#Check coupler
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/models/atm/cam/src/cpl
rc=`expr $? + $rc`

echo $rc

if [ $rc = 255 ]; then
   rc=1
fi

echo $rc

exit $rc

#!/bin/csh -f

# This script creates C header file gptl.h, Fortran header file gptl.inc, and
# Fortran module file gptlf.F90. DRY programming: define the header values
# in one place (this script)

# Set shell variables which will be substituted in the header files:

set GPTLsync_mpi        = 0
set GPTLwall            = 1
set GPTLcpu             = 2
set GPTLabort_on_error  = 3
set GPTLoverhead        = 4
set GPTLdepthlimit      = 5
set GPTLverbose         = 6
set GPTLnarrowprint     = 7
set GPTLpercent         = 9
set GPTLpersec          = 10
set GPTLmultiplex       = 11
set GPTLdopr_preamble   = 12
set GPTLdopr_threadsort = 13
set GPTLdopr_multparent = 14
set GPTLdopr_collision  = 15
set GPTLprint_method    = 16
set GPTLdopr_memusage   = 27
set GPTLtablesize       = 50
set GPTLmaxthreads      = 51

set GPTL_IPC           = 17
set GPTL_CI            = 18
set GPTL_FPC           = 19
set GPTL_FPI           = 20
set GPTL_LSTPI         = 21
set GPTL_DCMRT         = 22
set GPTL_LSTPDCM       = 23
set GPTL_L2MRT         = 24
set GPTL_LSTPL2M       = 25
set GPTL_L3MRT         = 26

set GPTLgettimeofday   = 1
set GPTLnanotime       = 2
set GPTLread_real_time = 3
set GPTLmpiwtime       = 4
set GPTLclockgettime   = 5
set GPTLpapitime       = 6
set GPTLplacebo        = 7

set GPTLfirst_parent  = 1
set GPTLlast_parent   = 2
set GPTLmost_frequent = 3
set GPTLfull_tree     = 4

# Create the sed script which will replace variables in the 3 header files
# with the above-defined values

cat >! sedscript <<EOF
s/#GPTLsync_mpi/$GPTLsync_mpi/1
s/#GPTLwall/$GPTLwall/1
s/#GPTLcpu/$GPTLcpu/1
s/#GPTLabort_on_error/$GPTLabort_on_error/1
s/#GPTLoverhead/$GPTLoverhead/1
s/#GPTLdepthlimit/$GPTLdepthlimit/1
s/#GPTLverbose/$GPTLverbose/1
s/#GPTLnarrowprint/$GPTLnarrowprint/1
s/#GPTLpercent/$GPTLpercent/1
s/#GPTLpersec/$GPTLpersec/1
s/#GPTLmultiplex/$GPTLmultiplex/1
s/#GPTLdopr_preamble/$GPTLdopr_preamble/1
s/#GPTLdopr_threadsort/$GPTLdopr_threadsort/1
s/#GPTLdopr_multparent/$GPTLdopr_multparent/1
s/#GPTLdopr_collision/$GPTLdopr_collision/1
s/#GPTLdopr_memusage/$GPTLdopr_memusage/1
s/#GPTLprint_method/$GPTLprint_method/1
s/#GPTLtablesize/$GPTLtablesize/1
s/#GPTLmaxthreads/$GPTLmaxthreads/1
s/#GPTL_IPC/$GPTL_IPC/1
s/#GPTL_CI/$GPTL_CI/1
s/#GPTL_FPC/$GPTL_FPC/1
s/#GPTL_FPI/$GPTL_FPI/1
s/#GPTL_LSTPI/$GPTL_LSTPI/1
s/#GPTL_DCMRT/$GPTL_DCMRT/1
s/#GPTL_LSTPDCM/$GPTL_LSTPDCM/1
s/#GPTL_L2MRT/$GPTL_L2MRT/1
s/#GPTL_LSTPL2M/$GPTL_LSTPL2M/1
s/#GPTL_L3MRT/$GPTL_L3MRT/1
s/#GPTLgettimeofday/$GPTLgettimeofday/1
s/#GPTLnanotime/$GPTLnanotime/1
s/#GPTLmpiwtime/$GPTLmpiwtime/1
s/#GPTLclockgettime/$GPTLclockgettime/1
s/#GPTLpapitime/$GPTLpapitime/1
s/#GPTLplacebo/$GPTLplacebo/1
s/#GPTLfirst_parent/$GPTLfirst_parent/1
s/#GPTLlast_parent/$GPTLlast_parent/1
s/#GPTLmost_frequent/$GPTLmost_frequent/1
s/#GPTLfull_tree/$GPTLfull_tree/1
s/#GPTLread_real_time/$GPTLread_real_time/1
EOF

# Run the sed script to create the 3 header files

sed -f sedscript gptl.h.template >! gptl.h       || (echo bad attempt to write to gptl.h && exit 1)
sed -f sedscript gptl.inc.template >! gptl.inc   || (echo bad attempt to write to gptl.inc && exit 1)
sed -f sedscript gptlf.F90.template >! gptlf.F90 || (echo bad attempt to write to gptlf.F90 && exit 1)

# Ensure that all variables got substituted by looking for "#" in the output files

grep -v '^#' gptl.h | grep '#'                   && (echo "found unsubstituted variables in gptl.h" && exit 1)
grep -v '^!' gptl.inc | grep '#'                 && (echo "found unsubstituted variables in gptl.inc" && exit 1)
grep -v '^!' gptlf.F90 | grep -v '^#' | grep '#' && (echo "found unsubstituted variables in gptlf.F90" && exit 1)

exit 0

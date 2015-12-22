#=======================================================
# Read the timebase register and return the result.
#
# C example:
#   long long tbvalue;
#   void timebase(long long *);
#   timebase(&tbvalue);
#
# Fortran example:
#   integer*8 tbvalue
#   call timebase(tbvalue)
#=======================================================

.globl	.timebase
.globl	.timebase_

.machine  "ppc64"

.csect	  timebase{PR}

.timebase:
.timebase_:
   mftb    6	   # move from timebase to gr6
   std	   6,0(3)  # store the value into the argument's address
   bclr    20,0    # return
.long	   0

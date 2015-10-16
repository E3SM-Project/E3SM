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
   mfspr   4, 269	     # move from timebase upper to gr4
   mfspr   5, 268	     # move from timebase lower to gr5
   mfspr   6, 269	     # move from timebase upper to gr6
   cmp	   0, 0, 4, 6	     # compare gr4 and gr6
   bc	   4, 2, .timebase_  # loop if upper register changed
   rldimi  5, 4, 32, 0	     # add lower and upper => 64-bit int
   std	   5, 0(3)	     # store into the argument's address
   bclr    20, 0	     # return
.long	   0

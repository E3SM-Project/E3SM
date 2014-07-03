
# *****************************COPYRIGHT****************************
# (c) British Crown Copyright 2009, the Met Office.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the
# following conditions are met:
# 
#     * Redistributions of source code must retain the above 
#       copyright  notice, this list of conditions and the following 
#       disclaimer.
#     * Redistributions in binary form must reproduce the above 
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials 
#       provided with the distribution.
#     * Neither the name of the Met Office nor the names of its 
#       contributors may be used to endorse or promote products
#       derived from this software without specific prior written 
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
# 
# *****************************COPYRIGHT*******************************
# *****************************COPYRIGHT*******************************
# *****************************COPYRIGHT*******************************

#$Id: Makefile,v 4.0 2009/02/13 08:19:14 hadmw Exp $

# sx-6 development
#F77=sxf90 -clear -ftrace -Wf"-L fmtlist transform -pvctl fullmsg" -Cvopt 

#Linux Lahey-Fujitsu compiler
#F77=lf95

#Linux gnu f77 compiler
F77=f77

DEBUG=-g -C

test_isccp_cloud_types: test_isccp_cloud_types.o isccp_cloud_types.o icarus.o scops.o Makefile
	$(F77) $(DEBUG) test_isccp_cloud_types.f isccp_cloud_types.o icarus.o scops.o -o test_isccp_cloud_types

test_congvec: test_congvec.o Makefile
	$(F77) $(DEBUG) test_congvec.o -o test_congvec

test_congvec.o: test_congvec.f congvec.f Makefile
	$(F77) $(DEBUG) -c test_congvec.f

clean:
	rm -f *.o test_isccp_cloud_types ftn* test_congvec

test:
	./test_isccp_cloud_types.ksh

test_isccp_cloud_types.o: test_isccp_cloud_types.f Makefile
	$(F77) $(DEBUG) -c test_isccp_cloud_types.f

isccp_cloud_types.o: isccp_cloud_types.f Makefile
	$(F77) $(DEBUG) -c isccp_cloud_types.f 

icarus.o: icarus.f Makefile
	$(F77) $(DEBUG) -c icarus.f 

scops.o: scops.f Makefile
	$(F77) $(DEBUG) -c scops.f 

tar: 
	./maketar


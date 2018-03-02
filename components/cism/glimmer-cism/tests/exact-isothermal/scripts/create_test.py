#!/usr/bin/env python
#
# Magnus Hagdorn
#
# Create test configuration based on Ed Bueler's exact solutions

import sys,string, Numeric, PyCF, datetime
from optparse import OptionParser
import exact_is

if __name__ == '__main__':

    # setup options
    parser = OptionParser(usage = "usage: %prog [options] name")
    parser.add_option("-m","--massbalance",metavar="MB",default=0.3,type="float",help="set accumulation [m/a] (default = 0.3)")
    parser.add_option("-n","--flow-law-exponent",metavar="N",default=3.,type="float",help="set Glen's flow law exponent (default = 3)")
    parser.add_option("-l","--margin-radius",metavar="L",default=750.,type="float",help="set radius of margin [km] (default = 750)")
    parser.add_option("-d","--delta",metavar="D",default=25.,type="float",help="set grid spacing [km] (default = 25)")

    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.error("no output file name given")

    # ice sheet
    print options.massbalance,options.flow_law_exponent,options.margin_radius
    ism = exact_is.ModelAE(options.massbalance,options.flow_law_exponent,options.margin_radius)

    # creating output netCDF file
    cffile = PyCF.CFcreatefile(args[0])
    # global attributes

    # history
    cffile.history = '%s: %s'%(datetime.datetime.today(),string.join(sys.argv))

    num = int(1.1*options.margin_radius/options.delta+0.5)
    num = 2*(num-1)+1
    # creating dimensions
    cffile.createDimension('x0',num-1)
    cffile.createDimension('x1',num)
    cffile.createDimension('y0',num-1)
    cffile.createDimension('y1',num)
    cffile.createDimension('level',1)
    cffile.createDimension('time',None)    
    #creating variables
    delta = options.delta*1000.
    varx=cffile.createVariable('x0')
    varx[:] = (delta/2.+delta*Numeric.arange(num-1)).astype(Numeric.Float32)
    varx=cffile.createVariable('x1')
    varx[:] = (delta*Numeric.arange(num)).astype(Numeric.Float32)

    vary=cffile.createVariable('y0')
    vary[:] = (delta/2.+delta*Numeric.arange(num-1)).astype(Numeric.Float32)
    vary=cffile.createVariable('y1')
    vary[:] = (delta*Numeric.arange(num)).astype(Numeric.Float32)

    varlevel=cffile.createVariable('level')
    varlevel[0] = 1

    vartime=cffile.createVariable('time')
    vartime[0] = 0

    varthck=cffile.createVariable('thk')
    varacab=cffile.createVariable('acab')
    
    centre = ((num-1)/2)*delta
    #print num*delta,centre
    for j in range(0,num):
        for i in range(0,num):
            r = ((varx[i]-centre)**2 + (vary[j]-centre)**2)**0.5
            varthck[0,j,i] = ism.getH(r)
            varacab[0,j,i] = ism.getMb(r)*ism.SperA
                        
    cffile.close()

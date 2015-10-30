#!/usr/bin/env python

import os
import sys
import cPickle as pick

try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

try:
    import matplotlib
except ImportError:
    print "Matplotlib was not found. "
matplotlib.use('Agg')

import uqtools as ut
from thisutils import *

## Main directory
thisdir=os.environ['ALMUQ_PATH']

## Checking input/output directories
filedir=thisdir+"/files/"
if(not os.path.isdir(filedir)):
    print "Directory of files %s does not exist. Exiting." % filedir
    sys.exit()
outfigdir=thisdir+"/outfig/"
if(not os.path.isdir(outfigdir)):
    os.makedirs(outfigdir)
    print "Directory of output figures %s does not exist. Creating one." % outfigdir


## Plotting beuatification
matplotlib.rc('legend',loc='upper left', fontsize=25)
matplotlib.rc('lines', linewidth=3, color='r')
matplotlib.rc('axes',linewidth=3,grid=True,labelsize=22)
matplotlib.rc('xtick',labelsize=18)
matplotlib.rc('ytick',labelsize=18)
###############################################################

## Argument parsing depends on run-mode
plid=int(sys.argv[1])
if(plid==0 or plid==1 or plid==2 or plid==3 or plid==4):
    varname=sys.argv[2]
if(plid==1 or plid==2 or plid==3):
    siteId=int(sys.argv[3])


## Read input
## HARDWIRED: use last 2000 samples removed one that had bogus output.
sampleInd=range(8000,10000)
del sampleInd[1529]

## Another scenario when the first 1000 samples were also used
#sampleInd=range(1000)
#sampleInd+=range(8000,10000)
#del sampleInd[2529]

## Get the input samples and parameter names
inputs_clm=np.loadtxt(filedir+'input_clm.dat')[sampleInd]
inputs_x=np.loadtxt(filedir+'input_x.dat')[sampleInd]
pnames_clm=np.genfromtxt(filedir+'pnames_clm.dat',dtype="string")
pnames_x=np.genfromtxt(filedir+'pnames_x.dat',dtype="string")


## Read netcdf file
## HARDWIRED dimensions, (years, sites, samples)=(20,96,2000), and remove a bogus sample
if (plid==0):
    ## Read netcdf into numpy file
    print "Reading Variable ", varname
    var_=read_nc(varname,filedir+os.sep+"FLUXNETUQ_8001-10000.nc",16, [20,96,2000], outfigdir)
    var_=removesam(outfigdir+os.sep+varname,1529,var_)
    
    ## Another scenario when the first 1000 samples are also taken
    #var_1=read_nc(varname,filedir+os.sep+"FLUXNETUQ_1-1000.nc",16, [20,96,1000], outfigdir)
    #var_2=read_nc(varname,filedir+os.sep+"FLUXNETUQ_8001-10000.nc",16, [20,96,2000], outfigdir)
    #var_=joinsam(outfigdir+os.sep+varname,var_1,var_2)
    #var_=removesam(outfigdir+os.sep+varname,2529,var_)
    sys.exit()

## Depending on run-mode, read site-specific information
if(plid==1 or plid==2 or plid==3):
    print "Reading Site # ", siteId
    var_site,out_ss,out_ave=readplot_qoi(varname,siteId,outfigdir=outfigdir)
    varname_siteId=varname+'_s'+str(siteId)

## Plot some diagnostic results
## HARDWIRED many inputs here - use as an example and tweak the inputs as necessary.
if (plid==1):

    ## Flag samples with vanishing output
    flag_nz=(out_ss==0.0) #(out<0.004)
    ## Plot parallel coordinates to highlight flagged samples
    plot_pcoord(inputs_x[:,:],flag_nz,11,6,outfigdir=outfigdir)
    ## Plot parameter 1(#11) vs parameter 2(#18) to highlight flagged samples
    plot_xxy(11,18,flag_nz,pnames_clm,inputs_clm,outfigdir=outfigdir)
    ## Read sensitivities if available
    #sens_all=np.loadtxt('sens_all.dat')
    
    ## Plot output wrt single input for a selected set of inputs
    for parid in [11,18]: #[1,2,8,14,27]:#[3,8,10,14,27]: #[1,2,8,14,27]:
        parid_clm=shift(parid)
        plot_output_x(varname_siteId,inputs_clm[:,parid_clm],out_ave, pnames_clm[parid_clm],outfigdir=outfigdir) #,tit=sens_all[siteId-1,parid])

## Read scaling info if computing surrogate or sensitivities
if(plid==2 or plid==3):
    logfl='unlog'
    out_scale=10./max(abs(out_ave))
    print "Scaling by a factor ",1./out_scale
    output=out_ave*out_scale

## Run surrogate construction via WBCS (Weighted Bayesian Compressive Sensing)
##     for a given output at a given site
if (plid==2):
    results=run_clwbcs(inputs_x,output,logfl) #,cl=[])
    pick.dump(results,open(outfigdir+os.sep+'results_'+varname_siteId+'.pk','wb'),-1)

## Compute sensitivities and visualize results for surrogate and sensitivity
if (plid==3):
    
    ## Load the results
    inputs_x,output,xtr,ytr,xval,yval,nval,yval_s,cfs_cur,mindex_cur,Sig,sigma2,erb=pick.load(open(outfigdir+os.sep+'results_'+varname_siteId+'.pk','rb'))
    
    ## Compute sensitivities
    msens,tsens,jsens= ut.pce_sens('LU',mindex_cur,cfs_cur)
    np.savetxt('pccf_'+varname_siteId+'.dat',cfs_cur)
    np.savetxt('mi_'+varname_siteId+'.dat', mindex_cur,fmt='%d')
    
    
    ## Plot surrogate versus model
    plot_results_xy(varname_siteId,yval,yval_s,erb,out_scale,logfl)
    ## Plot model/surrogate versus sample Id
    plot_results_surr(varname_siteId,yval,yval_s,erb,out_scale,logfl)
    ## Plot main and joint sensitivities, visualizeing on a circle
    plot_senscirc(varname_siteId,tsens,jsens,pnames_x)
    
    ## Copy figures to output folder and cleanup
    os.system('cp *.png '+outfigdir)
    os.system('rm -rf selected.dat regparams.dat Sig.dat sigma2* mi pccf varfrac.dat *sens.dat *mindex* args.in')

## Plot sensitivities across all sites for each output
## Run individually postp.py 4 <varname>
## Has to be done after uqmulti
if (plid==4):
    nsites=int(sys.argv[3])
    plotsites=range(1,nsites+1)
    
    npar=pnames_x.shape[0]
    allcolors=np.array(set_colors(npar))
    
    parrange=np.arange(npar)
    
    itask=0
    sens_all=np.zeros((nsites,npar))
    indpar=[]
    nlist=[]
    for siteId in plotsites:
        varname_siteId=varname+'_s'+str(siteId)
        
        inputs_x,output,xtr,ytr,xval,yval,nval,yval_s,cfs_cur,mindex_cur,Sig,sigma2,erb=pick.load(open(outfigdir+os.sep+'results_'+varname_siteId+'.pk','rb'))
        msens,tsens,jsens= ut.pce_sens('LU',mindex_cur,cfs_cur)
        sens_all[siteId-1,:]=tsens
        tmp=list(set(indpar) | set(parrange[tsens>0.2]))
        indpar=tmp
        nlist.append(str(siteId)) #('Site '+str(siteId))
    np.savetxt(outfigdir+os.sep+'sens_all_'+varname+'.dat',sens_all)
    plot_sensmat(sens_all,pnames_x,nlist,showplot=outfigdir+os.sep+'sensmat_'+varname+'.png')
    plot_sensbar(sens_all[:,indpar],range(len(indpar)),range(nsites),vis="bar",reverse=False,par_labels=pnames_x[indpar],case_labels=nlist,colors=allcolors[indpar,:],showplot=outfigdir+os.sep+'sensbar_'+varname+'.png')
    os.system('rm -rf mi pccf varfrac.dat *sens.dat *mindex*')

## Plot sensitivities across all outputs for each site
## Run individually postp.py 5 <siteId> GPP TLAI ...
## Has to be done after uqmulti
if (plid==5):
    siteId=int(sys.argv[2])
    plotvars=sys.argv[3:]
    nvars=len(plotvars)
    #print plotvars
    npar=pnames_x.shape[0]
    allcolors=np.array(set_colors(npar))
    
    parrange=np.arange(npar)
    
    itask=0
    sens_all=np.zeros((nvars,npar))
    indpar=[]
    nlist=[]
    ivar=0
    for varname in plotvars:
        varname_siteId=varname+'_s'+str(siteId)
        
        inputs_x,output,xtr,ytr,xval,yval,nval,yval_s,cfs_cur,mindex_cur,Sig,sigma2,erb=pick.load(open(outfigdir+os.sep+'results_'+varname_siteId+'.pk','rb'))
        msens,tsens,jsens= ut.pce_sens('LU',mindex_cur,cfs_cur)
        sens_all[ivar,:]=tsens
        tmp=list(set(indpar) | set(parrange[tsens>0.2]))
        indpar=tmp
        nlist.append(varname)
        ivar+=1
    
    np.savetxt(outfigdir+os.sep+'sens_all_s'+str(siteId)+'.dat',sens_all)
    plot_sensmat(sens_all,pnames_x,nlist,showplot=outfigdir+os.sep+'sensmat_s'+str(siteId)+'.png')
    plot_sensbar(sens_all[:,indpar],range(len(indpar)),range(nvars),vis="bar",reverse=False,par_labels=pnames_x[indpar],case_labels=nlist,colors=allcolors[indpar,:],showplot=outfigdir+os.sep+'sensbar_s'+str(siteId)+'.png')
    os.system('rm -rf mi pccf varfrac.dat *sens.dat *mindex*')






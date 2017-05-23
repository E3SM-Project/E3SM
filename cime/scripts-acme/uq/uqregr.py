#!/usr/bin/env python

import os
try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

import uqtools as ut

uqtkpath=os.environ['UQTK_PATH']

#############################################################
#############################################################
#############################################################
def pcfunc(xdata,param):
    
    sam=xdata.shape[0]
    dim=xdata.shape[1]
    np.savetxt('xdata.dat',xdata)
    

    mindex=param[0]
    pccf=param[1]
    pctype=param[2]
    
    npc=mindex.shape[0]
    
    #print "mindex=",mindex, pccf
    
    np.savetxt('pccf.dat',pccf)
    np.savetxt('mi.dat',mindex,fmt='%d')
    assert npc==pccf.shape[0]
    
    
    # Evaluate PC to get the training data
    cmd=uqtkpath+'/bin/pce_eval -x"PC_mi" -f"pccf.dat" -s'+pctype+' -r"mi.dat" > pceval.log'
    print "Running ", cmd
    os.system(cmd)
    ydata=np.loadtxt('ydata.dat')
  
    
    return ydata

#############################################################
def pce_lsq(xtr,ytr,xtarget,pars):
    basisparams=pars[0] #pctype,mindex
    regparams=pars[1] #'wbcs',np.ones((npc,))
    
    cfs_cur,mindex_cur,Sig=regression(xtr,ytr,basisparams,regparams) #pctype,'lsq',mindex,np.ones((npc,)))
    
    return pcfunc(xtarget,[mindex_cur,cfs_cur,basisparams[0]])

#############################################################

def pce_wbcs(xtr,ytr,xtarget,pars):
    
    basisparams=pars[0] #pctype,mindex
    regparams=pars[1] #'wbcs',np.ones((npc,))
    iterparams=pars[2] #niter,eps
    #print pars
    splitparams=pars[3] #split_method,KK,npt

    cfs_cur,mindex_cur,Sig,used=regression_splititer(xtr,ytr,basisparams,regparams,iterparams,splitparams)
    
    return (pcfunc(xtarget,[mindex_cur,cfs_cur,basisparams[0]]),cfs_cur,mindex_cur,Sig,used)

#############################################################
#############################################################

def regression_splititer(xdata,ydata,basisparams,regparams,iterparams,splitparams):
    
    pctype,mindex=basisparams
    #method,methodpars=regparams
    split_method,KK,npt=splitparams
    
    
    ns=xdata.shape[0]
    list_ind=ut.ind_split(ns,split_method,[KK,npt])
    nsets=len(list_ind)
    used_prev=np.arange(mindex.shape[0])
    for i in range(nsets):
        xdata_cur=xdata[list_ind[i],:]
        ydata_cur=ydata[list_ind[i]]
        print "Running regression for subset ", str(i+1), "/", str(nsets)
        cfs,mindex_new,Sig,used=regression_iter(xdata_cur,ydata_cur,basisparams,regparams,iterparams)
        print "The size of basis at subset   ", str(i+1), "/", str(nsets), ": ", used.shape[0]
        #print nzind,used
        #if (i==0):
        #    inter_ind=nzind
        #else:
        #    inter_ind=inter_ind*nzind
        used_prev=np.intersect1d(used,used_prev)
        #used_prev=list(set(used_prev).intersection(used))
        if (i>0):
            mindex_tmp=np.array([val for val in mindex_inter if ut.in_array(val,mindex_new)])
        else:
            mindex_tmp=mindex_new.copy()
            
        mindex_inter=mindex_tmp.copy()
    
    #mindex_inter=mindex[used_prev,:]
    nmi = mindex_inter.shape[0] #len(used_prev) #sum(inter_ind)
    print "The size of intersection basis:", nmi
#print used_prev
#npc_full=mindex.shape[0]
    #indx_orig=np.array(range(npc_full))[inter_ind]
    #mindex_inter=mindex[indx_orig]
#    mindex_inter=mindex[used_prev]
            
    # Final regression, if any
    
        #if KK>1:
    basisparams[1]=mindex_inter
    regparams[0]='lsq'
    print "Running %s with the intersection basis" % regparams[0]
    cfs_final,mindex_inter,Sig,_=regression(xdata,ydata,basisparams,regparams)
    
            #pctype,'lsq',mindex_inter,[]) # Do lsq or method?
        #else:
            #cfs_final=cfs

    
    return (cfs_final,mindex_inter,Sig,used_prev) #,indx_orig)

#############################################################
#############################################################

def regression_iter(xdata,ydata,basisparams,regparams,iterparams):
    
    pctype,mindex=basisparams
    method,methodpars=regparams

    basisparams_cur=[pctype,mindex]
    regparams_cur=[method,methodpars]
    
    niter,eps,update_weights,update_mindex=iterparams
    # update_weights=True
    #update_mindex=False
    nrange=np.arange(mindex.shape[0])
    cur_used=nrange
    npc=mindex.shape[0]
    for i in range(niter):
        print "Iteration %d / %d  " % (i+1,niter)
        #print regparams
        print "Initial mindex size ", basisparams_cur[1].shape[0]
        cfs_cur,mindex_cur,Sig,used=regression(xdata,ydata,basisparams_cur,regparams_cur)
        print "New mindex size     ", cfs_cur.shape[0], mindex_cur.shape[0], used.shape[0]
   
        #np.savetxt('mi.'+str(i+1)+'.dat',mindex_cur,fmt='%d')
        #np.savetxt('cf.'+str(i+1)+'.dat',cfs_cur)

        npc_cur=mindex_cur.shape[0]
        #nrange=np.array(range(npc_cur))
        #ind=[cfs_cur*cfs_cur>0]
        #indx=nrange[ind]
        #print cfs_cur
        #cfs_cur=cfs_cur[ind]
        if (update_weights==True):
            regparams_cur[1]=1./(abs(cfs_cur)+eps)
        else:
            tmp=regparams_cur[1]
            regparams_cur[1]=tmp[list(used)] #read used.dat and replace it here
        #mindex_cur=mindex_cur[indx,:]
        #np.savetxt('mipp.'+str(i+1)+'.dat',mindex_cur,fmt='%d')
        if (update_mindex==True and i<niter-1):
            mindex_new,mindex_add,mindex_f=ut.mi_addfront(mindex_cur)
            #np.savetxt('min.'+str(i+1)+'.dat',mindex_add,fmt='%d')

            mindex_cur=mindex_new.copy()
            basisparams_cur[1]=mindex_new
            regparams_new=np.ones(mindex_new.shape[0])/eps
            regparams_new[0:npc_cur]=regparams_cur[1]
            regparams_cur[1]=regparams_new
    return (cfs_cur,mindex_cur,Sig,cur_used)


#############################################################
#############################################################

def regression(xdata,ydata,basisparams,regparams):
    
    pctype,mindex=basisparams
    method,methodpars=regparams
    
    dim=mindex.shape[1]
    np.savetxt('xdata.dat',xdata)
    np.savetxt('ydata.dat',ydata)
    mindex_file='mindex.dat'
    np.savetxt(mindex_file,mindex,fmt='%d')
    params_file='regparams.dat'
    
    np.savetxt(params_file,np.array(methodpars).reshape(-1,1),fmt='%24.16f')
    #print "Current mindex ", mindex, mindex.shape
    #cmd='lin_reg -p'+pctype+' -m'+mindex_file+' -f'+params_file+' -t'+method+' > linreg.log'
    cmd=uqtkpath+'/bin/regression -x xdata.dat -y ydata.dat -b PC_MI -s '+pctype+' -p '+mindex_file+' -w '+params_file+' -m msc -r '+method #+' > regr.log'
    print "Running "+cmd
    os.system(cmd)
    cfs=np.loadtxt('coeff.dat')
    mindex=np.loadtxt('mindex.dat').reshape(-1,dim)
    used=np.loadtxt('selected.dat')
    #print "Selected bases ", used, used.shape
    ##npc=cfs.shape[0]
    ##used=np.array(range(npc))
    #phi=np.loadtxt('phi.dat')
    #mcoh=mut_coherence(phi)
    #print "Mutual coherence = ", mcoh
    Sig=np.loadtxt('Sig.dat')
    #print mindex,mindex.shape

    
    return (cfs,mindex,Sig,used)


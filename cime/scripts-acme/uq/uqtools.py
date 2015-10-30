#!/usr/bin/env python


import os
import sys

try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

uqtkpath=os.environ['UQTK_PATH']

###############################################################

def write_sp_mindex(mindex,filename='sp_mindex.dat'):
    sp_mindex=encodeMindex(mindex)
    thefile = open(filename, 'w')
    for item in sp_mindex:
        print>>thefile, item


##############################################################

def gen_mi(mi_type,params):

    if mi_type=='TO':
        nord=params[0]
        ndim=params[1]
        cmd=uqtkpath+'/bin/gen_mi -x' + mi_type + ' -p' + str(nord) + ' -q' + str(ndim)

    elif mi_type=='HDMR':
        hdmr_dims=params[0]
        dim=params[1]
        np.savetxt('hdmr_dims.dat',np.array(hdmr_dims),fmt='%d')
        print dim
        cmd=uqtkpath+'/bin/gen_mi -x' + mi_type + ' -f hdmr_dims.dat -q'+str(dim)

    else: 
        print "Multiindex type is not recognized. Exiting."
        sys.exit(1)

    os.system(cmd + ' > gen_mi.log')
    mindex=np.loadtxt('mindex.dat')
    
    return mindex
    
#############################################################

def mi_addfront(mindex):
    print "Adding multiindex front"
    npc=mindex.shape[0]
    ndim=mindex.shape[1]

    mindex_f=np.zeros((1,ndim),dtype=int)
    mindex_add=np.zeros((1,ndim),dtype=int)
    mindex_new=np.zeros((1,ndim),dtype=int)
    for i in range(npc):
        cur_mi=mindex[i,:]
        #print i
        #print "Current multiindex: ", cur_mi
        fflag=True
        for j in range(ndim):
            test_mi=np.copy(cur_mi)
            test_mi[j] += 1
            if not any(np.equal(mindex,test_mi).all(1)):
                if not any(np.equal(mindex_add,test_mi).all(1)):
                    mindex_add=np.vstack((mindex_add,test_mi))
                if fflag:
                    mindex_f=np.vstack((mindex_f,cur_mi))
                fflag=False

    mindex_f=mindex_f[1:]
    mindex_add=mindex_add[1:]
    mindex_new=np.vstack((mindex,mindex_add))


    print "Old MI size  : ", mindex.shape[0]
    print "New MI size  : ", mindex_new.shape[0]

    return [mindex_new,mindex_add,mindex_f]

#############################################################

def get_npc(ord,dim):

    npc = 1
    
    for i in range(ord):
        npc = npc * (dim+i+1)
    for i in range(ord):
        npc = npc / (i+1)

    return npc

#############################################################

def encodeMindex(mindex):

    npc=mindex.shape[0]
    ndim=mindex.shape[1]
    print "Multiindex has %d terms" % npc
    sp_mindex=[]
    for ipc in range(npc):
        nzs=np.nonzero(mindex[ipc,:])[0]
        if nzs.shape[0]==0:
            this_sp_mindex=np.zeros((1,2))
        else:
            this_sp_mindex=np.vstack((nzs+1,mindex[ipc,nzs])).T
            this_sp_mindex=this_sp_mindex.reshape(1,-1)

        #effdim=len(np.nonzero(mindex[ipc,:])[0])
        #print effdim
        sp_mindex.append(this_sp_mindex)
    return sp_mindex

#############################################################

def pce_sens(pctype,mi,pccf,mv=False):
    np.savetxt('mi',mi,fmt="%d")
    np.savetxt('pccf',pccf)

    cmd=uqtkpath+'/bin/pce_sens -m mi -f pccf -x '+pctype + ' > pcsens.log'
    os.system(cmd)

    mainsens=np.loadtxt('mainsens.dat')
    totsens=np.loadtxt('totsens.dat')
    jointsens=np.loadtxt('jointsens.dat')
    varfrac=np.loadtxt('varfrac.dat')
    
    if (mv):
        mean=pccf[0]
        var=mean**2 / varfrac[0]
    
        return mainsens,totsens,jointsens,mean,var

    else:
        return mainsens,totsens,jointsens

##############################################################

def in_array(arr1,arr2):
    assert(arr1.shape[0]==arr2.shape[1])
    for i in arr2:
        if np.array_equal(arr1,i):
            return True

    return False

##############################################################

def ind_split(ns,split_method,split_params):

    KK=split_params[0]
              
              
    if (split_method=='Kfold_small'):
        #ind=np.ones((ns))
        indp=np.random.permutation(ns)
        list_ind=np.array_split(indp,KK)
        #ind[i*sam_spl:(i+1)*sam_spl]=np.zeros((sam_spl))
        #intind=ind.astype(int)
            
    elif (split_method=='Kfold'):
        indp=np.random.permutation(ns)
        full_ind=range(ns)
        fold_ind=np.array_split(indp,KK)
        list_ind=[]
        for cur_ind in fold_ind:
            list_ind.append(np.array(list(set(full_ind)-set(cur_ind))))

    elif (split_method=='rand_fold'):
        npt=split_params[1]

        list_ind=[]
        for i in range(KK):
            list_ind.append(np.random.permutation(ns)[0:npt])



    return list_ind

    

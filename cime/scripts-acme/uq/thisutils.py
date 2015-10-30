#!/usr/bin/env python

import os
import sys
import cPickle as pick
from math import *

try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

try:
    import matplotlib
except ImportError:
    print "Matplotlib was not found. "
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm    as cm
from matplotlib.patches import Circle, Wedge, Polygon

import uqtools as ut
import uqregr as ur


uqtkpath=os.environ['UQTK_PATH']
ncdumppath=os.environ['NCDUMP_PATH']


############################################################
def read_nc(varname,ncfile, firstline,coords=[20,96,2000], outfigdir='.'):
    
    nyear,nsites, nsam= coords
    
    tmp=outfigdir.split(os.sep)
    if (len(tmp[-1])==0):
        del tmp[-1]
    outfigdir_short=tmp[-1]+os.sep
    
    cmd=ncdumppath+"ncdump -v "+varname+" "+ncfile+" | awk 'NR>"+str(firstline)+"{print}' | sed 's/,//g' | awk '{for ( i = 1 ; i <= NF ; i++ ) { printf ( \"%s\\n\",$i ) }}' | awk '$1!=\"\" && NR>=1{print}' | sed '$d' | sed '$d' > " + outfigdir+ varname+".out"
    print cmd
    os.system(cmd)
    print "Saved all data for %s in %s in a planar format" % (varname,outfigdir_short+varname+'.out')
    
    var=np.loadtxt(outfigdir+os.sep+varname+'.out')
    var_=var.reshape(nyear,nsam,nsites)
    
    np.save(outfigdir+os.sep+varname,var_)
    print "Saved all data for %s in %s in numpy format" % (varname,outfigdir_short+varname+'.py')
    
    return var_

###############################################################
def removesam(outfigdirvarname,remSamInd,var_):
    
    tmp=outfigdirvarname.split(os.sep)
    if (len(tmp[-2])==0):
        del tmp[-2]
    outfigdirvarname_short=os.sep.join(tmp[-2:])

    tmp=np.delete(var_,remSamInd,axis=1)


    var_=tmp.copy()
    np.save(outfigdirvarname,var_)
    print "Removed sample %d and stored in %s" % (remSamInd,outfigdirvarname_short+'.npy')
    return var_
###############################################################
def joinsam(outfigdirvarname,var_1,var_2):

    tmp=outfigdirvarname.split(os.sep)
    if (len(tmp[-2])==0):
        del tmp[-2]
    outfigdirvarname_short=os.sep.join(tmp[-2:])


    var_=np.concatenate((var_1,var_2),axis=1)
    np.save(outfigdirvarname,var_)
    print "Joined two arrays (%d,%d,%d)+(%d,%d,%d)=(%d,%d,%d) and stored in %s" % (var_1.shape[0],var_1.shape[1],var_1.shape[2],var_2.shape[0],var_2.shape[1],var_2.shape[2],var_.shape[0],var_.shape[1],var_.shape[2],outfigdirvarname_short+'.npy')
    return var_

############################################################
def readplot_qoi(varname,siteId,ind=[],outfigdir='.'):
    
    var_=np.load(outfigdir+os.sep+varname+'.npy')
    nyear,nsam,nsites=var_.shape
    print nyear,nsam,nsites
    if ind==[]:
        ind=range(nsam)
    
    var_site=var_[:,ind,siteId-1]
    print var_site.shape
    plt.plot(range(1,nyear+1),var_site,lw=1)
    plt.xlabel('Time')
    plt.ylabel(varname)
    plt.yscale('log')
    plt.savefig(outfigdir+os.sep+varname+'_s'+str(siteId)+'_t.png')
    #plt.show()
    plt.clf()
    np.savetxt(outfigdir+os.sep+varname+'_s'+str(siteId)+'.dat',var_site)
    
    out_ss=var_site[-1,:]
    out_ave=np.average(var_site,axis=0)
    
    #out_ave_sorted=np.sort(out_ave)
    
    np.savetxt(outfigdir+os.sep+varname+'_s'+str(siteId)+'_ss.dat',out_ss)
    np.savetxt(outfigdir+os.sep+varname+'_s'+str(siteId)+'_ave.dat',out_ave)
    
    plt.plot(out_ss,'bo',ms=3)
    plt.xlabel('Run Id')
    plt.ylabel(varname)
    plt.yscale('log')
    plt.savefig(outfigdir+os.sep+varname+'_s'+str(siteId)+'_ss.png')
    #plt.show()
    plt.clf()
    
    #out_ave_sorted=out_ave.sort()
    #plt.plot(out_ave_sorted,'bo',ms=3)
    plt.plot(out_ave,'bo',ms=3)
    plt.xlabel('Run Id')
    plt.ylabel(varname)
    plt.yscale('log')
    plt.savefig(outfigdir+os.sep+varname+'_s'+str(siteId)+'_ave.png')
    #plt.show()
    plt.clf()
    
    
    
    return var_site,out_ss,out_ave
###############################################################
def shift(parid):
    parid_clm=parid
    if parid>=23:
        parid_clm+=1
    if parid>=25:
        parid_clm+=1
    
    return parid_clm

###############################################################
def split_data(x,y,frac):
    
    nsam=x.shape[0]
    permuted=np.random.permutation(nsam)
    indtr=permuted[:int(frac*nsam)]
    indval=permuted[int(frac*nsam):]
    
    xtr=x[indtr,:]
    xval=x[indval,:]
    ytr=y[indtr]
    yval=y[indval]
    
    return xtr,xval,ytr,yval

####################################################
def run_clwbcs(xin_unsc,output,logfl='unlog'):
    
    xtr,xval,ytr,yval=split_data(xin_unsc,output,0.999)
    nval=xval.shape[0]
    nsam=xin_unsc.shape[0]
    ntr=nsam-nval
    nd=xin_unsc.shape[1]
    
    xval=np.vstack((xtr,xval))
    yval=np.append(ytr,yval)
    
    mindex=ut.gen_mi('TO',[1,nd])
    basisparams=['LU',mindex]
    regparams=['wbcs',np.ones((mindex.shape[0],))]
    iterparams=[2,1.e-3,True,True]
    splitparams=['rand_fold',3,int(0.9*ntr)]#['rand_fold',3,int(0.9*ntr)] #['Kfold',2,0]
    
    
    if (logfl=='log'):
        yval_s,cfs_cur,mindex_cur,Sig,used=ur.pce_wbcs(xtr,log(ytr),xval,(basisparams,regparams,iterparams,splitparams))
        yval_s=exp(yval_s)
    else:
        yval_s,cfs_cur,mindex_cur,Sig,used=ur.pce_wbcs(xtr,ytr,xval,(basisparams,regparams,iterparams,splitparams))

    sigma2=np.loadtxt('sigma2')
    #phi=get_phi(xval,mindex_cur)
    #erb=sqrt(np.diag(np.dot(phi,np.dot(Sig,phi.T)))+sigma2)
    print "Computing the errorbar"
    erb=get_erb(xval,mindex_cur, Sig,sigma2)
    
    
    
    relerr_tr=np.linalg.norm(yval_s[:nval]-yval[:nval])/np.linalg.norm(yval[:nval])
    print "Training error ",relerr_tr
    relerr_val=np.linalg.norm(yval_s[nval:]-yval[nval:])/np.linalg.norm(yval[nval:])
    print "Validation error ", relerr_val
    
    
    
    
    results=xin_unsc,output,xtr,ytr,xval,yval,nval,yval_s,cfs_cur,mindex_cur,Sig,sigma2,erb
    
    return results

###############################################################
def get_erb(xval,mindex, Sig,sigma2):
    npt=xval.shape[0]
    npc=mindex.shape[0]
    ndim=mindex.shape[1]
    erb=np.zeros((npt,))
    
    
    #for i in range(npt):
    #    for j in range(maxord):
    #        leg1d[i,j]=scipy.special.eval_legendre(j,xval[i,k])
    
    
    #This is commented out - it takes very long need to accelerate this somehow!
    """
        for i in range(npt):
        bas=np.zeros((npc,))
        for j in range(npc):
        mleg=1.0
        for k in range(ndim):
        mleg = mleg * scipy.special.eval_legendre(mindex[j,k],xval[i,k])
        bas[j]=mleg
        sum=0.0
        
        for j1 in range(npc):
        for j2 in range(npc):
        sum = sum + bas[j1]*bas[j2]*Sig[j1,j2]
        erb[i]=sqrt(sum+sigma2)
        """
    return sqrt(sigma2)*np.ones((npt,))#erb


############################################################
############################################################

def plot_output_x(varname,inp,out,pname,tit=[], outfigdir='.'):
    
    plt.plot(inp,out,'bo',ms=3)
    plt.xlabel(pname)
    plt.ylabel(varname)
    #plt.yscale('log')
    if tit!=[]:
        title('Sens '+str(tit))
    plt.savefig(outfigdir+os.sep+varname+'_'+pname+'.png')
    #plt.show()
    plt.clf()


############################################################

def plot_pcoord(inputs,labels,ndcut=9,ndg=6,outfigdir='.'):
    
    
    nd=inputs.shape[1]
    
    
    assert(ndg*ndcut==nd)
    for i in range(ndg):
        print "Plotting %d / %d " % (i+1,ndg)
        names=range(1+i*ndcut,1+(i+1)*ndcut)
        
        values=inputs[:,i*ndcut:(i+1)*ndcut].T
        #labels=(output[:nscut]<0.02) #values[0,:]<0.0#(np.dot(values.T, np.ones((ndcut,)))<0.0) #(output[:nscut]<0.02)
        
        labels_only=labels[labels==True]
        values_only=values[:,labels==True]
        print values_only.shape
        #pt.parallel_coordinates(names, values, labels)
        parallel_coordinates(names, values_only, labels_only, outfigdir+os.sep+'pcoord_'+str(names[0])+'_'+str(names[-1])+'.png')

###############################################################
def plot_xxy(x1,x2,flag,pnames,inputs,outfigdir='.'):
    print "Plotting inputs [%s x %s]" % (pnames[x1],pnames[x2])
    plt.figure(figsize=(10,10))
    
    intflag=[int(ifl) for ifl in flag]
    print sum(intflag)
    
    plt.plot(inputs[:,x1],inputs[:,x2],'o',ms=2)
    plt.plot(inputs[flag,x1],inputs[flag,x2],'ro',ms=7)
    plt.xlabel(pnames[x1])
    plt.ylabel(pnames[x2])
    plt.savefig(outfigdir+os.sep+'xxy__'+pnames[x1]+'__'+pnames[x2]+'.png')
    #plt.show()
    plt.clf()

###############################################################
def plot_results_xy(varname,yval,yval_s,erb,outsc,logfl="unlog"):
    
    yval=yval/outsc
    yval_s=yval_s/outsc
    erb=erb/outsc
    
    nval=yval_s.shape[0]
    if logfl=="log":
        errl=yval_s[:nval]*(1.-exp(-erb[:nval]))
        errh=yval_s[:nval]*(exp(erb[:nval])-1.)
    else:
        errl=erb[:nval]
        errh=erb[:nval]
    
    #plt.plot(yval[:nval],yval_s[:nval],'b*')
    plt.errorbar(yval[:nval],yval_s[:nval],yerr=[errl,errh],fmt='b*', label='Model')
    #fill_between(pr_grid,mmmr[inm,:]-mmsr[inm,:],mmmr[inm,:]+mmsr[inm,:],color='grey')
    plt.xlabel('Model')
    plt.ylabel('Surrogate')

    plt.savefig('res_xy_'+varname+'.png');
    plt.clf()


###############################################################
def plot_results_surr(varname,yval,yval_s,erb,outsc,logfl="unlog"):
    
    
    
    yval=yval/outsc
    yval_s=yval_s/outsc
    erb=erb/outsc
    
    nval=yval.shape[0]
    ppdat_uns=np.vstack((yval[:nval],yval_s[:nval])).T
    ppdat=ppdat_uns[ppdat_uns[:,1].argsort()]
    if logfl=="log":
        errl=ppdat[:,1]*(1.-exp(-erb[:nval]))
        errh=ppdat[:,1]*(exp(erb[:nval])-1.)
    else:
        errl=erb[:nval]
        errh=erb[:nval]
    runs=range(ppdat.shape[0]);
    fig = plt.figure(figsize=(8,6));
    ax=fig.add_axes([0.12, 0.12, 0.83, 0.83]);
    ax.errorbar(runs[::10], ppdat[::10,1], yerr=[errl[::10], errh[::10]], fmt='--',color='grey')
    pl1 =plt.plot(runs[::10],ppdat[::10,0],linestyle="none", marker=".", markersize=4, color="k",label='CLM data')
    pl2 =plt.plot(runs[::10],ppdat[::10,1],linestyle="none", marker=".", markersize=6, markeredgecolor='r', color="r",label='Surrogate model')
    """
        ax.annotate('CLM data', xy=(155, 1.5),  xycoords='data',
        xytext=(-50, 40.0), textcoords='offset points',
        size=20, color='k',
        bbox=dict(boxstyle="round4,pad=.5", fc="1.0", ec='k', lw=2),
        arrowprops=dict(arrowstyle="->",linewidth=2))
        ax.annotate('Surrogate Model', xy=(950, 5),  xycoords='data',
        xytext=(550, 15.0), textcoords='data',
        size=20, color='r',
        bbox=dict(boxstyle="round4,pad=.5", fc="1.0", ec='r', lw=2),
        arrowprops=dict(arrowstyle="->",color='r',linewidth=2))
        ax.annotate('5-95% Confidence Interval', xy=(900, 4.0),  xycoords='data',
        xytext=(-350, 80.0), textcoords='offset points',
        size=20, color='r',
        bbox=dict(boxstyle="round4,pad=.5", fc="0.9", ec='grey', lw=2),
        arrowprops=dict(arrowstyle="->",color='r',linewidth=2))
        """
    #ax.set_xlim([0,9900]);
    #ax.set_ylim([0,30]);
    #ax.set_ylim([0,60000]);
    
    #ax.set_xticks([200+200*i for i in range(4)]);
    #ax.set_yticks([5+5*i   for i in range(4)]);
    plt.xlabel("runID",fontsize=18)
    plt.ylabel(varname,  fontsize=18)
    #plt.yscale('log')
    plt.legend()
    #plt.show()
    plt.savefig('res_surr_'+varname+'.png');
    plt.clf()

###############################################################
def plot_senscirc(varname, msens,jsens,inpar_names):
    #varname,_=parse_varids(varid)
    vnames=inpar_names
    
    Nmain=6
    Nsec=Nmain-1
    lwMax=10
    lwCut=0.2
    radMain=50
    radOut=15
    lext=0.4
    verbose=0
    
    nx,ny=jsens.shape
    for i in range(nx):
        for j in range(ny):
            if ( i != j ):
                if jsens[i,j]<5.e-5:
                    jsens[i,j]=5.e-5;
    #jsens=np.log10(jsens);
    #print msens
    ind=msens.argsort()[::-1];
    msensShort=msens[ind[0:Nmain]]
    if verbose > 0:
        for i in range(Nmain):
            print "Variable ",ind[i],", main sensitivity ",msens[ind[i]]
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_axes([0.05, 0.05, 0.9, 0.9],aspect='equal')
    #circ=pylab.Circle((0,0),radius=0.5,color='r')
    circ=Wedge((0.0,0.0),1.01, 0, 360, width=0.02,color='r')
    ax.add_patch(circ)
    maxJfr=-1.e10;
    for i in range(Nmain):
        jfr_i=np.array(np.zeros(nx))
        iord=ind[i]
        for j in range(iord):
            jfr_i[j]=jsens[j,iord]
        for j in range(iord+1,nx):
            jfr_i[j]=jsens[iord,j]
        ind_j=jfr_i.argsort()[::-1];
        if jfr_i[ind_j[0]] > maxJfr: maxJfr = jfr_i[ind_j[0]];
        if verbose > 1:
            for j in range(Nsec):
                print iord," ",ind_j[j],jfr_i[ind_j[j]]
    if verbose > 1:
        print "Maximum joint sensitivity :",maxJfr
    gopar=[]
    for i in range(Nmain):
        jfr_i=np.array(np.zeros(nx))
        iord=ind[i]
        for j in range(iord):
            jfr_i[j]=jsens[j,iord]
        for j in range(iord+1,nx):
            jfr_i[j]=jsens[iord,j]
        ind_j=jfr_i.argsort()[::-1];
        elst=[]
        for j in range(Nsec):
            if jfr_i[ind_j[j]]/maxJfr >= lwCut:
                posj=[k for k,x in enumerate(ind[:Nmain]) if x == ind_j[j]]
                if verbose > 2:
                    print j," ",posj
                if len(posj) > 0 :
                    x1=np.cos(0.5*np.pi+(2.0*np.pi*posj[0])/Nmain)
                    x2=np.cos(0.5*np.pi+(2.0*np.pi*i      )/Nmain)
                    y1=np.sin(0.5*np.pi+(2.0*np.pi*posj[0])/Nmain)
                    y2=np.sin(0.5*np.pi+(2.0*np.pi*i      )/Nmain)
                    lw=lwMax*jfr_i[ind_j[j]]/maxJfr
                    plt.plot([x1,x2],[y1,y2],'g-',linewidth=lw)
                    if ( verbose > 2 ):
                        print iord," ",ind[posj[0]]
                else:
                    elst.append(j)
        if len(elst) > 0:
            asft=[0,-1,1]
            for k in range(min(len(elst),3)):
                ang=0.5*np.pi+(2.0*np.pi*i)/Nmain+2*np.pi/12*asft[k]
                x2=np.cos(0.5*np.pi+(2.0*np.pi*i)/Nmain)
                y2=np.sin(0.5*np.pi+(2.0*np.pi*i)/Nmain)
                x1=x2+lext*np.cos(ang)
                y1=y2+lext*np.sin(ang)
                lw=lwMax*jfr_i[ind_j[elst[k]]]/maxJfr
                plt.plot([x1,x2],[y1,y2],'g-',linewidth=lw)
                plt.plot([x1],[y1],"wo",markersize=radOut,markeredgecolor='k',
                         markeredgewidth=2)
                if ( ind_j[elst[k]] > 32 ):
                    ltext=str(ind_j[elst[k]]+3)
                elif ( ind_j[elst[k]] > 30 ):
                    ltext=str(ind_j[elst[k]]+2)
                else:
                    ltext=str(ind_j[elst[k]]+1)
                plt.text(x1+(0.15)*np.cos(ang),y1+(0.15)*np.sin(ang),ltext,
                            ha='center',va='center',fontsize=16)
                posj=[k1 for k1,x in enumerate(gopar) if x == ind_j[elst[k]]]
                if len(posj)==0:
                    gopar.append(ind_j[elst[k]])
        if ( verbose > 2 ):
            print "------------------------"
    for i in range(Nmain):
        angl=0.5*np.pi+(2.0*np.pi*i)/Nmain
        xc=np.cos(angl);
        yc=np.sin(angl);
        msize=radMain*msens[ind[i]]/msens[ind[0]]
        plt.plot([xc],[yc],"bo",markersize=msize,markeredgecolor='k',markeredgewidth=2)
        da=1.0
        lab=0.2
        llab=lab*msens[ind[i]]/msens[ind[0]]
        
        ltext=str(ind[i]+1)
        lleg=ltext+" - "+vnames[ind[i]]
        plt.text(xc+(0.08+llab)*np.cos(angl+da),yc+(0.08+llab)*np.sin(angl+da),ltext,
                 ha='center',va='center',fontsize=16)
        plt.text(1.6,1.2-0.15*i,lleg,fontsize=16)
    for k in range(len(gopar)):
        lleg=str(gopar[k]+1)+" - "+vnames[gopar[k]]
        plt.text(1.6,1.2-0.15*Nmain-0.15*k,lleg,fontsize=16)
    print gopar
    ax.set_xlim([-1-1.6*lext,1.8+1.6*lext])
    ax.set_ylim([-1-1.6*lext,1+1.6*lext])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig("res_senscirc_"+str(varname)+".png")


###############################################################
def plot_sensmat(sens_all,plist,nlist,showplot=[]):
    
    nobs=sens_all.shape[0]
    
    
    vlst=[]
    allSens=[]
    for iobs in range(nobs):
        vfr=sens_all[iobs,:] #np.loadtxt('si_'+str(nm)+'_0.dat')
        #print iobs,vfr
        allSens.append(vfr)
        vlst.append([ n for n,i in enumerate(vfr) if i>0.004 ])
    # Get union
    allV=[]
    for i in range(len(vlst)):
        allV=list(set(allV) | set(vlst[i]))
    allV=np.sort(allV)
    #print allV
    # Create matrix, populate, and rescale
    npar=len(allV);
    print "Number of observables = ", nobs
    print "Number of parameters  = ", npar
    jsens=np.array(np.zeros([nobs,npar]));
    for i in range(nobs):
        for j in range(npar):
            jsens[i,j]=max(allSens[i][allV[j]],0.0);
    for i in range(nobs):
        jsens[i]=jsens[i]/jsens[i].max();
    jsens[np.where(jsens==0)]=0.5*jsens[np.where(jsens>0)].min();
    for i in range(nobs):
        for j in range(npar):
            jsens[i,j]=np.log10(jsens[i,j]);


    # make fig
    fs1=9;
    fig = plt.figure(figsize=(10,3.9));
    ax=fig.add_axes([0.12, 0.27, 0.88, 0.68]);
    cp=get_cp()
    cs=ax.pcolor(jsens,cmap=cp);
    #cs=ax.pcolor(jsens,cmap=cm.jet)
    ax.set_xlim([0,npar]);
    ax.set_ylim([0,nobs]);
    ax.set_xticks([0.5+i for i in range(npar)]);
    ax.set_yticks([0.5+i for i in range(nobs)]);
    ax.set_yticklabels([nlist[i].upper() for i in range(nobs)],fontsize=fs1);
    ax.set_xticklabels([plist[allV[i]].upper() for i in range(npar)],rotation=90,fontsize=fs1);
    cbar=plt.colorbar(cs)
    cbar.set_ticks([-3, -2,-1,0])
    cbar.set_ticklabels(['$10^{-3}$', '$10^{-2}$','$10^{-1}$','$10^0$'])
    if (showplot==[]):
        plt.show()
    else:
        plt.savefig(showplot)
        plt.clf()
#############################################################

def plot_sensbar(sensdata,pars,cases,vis="bar",reverse=False,par_labels=[],case_labels=[],colors=[],showplot=[]):
    """Plots sensitivity for multiple observables"""
    
    ncases=sensdata.shape[0]
    npar=sensdata.shape[1]
    
    wd=0.6
    xlbl=''
    ylbl='Sensitivity'
    legend_show=True #False
    #print len(pars), npar
    assert set(pars) <= set(range(npar))
    assert set(cases) <= set(range(ncases))
    
    # Set up the figure
    # TODO need to scale figure size according to the expected amount of legends
    xticklabel_size=min(21,1000/ncases)
    fig = plt.figure(figsize=(20,12))
    #fig = plt.figure(figsize=(18,12))
    fig.add_axes([0.1,0.3,0.8,0.65])
    
    #########
    
    # Default parameter names
    if (par_labels==[]):
        for i in range(npar):
            par_labels.append(('par_'+str(i+1)))
    # Default case names
    if (case_labels==[]):
        for i in range(ncases):
            case_labels.append(('case_'+str(i+1)))


    if(reverse):
        tmp=par_labels
        par_labels=case_labels
        case_labels=tmp
        tmp=pars
        pars=cases
        cases=tmp
        sensdata=sensdata.transpose()
    ##############################################################################

    npar_=len(pars)
    ncases_=len(cases)

    # Create colors list
    if colors==[]:
        colors=ut.set_colors(npar_)
    
    
    case_labels_=[]
    for i in range(ncases_):
        case_labels_.append(case_labels[cases[i]])
    
    if (vis=="graph"):
        for i in range(npar_):
            plt.plot(np.array(range(1,ncases_+1)),sensdata[cases,i], '-o',color=colors[pars[i]], label=par_labels[pars[i]])
    elif (vis=="bar"):
        curr=np.zeros((ncases_))
        #print pars,colors
        for i in range(npar_):
            plt.bar(np.array(range(1,ncases_+1)),sensdata[cases,i], width=wd,color=colors[pars[i]], bottom=curr, label=par_labels[pars[i]])
            curr=sensdata[cases,i]+curr
        if ncases>20:
            plt.xticks(np.array(range(1,ncases_+1))+wd/2.,case_labels_,rotation='vertical')
        else:
            plt.xticks(np.array(range(1,ncases_+1))+wd/2.,case_labels_)
        plt.xlim(1-wd/2.,ncases_+1.5*wd)
    #plt.ylim(0.0,1.0)
    
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.ylim([0,1])
    if (legend_show):
        #plt.legend(bbox_to_anchor=(1.0, -0.05),fancybox=True, shadow=True,ncol=5,labelspacing=-0.1)
        plt.legend(bbox_to_anchor=(0.0, -0.05),fancybox=True, shadow=True,ncol=5,labelspacing=-0.1)

    zed = [tick.label.set_fontsize(xticklabel_size) for tick in plt.gca().xaxis.get_major_ticks()]

    if (showplot==[]):
        plt.show()
    else:
        plt.savefig(showplot)
        plt.clf()



##################################################
def get_cp():
    cdict = cm.jet._segmentdata.copy()
    cdict['red']=tuple([tuple([0.0,  1,   1  ]),
                        tuple([0.01, 0,   0  ]),
                        tuple([0.35, 0,   0  ]),
                        tuple([0.66, 1,   1  ]),
                        tuple([0.89, 1,   1  ]),
                        tuple([1,    0.5, 0.5])
                        ]
                       )
    cdict['green']=tuple([tuple([0.0,   1, 1]),
                         tuple([0.01,  0, 0]),
                         tuple([0.125, 0, 0]),
                         tuple([0.375, 1, 1]),
                         tuple([0.64,  1, 1]),
                         tuple([0.91,  0, 0]),
                         tuple([1,     0, 0])
                         ]
                        )
    cdict['blue']=tuple([tuple([0,    1.0,1.0]),
                        tuple([0.01, 0.5,0.5]),
                        tuple([0.11, 1,  1  ]),
                        tuple([0.34, 1,  1  ]),
                        tuple([0.65, 0,  0  ]),
                        tuple([1,    0,  0  ])
                        ]
                       )
   
    cp=matplotlib.colors.LinearSegmentedColormap('colormap',cdict,64)
    return cp

##################################################

def set_colors(npar):
    """ Sets a list of different colors of requested length"""
    colors = []
    pp=1+npar/6
    for i in range(npar):
        c=1-(float) (i/6)/pp
        b=np.empty((3))
        for jj in range(3):
            b[jj]=c*int(i%3==jj)
        a=int(i%6)/3
        colors.append(((1-a)*b[2]+a*(c-b[2]),(1-a)*b[1]+a*(c-b[1]),(1-a)*b[0]+a*(c-b[0])))
    
    return colors



#############################################################
def parallel_coordinates(coordinates, values, labels, savefig=[]):
    """Plot 2d array `values` using K parallel coordinates.
        see https://github.com/btel/matplotlib_2014/tree/zurich_2013/exercises
        and https://python.g-node.org/python-summerschool-2013/_media/wiki/datavis/exercises.html
        Arguments:
        
        coordinates -- list or array of K elements containg coordinate
        names,
        values -- (K,N)-shaped array of N data points with K
        coordinates,
        labels -- list or array of one string per data point
        describing its class membership (category)
        """
    
    # SOLUTION
    ax = plt.subplot(111)
    
    # find names and number of different classes
    ulabels = np.unique(labels)
    n_labels = len(ulabels)
    
    # for each select distinct colors from Accent pallette
    cmap = plt.get_cmap('Accent')
    colors = cmap(np.arange(n_labels)*cmap.N/(n_labels+1))
    
    # change the label strings to indices into class names array
    class_id = np.searchsorted(ulabels, labels)
    lines = plt.plot(values[:,:], 'ko-',ms=4,linewidth=0.3)
    [ l.set_color(colors[c]) for c,l in zip(class_id, lines) ]
    
    
    # add grid, configure labels and axes
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('none')
    
    plt.xticks(np.arange(len(coordinates)), coordinates)
    plt.grid(axis='x', ls='-')
    
    leg_handlers = [ lines[np.where(class_id==id)[0][0]]
                    for id in range(n_labels)]
    ax.legend(leg_handlers, ulabels, frameon=False, loc='upper left',
                    ncol=len(labels),
                    bbox_to_anchor=(0, -0.03, 1, 0))
    if (savefig==[]):
        plt.show()
    else:
        plt.savefig(savefig)
        plt.clf()



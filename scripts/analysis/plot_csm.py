#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
FOLD="_figs/plot_csm"

def main():
    os.popen("mkdir -p "+FOLD)
    if len(sys.argv)==1:
        use_folds=glob.glob("*/*/out_epw_2/*")
    else:
        use_folds=sys.argv[1:]
    
    for fold in use_folds:
        print "-"*80
        print "Using folder: ",fold
        print
        # read in loc.py and get data in there
        try:
            cur_loc=imp.load_source("loc",fold+"/loc.py")
        except:
            print "  Cant open file!"
            continue
        print " ecut: ",cur_loc.ECUT,"  kmesh: ",cur_loc.KMESH
        print
        for ln in cur_loc.EPW_SECOND_RUN.split("\n"):
            ln_use=ln.strip()
            if ln_use!="\n" and ln_use!="":
                print " "+ln_use
        print

        # read in lambda information on a regular mesh
        data_reshape=get_lambda(fold+"/epw.out",reshape=True)

        # find minimal and maximal omega over all smearing
        omega_0=np.min(data_reshape["omega"][:,:,:,:,:])
        #  rescale max omega a bit
        omega_1=1.3*np.max(data_reshape["omega"][:,:,:,:,:])

        # get number of q points in the mesh
        numq=float(data_reshape["lambda"].shape[0])*\
             float(data_reshape["lambda"].shape[1])*\
             float(data_reshape["lambda"].shape[2])
        # maximal lambda mean over q and summed over branches, for all smearings
        max_lmbd_all_smear=np.max(data_reshape["lambda"][:,:,:,:,:].sum(axis=3).sum(axis=2).sum(axis=1).sum(axis=0))
        max_lmbd_all_smear=max_lmbd_all_smear/numq

        # plot histograms of lambda
        for ism in range(data_reshape["num_smears"]):
            # get data_reshape first 
            lmbd=data_reshape["lambda"][:,:,:,:,ism]
            omega=data_reshape["omega"][:,:,:,:,ism]

            # get rid of all indices
            lmbd=lmbd.ravel()
            omega=omega.ravel()

            # doing mean over q points
            lmbd=lmbd/numq

            # make histogram of all contributions to lambda as a function of frequency
            hist,bin_edges=np.histogram(omega,bins=100,range=[omega_0,omega_1],weights=lmbd)
            bin_centers=(bin_edges[1:]+bin_edges[:-1])/2.0

            # make cumulative sum of all contributions
            csm=np.cumsum(hist)

            # plot cumulative sum
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.plot(bin_centers,csm,"ko-")
            ax.set_xlim(omega_0,omega_1)
            ax.set_ylim(0.0,max_lmbd_all_smear)
            ax.set_xlabel("Phonon energy (meV)")
            ax.set_ylabel("Cumulative $\lambda$")
            fig.savefig(FOLD+"/csm__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")

if __name__=="__main__":
    main()

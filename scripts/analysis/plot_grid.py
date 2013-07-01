#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
FOLD="_figs/plot_grid"

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

        # prepare text with data
        txt_info="ecut: "+str(cur_loc.ECUT)+"  kmesh: "+str(cur_loc.KMESH)+"\n"
        txt_info+="\n"
        for ln in cur_loc.EPW_SECOND_RUN.split("\n"):
            ln_use=ln.strip()
            if ln_use!="\n" and ln_use!="":
                txt_info+=ln_use+"\n"

        # read in lambda information on a regular mesh
        data_reshape=get_lambda(fold+"/epw.out",reshape=True)

        # plot phonon DOS
        fig=plt.figure()
        ax=fig.add_subplot(111)
        plt.hist(data_reshape["omega"][:,:,:,:,0].ravel(),50)
        fig.savefig(FOLD+"/phdos__"+fold.replace("/","__")+".pdf")
        
        # find maximal lambda over all smearing
        max_lmbd_all_smear=np.max(data_reshape["lambda"][:,:,:,:,:])

        # plot histograms of lambda
        for ism in range(data_reshape["num_smears"]):
            # get data_reshape first 
            lmbd=data_reshape["lambda"][:,:,:,:,ism]
            omega=data_reshape["omega"][:,:,:,:,ism]
            # sum over branches
            aver_lmbd=lmbd.sum(axis=3)
            
            # plot all contributions to lambda in a histogram
            fig=plt.figure()
            ax=fig.add_subplot(111)
            plt.hist(aver_lmbd.ravel(),50,range=[0.0,max_lmbd_all_smear*1.05])
            fig.savefig(FOLD+"/hist_lambda__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")
         
            # plot slices on the third direction of lambda
            #
            # first find range of lambda for this choice of smearing
            max_val=np.max(lmbd.sum(axis=3))
            if np.min(lmbd.sum(axis=3))<-1.0E-3:
                print "Negative lambda?"
                sys.exit(1)
            for ind_z in range(omega.shape[2]):
                lmbd_slice=lmbd.sum(axis=3)[:,:,ind_z]
                # plot stuff
                fig=plt.figure()
                ax=fig.add_subplot(111)
                ax.imshow(lmbd.sum(axis=3)[:,:,ind_z],interpolation='none',
                          vmin=0.0,vmax=max_val*1.01, origin='lower',
                          cmap=plt.cm.gray, extent=(0,1,0,1))
                fig.savefig(FOLD+"/slice_z__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+"__slice_"+"%03d" % ind_z+".pdf")

if __name__=="__main__":
    main()

#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
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

        # plot histograms of lambda
        for ism in range(data_reshape["num_smears"]):
            # get data_reshape first 
            lmbd=data_reshape["lambda"][:,:,:,:,ism]
            omega=data_reshape["omega"][:,:,:,:,ism]
            # sum over branches
            aver_lmbd=lmbd.sum(axis=3)
            
            # plot all contributions to lambda ina histogram
            fig=plt.figure()
            ax=fig.add_subplot(111)
            plt.hist(aver_lmbd.ravel(),50)
            fig.savefig("_figs/hist_lambda__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")

            # plot phonon DOS
            fig=plt.figure()
            ax=fig.add_subplot(111)
            plt.hist(omega.ravel(),50)
            fig.savefig("_figs/phdos__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")
          
            # plot slices on the third direction of lambda
            #
            # first find range of data_reshape
            max_val=np.max(lmbd.sum(axis=3))
            if np.min(lmbd.sum(axis=3))<-1.0E-3:
                print "Negative lambda?"
                sys.exit(1)
            for ind_z in range(omega.shape[2]):
                lmbd_slice=lmbd.sum(axis=3)[:,:,ind_z]
                # plot stuff
                fig=plt.figure()
                ax=fig.add_subplot(111)
                plt.imshow(lmbd.sum(axis=3)[:,:,ind_z],interpolation='nearest',
                           vmin=0.0,vmax=max_val*1.01)
                fig.savefig("_figs/slice_z__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+"__slice_"+"%03d" % ind_z+".pdf")

if __name__=="__main__":
    main()

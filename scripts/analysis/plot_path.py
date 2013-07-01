#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
FOLD="_figs/plot_path"

def main():
    os.popen("mkdir -p "+FOLD)
    #
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

        # read in lambda information in the order written in file
        data=get_lambda(fold+"/epw.out",reshape=False)
        
        # go over all smearings
        for ism in range(data["num_smears"]):
            # get data first
            if data.has_key("val_omega"):
                omega=data["val_omega"][:,:,ism]
            if data.has_key("val_nesting"):
                nesting=data["val_nesting"][:,ism]
            if data.has_key("val_lambda"):
                lmbd=data["val_lambda"][:,:,ism]

            if data.has_key("val_omega") and data.has_key("val_lambda"):
                nkpt=omega.shape[0]
                
                # plot freq, and lambda as size of the point
                fig=plt.figure()
                ax=fig.add_subplot(111)
                fig.subplots_adjust(0.15,0.15,0.8,0.9)
                for i in range(omega.shape[1]):
                    ax.plot(np.arange(omega.shape[0]),
                            omega[:,i],"k-",lw=0.75,zorder=5)
#                    ax.scatter(np.arange(omega.shape[0]),
#                               omega[:,i],c="k",s=0.2,lw=0.0,zorder=5)
                    ax.scatter(np.arange(omega.shape[0]),
                               omega[:,i],c="g",s=10.0*lmbd[:,i],lw=0.0,alpha=0.85,zorder=10)
#                ax.set_title("   ismear: "+str(ism+1))
                ax.set_xlabel("q path")
                ax.set_ylabel("Phonon frequency (meV)")
                ax.set_xlim(0.0,nkpt)
                ax.set_xticks([])
                ax.text(1.05, 1.00, txt_info,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=3)
                fig.savefig(FOLD+"/omla__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")
                
                # plot lambda
                fig=plt.figure()
                ax=fig.add_subplot(111)
                fig.subplots_adjust(0.15,0.15,0.8,0.9)
                ax.plot(np.arange(lmbd.shape[0]),lmbd[:,:].sum(axis=1),"go-",ms=3.0)
#                ax.set_title("   ismear: "+str(ism+1))
                ax.set_xlabel("q path")
                ax.set_ylabel("$\lambda$")
                ax.set_xticks([])
                ax.text(1.05, 1.00, txt_info,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=3)
                fig.savefig(FOLD+"/lambda__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")

            if data.has_key("val_nesting"):
                # plot nesting
                fig=plt.figure()
                ax=fig.add_subplot(111)
                fig.subplots_adjust(0.15,0.15,0.8,0.9)
                ax.plot(np.arange(nesting.shape[0]),nesting[:],"ro-",ms=3.0)
#                ax.set_title("   ismear: "+str(ism+1))
                ax.set_xlabel("q path")
                ax.set_ylabel("nesting")
                ax.set_xticks([])
                ax.text(1.05, 1.00, txt_info,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=3)
                fig.savefig(FOLD+"/nesting__"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".pdf")

if __name__=="__main__":
    main()

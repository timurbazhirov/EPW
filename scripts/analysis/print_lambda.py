#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys

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
        # read in lambda information
        data=get_lambda(fold+"/epw.out",reshape=False)
        print "Lambda summed over branches:"
        # go over all smearings
        for ism in range(data["num_smears"]):
            lmbd=data["val_lambda"][:,:,ism]
            # sum over branches
            lmbd=lmbd.sum(axis=1)
            print "  Smearing index: ",ism+1,
            print "  Average lambda: ",nice(np.mean(lmbd)),
            print "  Min,max lambda: ",nice(np.min(lmbd)),nice(np.max(lmbd))
        print

def nice(x): return str(round(x,4)).ljust(7)
        
if __name__=="__main__":
    main()

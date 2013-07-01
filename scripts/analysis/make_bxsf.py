#!/usr/bin/env python

from get_lmbd import *
import imp
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import xml.dom.minidom

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

        # get reciprocal vectors
        (rec_vec_0,rec_vec_1,rec_vec_2)=get_recip_vectors(fold+"/../../run/01/_work/pref.save/data-file.xml")
        
        # read in lambda information on a regular mesh
        data=get_lambda(fold+"/epw.out",reshape=True)
        # go over all smearings
        for ism in range(data["num_smears"]):
            # get data first 
            omega=data["omega"][:,:,:,:,ism]
            lmbd=data["lambda"][:,:,:,:,ism]
            # make bxsf file based on http://www.xcrysden.org/doc/XSF.html#bxsf
            f=open("_bxsf/"+fold.replace("/","__")+"__smear_"+"%02d" % (ism+1)+".bxsf","w")
            f.write("BEGIN_INFO\n")
            f.write("END_INFO\n")
            f.write("\nBEGIN_BLOCK_BANDGRID_3D\n")
            f.write("  --\n")
            f.write("  BEGIN_BANDGRID_3D_lambda_sum_over_branches\n")
            f.write("    1\n")
            f.write("   "+(" "+str(data["mesh_size"]+1))*3+"\n")
            f.write("    0.0 0.0 0.0\n")
            for i,vec in enumerate([rec_vec_0,rec_vec_1,rec_vec_2]):
                f.write("   ")
                for j in range(3):
                    f.write(" %.15f"%vec[j])
                f.write("\n")
            f.write("  BAND:  1\n")
            for i in range(data["mesh_size"]+1):
                for j in range(data["mesh_size"]+1):
                    f.write("    ")
                    for k in range(data["mesh_size"]+1):
                        f.write("  %.6f"%lmbd[i%data["mesh_size"],
                                              j%data["mesh_size"],
                                              k%data["mesh_size"],:].sum())
                    f.write("\n")
                if i<data["mesh_size"]:
                    f.write("\n")
            f.write("  END_BANDGRID_3D\n")
            f.write("END_BLOCK_BANDGRID_3D\n")
            f.close()
            
def get_recip_vectors(fname):
    # open xml file
    doc=xml.dom.minidom.parse(fname)
    # get reciprocal vectors
    cell=doc.getElementsByTagName("CELL")[0]
    recip=cell.getElementsByTagName("RECIPROCAL_LATTICE_VECTORS")[0]
    units=cell.getElementsByTagName("UNITS_FOR_RECIPROCAL_LATTICE_VECTORS")[0].getAttribute("UNITS")
    if units!="2 pi / a":
        print "Wrong units! ",units
        sys.exit(1)
    rec_vec_0=np.array(map(float,recip.getElementsByTagName("b1")[0].firstChild.data.split()))
    rec_vec_1=np.array(map(float,recip.getElementsByTagName("b2")[0].firstChild.data.split()))
    rec_vec_2=np.array(map(float,recip.getElementsByTagName("b3")[0].firstChild.data.split()))
    return (rec_vec_0,rec_vec_1,rec_vec_2)
    
if __name__=="__main__":
    main()

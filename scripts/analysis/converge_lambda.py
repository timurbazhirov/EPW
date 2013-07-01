#!/usr/bin/env python

from get_lmbd import *
import matplotlib.pyplot as plt
import imp
import glob
import sys

def main():
    if len(sys.argv)==1:
        use_folds=glob.glob("*/*")
    else:
        use_folds=sys.argv[1:]
    for sub_fold in use_folds:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        leg=[]
        leg_lab=[]
        for fold in glob.glob(sub_fold+"/out_epw_2/*"):
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
            
            # read in which smearings are used in the epw.in file
            f=open(fold+"/epw.in")
            ln=f.readlines()
            f.close()
            epw_params={"degaussw":None,
                        "nsmear":None,
                        "delta_smear":None,
                        "nkf1":None,
                        "nkf2":None,
                        "nkf3":None,
                        "nqf1":None,
                        "nqf2":None,
                        "nqf3":None,
                        }
            for l in ln:
                for k in epw_params:
                    if k in l and "!" not in l:
                        if epw_params[k]!=None:
                            print "Problem! More than one "+k+" in epw.in file."
                            sys.exit(1)
                        epw_params[k]=float(l.split("=")[1])
            # default values of params
            if epw_params["nsmear"]==None:
                epw_params["nsmear"]=1
            if epw_params["delta_smear"]==None:
                epw_params["delta_smear"]=0.0
            # consistency check
            if epw_params["nsmear"]!=data["num_smears"]:
                print "Wrong number of smearings!"
                sys.exit(1)
            # plot only if all variables are defined
            if None not in epw_params.values():
                # go over all smearings
                xdata=[]
                ydata=[]
                for ism in range(data["num_smears"]):
                    # compute value of smearing for this case
                    cur_smear=epw_params["degaussw"]+float(ism)*epw_params["delta_smear"]
                    # get lambda for this case
                    lmbd=data["val_lambda"][:,:,ism]
                    # sum over branches
                    lmbd=lmbd.sum(axis=1)
                    # average lambda value
                    av_lmbd=np.mean(lmbd)
                    # store values
                    xdata.append(cur_smear)
                    ydata.append(av_lmbd)
                # get legend for this case
                txt_leg="q: "
                for k in ["nqf1","nqf2","nqf3"]:
                    txt_leg+=str(int(epw_params[k]))+" "
                txt_leg+=" k: "
                for k in ["nkf1","nkf2","nkf3"]:
                    txt_leg+=str(int(epw_params[k]))+" "
                # plot and store this case if had at least one data point
                if len(xdata)>0:
                    leg.append(ax.plot(xdata,ydata,"o-",ms=4)[0])
                    leg_lab.append(txt_leg)
            else:
                print "Skipping this case."
            print
        # make figure if had at least one datapoint
        if len(leg)>0:
            ax.legend(leg,leg_lab,prop={'size':4})
            ax.set_xlabel("degaussw")
            ax.set_ylabel("lambda")
            ax.set_xlim(0.0)#,0.03)
            ax.set_ylim(0.0,2.0)
            fig.savefig("_figs/converge__"+sub_fold.replace("/","_")+".pdf")
    
if __name__=="__main__":
    main()

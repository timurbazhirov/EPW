#!/usr/bin/env python

import numpy as np
import sys

def get_lambda(fname,reshape=True):
    # open entire file
    f=open(fname,"r")
    ln=f.readlines()
    f.close()
    
    # get number of q points in the fine mesh
    for i,l in enumerate(ln):
        if len(l)>44:
            if l[:44]=="     Size of q point mesh for interpolation:":
                num_qf=int(l.split()[-1])
#    print "Number of phonons in the fine mesh: ",num_qf
    
    # find number of phonon branches
    num_branch=None
    for i,l in enumerate(ln):
        if len(l)>20:
            if l[:20]=="     lambda( sum )= ":
                num_branch=int(ln[i-1].split()[1])
#    print "Number of phonon branches: ",num_branch

    # find how many kinds of smearings are there
    all_smears=[]
    for i,l in enumerate(ln):
        if "ismear:" in l:
            smear=int(l.split("ismear:")[1])
            if smear not in all_smears:
                all_smears.append(smear)
    # in case lambda was not run but only nesting then need to count differently
    if len(all_smears)==0:
        tmp_sm=1
        for i,l in enumerate(ln):
            if "Gaussian Broadening:" in l:
                all_smears.append(tmp_sm)
                tmp_sm+=1
    all_smears.sort()
    num_smears=len(all_smears)
    if all_smears!=range(1,num_smears+1):
        print "Inconsistent smears! ",all_smears
        sys.exit(1)
    
    # get lambda gamma and omega and q vectors coordinates and weights
    if num_branch!=None:
        val_lambda=np.zeros((num_qf,num_branch,num_smears),dtype=float)
        val_gamma=np.zeros((num_qf,num_branch,num_smears),dtype=float)
        val_omega=np.zeros((num_qf,num_branch,num_smears),dtype=float)
    val_nesting=np.zeros((num_qf,num_smears),dtype=float)
    q_vectors=np.zeros((num_qf,3),dtype=float)
    q_vectors_weight=np.zeros((num_qf),dtype=float)
    #
    for i,l in enumerate(ln):
        if len(l)>76:
            if l[:76]=="     Phonon (Imaginary) Self-Energy in the Migdal Approximation (on the fly)":
                start_ind=i
                break
        if len(l)>48:
            if l[:48]=="     Nesting Function in the double delta approx":
                start_ind=i
                break
    i=start_ind
    cur_smear=0 # smearings for Nesting function case
    while i<len(ln):
        # go over all q points
        if len(ln[i])>9:
            if ln[i][:9]=="     iq =":
                # make sure this is not Nesting function part
                if ln[i+2][:22]!="      Nesting function":
                    # index of the q vector
                    iq=int(ln[i].split("=")[1].split()[0])
                    # get first q-vector coordinates in the crystal coordinates units from 0 to 1
                    q_vectors[iq-1,:]=map(float,ln[i].split(":")[1].split()[:3])
                    # get weight of this q-vector
                    q_vectors_weight[iq-1]=float(ln[i].split(":")[2].split()[0])
                    # get which smearing is this
                    smear=int(ln[i].split("ismear:")[1])
                    # now get lambda gamma and omega
                    for j in range(num_branch):
                        val_lambda[iq-1,j,smear-1]=float(ln[i+2+j].split("=")[1].split()[0])
                        val_gamma[iq-1,j,smear-1] =float(ln[i+2+j].split("=")[2].split()[0])
                        val_omega[iq-1,j,smear-1] =float(ln[i+2+j].split("=")[3].split()[0])
                    i=i+num_branch+3
                    # break if this was the last point
                    if iq==num_qf:
                        break
                # get nesting part
                else:
                    # index of the q vector
                    iq=int(ln[i].split("=")[1].split()[0])
                    # store nesting value
                    vl=ln[i+2].split("=")[1].split()[0]
                    if vl=="*********" or vl=="*************************":
                        val_nesting[iq-1,cur_smear]=0.0
                    else:
                        val_nesting[iq-1,cur_smear]=float(vl)
                    # increase smearing index
                    cur_smear=(cur_smear+1)%num_smears
        i=i+1
    # reshape matrices if needed
    if reshape==True:
        # get mesh size, only works for uniform meshes at the moment 10x10x10 or so
        mesh_size=np.power(num_qf,1.0/3.0)
        if np.abs(round(mesh_size)-mesh_size)>1.0E-6:
            print "Not uniform mesh! ",mesh_size
            sys.exit(1)
        mesh_size=int(round(mesh_size))
    
        # reshape matrices
        if num_branch!=None:
            shp_lambda=np.zeros((mesh_size,mesh_size,mesh_size,num_branch,num_smears),dtype=float)
            shp_gamma=np.zeros((mesh_size,mesh_size,mesh_size,num_branch,num_smears),dtype=float)
            shp_omega=np.zeros((mesh_size,mesh_size,mesh_size,num_branch,num_smears),dtype=float)
        shp_nesting=np.zeros((mesh_size,mesh_size,mesh_size,num_smears),dtype=float)
        shp_weight=np.zeros((mesh_size,mesh_size,mesh_size),dtype=float)
        for j in range(num_smears):
            if num_branch!=None:
                for i in range(num_branch):
                    shp_lambda[:,:,:,i,j]=reshape_matrix(val_lambda[:,i,j],q_vectors,mesh_size)
                    shp_gamma[:,:,:,i,j] =reshape_matrix(val_gamma[:,i,j],q_vectors,mesh_size)
                    shp_omega[:,:,:,i,j] =reshape_matrix(val_omega[:,i,j],q_vectors,mesh_size)
            shp_nesting[:,:,:,j]=reshape_matrix(val_nesting[:,j],q_vectors,mesh_size)
        shp_weight[:,:,:]=reshape_matrix(q_vectors_weight[:],q_vectors,mesh_size)
        
        # return data
        ret={}
        ret["num_smears"]=num_smears
        ret["mesh_size"]=mesh_size
        ret["num_branch"]=num_branch
        if num_branch!=None:
            ret["lambda"]=shp_lambda
            ret["gamma"]=shp_gamma
            ret["omega"]=shp_omega
        ret["nesting"]=shp_nesting
        ret["weight"]=shp_weight
        return ret
    else:
        # return data
        ret={}
        ret["num_smears"]=num_smears
        ret["num_qf"]=num_qf
        ret["num_branch"]=num_branch
        if num_branch!=None:
            ret["val_lambda"]=val_lambda
            ret["val_gamma"]=val_gamma
            ret["val_omega"]=val_omega
        ret["val_nesting"]=val_nesting
        ret["q_vectors"]=q_vectors
        ret["q_vectors_weight"]=q_vectors_weight
        return ret
    
       
def reshape_matrix(mat,q_vectors,mesh_size):
    # reshape matrices, works only on fine meshes that are
    # equally dense in reduced coordinates in each direction, i.e. 10x10x10 or so
    out=np.zeros((mesh_size,mesh_size,mesh_size),dtype=float)
    idx=np.zeros((3),dtype=float)
    for i in range(q_vectors.shape[0]):
        # get index from reduced coordinate
        for j in range(3):
            idx[j]=q_vectors[i,j]*float(mesh_size)
            # check that it is really uniform mesh
            if np.abs(round(idx[j])-idx[j])>5.0E-4:
                print "Wrong mesh!! ",idx
                sys.exit(1)
        out[int(round(idx[0])),int(round(idx[1])),int(round(idx[2]))]=mat[i]
    return out

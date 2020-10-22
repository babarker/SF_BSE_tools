#!/usr/bin/env python

###########################################################################################
#                                                                                         #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np
import constants

def get_eks_from_wfn(wfnfile):

    # INPUT: wfnfile :: type(string) < string for wfn.h5 file, including path

    #wfn = "../04-wfn_flipped_spin/wfn_old.h5"
    wfn = wfnfile

    f_in = h5py.File(wfn)

    el = f_in['/mf_header/kpoints/el'][()] # Dims(1): mnband Dims(2): nrk Dims(3): nspin
    ifmin = f_in['/mf_header/kpoints/ifmin'][()] # Dims(1): nrk, Dims(2): nspin
    ifmax = f_in['/mf_header/kpoints/ifmax'][()] # Dims(1): nrk, Dims(2): nspin

    f_in.close()
    
    # hard-code which bands contribute to valence...
    #nvbb = 5
    #ncbb = 16
    nvbb = 12
    ncbb = 12
    
    ifmax1 = ifmax[0,0]
    ifmax2 = ifmax[1,0]
    
    #print(" ")
    #print("ifmax1")
    #print(ifmax1)
    #print("ifmax2")
    #print(ifmax2)
    #print(" ")

    ekv = el[0,0,ifmax1-nvbb:ifmax1]
    ekc = el[1,0,ifmax2:ifmax2+ncbb]

    #print("ekv")
    #print(ekv)
    #print("ekc")
    #print(ekc)
    #print(" ")

    return (ekv[::-1], ekc)

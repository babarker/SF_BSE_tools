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

def get_enk_from_wfnh5(wfnfile):

    # INPUT: wfnfile :: type(string) < string for wfn.h5 file, including path

    #wfn = "../04-wfn_flipped_spin/wfn_old.h5"
    wfn = wfnfile

    f_in = h5py.File(wfn)

    el = f_in['/mf_header/kpoints/el'][()] # Dims(1): mnband Dims(2): nrk Dims(3): nspin
    occ = f_in['/mf_header/kpoints/occ'][()] 
    ifmin = f_in['/mf_header/kpoints/ifmin'][()] # Dims(1): nrk, Dims(2): nspin
    ifmax = f_in['/mf_header/kpoints/ifmax'][()] # Dims(1): nrk, Dims(2): nspin

    f_in.close()
    
    return (el,occ,ifmin,ifmax)

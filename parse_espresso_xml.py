#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#                                                                                         #
# Author:                                                                       #
# Date:                                                                         #
#                                                                                         #
# Thanks to Diana Qiu, Ting Cao, Felipe da Jornada, Meng Wu                               #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np
import constants

def get_eks(nvb,ncb):

    # parse an eqp1.dat file to read in Kohn-Sham eigenvalues
    # (if needed, this can be lightly adapted to read in Eqp1 eigenvalues)
    
    # Dim 0: spin index (float; force to be int)
    # Dim 1: band index (float; force to be int)
    # Dim 2: KS eigenvalue (float), eV
    # Dim 3: Real part Eqp1, eV
    # Dim 4: Imag part Eqp1, eV

    infile = 'eqp1.dat'
    raw_data = np.loadtxt(infile,skiprows=1)
    #print raw_data

    # hard-code the ifmax values...
    ifmax1 = 128
    ifmax2 = 126

    # hard-code which bands contribute to valence...
    #nvbb = 5
    #ncbb = 16
    nvbb = 12
    ncbb = 12
    
    #print raw_data.shape
    ekv = []
    ekc = []
    for ii in range(raw_data.shape[0]):
        for iv in range(nvbb):
            if int(raw_data[ii,0]) == 1 and int(raw_data[ii,1]) == ifmax1 - nvbb + 1 + iv:
                ekv.append(raw_data[ii,2]/RYD)
        
        for ic in range(ncbb):
            if int(raw_data[ii,0]) == 2 and int(raw_data[ii,1]) == ifmax2 + ic + 1:
                ekc.append(raw_data[ii,2]/RYD)

    # reverse ekv since we have to count backward from Ef

    return (ekv[::-1],ekc)

def get_eks_from_wfn(nvb,ncb):

    wfn = "../04-wfn_flipped_spin/wfn_old.h5"

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
    
    print(" ")
    print("ifmax1")
    print(ifmax1)
    print("ifmax2")
    print(ifmax2)
    print(" ")

    ekv = el[0,0,ifmax1-nvbb:ifmax1]
    ekc = el[1,0,ifmax2:ifmax2+ncbb]

    #print("ekv")
    #print(ekv)
    #print("ekc")
    #print(ekc)
    #print(" ")

    return (ekv[::-1], ekc)

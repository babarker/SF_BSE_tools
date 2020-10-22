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
import matplotlib.pyplot as plt
import constants

def get_singlet_energy():

    # open info file in 01/static directory
    # get total energy of groundstate singlet system, in eV

    f=open('../no_spin_excitations/01-scf/static/info')
    ln=f.readlines()
    f.close()

    for ll in ln:
        if "Total       =" in ll:
            Egs = float(ll.split()[2])

    return Egs
    
def get_singlet_bse_singlet():

    # get singlet-ground-state-at-zero-degrees
    # for plotting, along with DFT energy

    fname = '../no_spin_excitations/14-absorption_singlet/eigenvalues.dat'
    data = np.loadtxt(fname)
    w_s_bse = data[:,0]

    return w_s_bse

def get_singlet_bse_triplet():

    fname = '../no_spin_excitations/14-absorption_triplet/eigenvalues.dat'
    data = np.loadtxt(fname)
    w_t_bse = data[:,0]

    return w_t_bse

#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# Execute this script in a directory with FILE: eigvec.dat
#
# run via command: python print_eigvec.py --nv_in NV_IN
# in directory with eigvec.dat pickle file.
#
# Author: Bradford A. Barker,
# Last Modified: October 29, 2020.
#
#                                                                                         #
###########################################################################################

#
#
import sys
import h5py
import numpy as np
#import constants
import pickle, os
#import parse_overlaps
#import parse_bsemat

def get_val():

    from argparse import ArgumentParser

    desc = ('Read in user-specified valence bands.')

    parser = ArgumentParser(description=desc)

    group = parser.add_argument_group('bands')
    group.add_argument('--nv_in', type=int, default=1, nargs=1,
        metavar=('nv_in'),
        help='Number of valence bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')

    args = parser.parse_args()

    return args.nv_in

def read_eigvec_data():

    # Read in Avec from pickle file
    pfile = open('eigvec.dat','rb')
    Avec = pickle.load(pfile)
    pfile.close()

    return Avec

def size_of_Avec(Avec):

    Nbse = np.shape(Avec)[0]
    Nexc = np.shape(Avec)[1]

    return Nbse, Nexc

def orb_indices(ibse,nvb):

    ic = int(ibse/nvb)
    iv = ibse%nvb

    return iv,ic

def main():

    nvb = get_val()
    nvb = nvb[0]

    Avec = read_eigvec_data()
    Nbse, Nexc = size_of_Avec(Avec)

    for iexc in range(Nexc):

        print("# iv  ic  Re(Avec)  Im(Avec)", file=open("Avec"+str(iexc)+".txt", "a"))
        
        for ibse in range(Nbse):

            iv,ic = orb_indices(ibse,nvb)
            print('{0:4d} {1:3d} {2: 8f} {3: 8f}'.format(iv,ic,np.real(Avec[ibse,iexc]),np.imag(Avec[ibse,iexc])), file=open("Avec"+str(iexc)+".txt", "a"))

    return
    
main()

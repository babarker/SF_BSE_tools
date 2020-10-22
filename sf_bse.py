#!/usr/bin/env python

###########################################################################################
#                                                                                         #
#   This python script is used for spin-flip method.                                      #
#                                                                                         #
#   Run this script in a directory with a bsemat.h5 file that was computed
#   via SF-BSE compatible BerkeleyGW.
#
#   Run by command
#                     python sf_bse.py --nv_in V --nc_in C > outfile
#   (with desired integers V and C)
#
#
#
#   Author: Bradford A. Barker                                                            #
#   Last edited: July 20, 2020                                                            #
#                                                                                         #
#   Thanks to Diana Qiu, Ting Cao, Felipe da Jornada, Zhenglu Li                          #
#                                                                                         #
###########################################################################################

#
#
import sys
import pickle, os
import h5py
import numpy as np
import constants
import parse_bsemat
# get_bsemat
#   in:  none
#   out: head, wing, body
# get_enk_from_bsemat
#   in:  none
#   out: el, occ, ifmin, ifmax
import make_hsfbse
# make_hbse
#   in:  ekv, ekc, head, wing, body
#   out: hbse
import diagonalize_hsfbse
# eigenstuff
#   in:  hbse 
#   out: weig*RYD,Avec
# order_excitations
#   in:  weig, Avec
#   out: weig_reorder, Avec_reorder
#import read_input
# get_val_and_cond
#   in: none (from command line: nv_in, nc_in)
#   out: nv_in, nc_out

def get_val_and_cond():

    from argparse import ArgumentParser

    desc = ('Read in user-specified valence and conduction bands.')
    # To-do: make optional and use whatever is in bsemat

    parser = ArgumentParser(description=desc)

    group = parser.add_argument_group('bands')
    group.add_argument('--nv_in', type=int, default=1, nargs=1,
        metavar=('nv_in'),
        help='Number of valence bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')
    group.add_argument('--nc_in', type=int, default=1, nargs=1,
        metavar=('nc_in'),
        help='Number of conduction bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')

    args = parser.parse_args()
    
    return args.nv_in, args.nc_in

def sf_bse(nv,nc):

    # read in number of valence and conduction bands from input:
    #nv, nc = read_input.get_val_and_cond()
    #nv = nv[0] ; nc = nc[0]
    
    # Parse bsemat.h5 for bsemat information as well as WFN data
    head, wing, body = parse_bsemat.get_bsemat()
    el, occ, ifmin, ifmax = parse_bsemat.get_enk_from_bsemat()

        
    # Construct the SF-BSE Hamiltonian
    # (Must first pick out appropriate sub-sections of el, head, and body)
    # (while noting that the spin ordering is in fact reversed!)
    
    # Get the valence band energies from "spin down" (the spin was reversed)
    ekv = el[1,0,ifmax[1,0]-nv:ifmax[1,0]]
    print("ekv before reversing order: ")
    print(ekv)
    print()
    # BAB_debug: reverse array?
    ekv = ekv[::-1]
    print("ekv after reversing order: ")
    print(ekv)
    print()
    # Get the conduction band energies from "spin up" (the spin was reversed)
    ekc = el[0,0,ifmax[0,0]:ifmax[0,0]+nc]
    # Pick out the subsections of BSE Kernel head and body:
    # N.B.: these sub-matrices now only have shape [nc,nc,nv,nv]
    sub_head = head[0,0,0:nc,0:nc,0:nv,0:nv,0]
    sub_wing = wing[0,0,0:nc,0:nc,0:nv,0:nv,0]
    sub_body = body[0,0,0:nc,0:nc,0:nv,0:nv,0]

    # Now we finally construct the SF-BSE Hamiltonian:    
    hbse = make_hsfbse.make_hsfbse(ekv, ekc, sub_head, sub_wing, sub_body)

    # Diagonalize:
    weig,Avec = diagonalize_hsfbse.eigenstuff(hbse)
    # Reorder:
    # BAB edit: let us not reorder, at first...
    weig_ro, Avec_ro = diagonalize_hsfbse.order_excitations(weig,Avec)
    #weig_ro = weig
    #Avec_ro = Avec
    # end BAB edit
    
    return weig_ro, Avec_ro

def main():

    # BAB note: erase this block
    nv_in, nc_in = get_val_and_cond()
    nv_in = nv_in[0] ; nc_in = nc_in[0]
    # Debug
    #print(nv_in,nc_in)
    
    #weig_ro, Avec_ro = sf_bse(nv_in[0], nc_in[0])
    weig_ro, Avec_ro = sf_bse(nv_in,nc_in)
    # Debug
    print("weig_ro")
    print(weig_ro)
    print()
    print("Avec_ro")
    print(Avec_ro)
    print()
    
    # Pickle the energy eigenvalues
    energy_file = open('weig.dat','wb')
    pickle.dump(weig_ro,energy_file)
    energy_file.close()    

    # Pickle the eigenvectors
    
    #Lets also print out each eigenvector, so I know better how to read them in...
    print("Avecs, by column: ")
    print()
    col = np.shape(Avec_ro)[1]
    for icol in range(col):
        print(Avec_ro[:,icol])
        print()
    print()
    print("Avecs, as matrix: ")
    print(Avec_ro)
    print()
    print()
    # end debug
    
    eigvec_file = open('eigvec.dat','wb')
    pickle.dump(Avec_ro,eigvec_file)
    eigvec_file.close()    

    return

main()

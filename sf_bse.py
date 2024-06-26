#!/usr/bin/env python

###########################################################################################
#                                                                                         #
#   This python script is used for spin-flip method.                                      #
#                                                                                         #
#   Run this script in a directory with a bsemat.h5 file that was computed
#   via SF-BSE compatible BerkeleyGW.
#
#   Run by command
#       python sf_bse.py --wfn wfn.h5 --nv_in V --nc_in C (--eig_type qp or ks --kernel on or off) > outfile
#   (with wavefunction file name wfn.h5 in hdf5 format, desired integers V and C,
#      and optional choice of eigenvalue type: quasiparticle or Kohn-Sham)
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
from constants import Ry2eV as RYD
import parse_bsemat
# get_bsemat
#   in:  none
#   out: head, wing, body
# get_enk_from_bsemat
#   in:  none
#   out: el, occ, ifmin, ifmax
import parse_wfnh5
# get_enk_from_wfn
#   in:  wfn file name (string)
#   out: el, occ, ifmin, ifmax
import make_hsfbse
# make_hbse
#   in:  ekv, ekc, head, wing, body
#   out: hbse
# make_hbse_no_kernel
#   in:  ekv, ekc
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

def get_input():

    from argparse import ArgumentParser

    desc = ('Read in user-specified valence and conduction bands.')
    # To-do: make optional and use whatever is in bsemat

    parser = ArgumentParser(description=desc)

    group = parser.add_argument_group('bands')
    group.add_argument('--wfn', type=str, default='wfn.h5', nargs=1,
        metavar=('wfn'),
        help='Filename for h5 format wavefunction file. Defaults to wfn.h5.')
    group.add_argument('--nv_in', type=int, default=1, nargs=1,
        metavar=('nv_in'),
        help='Number of valence bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')
    group.add_argument('--nc_in', type=int, default=1, nargs=1,
        metavar=('nc_in'),
        help='Number of conduction bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')
    group.add_argument('--eig_type', type=str, default='ks', nargs=1,
        metavar=('eig_type'),
        help='Specify ks for Kohn-Sham eigenvalues or qp for Quasiparticle eigenvalues. Defaults to ks.')
    group.add_argument('--kernel', type=str, default='on', nargs=1,
        metavar=('kernel'),
        help='Specify off to ignore BSE Kernel. Defaults to on.')

    args = parser.parse_args()
    
    return args.wfn, args.nv_in, args.nc_in, args.eig_type, args.kernel

def sf_bse(wfn,nv,nc,eig_type,kernel_off):

    # read in number of valence and conduction bands from input:
    #nv, nc = read_input.get_val_and_cond()
    #nv = nv[0] ; nc = nc[0]
    
    # Parse bsemat.h5 for bsemat information as well as WFN data
    head, wing, body = parse_bsemat.get_bsemat()
    #el, occ, ifmin, ifmax = parse_bsemat.get_enk_from_bsemat()
    # We will now parse wfn.h5 files for WFN data
    el, occ, ifmin, ifmax = parse_wfnh5.get_enk_from_wfnh5(wfn)

    if eig_type=='qp' :
        el = read_eqp()
        
    # Construct the SF-BSE Hamiltonian
    # (Must first pick out appropriate sub-sections of el, head, and body)
    # (while noting that the spin ordering is in fact reversed!)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Get the valence band energies from "spin down" (the spin was reversed)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #ekv = el[1,0,ifmax[1,0]-nv:ifmax[1,0]]
    ekv = el[0,0,int(ifmax[0,0])-nv:int(ifmax[0,0])]
    print("ekv before reversing order: ")
    print(ekv)
    print()
    # BAB_debug: reverse array?
    ekv = ekv[::-1]
    print("ekv after reversing order: ")
    print(ekv)
    print()
    # Get the conduction band energies from "spin up" (the spin was reversed)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #ekc = el[0,0,ifmax[0,0]:ifmax[0,0]+nc]
    ekc = el[1,0,int(ifmax[1,0]):int(ifmax[1,0])+nc]
    # Pick out the subsections of BSE Kernel head and body:
    # N.B.: these sub-matrices now only have shape [nc,nc,nv,nv]
    sub_head = head[0,0,0:nc,0:nc,0:nv,0:nv,0]
    sub_wing = wing[0,0,0:nc,0:nc,0:nv,0:nv,0]
    sub_body = body[0,0,0:nc,0:nc,0:nv,0:nv,0]

    # Now we finally construct the SF-BSE Hamiltonian:    
    if (kernel_off=='off'):
        hbse = make_hsfbse.make_hsfbse_no_kernel(ekv, ekc)
    else:
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

def read_eqp():

    # When the --eig_type qp flag is used,
    # need to have the eqp1.dat file in directory where this script is executed
    # and a number of Eqp matrix elements that exceed
    # ifmax_spin_up - nv + 1 for the occupied, spin-up channel, and
    # nc for the unoccupied, spin-down channel
    
    
    ff = open('eqp1.dat','r')
    fline = ff.readlines()[1:]

    ibmin = int(fline[0].split()[1])
    ibmax = int(fline[-1].split()[1])

    # We pad the array with zeros from (0,ibmax-1),
    # and assume that the user will not ask for more valence bands
    # than available in the eqp1.dat file.
    el = np.zeros((2,1,ibmax))
    #el = np.zeros((2,1,ibmax-ibmin+1))
    
    for ii, line in enumerate(fline):
        spin = int(line.split()[0])
        # For some confusing reason, I have spins up and down reversed.
        # This makes spin-up: 1 --> 1; spin-down: 2 --> 0
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #spin=spin%2
        spin=spin-1
        band = int(line.split()[1])
        eks = float(line.split()[2])
        eqp = float(line.split()[3])
        print(spin)
        print(band)
        print(eqp)
        #el[spin,0,ii+ibmin-1]
        #ind = ii%(ibmax-ibmin+1)
        #print(ind)
        #el[spin,0,ind+ibmin-1] = eqp/RYD
        #el[spin,0,band-ibmin] = eqp/RYD
        el[spin,0,band-1] = eqp/RYD

    print(el[0,0,:])
    print(el[1,0,:])
    
    ff.close()
    
    return el

def main():

    # BAB note: erase this block
    wfn, nv_in, nc_in, eig_type, kernel = get_input()
    wfn = wfn[0] ; nv_in = nv_in[0] ; nc_in = nc_in[0] ; eig_type = eig_type[0] ; kernel = kernel[0]
    # Debug
    #print(nv_in,nc_in)
    
    #weig_ro, Avec_ro = sf_bse(nv_in[0], nc_in[0])
    weig_ro, Avec_ro = sf_bse(wfn,nv_in,nc_in,eig_type,kernel)
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

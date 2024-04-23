#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# Execute this script in a directory with
# 1) Data Files: eigvec.dat, overlaps.dat, and bsemat.h5
# 2) Scripts: constants.py, parse_overlaps.py, parse_bsemat.py
#
# New Method based on Ipatov, original work
#
# To run, execute python ref_s_sq.py --wfn wfn.h5 --nv_in V --nc_in C
#
#
# Author: Bradford A. Barker,
# Last Modified: June 26, 2022.
#
#                                                                                         #
###########################################################################################

#
#
import sys
import h5py
import numpy as np
import constants
import pickle, os
import parse_overlaps
import parse_bsemat
import parse_wfnh5
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
    group.add_argument('--wfn', type=str, default='wfn.h5', nargs=1,
        metavar=('wfn'),
        help='Filename for h5 format wavefunction file. Defaults to wfn.h5.')
    group.add_argument('--nv_in', type=int, default=1, nargs=1,
        metavar=('nv_in'),
        help='Number of valence bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')
    group.add_argument('--nc_in', type=int, default=1, nargs=1,
        metavar=('nc_in'),
        help='Number of conduction bands requested for construction of SF-BSE Hamiltonian. Defaults to 1.')

    args = parser.parse_args()

    return args.wfn, args.nv_in, args.nc_in

def read_data(wfn,nvb,ncb):

    # This will be called in function calculate_s2
    
    # First read in ifmin and ifmax via parse_bsemat,
    # for determining set of occupied orbitals for each excitation.
    # We need ifmin and ifmax for parse_overlaps, as well.
    #

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #el, occ, ifmin, ifmax = parse_bsemat.get_enk_from_bsemat()
    el, occ, ifmin, ifmax = parse_wfnh5.get_enk_from_wfnh5(wfn)

    # Read in the overlap matrix information by calling the appropriate function
    reoarray, imoarray, osqarray = parse_overlaps.parse_overlaps(ifmin,ifmax,nvb,ncb)

    # Read in Avec from pickle file
    pfile = open('eigvec.dat','rb')
    Avec = pickle.load(pfile)
    pfile.close()

    return el, ifmin, ifmax, reoarray, imoarray, osqarray, Avec

def determine_case(iv,ic,ivp,icp):

    # This will be called in function calculate_s2
    
    # for a pair of excitations, determine whether they are the same, differ by one set of orbitals, or differ by two
    # call case_one, case_two_a, case_two_b, or case_three, as appropriate

    if iv == ivp and ic == icp:
        cflag = "one"
    elif iv != ivp and ic == icp:
        cflag = "two_a"
    elif iv == ivp and ic != icp:
        cflag = "two_b"
    elif iv != ivp and ic != icp:
        cflag = "three"
        
    return cflag

def occ_ref(ifmin,ifmax,nvb):

    # This will be called in function calculate_s2

    # initialize list of occupied bands, given the occupations of the high-spin ref. state

    # The spin order from BSE Kernel is actually swapped    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #ifmina = ifmin[1,0] # no need to worry about k-points for SF-BSE
    #ifminb = ifmin[0,0]
    #ifmaxa = ifmax[1,0]
    #ifmaxb = ifmax[0,0]
    ifmina = int(ifmin[0,0]) # no need to worry about k-points for SF-BSE
    ifminb = int(ifmin[1,0])
    ifmaxa = int(ifmax[0,0])
    ifmaxb = int(ifmax[1,0])

    # What we actually need is to create a list with the user-specified number_valence_bands,
    # mapped to low-to-high indexing (from BSE indexing),
    # and then perform the following check to determine the number of occupied spin-down bands:
    # The difference (ifmaxa - ifmina) - (ifmaxb - ifminb) then gets subtracted from Nv.
    # E.g., Ethylene has ifmaxa = 7, ifmina = 1, ifmaxb =5, ifminb = 1:
    #       (7-1) - (5-1) = 2. For Nv > 2, we compute Nv - 2 (when Nv = 2, we return an null list; Nv<2 is an error);
    # This is the number of entries in the list of occupied down-spin states, with band index up to and including ifmaxb.

    occ_ref_alpha = np.asarray(list(range(ifmaxa+1-nvb,ifmaxa+1)))

    excess_spin = (ifmaxa - ifmina) - (ifmaxb - ifminb)

    n_occ_beta = nvb - excess_spin
    
    if n_occ_beta < 0:
        quit('Need more valence bands!')
    elif n_occ_beta == 0:
        occ_ref_beta = []
    elif n_occ_beta > 0:
        occ_ref_beta = list(range(ifmaxb+1-n_occ_beta,ifmaxb+1))

    return occ_ref_alpha, occ_ref_beta

def occupied_alpha(occ_ref_alpha,iv):

    # This will be called in a loop over row/col excitation vectors in the function calculate_s2

    # for a given excitation vector, determine the full set of occupied up-spin states
    # return  a list of the occupied states for this excitation
    
    # N.B.: iv should *not* be referenced from E_{Fermi} for this application:
    #       If ifmaxa = 10 and iv=1 --> new_iv = 10.
    #       If ifmaxa = 10 and iv=2 --> new_iv = 9, etc.
    #       I.e., new_iv = ifmaxa - iv + 1

    # The above comment is handled in `calculate...` function    
    #new_iv = occ_ref_alpha[-1] - iv + 1
    
    occ_alpha = []
    for aa in occ_ref_alpha:
        if aa != iv:
            occ_alpha.append(aa)
    
    return np.asarray(occ_alpha)

def occupied_beta(occ_ref_beta,ic):

    # This will be called in a loop over row/col excitation vectors in the function calculate_s2

    # for a given excitation vector, determine the full set of occupied down-spin states
    #return  a list of the occupied states for this excitation

    # N.B.: ic should *not* be referenced from E_{Fermi} for this application:
    #       If ifmaxb = 10 and ic=1 --> new_ic = 11.
    #       If ifmaxb = 10 and ic=2 --> new_ic = 12, etc.
    #       I.e., new_ic = ifmaxb + ic
    
    # The above comment is handled in `calculate...` function    
    #new_ic = occ_ref_beta[-1] + ic
    
    occ_beta = occ_ref_beta[:]
    occ_beta.append(ic)
    
    return np.asarray(occ_beta)

def compute_overlaps(occ_ref_alpha,occ_ref_beta,osqarray):

    # Loewdin formula is M_S ( M_S + 1 ) - sum_{m up, n down} | < m | n > |^2
    # so we need to evaluate the rightward portion (the left, static factor)
    
    temp = 0.0
    for mm in occ_ref_alpha:
        for nn in occ_ref_beta:
            temp += osqarray[mm-1,nn-1]
            print(temp)
                
    return -temp

def make_slater(el,occ_alpha,occ_beta):

    # order states in slater determinant list from low energy to high energy
    # Careful: The spin order from BSE Kernel is actually swapped,
    # for using spin in el
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Start with spin-up states
    slater = []    
    for ii in occ_alpha:
        slater.append([ii,0])

    # Now populate with down-spin states, ordered by energy
    for jj in occ_beta:
        for ii in range(len(slater)):
            #if el[0,0,jj-1] < el[1,0,slater[ii][0]-1]:
            if el[1,0,jj-1] < el[0,0,slater[ii][0]-1]:
                new_slater = slater[:ii]+[[jj,1]]+slater[ii:]
                slater = new_slater
                break
            if ii == len(slater)-1 :
                new_slater = slater+[[jj,1]]
                slater = new_slater
                
    return slater    

def slater_phase(el,sl1,sl2,alpha_orb_idx,beta_orb_idx):

    # each slater determinant has form
    # [ [ib1, is1 ] , [ib2, is2] , ...  ]  
    
    slen = len(sl1)
    if len(sl2) != slen:
        quit('Slater determinants are of unequal size')
    
    # If beta_orb_idx is a null list, we are in Case 2A,
    #  meaning, sl1 has an orbital [ alpha_orb_idx[0] , 0 ]
    #  and sl2 has an orbital [ alpha_orb_idx[1] , 0 ]
    #  that need to be placed in the first position,
    #  with a phase contribution of (-1)**(position_of_orbital_in_det - 1)
    #  N.B.: the index of these orbitals is position_of_orbital_in_det - 1

    if (len(alpha_orb_idx) != 0 and len(beta_orb_idx) == 0):
     
        pos1a = sl1.index([alpha_orb_idx[0] , 0])
        pos2a = sl2.index([alpha_orb_idx[1] , 0])

        sl_phase = (-1)**(pos1a - pos2a)
    
    # If alpha_orb_idx is a null list, we are in Case 2B,
    #  meaning, sl1 has an orbital [ beta_orb_idx[0] , 1 ]
    #  and sl2 has an orbital [ beta_orb_idx[1] , 1 ]
    #  that need to be placed in the first position,
    #  with a phase contribution of (-1)**(position_of_orbital_in_det - 1)

    elif ( len(beta_orb_idx) != 0 and len(alpha_orb_idx) == 0 ):
    
        pos1b = sl1.index([beta_orb_idx[0] , 1])
        pos2b = sl2.index([beta_orb_idx[1] , 1])

        sl_phase = (-1)**(pos1b - pos2b)
    
    # If neither alpha_ nor beta_orb_idx are null lists, we are in Case 3,
    #  meaning, sl1 has orbital [alpha_orb_idx[0] , 0 ]
    #  and sl2 has orbital [alpha_orb_idx[1], 0 ]
    #  that need to be placed in the first position, 
    #  with a phase contribution of (-1)**(position_of_orbital_in_det - 1)
    #  AND THEN
    #  sl1 has orbital [ beta_orb_idx[0] , 1 ]
    #  and sl2 has an orbital [ beta_orb_idx[1] , 1 ]
    #  that need to be placed in the second position,
    #  with a phase contribution of (-1)**(position_of_orbital_in_det - 2)
    
    elif ( len(beta_orb_idx) != 0 and len(alpha_orb_idx) != 0 ):

        # code here    
        pos1a = sl1.index([alpha_orb_idx[0] , 0])
        pos2a = sl2.index([alpha_orb_idx[1] , 0])

        sl_phase = (-1)**(pos1a - pos2a)

        # two cases for the down-spins:
        # if pos1b < pos1a, then the beta orbital does not see the alpha orbital
        # if pos1b > pos1a, there is one fewer orbital for the beta orbital to permute with
        
        pos1b = sl1.index([beta_orb_idx[0] , 1])
        pos2b = sl2.index([beta_orb_idx[1] , 1])

        if (pos1b > pos1a):
            pos1b -= 1
        if (pos2b > pos2a):
            pos2b -= 1        

        sl_phase = sl_phase*(-1)**(pos1b - pos2b)

        # BAB: Oct. 4, 2021
        # We need to impose a check on the degeneracy of orbitals, across spin channels.
        # Relevant especially for NV-minus center.
        # Consider Restricted Kohn-Sham case, with Slater 1 | e_x up, e_y down > and
        # Slater 2 | e_y up, e_x down > .
        # Without the below correction, <S2> Ipatov formula will give an overall contribution of -1.
        # Need reordering in Slater 2 to change order, | e_x down, e_y up >.
        # This also holds even in Unrestricted case, despite energy being higher for e_x down.
        
        # For each determinant, locate degenerate space of beta_orb_idx,
        # check if alpha_orb_idx is in that space,
        # and multiply sl_phase by -1 if beta_orb_idx is smaller than alpha_orb_idx.

        for jj in [0,1]:
            deg = []    
            for ii in list(range(len(el[1,0,:]))):
                if (abs(el[1,0,ii] - el[1,0,beta_orb_idx[jj]-1]) < 0.000001):
                    deg.append(ii+1) # since the list above for el is indexed from one; see correction above, too.
            print('deg')
            print(deg)
            if alpha_orb_idx[jj] in deg and alpha_orb_idx[jj] > beta_orb_idx[jj]:
                sl_phase *= -1

    return sl_phase
            
def static_factor(occ_ref_alpha,occ_ref_beta):

    # M_S(High_Spin Ref State) * (M_S + 1) 
    
    N_a = int(np.shape(occ_ref_alpha)[0])
    N_b = int(np.shape(occ_ref_beta)[0])
    M_s = float(N_a - N_b)/2.0
    
    return M_s*(M_s + 1.0) + N_b

def calculate_s2(wfn,nvb,ncb):

    # call read_data
    el, ifmin, ifmax, reoarray, imoarray, osqarray, Avec = read_data(wfn,nvb,ncb)
    #ifmax_v = ifmax[1,0]
    #ifmax_c = ifmax[0,0]
    ifmax_v = int(ifmax[0,0])
    ifmax_c = int(ifmax[1,0])

    # create list of occupied orbitals in the high-spin reference state
    occ_ref_alpha, occ_ref_beta = occ_ref(ifmin,ifmax,nvb)

    # BAB debug:
    print('occ_ref_alpha')
    print(occ_ref_alpha)
    print('occ_ref_beta')
    print(occ_ref_beta)
    print(' ')
    print(' ')


    # initialize <S^2>:
    s2 = 0.0
    
    # compute sum of overlaps between occupied up and down channels:
    s2 = compute_overlaps(occ_ref_alpha,occ_ref_beta,osqarray)
    
    # now tack on the static factor
    sfac = static_factor(occ_ref_alpha,occ_ref_beta)
    print("static factor: ")
    print(sfac)
    print()
    s2 += sfac
    
    return s2

def main():

    wfn, nvb, ncb = get_val_and_cond()
    wfn = wfn[0] ; nvb = nvb[0] ; ncb = ncb[0]

    s2 = calculate_s2(wfn,nvb,ncb)

    # Print the final result!
    print("Estimate for <S^2> for excitations, ordered by energy: ")
    print(s2)
    print()

    return
    
main()

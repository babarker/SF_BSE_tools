#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# Execute this script in a directory with
# 1) Data Files: eigvec.dat, overlaps.dat, and bsemat.h5
# 2) Scripts: constants.py, parse_overlaps.py, parse_bsemat.py
#
# Method based on J. Chem. Phys. 102, 3477 (1995); https://doi.org/10.1063/1.468585
#
# To run, execute python calculate_s_sq.py --wfn wfn.h5 --nv_in V --nc_in C
#
#
# Author: Bradford A. Barker,
# Last Modified: October 26, 2020.
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

def case_one(occ_ref_alpha,occ_ref_beta,iv,ic,osqarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - sum_{m in occ_alpha} sum_{n in occ_beta} |\Delta_{m,n}|^2

    occ_alpha = occupied_alpha(occ_ref_alpha,iv)
    occ_beta = occupied_beta(occ_ref_beta,ic)

    print('case one')
    print('occ_alpha')
    print(occ_alpha)
    print('occ_beta')
    print(occ_beta)

    temp = complex(0.0)
    for mm in occ_alpha:
        for nn in occ_beta:
            print('osqarray')
            print(osqarray[mm-1,nn-1])
            temp = temp + osqarray[mm-1,nn-1]

    # Need the negative of the above
    temp = -temp

    # now tack on the static factor
    sfac = static_factor(occ_ref_alpha,occ_ref_beta)
    print('static factor')
    print(sfac)
    print()
    temp += sfac
    print()

    return temp

def case_two_a(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{iv-ivp} sum_{n in occ_beta} \Delta_{iv,n}\Delta^*_{ivp,n}

    occ_beta = occupied_beta(occ_ref_beta,ic)
    occ_alpha1 = occupied_alpha(occ_ref_alpha,iv)
    occ_alpha2 = occupied_alpha(occ_ref_alpha,ivp)
    alpha_orb_idx = [x for x in occ_alpha1 if x not in set(occ_alpha2)]
    alpha_orb_idx = alpha_orb_idx + [x for x in occ_alpha2 if x not in set(occ_alpha1)]
    # define beta_orb_idx as null list, for use in slater determinant phase calculation
    beta_orb_idx = []

    print('case two a')
    print('occ_alpha1')
    print(occ_alpha1)
    print('occ_alpha2')
    print(occ_alpha2)
    print('alpha_orb_idx')
    print(alpha_orb_idx)
    print()
    print('occ_beta')
    print(occ_beta)
    print()

    temp = complex(0.0)
    for nn in occ_beta:
        # BAB debug
        print('nn, iv, ivp')
        print(nn, iv, ivp)
        print('Real part of second overlap matrix:')
        print(reoarray[alpha_orb_idx[0]-1,nn-1])
        print('Real part of first overlap matrix:')
        print(reoarray[alpha_orb_idx[1]-1,nn-1])
        print()
        # end BAB debug
        temp = temp + reoarray[alpha_orb_idx[0]-1,nn-1]*reoarray[alpha_orb_idx[1]-1,nn-1] 
        temp = temp + imoarray[alpha_orb_idx[0]-1,nn-1]*imoarray[alpha_orb_idx[1]-1,nn-1]
        temp = temp - 1j*reoarray[alpha_orb_idx[0]-1,nn-1]*imoarray[alpha_orb_idx[1]-1,nn-1] 
        temp = temp + 1j*imoarray[alpha_orb_idx[0]-1,nn-1]*reoarray[alpha_orb_idx[1]-1,nn-1]

    sl1 = make_slater(el,occ_alpha1,occ_beta)
    sl2 = make_slater(el,occ_alpha2,occ_beta)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(el,sl1,sl2,alpha_orb_idx,beta_orb_idx)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp

    return -temp

def case_two_b(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{ic-icp} sum_{m in occ_alpha} \Delta_{m,icp}\Delta^*_{m,ic}

    occ_alpha = occupied_alpha(occ_ref_alpha,iv)
    occ_beta1 = occupied_beta(occ_ref_beta,ic)
    occ_beta2 = occupied_beta(occ_ref_beta,icp)
    beta_orb_idx = [x for x in occ_beta1 if x not in set(occ_beta2)]
    beta_orb_idx = beta_orb_idx + [x for x in occ_beta2 if x not in set(occ_beta1)]
    # initialize alpha_orb_idx as null list for slater determinant phase calculation
    alpha_orb_idx = []

    print('case two b')
    print('occ_beta1')
    print(occ_beta1)
    print('occ_beta2')
    print(occ_beta2)
    print()
    print('occ_alpha')
    print(occ_alpha)
    print()

    temp = complex(0.0)
    for mm in occ_alpha:
        # BAB debug
        print('mm,ic,icp')
        print(mm, ic, icp)
        print('Real part of second overlap matrix:')
        print(reoarray[mm-1,ic-1])
        print('Real part of first overlap matrix:')
        print(reoarray[mm-1,icp-1])
        print()
        # end BAB debug
        temp = temp + reoarray[mm-1,icp-1]*reoarray[mm-1,ic-1] + imoarray[mm-1,icp-1]*imoarray[mm-1,ic-1]
        temp = temp - 1j*reoarray[mm-1,icp-1]*imoarray[mm-1,ic-1] + 1j*imoarray[mm-1,icp-1]*reoarray[mm-1,ic-1]

    sl1 = make_slater(el,occ_alpha,occ_beta1)
    sl2 = make_slater(el,occ_alpha,occ_beta2)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(el,sl1,sl2,alpha_orb_idx,beta_orb_idx)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp

    return -temp

def case_three(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):
    
    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{iv-ivp} (-1)^{ic-icp} \Delta_{iv,icp}\Delta^*_{ivp,ic}

    # Solution to problem of comparing two lists and determining the different elements,
    # used in this case to create a list of the different orbital indices for the pair
    # of Slater determinants, found from
    # https://stackoverflow.com/questions/3462143/get-difference-between-two-lists

    occ_alpha1 = occupied_alpha(occ_ref_alpha,iv)
    occ_alpha2 = occupied_alpha(occ_ref_alpha,ivp)
    occ_beta1 = occupied_beta(occ_ref_beta,ic)
    occ_beta2 = occupied_beta(occ_ref_beta,icp)

    print('case three')
    print('occ_alpha1')
    print(occ_alpha1)
    print('occ_alpha2')
    print(occ_alpha2)
    print('occ_beta1')
    print(occ_beta1)
    print('occ_beta2')
    print(occ_beta2)

    # BAB debug: output list of the different up-spin orbital, in order [orbital_alpha_no_prime,orbital_alpha_prime]
    alpha_orb_idx = [x for x in occ_alpha1 if x not in set(occ_alpha2)]
    alpha_orb_idx = alpha_orb_idx + [x for x in occ_alpha2 if x not in set(occ_alpha1)]

    # BAB debug: output list of the different down-spin orbital, in order [orbital_beta_no_prime,orbital_beta_prime]
    beta_orb_idx = [x for x in occ_beta1 if x not in set(occ_beta2)]
    beta_orb_idx = beta_orb_idx + [x for x in occ_beta2 if x not in set(occ_beta1)]

    # BAB debug: Use the above lists
    print('alpha_orb_idx')
    print(alpha_orb_idx)
    print('beta_orb_idx')
    print(beta_orb_idx)
    print('Real part of first overlap matrix:')
    print(reoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1])
    print('Real part of second overlap matrix:')
    print(reoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1])
    print()

    # BAB debug: do not use iv, ivp directly! Use the alpha_orb_idx, beta_orb_idx lists to keep track
    #          of the different sets of orbitals from the two different determinants
    #temp = reoarray[iv-1,icp-1]*reoarray[ivp-1,ic-1] + imoarray[iv-1,icp-1]*imoarray[ivp-1,ic-1]
    #temp = temp - 1j*reoarray[iv-1,icp-1]*imoarray[ivp-1,ic-1] + 1j*imoarray[iv-1,icp-1]*reoarray[ivp-1,ic-1]
    #temp = (-1.0)**(iv-ivp)*(-1.0)**(ic-icp)*temp
    temp = complex(reoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*reoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1]) 
    temp = complex(temp + imoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*imoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1])
    temp = complex(temp - 1j*reoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*imoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1]) 
    temp = complex(temp + 1j*imoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*reoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1])

    # BAB debug: use slater determinant method...
    sl1 = make_slater(el,occ_alpha1,occ_beta1)
    sl2 = make_slater(el,occ_alpha2,occ_beta2)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(el,sl1,sl2,alpha_orb_idx,beta_orb_idx)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp

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

    # This computes the first few terms for <S^2> according to Wang 1995:
    # static_factor = (N_a - N_b)/2 * ((N_a - N_b/2)+1) + N_b,
    # and <S^2> = static_factor + 2\int \Gamma^{abba} (r_1,r_2|r_1,r_2) d r_1 d r_2
    
    # Make sure that N_a, N_b is WITH TRANSITION
    N_a = int(np.shape(occ_ref_alpha)[0] - 1)
    N_b = int(np.shape(occ_ref_beta)[0] + 1)

    diff = float(N_a - N_b)
    
    return diff/2.00*(diff/2.00 + 1) + N_b

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
    s2 = np.zeros((nvb*ncb),dtype="complex")
    
    # loop through excitation number,
    # then loop through valence and conduction bands in row and col dimensions,
    # determine the CASE for each pair of excitations:    
    # N.B.: the band indices in the loop are in `BSE` index format.

    for iexc in range(nvb*ncb):
    #for iexc in range(40):
        for iv in range(nvb):
            for ic in range(ncb):
                for ivp in range(nvb):
                    for icp in range(ncb):

                        # Need to convert BSE band index read in from overlaps.dat file into RAW format band index,
                        # with ib1_raw = ifmax_v - ib1_bse (counting from 0 to Nv-1), for valence, and
                        #      ib2_raw = ib2_bse + ifmax_c + 1 (counting from 0 to Nc-1), for conduction.

                        ivr = int(ifmax_v - iv)
                        ivpr = int(ifmax_v - ivp)
                        icr = int(ic + ifmax_c + 1)
                        icpr = int(icp + ifmax_c + 1)
                 
                        cflag = determine_case(ivr,icr,ivpr,icpr)

                        ibse = ic*nvb + iv
                        ibsep = icp*nvb + ivp
                    
                        if cflag == "one":
                            temp = case_one(occ_ref_alpha,occ_ref_beta,ivr,icr,osqarray)
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                        elif cflag == "two_a":
                            temp = case_two_a(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                        elif cflag == "two_b":
                            temp = case_two_b(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                        elif cflag == "three":
                            temp = case_three(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                         
                        print(' ')
                        print('iexc')
                        print(iexc)
                        print('Avec[ibse,iexc]')
                        print(Avec[ibse,iexc])
                        print('Avec[ibsep,iexc]')
                        print(Avec[ibsep,iexc])
                        s2[iexc] += np.conj(Avec[ibsep,iexc])*temp*Avec[ibse,iexc]
                        print('s2[iexc]')
                        print("{:.6f}".format(s2[iexc]))
                        print(' ')

        # now tack on the static factor
        #sfac = static_factor(occ_ref_alpha,occ_ref_beta)
        #print('static factor')
        #print(sfac)
        #print()
        #temp += sfac
        #s2[iexc] += sfac
        #print('s2[iexc] plus static factor')
        #print(s2[iexc])
        #print()

    return s2

def main():

    wfn, nvb, ncb = get_val_and_cond()
    wfn = wfn[0] ; nvb = nvb[0] ; ncb = ncb[0]

    s2 = calculate_s2(wfn,nvb,ncb)

    # Print the final result!
    print("Estimate for <S^2> for excitations, ordered by energy: ")
    print(np.absolute(s2))
    print()

    return
    
main()

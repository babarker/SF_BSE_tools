#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# Execute this script in a directory with
# 1) Data Files: eigvec.dat, overlaps.dat, and bsemat.h5
# 2) Scripts: constants.py, parse_overlaps.py, parse_bsemat.py
#
# Method based on J. Chem. Phys. 102, 3477 (1995); https://doi.org/10.1063/1.468585
#
# To run, execute python calculate_s_sq.py --nv_in V --nc_in C
#
#
# Author: Bradford A. Barker,
# Last Modified: August 4, 2020.
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

def read_data(nvb,ncb):

    # This will be called in function calculate_s2
    
    # First read in ifmin and ifmax via parse_bsemat,
    # for determining set of occupied orbitals for each excitation.
    # We need ifmin and ifmax for parse_overlaps, as well.
    #

    el, occ, ifmin, ifmax = parse_bsemat.get_enk_from_bsemat()

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
    ifmina = ifmin[1,0] # no need to worry about k-points for SF-BSE
    ifminb = ifmin[0,0]
    ifmaxa = ifmax[1,0]
    ifmaxb = ifmax[0,0]

    #occ_ref_alpha = np.asarray(range(ifmina,ifmaxa+1))
    #occ_ref_beta = np.asarray(range(ifminb,ifmaxb+1))

    # What we actually need is to create a list with the user-specified number_valence_bands,
    # mapped to low-to-high indexing (from BSE indexing),
    # and then perform the following check to determine the number of occupied spin-down bands:
    # The difference (ifmaxa - ifmina) - (ifmaxb - ifminb) then gets subtracted from Nv.
    # E.g., Ethylene has ifmaxa = 7, ifmina = 1, ifmaxb =5, ifminb = 1:
    #       (7-1) - (5-1) = 2. For Nv > 2, we compute Nv - 2 (when Nv = 2, we return an null list; Nv<2 is an error);
    # This is the number of entries in the list of occupied down-spin states, with band index up to and including ifmaxb.

    occ_ref_alpha = np.asarray(range(ifmaxa+1-nvb,ifmaxa+1))

    excess_spin = (ifmaxa - ifmina) - (ifmaxb - ifminb)

    n_occ_beta = nvb - excess_spin
    
    if n_occ_beta < 0:
        quit('Need more valence bands!')
    elif n_occ_beta == 0:
        #occ_ref_beta = np.asarray([])
        occ_ref_beta = []
    elif n_occ_beta > 0:
        #occ_ref_beta = range(ifmaxb+1-n_occ_beta,ifmaxb+1)
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
    #    if aa != new_iv:
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
    
    #occ_beta = np.concatenate((occ_ref_beta,np.asarray([int(ic)])))
    occ_beta = occ_ref_beta[:]
    #print(occ_beta)
    #occ_beta.tolist()
    #print(occ_beta)
    occ_beta.append(ic)
    #print(occ_beta)
    #print(occ_ref_beta)
    #quit()
    #occ_beta = [ int(occ) for occ in occ_beta ]
    
    return np.asarray(occ_beta)

def case_one(occ_ref_alpha,occ_ref_beta,iv,ic,osqarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - sum_{m in occ_alpha} sum_{n in occ_beta} |\Delta_{m,n}|^2

    occ_alpha = occupied_alpha(occ_ref_alpha,iv)
    occ_beta = occupied_beta(occ_ref_beta,ic)

    print('case one')
    print('iv')
    print(iv)
    print('ic')
    print(ic)
    print()
    print('occ_alpha')
    print(occ_alpha)
    print('occ_beta')
    print(occ_beta)
    print()
    #print('osqarray')
    #print(osqarray)
    #print()

    temp = 0.0
    for mm in occ_alpha:
        for nn in occ_beta:
            # temp = temp + osqarray[mm,nn]
            # BAB debug
            print('mm,nn')
            print(mm,nn)
            print('Square of overlap matrix:')
            print(osqarray[mm-1,nn-1])
            print()
            # end BAB debug
            temp = temp + osqarray[mm-1,nn-1]
    # when the excitations are identical, use formula for identical.
    # needed: overlaps-squared for all occupied orbitals
    # loop through occupied up-spins (m), then down-spins (n), use |Delta_(mn)|**2

    # Need the negative of the above
    temp = -temp

    # now tack on the static factor

    sfac = static_factor(occ_ref_alpha,occ_ref_beta)
    print('static factor')
    print(sfac)
    print()
    temp += sfac
    #s2[iexc] += sfac
    #print('s2[iexc] plus static factor')
    #print(s2[iexc])
    print()

    #return -temp
    return temp

def case_two_a(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{iv-ivp} sum_{n in occ_beta} \Delta_{iv,n}\Delta^*_{ivp,n}

    occ_beta = occupied_beta(occ_ref_beta,ic)
    # BAB debug: need to compare list of occupied spin-up channels from each determinant:
    occ_alpha1 = occupied_alpha(occ_ref_alpha,iv)
    occ_alpha2 = occupied_alpha(occ_ref_alpha,ivp)
    # BAB debug: output list of the different up-spin orbital, in order [orbital_alpha_no_prime,orbital_alpha_prime]
    alpha_orb_idx = [x for x in occ_alpha1 if x not in set(occ_alpha2)]
    alpha_orb_idx = alpha_orb_idx + [x for x in occ_alpha2 if x not in set(occ_alpha1)]
    #alpha_orb_idx = list(set(occ_alpha1).symmetric_difference(occ_alpha2))
    # end BAB

    print('case two a')
    print('iv')
    print(iv)
    print('ivp')
    print(ivp)
    print('occ_alpha1')
    print(occ_alpha1)
    print('occ_alpha2')
    print(occ_alpha2)
    print('alpha_orb_idx')
    print(alpha_orb_idx)
    print('icp')
    print(icp)
    print('ic')
    print(ic)
    print()
    print('occ_beta')
    print(occ_beta)
    print()

    temp = 0.0
    for nn in occ_beta:
        # BAB debug
        print('nn, iv, ivp')
        print(nn, iv, ivp)
        print('Real part of second overlap matrix:')
        #print(reoarray[iv-1,nn-1])
        print(reoarray[alpha_orb_idx[0]-1,nn-1])
        print('Real part of first overlap matrix:')
        #print(reoarray[ivp-1,nn-1])
        print(reoarray[alpha_orb_idx[1]-1,nn-1])
        print()
        # end BAB debug
        #temp = temp + reoarray[iv-1,nn-1]*reoarray[ivp-1,nn-1] 
        #temp = temp + imoarray[iv-1,nn-1]*imoarray[ivp-1,nn-1]
        #temp = temp - 1j*reoarray[iv-1,nn-1]*imoarray[ivp-1,nn-1] 
        #temp = temp + 1j*imoarray[iv-1,nn-1]*reoarray[ivp-1,nn-1]
        temp = temp + reoarray[alpha_orb_idx[0]-1,nn-1]*reoarray[alpha_orb_idx[1]-1,nn-1] 
        temp = temp + imoarray[alpha_orb_idx[0]-1,nn-1]*imoarray[alpha_orb_idx[1]-1,nn-1]
        temp = temp - 1j*reoarray[alpha_orb_idx[0]-1,nn-1]*imoarray[alpha_orb_idx[1]-1,nn-1] 
        temp = temp + 1j*imoarray[alpha_orb_idx[0]-1,nn-1]*reoarray[alpha_orb_idx[1]-1,nn-1]

    # BAB debug: use slater determinant method...
    sl1 = make_slater(el,occ_alpha1,occ_beta)
    sl2 = make_slater(el,occ_alpha2,occ_beta)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(sl1,sl2)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp

    # Using phases this way is actually too difficult
    #ph1 = find_phase(el,occ_alpha1,occ_beta,alpha_orb_idx[0],0)
    #ph2 = find_phase(el,occ_alpha2,occ_beta,alpha_orb_idx[1],0)
    #print('ph1')
    #print(ph1)
    #print('ph2')
    #print(ph2)
    #print()
    #temp = (-1.0)**(ph1-ph2)*temp

    #temp = (-1.0)**(iv-ivp)*temp
    # when excitations are off by one up-spin orbital.
    # needed: the overlaps of the pair of the relevant up-spin orbitals (k,l) and the occupied down-spin orbitals
    # also needed: the sum of the `order` of the orbitals k and l
    # loop through occupied down-spins (n), use Delta_(kn), Delta(ln), according to Wang 1995.
    return -temp

def case_two_b(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):

    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{ic-icp} sum_{m in occ_alpha} \Delta_{m,icp}\Delta^*_{m,ic}

    occ_alpha = occupied_alpha(occ_ref_alpha,iv)
    # BAB debug: need to compare list of occupied spin-down channels from each determinant:
    occ_beta1 = occupied_beta(occ_ref_beta,ic)
    occ_beta2 = occupied_beta(occ_ref_beta,icp)
    beta_orb_idx = [x for x in occ_beta1 if x not in set(occ_beta2)]
    beta_orb_idx = beta_orb_idx + [x for x in occ_beta2 if x not in set(occ_beta1)]
    # end BAB

    print('case two b')
    print('iv')
    print(iv)
    print('ivp')
    print(ivp)
    print('ic')
    print(ic)
    print('icp')
    print(icp)
    print('occ_beta1')
    print(occ_beta1)
    print('occ_beta2')
    print(occ_beta2)
    print()
    print('occ_alpha')
    print(occ_alpha)
    print()

    temp = 0.0
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

    # BAB debug: use slater determinant method...
    sl1 = make_slater(el,occ_alpha,occ_beta1)
    sl2 = make_slater(el,occ_alpha,occ_beta2)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(sl1,sl2)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp

    #ph1 = find_phase(el,occ_alpha,occ_beta1,beta_orb_idx[0],1)
    #ph2 = find_phase(el,occ_alpha,occ_beta2,beta_orb_idx[1],1)
    #print('ph1')
    #print(ph1)
    #print('ph2')
    #print(ph2)
    #print()
    #temp = (-1.0)**(ph1-ph2)*temp

    #temp = (-1.0)**(ic-icp)*temp
    # when excitations are off by one down-spin orbital.
    # needed: the overlaps of the pair of the relevant down-spin orbitals (k,l) and the occupied up-spin orbitals
    # also needed: the sum of the `order` of the orbitals k and l
    # loop through occupied up-spins (m), use Delta_(mk), Delta(ml), according to Wang 1995.
    return -temp

def case_three(el,occ_ref_alpha,occ_ref_beta,iv,ic,ivp,icp,reoarray,imoarray):
    
    # Calculation of 2\int \Gamma^{\alpha\beta\beta\alpha}(r_1,r_2|r_1,r_2) dr_1 dr_2
    # = - (-1)^{iv-ivp} (-1)^{ic-icp} \Delta_{iv,icp}\Delta^*_{ivp,ic}

    # Solution to problem of comparing two lists and determining the different elements,
    # used in this case to create a list of the different orbital indices for the pair
    # of Slater determinants, found from
    # https://stackoverflow.com/questions/3462143/get-difference-between-two-lists

    # BAB debug: need to compare list of occupied spin-up channels from each determinant:
    occ_alpha1 = occupied_alpha(occ_ref_alpha,iv)
    occ_alpha2 = occupied_alpha(occ_ref_alpha,ivp)
    # end BAB
    # BAB debug: need to compare list of occupied spin-down channels from each determinant:
    occ_beta1 = occupied_beta(occ_ref_beta,ic)
    occ_beta2 = occupied_beta(occ_ref_beta,icp)
    # end BAB


    print('case three')
    print('iv')
    print(iv)
    print('ivp')
    print(ivp)
    print('occ_alpha1')
    print(occ_alpha1)
    print('occ_alpha2')
    print(occ_alpha2)
    print('ic')
    print(ic)
    print('icp')
    print(icp)
    print('occ_beta1')
    print(occ_beta1)
    print('occ_beta2')
    print(occ_beta2)

    # BAB debug: output list of the different up-spin orbital, in order [orbital_alpha_no_prime,orbital_alpha_prime]
    #alpha_orb_idx = list(set(occ_alpha1).symmetric_difference(occ_alpha2))
    alpha_orb_idx = [x for x in occ_alpha1 if x not in set(occ_alpha2)]
    alpha_orb_idx = alpha_orb_idx + [x for x in occ_alpha2 if x not in set(occ_alpha1)]

    # BAB debug: output list of the different down-spin orbital, in order [orbital_beta_no_prime,orbital_beta_prime]
    #beta_orb_idx = list(set(occ_beta1).symmetric_difference(occ_beta2))
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
    #print('Real part of first overlap matrix:')
    #print(reoarray[iv-1,icp-1])
    #print('Real part of second overlap matrix:')
    #print(reoarray[ivp-1,ic-1])
    #print()

    # BAB debug: do not use iv, ivp directly! Use the alpha_orb_idx, beta_orb_idx lists to keep track
    #          of the different sets of orbitals from the two different determinants
    #temp = reoarray[iv-1,icp-1]*reoarray[ivp-1,ic-1] + imoarray[iv-1,icp-1]*imoarray[ivp-1,ic-1]
    #temp = temp - 1j*reoarray[iv-1,icp-1]*imoarray[ivp-1,ic-1] + 1j*imoarray[iv-1,icp-1]*reoarray[ivp-1,ic-1]
    #temp = (-1.0)**(iv-ivp)*(-1.0)**(ic-icp)*temp
    temp = reoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*reoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1] 
    temp = temp + imoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*imoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1]
    temp = temp - 1j*reoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*imoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1] 
    temp = temp + 1j*imoarray[alpha_orb_idx[0]-1,beta_orb_idx[1]-1]*reoarray[alpha_orb_idx[1]-1,beta_orb_idx[0]-1]


    # BAB debug: use slater determinant method...
    sl1 = make_slater(el,occ_alpha1,occ_beta1)
    sl2 = make_slater(el,occ_alpha2,occ_beta2)
    print('Slater 1')
    print(sl1)
    print('Slater 2')
    print(sl2)
    sl_phase = slater_phase(sl1,sl2)
    print('sl_phase')
    print(sl_phase)
    print()
    temp = sl_phase*temp


    #ph1 = find_phase(el,occ_alpha1,occ_beta1,alpha_orb_idx[0],0)
    #ph2 = find_phase(el,occ_alpha2,occ_beta2,alpha_orb_idx[1],0)
    #ph3 = find_phase(el,occ_alpha1,occ_beta1,beta_orb_idx[0],1)
    #ph4 = find_phase(el,occ_alpha2,occ_beta2,beta_orb_idx[1],1)
    #print('ph1')
    #print(ph1)
    #print('ph2')
    #print(ph2)
    #print('ph3')
    #print(ph3)
    #print('ph4')
    #print(ph4)
    #print()
    #temp = (-1.0)**(ph1-ph2)*(-1.0)**(ph3-ph4)*temp

    #temp = (-1.0)**(alpha_orb_idx[0]-alpha_orb_idx[1])*(-1.0)**(beta_orb_idx[0]-beta_orb_idx[1])*temp
    # when excitations are off for both the up- and down-spin orbitals
    # needed: the overlaps of the pairs of orbitals
    # also needed: the sum of the `order` of the orbitals k and l, d and f    
    # use Delta_(mk), Delta(ml), according to Wang 1995.
    return -temp

def find_phase(el,occ_alpha,occ_beta,idx,idx_spin):

    # THIS FUNCTION SHOULD NOT BE USED; IT DOES NOT WORK

    # idx is index from one of alpha_orb_idx[0], alpha_orb_idx[1], beta_orb_idx[0], beta_orb_idx[1]
    # if idx is from alpha_orb_idx, set idx_spin == 0
    # if idx is from beta_orb_idx, set idx_spin == 1

    ph = 0    
    for ii in occ_alpha:
        if el[0,0,ii-1] < el[idx_spin,0,idx-1]:
            ph += 1

    for jj in occ_beta:
        if el[1,0,ii-1] < el[idx_spin,0,idx-1]:
            ph += 1

    # A down-spin does not go all the way to the left in the Slater determinant; that spot is for the up-spin
    # Subtract one from the down-spin position in the determinant
    if idx_spin == 1:
        ph -= 1
    if ph < 0:
        ph = 0

    return ph

def make_slater(el,occ_alpha,occ_beta):
    # order states in slater determinant list from low energy to high energy
    # Careful: The spin order from BSE Kernel is actually swapped,
    # for using spin in el

    # Start with spin-up states
    slater = []    
    for ii in occ_alpha:
        slater.append([ii,0])

    # Now check 
    for jj in occ_beta:
        for ii in range(len(slater)):
            #if el[1,0,jj-1] < el[slater[ii][1],0,slater[ii][0]-1]:
            if el[0,0,jj-1] < el[1,0,slater[ii][0]-1]:
                print('el[0,0,jj-1] in eV')
                print(el[0,0,jj-1]*13.6)
                print('el[1,0,slater[ii][0]-1] in eV')
                print(el[1,0,slater[ii][0]-1]*13.6)
                new_slater = slater[:ii]+[[jj,1]]+slater[ii:]
                slater = new_slater
                break
            if ii == len(slater)-1 :
                new_slater = slater+[[jj,1]]
                slater = new_slater
                
    return slater    

def slater_phase(sl1,sl2):
    # each slater determinant has form
    # [ [ib1, is1 ] , [ib2, is2] , ...  ]
    
    # we count permutations by tracking swapped spins, between pairs of slater determinants
    
    slen = len(sl1)
    if len(sl2) != slen:
        quit('Slater determinants are of unequal size')
    
    count = 0
    for ii in range(slen):
        if sl1[ii][1] != sl2[ii][1]:
            count += 1
    nswap = count/2
    
    sl_phase = (-1)**(nswap)
    
    return sl_phase

            
def static_factor(occ_ref_alpha,occ_ref_beta):

    # BAB note:
    # Make sure that N_a, N_b is WITH TRANSITION

    # This will be called in function calculate_s2

    # This computes the first few terms for <S^2> according to Wang 1995:
    # static_factor = (N_a - N_b)/2 * ((N_a - N_b/2)+1) + N_b,
    # and <S^2> = static_factor + 2\int \Gamma^{abba} (r_1,r_2|r_1,r_2) d r_1 d r_2
    
    N_a = int(np.shape(occ_ref_alpha)[0] - 1)
    N_b = int(np.shape(occ_ref_beta)[0] + 1)

    diff = float(N_a - N_b)
    print('N_a')
    print(N_a)
    print('N_b')
    print(N_b)
    #quit()
    
    return diff/2.00*(diff/2.00 + 1) + N_b

def calculate_s2(nvb,ncb):

    #nvb, ncb = read_input.get_val_and_cond()
    #nvb = nvb[0] ; ncb = ncb[0]

    # Find nvmax and ncmax from Avec or bsemat.h5
    #nvb, ncb = get_nvb_ncb_from_bsemat.parse_bsemat()
    #nvb, ncb = parse_bsemat.get_nvb_ncb_from_bsemat()
    #print('nvb')
    #print(nvb)
    #print('ncb')
    #print(ncb)

    # call read_data
    el, ifmin, ifmax, reoarray, imoarray, osqarray, Avec = read_data(nvb,ncb)
    #print('read data correctly')
    #quit()
    #nbmax = np.shape(reoarray)[0]
    ifmax_v = ifmax[1,0]
    ifmax_c = ifmax[0,0]

    # create list of occupied orbitals in the high-spin reference state
    occ_ref_alpha, occ_ref_beta = occ_ref(ifmin,ifmax,nvb)

    #print('occ_ref_alpha')
    #print(occ_ref_alpha)
    #print('occ_ref_beta')
    #print(occ_ref_beta)
    #print()
    
    # initialize <S^2>:
    s2 = np.zeros((nvb*ncb))
    
    # loop through excitation number,
    # then loop through valence and conduction bands in row and col dimensions,
    # determine the CASE for each pair of excitations:    
    # N.B.: the band indices in the loop are in `BSE` index format.

    for iexc in range(nvb*ncb):
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

                        #ibse = ic*len(nvb) + iv
                        #ibsep = icp*len(nvb) + ivp
                    
                        ibse = ic*nvb + iv
                        ibsep = icp*nvb + ivp
                    
                        if cflag == "one":
                            temp = case_one(occ_ref_alpha,occ_ref_beta,ivr,icr,osqarray)
                            # replace this with the correct version....
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                            #print('Avec[ibse,iexc]')
                            #print(Avec[ibse,iexc])
                            #print(' ')
                            #s2[iexc] += np.conj(Avec[ibsep,iexc])*temp*Avec[ibse,iexc]
                        elif cflag == "two_a":
                            temp = case_two_a(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            # replace this with the correct version....
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                            #print('Avec[ibse,iexc]')
                            #print(Avec[ibse,iexc])
                            #print(' ')
                            #s2[iexc] += np.conj(Avec[ibsep,iexc])*temp*Avec[ibse,iexc]
                        elif cflag == "two_b":
                            temp = case_two_b(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            # replace this with the correct version....
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                            #print('Avec[ibse,iexc]')
                            #print(Avec[ibse,iexc])
                            #print(' ')
                            #s2[iexc] += np.conj(Avec[ibsep,iexc])*temp*Avec[ibse,iexc]
                        elif cflag == "three":
                            temp = case_three(el,occ_ref_alpha,occ_ref_beta,ivr,icr,ivpr,icpr,reoarray,imoarray)
                            # replace this with the correct version....
                            print(' ')
                            print('temp')
                            print(temp)
                            print(' ')
                            #print('Avec[ibse,iexc]')
                            #print(Avec[ibse,iexc])
                            #print(' ')
                            #s2[iexc] += np.conj(Avec[ibsep,iexc])*temp*Avec[ibse,iexc]
                         
                        print(' ')
                        print('iexc')
                        print(iexc)
                        #print(s2[iexc])
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

    nvb, ncb = get_val_and_cond()
    nvb = nvb[0] ; ncb = ncb[0]

    s2 = calculate_s2(nvb,ncb)

    # Print the final result!
    print("Estimate for <S^2> for excitations, ordered by energy: ")
    print(s2)
    print()

    return
    
main()

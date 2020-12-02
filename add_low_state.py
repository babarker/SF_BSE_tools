#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#                                                                                         #
# Thanks to Diana Qiu, Ting Cao, Felipe da Jornada                               #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np

# This is from making wavefunction smaller by reducing number of bands
#nbsmall = 30


def swap_spins():
    r"""

    Read in old wavefunction, spit out new wavefunction with additional garbage state
    in lowest position, with negative-infinity energy, for the spin channel with
    fewer occupied states.

    This function is modified based on Felipe's pseudobands.py script
    """

    wfnold = "wfn_old.h5"
    wfnnew = "wfn_new.h5"

    f_in = h5py.File(wfnold)
    f_out = h5py.File(wfnnew, "w")


    f_out.create_group('mf_header')
    f_out.copy(f_in['mf_header/versionnumber'], 'mf_header/versionnumber')
    f_out.copy(f_in['mf_header/flavor'], 'mf_header/flavor')

    f_out.copy(f_in['mf_header/kpoints/nspin'], 'mf_header/kpoints/nspin')
    f_out.copy(f_in['mf_header/kpoints/nspinor'], 'mf_header/kpoints/nspinor')

    # This is from making wavefunction smaller
    #bshape = list(f_in['mf_header/kpoints/mnband'].shape)
    #f_out.create_dataset('mf_header/kpoints/mnband', bshape, 'd')
    #f_out['mf_header/kpoints/mnband'][()] = nbsmall
    
    f_out.copy(f_in['mf_header/kpoints/mnband'], 'mf_header/kpoints/mnband')

    f_out.copy(f_in['mf_header/kpoints/nrk'], 'mf_header/kpoints/nrk')
    f_out.copy(f_in['mf_header/kpoints/ngkmax'], 'mf_header/kpoints/ngkmax')
    f_out.copy(f_in['mf_header/kpoints/ecutwfc'], 'mf_header/kpoints/ecutwfc')
    f_out.copy(f_in['mf_header/kpoints/kgrid'], 'mf_header/kpoints/kgrid')
    f_out.copy(f_in['mf_header/kpoints/shift'], 'mf_header/kpoints/shift')
    f_out.copy(f_in['mf_header/kpoints/ngk'], 'mf_header/kpoints/ngk')
    f_out.copy(f_in['mf_header/kpoints/ifmin'], 'mf_header/kpoints/ifmin')

    ######################################################################
    # This appears to work correctly, 12/1/2020
    #f_out.copy(f_in['mf_header/kpoints/ifmax'], 'mf_header/kpoints/ifmax')
    # Figure out which channel has fewer occupied bands; call that by some name.
    # Then set the number of occupied bands to be the larger of the two.

    ifmaxin = f_in['mf_header/kpoints/ifmax'][()] #ifmaxin for ifmax from f_in
    #print(ifmaxin)
    #print(list(ifmaxin.shape))
    if list(ifmaxin.shape)[1] != 1 :
        quit('Can not do SF-BSE for WFN files with more than one k-point yet')
    ifmaxin0 = ifmaxin[0,0]
    ifmaxin1 = ifmaxin[1,0]
    #print(ifmaxin0)
    #print(ifmaxin1)
    
    delta_ifmaxin = ifmaxin0 - ifmaxin1
    
    if delta_ifmaxin > 0 :
        is_large = 0 # up-spin channel; usual case, until multi-configuration SF is implemented
        is_small = 1
    elif delta_ifmaxin < 0 :
        is_large = 1 # down-spin channel has more occupied states; not expected to be common, yet
        is_small = 0
    elif delta_ifmaxin == 0 :
        quit('This is not a high-spin reference state WFN')
    
    f_out.create_dataset('mf_header/kpoints/ifmax', list(ifmaxin.shape), 'int')
    # FIXME: hard-coding that there is only one k-point
    f_out['mf_header/kpoints/ifmax'][0,0] = f_in['mf_header/kpoints/ifmax'][is_large,0]
    f_out['mf_header/kpoints/ifmax'][1,0] = f_in['mf_header/kpoints/ifmax'][is_large,0]

    #print(is_large)
    #print(f_out['mf_header/kpoints/ifmax'][0,0])
    #print(f_out['mf_header/kpoints/ifmax'][1,0])
    #print(f_in['mf_header/kpoints/ifmax'][0,0])
    #print(f_in['mf_header/kpoints/ifmax'][1,0])

    #Dataset: /mf_header/kpoints/ifmax
    #Type: integer
    #Rank: 2
    #Dims(1): nrk
    #Dims(2): nspin
    #Value: highest occupied band at each k-point. Numbering starts with 1 (not 0).
    ######################################################################

    f_out.copy(f_in['mf_header/kpoints/w'], 'mf_header/kpoints/w')
    f_out.copy(f_in['mf_header/kpoints/rk'], 'mf_header/kpoints/rk')

    ######################################################################    
    # This appears to work correctly, 12/1/2020

    #shape = list(f_in['mf_header/kpoints/el'].shape)
    #shape[-1] = nbsmall
    #f_out.create_dataset('mf_header/kpoints/el', shape, 'f')
    #f_out['mf_header/kpoints/el'][:,:,:] = f_in['mf_header/kpoints/el'][:,:,:nbsmall]

    #shape = list(f_in['mf_header/kpoints/occ'].shape)
    #shape[-1] = nbsmall
    #f_out.create_dataset('mf_header/kpoints/occ', shape, 'f')
    #f_out['mf_header/kpoints/occ'][:,:,:] = f_in['mf_header/kpoints/occ'][:,:,:nbsmall]

    # For the channel of interest, move band 0 to position (delta_ifmax),
    # all the way on up to (Nmax - delta_ifmax) to (Nmax), wiping out
    # the old bands (Nmax - delta_ifmax + 1) through (Nmax).
    # Do it this way, via StackOverflow question with title
    #    Shift list elements to the right and shift list element at the end :
    # for ii in range(abs(delta_ifmaxin)):
    #     bb.insert(0,bb.pop())
  
    
    # Then, for the new bands 0 to (delta_ifmax - 1), el = -1.0E10 + something;
    # and occupations are unity, of course.
    
    # For the channel not of interest, identify by (spin_index_of_interest+1)%2,
    # then just copy data.

    shape = list(f_in['mf_header/kpoints/el'].shape)
    f_out.create_dataset('mf_header/kpoints/el', shape, 'f')
    f_out['mf_header/kpoints/el'][is_large,:,:] = f_in['mf_header/kpoints/el'][is_large,:,:]
    # FIXME: no k-points allowed; in principle, just add a loop over kpoints
    eltemp = list(f_in['mf_header/kpoints/el'][is_small,0,:])
    for ii in range(abs(delta_ifmaxin)):
        eltemp.insert(0,eltemp.pop())
    f_out['mf_header/kpoints/el'][is_small,0,:] = np.asarray(eltemp)
    for ii in range(abs(delta_ifmaxin)):
        f_out['mf_header/kpoints/el'][is_small,0,ii] = -1.0e10 + ii    
    print(f_in['mf_header/kpoints/el'][is_small,0,:])
    print(f_out['mf_header/kpoints/el'][is_small,0,:])
    
    shape = list(f_in['mf_header/kpoints/occ'].shape)
    f_out.create_dataset('mf_header/kpoints/occ', shape, 'f')
    f_out['mf_header/kpoints/occ'][is_large,:,:] = f_in['mf_header/kpoints/occ'][is_large,:,:]
    # FIXME: no k-points allowed; in principle, just add a loop over kpoints
    occtemp = list(f_in['mf_header/kpoints/occ'][is_small,0,:])
    for ii in range(abs(delta_ifmaxin)):
        occtemp.insert(0,occtemp.pop())
    f_out['mf_header/kpoints/occ'][is_small,0,:] = np.asarray(occtemp)
    for ii in range(abs(delta_ifmaxin)):
        f_out['mf_header/kpoints/occ'][is_small,0,ii] = 1.
    print(f_in['mf_header/kpoints/occ'][is_small,0,:])
    print(f_out['mf_header/kpoints/occ'][is_small,0,:])
    
    #Dataset: /mf_header/kpoints/el
    #Type: double
    #Rank: 3
    #Dims(1): mnband
    #Dims(2): nrk
    #Dims(3): nspin
    #Value: mean-field energies. (In Ry)
    #Dataset: /mf_header/kpoints/occ

    #Type: double
    #Rank: 3
    #Dims(1): mnband
    #Dims(2): nrk
    #Dims(3): nspin
    #Value: occupations, between 0 and 1

    ######################################################################    

    f_out.copy(f_in['mf_header/gspace'], 'mf_header/gspace')
    f_out.copy(f_in['mf_header/symmetry'], 'mf_header/symmetry')
    f_out.copy(f_in['mf_header/crystal'], 'mf_header/crystal')

    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')


    ######################################################################    
    # This needs to be edited, the info for wfns/coeffs

    shape = list(f_in['wfns/coeffs'].shape)
    nbmax = shape[0]
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    #f_out['wfns/coeffs'][:,:,:] = f_in['wfns/coeffs'][:nbsmall,:,:]
    f_out['wfns/coeffs'][:,is_large,:,:] = f_in['wfns/coeffs'][:,is_large,:,:]
    # FIXME: no k-points allowed; in principle, just add a loop over kpoints
    for ib in range(nbmax):
        realcgtemp = list(f_in['wfns/coeffs'][ib,is_small,:,0])
        imagcgtemp = list(f_in['wfns/coeffs'][ib,is_small,:,1])

        for ii in range(abs(delta_ifmaxin)):
            realcgtemp.insert(0,realcgtemp.pop())
            imagcgtemp.insert(0,imagcgtemp.pop())

        f_out['wfns/coeffs'][ib,is_small,:,0] = np.asarray(realcgtemp)
        f_out['wfns/coeffs'][ib,is_small,:,1] = np.asarray(imagcgtemp)

        for ii in range(abs(delta_ifmaxin)):
            f_out['wfns/coeffs'][ii,is_small,:ii,0] = 0.0
            f_out['wfns/coeffs'][ii,is_small,:ii,1] = 0.0
            f_out['wfns/coeffs'][ii,is_small,ii,0] = 1.0
            f_out['wfns/coeffs'][ii,is_small,ii,1] = 0.0
            f_out['wfns/coeffs'][ii,is_small,ii+1:,0] = 0.0
            f_out['wfns/coeffs'][ii,is_small,ii+1:,1] = 0.0

    print(f_out['wfns/coeffs'][0,is_small,:,0])
    print(f_out['wfns/coeffs'][1,is_small,:,0])
    print(f_out['wfns/coeffs'][2,is_small,:,0])

    # For the channel of interest, move band 0 to position (delta_ifmax),
    # all the way on up to (Nmax - delta_ifmax) to (Nmax), wiping out
    # the old bands (Nmax - delta_ifmax + 1) through (Nmax).
    
    # For new bands 0 to (delta_ifmax - 1), coeffs = 1.0 for some gvector, 0.0 otherwise.
    #
    # Dataset: /wfns/coeffs
    # Type: double 
    # Rank: 4
    # Dims(1): 1 for REAL, 2 for CPLX
    # Dims(2): ngktot
    # Dims(3): nspin*nspinor
    # Dims(4): mnband
    
    # For the channel not of interest, identify by (spin_index_of_interest+1)%2,
    # then just copy data.
    ######################################################################    

# debug
#    print f_in['wfns/coeffs'][0,0,:]
#    print
#    print f_in['wfns/coeffs'][0,1,:]
#    print
#    print f_out['wfns/coeffs'][0,0,:]
#    print
#    print f_out['wfns/coeffs'][0,1,:]
# end debug


    f_in.close()
    f_out.close()


    return

def scgwtool():
    r"""
    This is the major function that controls the flow
    """

    swap_spins()
 
    return





def main():

    scgwtool()

    return





if __name__ == "__main__":
    main()




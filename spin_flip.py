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



def swap_spins():
    r"""

    Read in old wavefunction, spit out new wavefunction with spins reversed

    This function is modified based on Felipe's pseudobands.py script
    """

# things to change:
# mf_header/kpoints/ifmin
# mf_header/kpoints/ifmax [nspin, nrk]
# mf_header/kpoints/el [nspin,nrk,mnband]
# mf_header/kpoints/occ [nspin,nrk,mnband]
# wfns/coeffs [mnband,nspin*nspinor,ngktot,1 real 2 cplx]


    wfnold = "wfn_new.h5"
    wfnnew = "wfn_flip.h5"

    f_in = h5py.File(wfnold)
    f_out = h5py.File(wfnnew, "w")

    nbtot = f_in['mf_header/kpoints/mnband'][()]

    f_out.copy(f_in['mf_header'], 'mf_header')
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')

    f_out['mf_header/kpoints/mnband'][()] = nbtot


# debug
#    ifmax = f_out['mf_header/kpoints/ifmax'][0,:]
#    print "ifmax 1", ifmax
#    ifmax = f_out['mf_header/kpoints/ifmax'][1,:]
#    print "ifmax 2", ifmax
# end debug

    f_out['mf_header/kpoints/ifmin'][0,:] = f_in['mf_header/kpoints/ifmin'][1,:]
    f_out['mf_header/kpoints/ifmin'][1,:] = f_in['mf_header/kpoints/ifmin'][0,:]
    f_out['mf_header/kpoints/ifmax'][0,:] = f_in['mf_header/kpoints/ifmax'][1,:]
    f_out['mf_header/kpoints/ifmax'][1,:] = f_in['mf_header/kpoints/ifmax'][0,:]

# debug
#    ifmax = f_out['mf_header/kpoints/ifmax'][0,:]
#    print "ifmax 1", ifmax
#    ifmax = f_out['mf_header/kpoints/ifmax'][1,:]
#    print "ifmax 2", ifmax
# end debug

# debug
#    temp = f_in['mf_header/kpoints/el'][1,:,:]
#    print "eigenvalues"
#    print temp
#    f_out['mf_header/kpoints/el'][0,:,:] = temp
#    print "eigenvalues"
#    print f_out['mf_header/kpoints/el'][0,:,:]
# end debug


    f_out['mf_header/kpoints/el'][0,:,:] = f_in['mf_header/kpoints/el'][1,:,:]
    f_out['mf_header/kpoints/el'][1,:,:] = f_in['mf_header/kpoints/el'][0,:,:]
    f_out['mf_header/kpoints/occ'][0,:,:] = f_in['mf_header/kpoints/occ'][1,:,:]
    f_out['mf_header/kpoints/occ'][1,:,:] = f_in['mf_header/kpoints/occ'][0,:,:]

# debug
#    print f_in['mf_header/kpoints/el'][0,:,:]
#    print f_out['mf_header/kpoints/el'][0,:,:]
#    print f_in['mf_header/kpoints/el'][1,:,:]
#    print f_out['mf_header/kpoints/el'][1,:,:]
#    print f_in['mf_header/kpoints/occ'][0,:,:]
#    print f_out['mf_header/kpoints/occ'][0,:,:]
#    print f_in['mf_header/kpoints/occ'][1,:,:]
#    print f_out['mf_header/kpoints/occ'][1,:,:]
#    return
# end debug


    shape = list(f_in['wfns/coeffs'].shape)
    shape[0] = nbtot
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    f_out['wfns/coeffs'][:nbtot,0,:,:] = f_in['wfns/coeffs'][:nbtot,1,:,:]
    f_out['wfns/coeffs'][:nbtot,1,:,:] = f_in['wfns/coeffs'][:nbtot,0,:,:]

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




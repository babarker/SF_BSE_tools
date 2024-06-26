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

nbsmall = 30


def swap_spins():
    r"""

    Read in old wavefunction, spit out new wavefunction with reduced empty states

    This function is modified based on Felipe's pseudobands.py script
    """

# things to change:
# mf_header/kpoints/ifmin
# mf_header/kpoints/ifmax [nspin, nrk]
# mf_header/kpoints/el [nspin,nrk,mnband]
# mf_header/kpoints/occ [nspin,nrk,mnband]
# wfns/coeffs [mnband,nspin*nspinor,ngktot,1 real 2 cplx]


    wfnold = "wfn_old.h5"
    wfnnew = "wfn_flip.h5"

    f_in = h5py.File(wfnold)
    f_out = h5py.File(wfnnew, "w")


    f_out.create_group('mf_header')
    f_out.copy(f_in['mf_header/versionnumber'], 'mf_header/versionnumber')
    f_out.copy(f_in['mf_header/flavor'], 'mf_header/flavor')

    f_out.copy(f_in['mf_header/kpoints/nspin'], 'mf_header/kpoints/nspin')
    f_out.copy(f_in['mf_header/kpoints/nspinor'], 'mf_header/kpoints/nspinor')

    bshape = list(f_in['mf_header/kpoints/mnband'].shape)
    f_out.create_dataset('mf_header/kpoints/mnband', bshape, 'd')
    f_out['mf_header/kpoints/mnband'][()] = nbsmall

    f_out.copy(f_in['mf_header/kpoints/nrk'], 'mf_header/kpoints/nrk')
    f_out.copy(f_in['mf_header/kpoints/ngkmax'], 'mf_header/kpoints/ngkmax')
    f_out.copy(f_in['mf_header/kpoints/ecutwfc'], 'mf_header/kpoints/ecutwfc')
    f_out.copy(f_in['mf_header/kpoints/kgrid'], 'mf_header/kpoints/kgrid')
    f_out.copy(f_in['mf_header/kpoints/shift'], 'mf_header/kpoints/shift')
    f_out.copy(f_in['mf_header/kpoints/ngk'], 'mf_header/kpoints/ngk')

    f_out.copy(f_in['mf_header/kpoints/ifmin'], 'mf_header/kpoints/ifmin')
    f_out.copy(f_in['mf_header/kpoints/ifmax'], 'mf_header/kpoints/ifmax')
    f_out['mf_header/kpoints/ifmin'][0,:] = f_in['mf_header/kpoints/ifmin'][1,:]
    f_out['mf_header/kpoints/ifmin'][1,:] = f_in['mf_header/kpoints/ifmin'][0,:]
    f_out['mf_header/kpoints/ifmax'][0,:] = f_in['mf_header/kpoints/ifmax'][1,:]
    f_out['mf_header/kpoints/ifmax'][1,:] = f_in['mf_header/kpoints/ifmax'][0,:]


    f_out.copy(f_in['mf_header/kpoints/w'], 'mf_header/kpoints/w')
    f_out.copy(f_in['mf_header/kpoints/rk'], 'mf_header/kpoints/rk')

    shape = list(f_in['mf_header/kpoints/el'].shape)
    shape[-1] = nbsmall
    f_out.create_dataset('mf_header/kpoints/el', shape, 'f')
    f_out['mf_header/kpoints/el'][0,:,:] = f_in['mf_header/kpoints/el'][1,:,:nbsmall]
    f_out['mf_header/kpoints/el'][1,:,:] = f_in['mf_header/kpoints/el'][0,:,:nbsmall]

    shape = list(f_in['mf_header/kpoints/occ'].shape)
    shape[-1] = nbsmall
    f_out.create_dataset('mf_header/kpoints/occ', shape, 'f')
    f_out['mf_header/kpoints/occ'][0,:,:] = f_in['mf_header/kpoints/occ'][1,:,:nbsmall]
    f_out['mf_header/kpoints/occ'][1,:,:] = f_in['mf_header/kpoints/occ'][0,:,:nbsmall]

    f_out.copy(f_in['mf_header/gspace'], 'mf_header/gspace')
    f_out.copy(f_in['mf_header/symmetry'], 'mf_header/symmetry')
    f_out.copy(f_in['mf_header/crystal'], 'mf_header/crystal')

    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')


    shape = list(f_in['wfns/coeffs'].shape)
    shape[0] = nbsmall
    f_out.create_dataset('wfns/coeffs', shape, 'd')
    f_out['wfns/coeffs'][:,0,:] = f_in['wfns/coeffs'][:nbsmall,1,:,:]
    f_out['wfns/coeffs'][:,1,:] = f_in['wfns/coeffs'][:nbsmall,0,:,:]

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




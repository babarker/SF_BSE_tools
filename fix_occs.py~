#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#                                                                                         #
# Author: Zhenglu Li                                                                      #
# Date: 08/10/2015                                                                        #
#                                                                                         #
# Thanks to Diana Qiu, Ting Cao, Felipe da Jornada, Meng Wu                               #
#                                                                                         #
###########################################################################################

#
# Present version only supports complex spinless, with setting the following line in sigma.inp:
#      sigma_matrix   -1/-2   0
# we are averaging Sigma evaluated at two energies, from colomn and row
# However they are usually very close to each other, so using one might be enough
#
#
# To keep consistency with BerkeleyGW, we now force WFN_outer and WFN_iner use same number
# of kpoints in the same order, although not necessarily.
# kmap is defined under this assumption, or errors will occur.
#


import sys
import h5py
import numpy as np



def fix_occs():
    r"""

    Read in old wavefunction, spit out new wavefunction with occupations fixed

    This function is modified based on Felipe's pseudobands.py script
    """

    wfnold = "wfn_old.h5"
    wfnnew = "wfn_new.h5"

    f_in = h5py.File(wfnold)
    f_out = h5py.File(wfnnew, "w")

    occ = f_in['mf_header/kpoints/occ'][()]

    f_out.copy(f_in['mf_header'], 'mf_header')
    f_out.create_group('wfns')
    f_out.copy(f_in['wfns/gvecs'], 'wfns/gvecs')
    f_out.copy(f_in['wfns/coeffs'], 'wfns/coeffs')

    f_out['mf_header/kpoints/ifmax'][0,0] = int(np.sum(occ[0,0,:]))
    f_out['mf_header/kpoints/ifmax'][1,0] = int(np.sum(occ[1,0,:]))

    f_in.close()
    f_out.close()

    return

def scgwtool():
    r"""
    This is the major function that controls the flow
    """

    # wavefunctions
    print
    print "============================== Step 1/1 =============================="
    print "==                                                                  =="
    print "==                     Updating wavefunctions                       =="
    print "==                                                                  =="
    print "======================================================================"
    print

    fix_occs()
 
    return





def main():

    scgwtool()

    return





if __name__ == "__main__":
    main()




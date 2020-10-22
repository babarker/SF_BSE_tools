#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#
# Read bsemat.h5 file and return bsemat (direct_head, direct_body)
# as well as DFT eignvalues and occupations
#                                                                                         #
# Author:                                                                       #
# Date:                                                                         #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np
import constants

def get_bsemat():
    r"""
    Read in Direct_Head and Direct_Body matrix elements from bsemat.h5 file
    """

    bsemat = "bsemat.h5"

    f_in = h5py.File(bsemat)

    #nvb = f_in['bse_header/bands/nvb'][()]
    #ncb = f_in['bse_header/bands/ncb'][()]

    celvol = f_in['mf_header/crystal/celvol'][()]
    fac = bse_fac(celvol)

   #     Rank: 6
   #     Dims(1): flavor
   #     Dims(2:3): n1b, aka nvb
   #     Dims(4:5): n2b, aka ncb
   #     Dims(6:7): nk*ns
   #     Value: Kernel matrix elements.
                                
    head = f_in['mats/head'][()]
    wing = f_in['mats/wing'][()]
    body = f_in['mats/body'][()]

    #print()
    #print( "Number valence bands ...")
    #print( nvb)

    #print()
    #print( "Number conduction bands ...")
    #print( ncb)

    f_in.close()

    #return (nvb,ncb,fac*head,fac*body)
    #return (head,wing,body)
    return (fac*head,fac*wing,fac*body)

def get_enk_from_bsemat():
    r"""
    Read DFT energy eigenvalues from bsemat.h5 file
    """

    # For el and occ: Dims(3): nspin; Dims(2): nrk; Dims(1): mnband
    
    bsemat = "bsemat.h5"

    f_in = h5py.File(bsemat)

    el = f_in['mf_header/kpoints/el'][()]
    occ = f_in['mf_header/kpoints/occ'][()]
    ifmin = f_in['mf_header/kpoints/ifmin'][()] # Dims(1): nrk, Dims(2): nspin
    ifmax = f_in['mf_header/kpoints/ifmax'][()] # Dims(1): nrk, Dims(2): nspin

    f_in.close()

    #print()
    #print( "Number conduction bands ...")
    #print( ncb)

    return (el,occ,ifmin,ifmax)

def get_nvb_ncb_from_bsemat():

    r"""
    Read in number of valence and conduction bands from bsemat.h5 file
    """

    bsemat = "bsemat.h5"

    f_in = h5py.File(bsemat)

    nvb = f_in['bse_header/bands/nvb'][()]
    ncb = f_in['bse_header/bands/ncb'][()]

    f_in.close()

    return (nvb,ncb)

def bse_fac(celvol):

    # scaling and Nk are both 1.0
    fac = -8.0*np.pi/celvol

    return fac

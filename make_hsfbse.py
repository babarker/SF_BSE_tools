#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#                                                                                         #
# With bsemat.h5 contents already read, we now construct the Hamiltonian for the SF-BSE
#
#
# Author:                                                                       #
# Date:                                                                         #
#                                                                                         #
# Thanks to Diana Qiu, Ting Cao, Felipe da Jornada, Meng Wu                               #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np
from constants import Ry2eV as RYD

def read_weff():

    fin = open("w_eff.txt","r")
    lines = fin.readlines()
    weff = float(lines[1])
    fin.close()

    return weff

def make_hsfbse(ekv,ekc,head,wing,body):

    #print( [vv*RYD for vv in ekv])
    #print( [cc*RYD for cc in ekc])

    weff = read_weff()

    # initialize
    hbse = np.zeros((len(ekv)*len(ekc),len(ekv)*len(ekc)),dtype=float)

    print("From make_hsfbse.py ...")
    print()
    
    for iv in range(len(ekv)):
        for ic in range(len(ekc)):
            # BAB note: make sure this is consistent with below!!!
            #ibse = iv*len(ekc)+ic  
            ibse = ic*len(ekv)+iv  
            hbse[ibse,ibse] = ekc[ic] - ekv[iv]

            print("hbse["+str(ibse)+","+str(ibse)+"] = ekc["+str(ic)+"] - ekv["+str(iv)+"]")
            
            #print( hbse[ibse,ibse]*RYD)
            #print( "iv")
            #print( iv)
            #print( "ic")
            #print( ic)
            #print( " ")

    print("ekv[0]")
    print(ekv[0])
    print("ekv[1]")
    print(ekv[1])
    print("ekc[0]")
    print(ekc[0])
    print("ekc[1]")
    print(ekc[1])

    print()
    print()
    print("Energy differences: ")
    print(hbse)
    print()
    print()
    print("hbse[0,0]: ")
    print(hbse[0,0])
    print("hbse[1,1]: ")
    print(hbse[1,1])
    print("hbse[2,2]: ")
    print(hbse[2,2])
    print("hbse[3,3]: ")
    print(hbse[3,3])
    print()
    print()

    #print( "IPA:")
    #print( hbse*RYD)
    #print( " ")

    # then put in bsemat (kernel, direct-only) along diagonals...
    # recall that the kernel has "occ,up; unocc,down" in its "spin up" component;
    # no need to use "spin down" components of BSE kernel at all.

    #     Dims(1): flavor
    #     Dims(2:3): n1b, aka nvb
    #     Dims(4:5): n2b, aka ncb
    #     Dims(6:7): nk*ns

    # BIG BAB EDIT:
    # GETTING THE ORDERING RIGHT ON THE LOOPS IS CRITICAL

    # Via DAS (BerkeleyGW/Common/misc.f90):
    # bse_index = is + (iv - 1 + (ic - 1 + (ik - 1)*ncband_)*nvband_)*xct%nspin
    # N. B.: Fortran indexes from 1; Python indexes from 0.
    
    # In our case, is=1, ik=1, xct%nspin=1 for the above equation (Fortran indexing):
    # bse_index = 1 + (iv -1 + (ic - 1)*nvband); FORTRAN
    # bse_index = iv + ic*Nv; Python
    
    #for iv in range(len(ekv)):
    #    for ivp in range(len(ekv)):
    #        for ic in range(len(ekc)):
    #            for icp in range(len(ekc)):
                
                    #ii = iv*len(ekv) + ic
                    #ip = ivp*len(ekv) + icp
    #                ii = iv*len(ekc) + ic
    #                ip = ivp*len(ekc) + icp

                    #hbse[ii,ip] += head[0,0,ic,icp,iv,ivp,0] + body[0,0,ic,icp,iv,ivp,0]
    #                hbse[ii,ip] += head[ic,icp,iv,ivp] + wing[ic,icp,iv,ivp] + body[ic,icp,iv,ivp]

    print("From make_hsfbse.py ...")
    print()
    
    for iv in range(len(ekv)):
        for ic in range(len(ekc)):

            #ii = iv*len(ekc) + ic
            #ii = ic*len(ekc) + iv
            ii = ic*len(ekv) + iv

            for ivp in range(len(ekv)):
                for icp in range(len(ekc)):
                
                    #ip = ivp*len(ekc) + icp
                    #ip = icp*len(ekc) + ivp
                    ip = icp*len(ekv) + ivp

                    #hbse[ii,ip] += head[0,0,ic,icp,iv,ivp,0] + body[0,0,ic,icp,iv,ivp,0]
                    #hbse[ii,ip] += head[ic,icp,iv,ivp] + wing[ic,icp,iv,ivp] + body[ic,icp,iv,ivp]
                    hbse[ii,ip] += head[ic,icp,iv,ivp]*weff + wing[ic,icp,iv,ivp] + body[ic,icp,iv,ivp]
                    print("hbse["+str(ii)+","+str(ip)+"] += kernel["+str(ic)+","+str(icp)+","+str(iv)+","+str(ivp)+"]")

    print()
    print()
    print("H_BSE: ")
    print(hbse)
    print()
    print()

    #print( "BSE:")
    #print( hbse*RYD)
    #print( " ")
    return hbse

def make_hsfbse_no_kernel(ekv,ekc):

    #print( [vv*RYD for vv in ekv])
    #print( [cc*RYD for cc in ekc])


    # initialize
    hbse = np.zeros((len(ekv)*len(ekc),len(ekv)*len(ekc)),dtype=float)

    print("From make_hsfbse.py ...")
    print()
    
    for iv in range(len(ekv)):
        for ic in range(len(ekc)):
            # BAB note: make sure this is consistent with below!!!
            #ibse = iv*len(ekc)+ic  
            ibse = ic*len(ekv)+iv  
            hbse[ibse,ibse] = ekc[ic] - ekv[iv]

            print("hbse["+str(ibse)+","+str(ibse)+"] = ekc["+str(ic)+"] - ekv["+str(iv)+"]")
            
            #print( hbse[ibse,ibse]*RYD)
            #print( "iv")
            #print( iv)
            #print( "ic")
            #print( ic)
            #print( " ")

    print("ekv[0]")
    print(ekv[0])
    print("ekv[1]")
    print(ekv[1])
    print("ekc[0]")
    print(ekc[0])
    print("ekc[1]")
    print(ekc[1])

    print()
    print()
    print("Energy differences: ")
    print(hbse)
    print()
    print()
    print("hbse[0,0]: ")
    print(hbse[0,0])
    print("hbse[1,1]: ")
    print(hbse[1,1])
    print("hbse[2,2]: ")
    print(hbse[2,2])
    print("hbse[3,3]: ")
    print(hbse[3,3])
    print()
    print()

    return hbse

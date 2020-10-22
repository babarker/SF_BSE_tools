#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# This python script is used for spin-flip method.   #
#                                                                                         #
# After constructing the Hamiltonian for the SF-BSE,
# we diagonalize (eigenstuff) and re-order the columns by energy (order_excitations)
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

def eigenstuff(hbse):

    weig, Avec = np.linalg.eig(hbse)
    return (weig*RYD,Avec)

def order_excitations(weig,Avec):

    # The eigenvalues are not ordered lowest to highest automatically;
    # we find the order that would do so, then apply that re-ordering
    # to the excited-state vectors themselves.

    order = np.argsort(weig)
    tempw = np.array(weig)[order]

    # BAB edit: reorder rows and columns?
    temp = np.array(Avec)[:,order]
    # Attempt1
    #temp = np.array(Avec)[order,order]
    # Attempt2
    #temp = np.array(Avec)[:,order]
    #temp = temp[order,:]
    # Attempt3
    #temp = np.array(Avec)[:,order]
    #temp = temp[np.transpose(order),:]
    # Attempt4
    #temp = np.array(Avec)[:,order]
    # Attempt5
    #nrow = np.shape(Avec)[0]
    #for ii in range(nrow):
    #    temp = np.array(Avec[ii,:])[order]
    #    Avec[ii,:] = temp
    # Attempt6
    #nrow = np.shape(Avec)[0]
    #for ii in range(nrow):
    #    temp = np.array(Avec[ii,:])[:,order]
    #    Avec[ii,:] = temp
    # end BAB

    return (tempw,temp)

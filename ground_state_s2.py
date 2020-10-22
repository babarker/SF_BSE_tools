#!/usr/bin/env python

###########################################################################################
#                                                                                         #
#                                                                                         #
###########################################################################################

#
#


import sys
import h5py
import numpy as np

RYD = 13.60569253

def read_overlaps():       

    fname="overlap.out" # this is standard
    #rawtxt = np.loadtxt(fname,skiprows=9)
    
    f=open(fname)
    ln=f.readlines()[9:]
    f.close()

    # ordering of content of rawtxt:
    # band 1, energy 1, band 2, energy 2, spin, Re overlap, Im overlap, |overlap|^2

    reoarray=np.zeros((7,7))
    imoarray=np.zeros((7,7))
    osqarray=np.zeros((7,7))

    for l in ln:
        ib1=int(l.split()[0])
        ib2=int(l.split()[2])
        isp=int(l.split()[4])
        reo=float(l.split()[5])
        imo=float(l.split()[6])
        osq=float(l.split()[7])
                        
        #if ib1 >= 6 and ib1 <= 7 and ib2 >= 6 and ib2 <= 7 and isp == 1:
        if ib1 <= 7 and ib2 <= 7 and isp == 1:
#            reoarray[ib1-6][ib2-6]=reo
#            imoarray[ib1-6][ib2-6]=imo
#            osqarray[ib1-6][ib2-6]=osq
            reoarray[ib1-1][ib2-1]=reo
            imoarray[ib1-1][ib2-1]=imo
            osqarray[ib1-1][ib2-1]=osq

    #print reoarray
    #print imoarray
    #print osqarray

    return (reoarray,imoarray,osqarray)
    
def ground_state_spin(osqarray):

    # Term1 + Term2 + Term3
    # Term1: -1.0 * sum_{ij} |overlap(i jbar)|**2
    # Term2: (nup + ndn)/2 
    # Term3: ((nup - ndn)/2)**2
    
    # Term1 is np.sum(array)
    # Term2 fixed at 6.0
    # Term3 fixed at 1.0

    #summ = np.sum(osqarray[:5-1][:5-1])
    summ = np.sum(osqarray[:5][:5])
#    print summ
    s2 = -1.0*summ+6.0+1.0
    print s2

    return
    
(reoarray,imoarray,osqarray) = read_overlaps()
print " "
ground_state_spin(osqarray)

#!/usr/bin/env python

###########################################################################################
#                                                                                         #
# What we actually need is to read in all data for all potentially-occupied orbitals
# for any transition within the (nvb alpha, ncb beta) space
# This includes: all occupied up-spin orbitals (with particular orbital involved in transition being removed later)
#  all occupied down-spin orbitals,
#  all unoccupied down-spin orbitals potentially used in transitions.
#                
#                                                                                         #
###########################################################################################

#
#


import numpy as np

def parse_overlaps(ifmin,ifmax,nvb,ncb):       

    # N.B.: Since ifmin and ifmax were read from the bse_kernel file,
    #   and the kernel was calculated with the spins of WFNq_co reversed,
    #   ifmax[1,0] is used for the valence bands,
    #   ifmax[0,0] is used for the conduction bands.
    # This is confusing and perhaps un-satisfactory.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # March 2021: using WFN file data and not BSEMAT, spin no longer reversed!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #ifmax_v = ifmax[1,0]
    #ifmax_c = ifmax[0,0]
    ifmax_v = ifmax[0,0]
    ifmax_c = ifmax[1,0]

    fname="overlaps.dat" # this is standard
    
    # skip the first nine lines of text in the file
    f=open(fname)
    ln=f.readlines()[9:]
    f.close()

    # ordering of content of rawtxt:
    # band 1, energy 1, band 2, energy 2, spin, Re overlap, Im overlap, |overlap|^2

    # size of arrays needs to properly include all potential occupied orbitals

    #reoarray=np.zeros((nvb*ncb,nvb*ncb))
    #imoarray=np.zeros((nvb*ncb,nvb*ncb))
    #osqarray=np.zeros((nvb*ncb,nvb*ncb))
    #reoarray=np.zeros((nvb,ncb))
    #imoarray=np.zeros((nvb,ncb))
    #osqarray=np.zeros((nvb,ncb))
    reoarray=np.zeros((ifmax_v,ifmax_c+ncb))
    imoarray=np.zeros((ifmax_v,ifmax_c+ncb))
    osqarray=np.zeros((ifmax_v,ifmax_c+ncb))


    for l in ln:
        ib1=int(l.split()[0])
        ib2=int(l.split()[2])
        isp=int(l.split()[4])
        reo=float(l.split()[5])
        imo=float(l.split()[6])
        osq=float(l.split()[7])

        if isp == 1:
            #if ( (ib1 > ifmax_v-nvb and ib1 <= ifmax_v) and (ib2 > ifmax_c and ib2 <= ifmax_c+ncb) ):
            if ( (ib1 <= ifmax_v) and (ib2 <= ifmax_c+ncb) ):

                #print("ib1")
                #print(ib1)
                #print("ib2")
                #print(ib2)



                # Useful comment about converting raw index to BSE-compatible index, but not useful for this application:
                # Need to convert RAW band index read in from overlaps.dat file into BSE-format band index,
                # with ib1_raw --> ib1_bse = ifmax_v - ib1_raw (counting from 0 to Nv-1), for valence, and
                #      ib2_raw --> ib2_bse = ib2_raw - ifmax_c - 1 (counting from 0 to Nc-1), for conduction.
                
                #ib1 = ifmax_v - ib1
                #ib2 = ib2 - ifmax_c - 1
                
                #print("ib1_bse")
                #print(ib1) 
                #print("ib2_bse")
                #print(ib2) 
                #print()

                reoarray[ib1-1][ib2-1]=reo
                imoarray[ib1-1][ib2-1]=imo
                osqarray[ib1-1][ib2-1]=osq

    return (reoarray,imoarray,osqarray)

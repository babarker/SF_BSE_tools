#! /usr/bin/env python

#########
## This script is for use with artificial ROKS methods, helpful for SF-BSE,
## where a no-spin EQP calculation is performed,
## and only quasiparticle energies have only one spin index.

## The task is to more or less just double the file, with the second copy having the appropriate spin index.

import numpy as np


    # Note: This code assumes Gamma-Point Only kpoint sampling, for molecules or supercells
    
def main():

    # read eqp1.dat file
    file1 = open('eqp1.dat', 'r')
    Lines = file1.readlines()
    file1.close()
    
    # open processed file, write spin-up data:
    file2 = open('eqp1_roks.dat','w')
    file2.writelines(Lines)
    file2.close()

    # re-read eqp1.dat data, but change spin index from 1 to 2
    data = np.loadtxt('eqp1.dat',skiprows=1)
    #print(data[0])
    for dline in data:
        nspin = int(dline[0])
        #print(nspin)
        nband = int(dline[1])
        #print(nband)
        eks = float(dline[2])
        #print(eks)
        eqp = float(dline[3])
        #print(eqp)

        # format string for writing line in same style
        # 8 spaces for nspin
        # 8 spaces for nband
        # 5 spaces for sign and before decimal; 9 spaces after decimal
        # and same for eks and eqp
        line = '{0:8d} {1:7d} {2:>-5.9f} {3:>-5.9f}\n'.format(nspin+1,nband,eks,eqp)

        file_object = open('eqp1_roks.dat', 'a')
        file_object.write(line)
        file_object.close()


    
    return
    
main()

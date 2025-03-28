#!/usr/bin/env python

"""

Chemical shift referencing using DSS


Erik Walinda
Kyoto University
Graduate School of Medicine

Last change: 2015/03/02

    1. Record 1D spectrum of DSS in the same buffer (external)
       or the actual protein sample (internal) standard
       using 1D pulse with water suppression (e.g. p3919fpgp)
    2. Process the 1D in NMR draw and note the frequency of
       the DSS methyl resonance line (sharp peak around 0 ppm)
    3. Run this script in the folder containing the nmr data
       (acqu files etc.). Plug in the DSS-methyl shift in ppm.
    4. The script will give you the new, referenced spectral
       centers in ppm. Use these in your fid.com and you will
       obtain properly referenced spectra.

    Note: This script assumes 15N is on O3 and 13C is on O2.

"""


import os.path
import sys
from decimal import *


def main():
    fileCheck()
    SFO1, O1, SFO2, O2, SFO3, O3, DSS_ppm = getPars()
    calc(SFO1, O1, SFO2, O2, SFO3, O3, DSS_ppm)


def fileCheck():

    # Dimensionality and file check

    dimensions = 0

    if os.path.exists("./acqu") is True:
        dimensions += 1
    else:
        print "NMR file [./acqu] not found!"
        print "Please make sure you are in the correct directory."
        sys.exit()

    if os.path.exists("./acqu2") is True:
        dimensions += 1

    if os.path.exists("./acqu3") is True:
        dimensions += 1

    print str(dimensions) + "-dimensional experiment detected!"


def getPars():

    # Extract parameters

    acqu = open('acqus')

    for line in acqu:
        if "SFO1" in line:
            aa = line.split()
            SFO1 = float(aa[1])
        elif "$O1" in line:
            ab = line.split()
            O1 = float(ab[1])
        elif "SFO3" in line:
            ac = line.split()
            SFO3 = float(ac[1])
        elif "$O3" in line:
            ad = line.split()
            O3 = float(ad[1])
        elif "SFO2" in line:
            ae = line.split()
            SFO2 = float(ae[1])
        elif "$O2" in line:
            af = line.split()
            O2 = float(af[1])

    DSS_ppm = input('Chemical shift of DSS [ppm] ....  ')

    return SFO1, O1, SFO2, O2, SFO3, O3, DSS_ppm


def calc(SFO1, O1, SFO2, O2, SFO3, O3, DSS_ppm):

    getcontext().prec = 10

    carrier = (((Decimal(SFO1) * Decimal(1E6)) - Decimal(O1)) / Decimal(1E6))
    DSS_zure_Hz = carrier * Decimal((DSS_ppm * 1E-6)) * Decimal(1E6)
    DSS_shift_Hz = carrier + (DSS_zure_Hz * Decimal(1E-6))

    # Proton

    conversion_factor_hydrogen = Decimal(1)
    zero_frequency_hydrogen = DSS_shift_Hz * conversion_factor_hydrogen
    center_hydrogen = ((Decimal(SFO1) - zero_frequency_hydrogen) /
                       zero_frequency_hydrogen) * Decimal(1E6)

    # Nitrogen

    conversion_factor_nitrogen = Decimal(0.101329118)
    zero_frequency_nitrogen = DSS_shift_Hz * conversion_factor_nitrogen
    center_nitrogen = ((Decimal(SFO3) - zero_frequency_nitrogen) /
                       zero_frequency_nitrogen) * Decimal(1E6)

    # Carbon

    conversion_factor_carbon = Decimal(0.251449530)
    zero_frequency_carbon = DSS_shift_Hz * conversion_factor_carbon
    center_carbon = ((Decimal(SFO2) - zero_frequency_carbon) /
                     zero_frequency_carbon) * Decimal(1E6)

    print
    print "Parameters"
    print
    print "SFO1 [MHz] ..................... ", SFO1
    print "O1 [Hz] ........................ ", O1
    print "Carrier [MHz] .................. ", carrier
    print "DSS drift from 0.00 ppm [Hz] ... ", DSS_zure_Hz
    print "Frequency of DSS at 0 ppm [Hz] . ", DSS_shift_Hz

    print
    print
    print "DSS-referenced parameters"
    print
    print "1H new center [ppm] ............ ", center_hydrogen
    print "15N new center [ppm] ........... ", center_nitrogen
    print "13C new center [ppm] ........... ", center_carbon
    print


main()

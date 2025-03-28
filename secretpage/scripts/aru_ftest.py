#!/usr/bin/env python

"""

F-test script for glove3


Erik Walinda
Kyoto University
Graduate School of Medicine

Last change: 2016/03/02

"""

import argparse
import numpy
from scipy.stats import f


def readGloveout(filename):

    c = open(filename,'r')

    dof = []
    x2dof = []
    residues = []

    for line in c:
        columns = line.split()
        if len(columns) == 2:
            if columns[0] == 'DoF':
                dof.append(float(columns[1]))
            elif columns[0] == 'X2/DoF':
                x2dof.append(float(columns[1]))
            elif columns[0] == 'SET':
                residues.append(columns[1])
    c.close()

    x2dof = numpy.array(x2dof)
    dof = numpy.array(dof)

    x2 = x2dof * dof   # convert reduced chi-squares to chi-squares

    return x2, dof, residues


def computeFRatio(x2_const, x2_richards, dof_const, dof_richards):

    F = (x2_const - x2_richards) / x2_richards
    G = (dof_const - dof_richards) / dof_richards

    fRatio = F/G

    return fRatio 


def doFtest(constantFilename, richardsFilename, alpha):

    x2_const, dof_const, residues_const = readGloveout(constantFilename)
    x2_richards, dof_richards, residues_richards = readGloveout(richardsFilename)

    fRatio = computeFRatio(x2_const, x2_richards, dof_const, dof_richards)

    # compute pValues

    p = []

    for i in range(len(fRatio)):
        pValue = f.cdf(fRatio[i], dof_const[i], dof_richards[i])
        p.append(pValue)
        if pValue > alpha:
            reportPassedResidue(residues_const[i], residues_richards[i], pValue)



def reportPassedResidue(i, j, pValue):
    if i == j:
        print "SET ", i, " passed the F-test! (p: ", str(pValue), ")"
    if i != j:
        print "Mismatch in input files. Cannot compare datasets ", i, " and ", j, ". Please make sure both input files contain the same residues."
    # i may not be actual set number! keep assignment!



def main():
    parser = argparse.ArgumentParser(description='Perform a F-test on relaxation dispersion output from glove3.')

    parser.add_argument('constant_file', help='glove output file fitted to a constant model')
    parser.add_argument('dispersion_file', help='glove output file fitted to a relaxation dispersion model')
    parser.add_argument('alpha', nargs='?', help='alpha threshold for the F-test [default: 0.99]', type=float, default=0.99)

    args = parser.parse_args()

    doFtest(args.constant_file, args.dispersion_file, args.alpha)


main()

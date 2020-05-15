#!/usr/bin/python

import sys, string
import numpy as N

validkeywords = ['-prtvol','-figure']
arguments = sys.argv[1:]

#Checks usage and arguments
if len(arguments)>4:
    print 'Usage: python set_ecut_etotal.py abinitoutfile.out[X] [etotal tolerance in hartrees] [-prtvol] [-figure]'
    sys.exit()

#Taking care of keyword options
for arg in arguments:
    if arg[0] == '-':
        if arg not in validkeywords:
            print 'ERROR: %s is not a valid keyword ... exit' %arg
            sys.exit()
for keyword in validkeywords:
    if arguments.count(keyword) > 1:
        print 'ERROR: keyword "%s" is repeated %s times ... exit' %(keyword,arguments.count(keyword))
        sys.exit()

#Get keywords
prtvol = 0
etotalfigures = 0
if '-prtvol' in arguments:
    arguments.pop(arguments.index('-prtvol'))
    prtvol = 1
if '-figure' in arguments:
    arguments.pop(arguments.index('-figure'))
    etotalfigures = 1

#Optional import
if (etotalfigures == 1):
    import matplotlib.pyplot as plt

#Function definitions
def istointcastable(string):
    try:
        N.int(string)
        return 1
    except:
        return 0

def istofloatcastable(string):
    try:
        N.float(string)
        return 1
    except:
        return 0

# Parse abinit file 1, setup etotal tolerance
etotal_tol_ha = N.float(0.005)
input_file1_name = str(arguments[0]) # name of first input file (first argument)
if (len(arguments) == 2):
    if istofloatcastable(arguments[1]):
        etotal_tol_ha = N.float(arguments[1])
    else:
        print 'ERROR: the second argument "%s" cannot be casted to a float value ... exit' %arguments[1]
        sys.exit()
input_file1_r = open(input_file1_name,'r')  # open it as read file       
abinit_file1_data = input_file1_r.readlines()
input_file1_r.close()
ndtset = -1
ecuts_jdtset = list()
etotals_jdtset = list()
add_ecut_values = 1
for iline in range(len(abinit_file1_data)):
    if abinit_file1_data[iline].find(' ecut') > -1 and (add_ecut_values == 1) and not (abinit_file1_data[iline].find(' ecutsm') > -1):
        ecuts_jdtset.append(N.float(abinit_file1_data[iline].split()[1]))
    if abinit_file1_data[iline].find('ndtset') > -1:
        ndtset = N.int(abinit_file1_data[iline].split()[1])
        jdtset = N.zeros(ndtset,N.int)
    if abinit_file1_data[iline].find('       ionmov    ') > -1:
        if abinit_file1_data[iline][1] != "0":
            geometryconvergence = 1
    if abinit_file1_data[iline].find(' etotal') > -1:
        etotals_jdtset.append(N.float(abinit_file1_data[iline].split()[1]))
    if abinit_file1_data[iline].find('       znucl    ') > -1:
        add_ecut_values = 0


iconverged_cutoff = -1
converged_cutoff = N.float(-1)

if prtvol > 0:
    print "/------------------------------------\\"
    print "===>>> Total energy convergence <<<==="
    print "\\------------------------------------/"
    print ""
    print "Total energy tolerance for convergence (in Hartrees) :",etotal_tol_ha
    print "Highest (supposed converged) cutoff energy (in Hartrees) :",ecuts_jdtset[-1]
    print "Second highest cutoff energy (in Hartrees) :",ecuts_jdtset[-2]

for iecut in range(len(ecuts_jdtset)-1):
    if (iconverged_cutoff >= 0):
        if (N.abs(etotals_jdtset[iecut] - etotals_jdtset[-1]) > etotal_tol_ha):
            iconverged_cutoff = -1
            converged_cutoff = N.float(-1)
            if prtvol > 0:
                print "... but then not converged for the higher ecut=",ecuts_jdtset[iecut]
    else:
        if (N.abs(etotals_jdtset[iecut] - etotals_jdtset[-1]) < etotal_tol_ha):
            iconverged_cutoff = iecut
            converged_cutoff = N.float(ecuts_jdtset[iecut])
            if prtvol > 0:
                print "Converged cutoff found :", converged_cutoff

if (iconverged_cutoff >= 0):
    if prtvol == 0:
        print converged_cutoff
    elif prtvol > 0:
        print "Final converged cutoff : ", converged_cutoff
        print "Difference in total energy Etot(converged_cutoff) - Etot(max_cutoff) :",etotals_jdtset[iconverged_cutoff]-etotals_jdtset[-1]
    else:
        print "ERROR: prtvol =",prtvol,"is not allowed ... exit"
        sys.exit()
else:
    if prtvol == 0:
        print "-1"
    elif prtvol > 0:
        print "WARNING : no converged cutoff was found !"
    else:
        print "ERROR: prtvol =",prtvol,"is not allowed ... exit"
        sys.exit()

# Optional figures
if etotalfigures == 1:
#    ylim_expand_tol = 1.25
#    ylim_expand_max = 5.0
#    ylim_min = etotals_jdtset[-1] - N.min([ylim_expand_max*etotal_tol_ha,etotals_jdtset[-1]-N.min(etotals_jdtset)])
#    ylim_max = etotals_jdtset[-1] + N.min([ylim_expand_max*etotal_tol_ha,N.max(etotals_jdtset)-etotals_jdtset[-1]])
    etotals_plot = N.array(etotals_jdtset[:-1],N.float)-etotals_jdtset[-1]
    plt.figure(1)
    plt.hold('on')
    plt.plot([ecuts_jdtset[0],ecuts_jdtset[-1]],[etotals_jdtset[-1]+etotal_tol_ha,etotals_jdtset[-1]+etotal_tol_ha],'-',linewidth=1.5)
    plt.plot([ecuts_jdtset[0],ecuts_jdtset[-1]],[etotals_jdtset[-1]-etotal_tol_ha,etotals_jdtset[-1]-etotal_tol_ha],'-',linewidth=1.5)
    plt.plot(ecuts_jdtset,etotals_jdtset,'x-',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
#    plt.ylim([etotals_jdtset[-1]-etotal_tol_ha*ylim_expand_tol,etotals_jdtset[-1]+etotal_tol_ha*ylim_expand_tol])
    plt.grid('on')
    plt.title('Total energy convergence')
    plt.xlabel('Energy cutoff')
    plt.show()

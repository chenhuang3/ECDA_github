#!/usr/bin/python

import os, sys, string
import numpy as N
import netCDF4
import matplotlib.pyplot as plt
import etsfIO
import basic_utils

#############
# Variables #
#############

#Container for any group of variables
class VariableContainer:pass

#Control variables
ctrl = VariableContainer()

ctrl.number_of_etsf_wfk_files_to_read    = N.int(1)
ctrl.etsf_wfk_files_reading_option       = N.int(1)
ctrl.etsf_wfk_files_to_read              = list()

ctrl.geometry_primitive_vectors          = N.int(-1)
ctrl.geometry_red_coord                  = N.int(-1)
ctrl.geometry_cart_coord                 = N.int(-1)
ctrl.geometry                            = N.int(-1)
ctrl.kinetic_energy_cutoff               = N.int(-1)
ctrl.kpoints                             = N.int(-1)

ctrl.tests                               = dict()
ctrl.wfk_conv_indices                    = N.array([],N.int)

ctrl.prtvol                              = N.int(0)

#List containing all the wavefunction containers
wavefunctions_list                       = list()

constants = etsfIO.constants

# Beginning of script
print '+------------------------------------+'
print '| Reading of ETSF wavefunction files |'
print '+------------------------------------+\n'

#Get filenames and check they exist and that they are in ETSF file format
ctrl.number_of_etsf_wfk_files_to_read = basic_utils.raw_input_int('\nEnter the number of ETSF wavefunction files to read : \n')
print ctrl.number_of_etsf_wfk_files_to_read
if ctrl.number_of_etsf_wfk_files_to_read <= 0:
    basic_utils.error_exit('The number of ETSF wavefunction files to read is lower than or equal to 0')
else:
    for iread in range(ctrl.number_of_etsf_wfk_files_to_read):
        current_file = raw_input('\nEnter the name of file number %s : ' %(iread+1))
        print current_file
        etsfIO.check_file_etsf_wavefunction(current_file)
        if (current_file in ctrl.etsf_wfk_files_to_read):
            error_message = 'File "%s" is already in the list of files to be read' %current_file
            basic_utils.error_exit(error_message)
        else:
            ctrl.etsf_wfk_files_to_read.append(current_file)

#Get reading/checking choice
print '\nReading option'
print '   "1" => Read geometry, k-points and basic input/output variables'
print '   "2" => "1" + Read eigenvalues and occupation numbers'
print '   "3" => "2" + Read coefficients'
ctrl.etsf_wfk_files_reading_option = basic_utils.raw_input_int('\nEnter option for the reading of ETSF wavefunction files : \n')

if (ctrl.etsf_wfk_files_reading_option not in [1,2,3]):
    error_message = 'Wrong choice for reading options'
    basic_utils.error_exit(error_message)
else:
    print ctrl.etsf_wfk_files_reading_option

#Setup reading/checking options according to option chosen
etsfIO.set_read_basic()
if (ctrl.etsf_wfk_files_reading_option == 2):
    etsfIO.set_read_eigenvalues()
elif (ctrl.etsf_wfk_files_reading_option == 3):
    etsfIO.set_read_eigenvalues()
    etsfIO.set_read_coefficients()

#Read attributes, dimensions and variables asked
for ifile in range(ctrl.number_of_etsf_wfk_files_to_read):
    wfkdata = etsfIO.WaveFunctionData(ctrl.etsf_wfk_files_to_read[ifile])
    wavefunctions_list.append(wfkdata)

#Check some variables
#####################
#Units
units = etsfIO.get_units(wavefunctions_list)
comment_message = 'Units are "%s" throughout all the wavefunction files' %units
basic_utils.print_comment(comment_message)
#System
etsfIO.check_system(wavefunctions_list)
#Kinetic energy cutoff
if etsfIO.check_kinetic_energy_cutoff(wavefunctions_list) == constants.VARIABLE:
    ctrl.kinetic_energy_cutoff = constants.VARIABLE
    comment_message = 'The kinetic energy cutoff for the plane waves is variable from one wavefunction to another'
    basic_utils.print_comment(comment_message)
else:
    ctrl.kinetic_energy_cutoff = constants.CONSTANT
    if ctrl.prtvol > 0:
        comment_message = 'The kinetic energy cutoff for the plane waves is the same for all wavefunctions'
        basic_utils.print_comment(comment_message)
#Kpoints
if etsfIO.check_kpoints(wavefunctions_list) == constants.VARIABLE:
    ctrl.kpoints = constants.VARIABLE
    comment_message = 'The k-points used are changing from one wavefunction to another'
    basic_utils.print_comment(comment_message)
else:
    ctrl.kpoints = constants.CONSTANT
    if ctrl.prtvol > 0:
        comment_message = 'The k-points used are the same for all wavefunctions'
        basic_utils.print_comment(comment_message)
#Primitive translation vectors
if etsfIO.check_primitive_vectors(wavefunctions_list) == constants.VARIABLE:
    ctrl.geometry_primitive_vectors = constants.VARIABLE
    comment_message = 'The unit cell is changing from one wavefunction to another'
    basic_utils.print_comment(comment_message)
else:
    ctrl.geometry_primitive_vectors = constants.CONSTANT
    if ctrl.prtvol > 0:
        comment_message = 'The unit cell is the same for all wavefunctions'
        basic_utils.print_comment(comment_message)
#Reduced coordinates
if etsfIO.check_red_coord(wavefunctions_list) == constants.VARIABLE:
    ctrl.geometry_red_coord = constants.VARIABLE
    comment_message = 'The reduced coordinates of the atoms are changing from one wavefunction to another'
    basic_utils.print_comment(comment_message)
else:
    ctrl.geometry_red_coord = constants.CONSTANT
    if ctrl.prtvol > 0:
        comment_message = 'The reduced coordinates of the atoms are the same for all wavefunctions'
        basic_utils.print_comment(comment_message)
#Cartesian coordinates
if etsfIO.check_cart_coord(wavefunctions_list) == constants.VARIABLE:
    ctrl.geometry_cart_coord = constants.VARIABLE
    comment_message = 'The cartesian coordinates of the atoms are changing from one wavefunction to another'
    basic_utils.print_comment(comment_message)
else:
    ctrl.geometry_cart_coord = constants.CONSTANT
    if ctrl.prtvol > 0:
        comment_message = 'The cartesian coordinates of the atoms are the same for all wavefunctions'
        basic_utils.print_comment(comment_message)
#Geometry
if ctrl.geometry_primitive_vectors == constants.CONSTANT and \
   ctrl.geometry_red_coord == constants.CONSTANT and \
   ctrl.geometry_cart_coord == constants.CONSTANT:
    ctrl.geometry = constants.CONSTANT
else:
    ctrl.geometry = constants.VARIABLE
    error_message = 'Geometry changes in wavefunctions is not yet implemented'
    basic_utils.error_exit(error_message)

#Choices for case where only Kinetic energy cutoff is changing
print 'Enter option choice : '
print ' 1 => eigenvalue analysis'
print ' 2 => total energy analysis'
option_choice = raw_input(' ... ')

if option_choice == '1':
    print 'TODO ...'
elif option_choice == '2':
    error_message = 'total energy analysis not yet implemented'
    basic_utils.error_exit(error_message)
else:
    error_message = 'Wrong choice'
    basic_utils.error_exit(error_message)


#if ctrl.geometry == constants.CONSTANT:
#    if ctrl.kinetic_energy_cutoff == constants.VARIABLE:
#        if ctrl.kpoints == constants.VARIABLE:
#            error_message = 'Kinetic energy cutoff and k-points are both changing from one wavefunction to another'
#            basic_utils.error_exit(error_message)
#        ctrl.tests['kinetic_energy_cutoff_convergence_fixed_geometry'] = constants.ACTIVATED
#    elif ctrl.kinetic_energy_cutoff == constants.CONSTANT:
#        if ctrl.kpoints == constants.CONSTANT:
#            error_message = 'Kinetic energy cutoff and k-points are constant in all wavefunctions for a fixed geometry (same parameters everywhere ?)'
#            basic_utils.error_exit(error_message)
#        ctrl.tests['kpoints_convergence_fixed_geometry'] = constants.ACTIVATED
#
#for test,test_active in ctrl.tests.iteritems():
#    if test_active:
#        if test == 'kinetic_energy_cutoff_convergence_fixed_geometry':
ibandmin = 0
mynband = 8
ibandmax = ibandmin + mynband
band_block = N.arange(ibandmin,ibandmax)
ctrl.wfk_conv_indices,ecuts = etsfIO.get_indices_kinetic_energy(wavefunctions_list)
minvalues,maxvalues,meanvalues,stdvalues,varvalues = etsfIO.eigenvalues_analysis(wavefunctions_list,sort_indices=ctrl.wfk_conv_indices,\
      bands=band_block)
#            plt.figure(2)
#            plt.hold('on')
#            plt.grid('on')
##            plt.plot(ecuts,minvalues,'x-',label='Minimum difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
##            plt.plot(ecuts,maxvalues,'x-',label='Maximum difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
##            plt.plot(ecuts,meanvalues,'x-',label='Mean absolute difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
#            plt.plot(ecuts,27.211*minvalues,'x-',label='Minimum difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
#            plt.plot(ecuts,27.211*maxvalues,'x-',label='Maximum difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
#            plt.plot(ecuts,27.211*meanvalues,'x-',label='Mean absolute difference in eigenvalues',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
#            plt.title('Convergence of eigenvalues')
#            plt.xlabel('Energy cutoff (in Ha)')
#            plt.ylabel('Eigenvalue convergence (in eV)')
#            plt.legend(loc='best')
#            plt.show()
#            if etsfIO.have_keys(wavefunctions_list,'etot'):
#                values = etsfIO.get_values_from_key(wavefunctions_list,'etot')
#                etotals = values[ctrl.wfk_conv_indices]
#                etotal_tol_ha = 0.001
#                plt.figure(1)
#                plt.hold('on')
#                plt.plot([ecuts[0],ecuts[-1]],[etotals[-1]+etotal_tol_ha,etotals[-1]+etotal_tol_ha],'-',linewidth=1.5)
#                plt.plot([ecuts[0],ecuts[-1]],[etotals[-1]-etotal_tol_ha,etotals[-1]-etotal_tol_ha],'-',linewidth=1.5)
#                plt.plot(ecuts,etotals,'x-',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
##                plt.ylim([etotals_jdtset[-1]-etotal_tol_ha*ylim_expand_tol,etotals_jdtset[-1]+etotal_tol_ha*ylim_expand_tol])
#                plt.grid('on')
#                plt.title('Total energy convergence')
#                plt.xlabel('Energy cutoff')
#                plt.show()
#        else:
#            comment_message = 'Test "%s" is not yet implemented' %test
#            basic_utils.print_comment(comment_message)
#
##for wfk in wavefunctions_list[ctrl.wfk_conv_indices]:
#conveig = wavefunctions_list[ctrl.wfk_conv_indices[-1]].align()
#for iwfk in ctrl.wfk_conv_indices:
#    wfk = wavefunctions_list[iwfk]
#    neweig=wfk.align()
#    if ctrl.prtvol > 0:
#        print 'For wavefunction number %s, these are the results' %iwfk
#        print '  maximum difference (eigenvalue bigger than converged one) :',N.max(neweig-conveig)
#        print '  minimum difference (eigenvalue smaller than converged one) :',N.min(neweig-conveig)
#        print '  mean absolute difference:',N.mean(N.abs(neweig-conveig))
#        print '  std:',N.std(neweig-conveig)
#        print '  var:',N.var(neweig-conveig)
#
eig_ev_convtol=0.01
#
myparam=basic_utils.set_convergence_parameter_scalar(ecuts,etotals,0.001,parameter_type='Kinetic energy cutoff',\
           property_type='Total energy',prtvol=1,plot_figures=1)
#myparam=basic_utils.set_convergence_parameter_multiple_scalar(ecuts,minvalues,eig_ev_convtol/27.211,parameter_type='Kinetic energy cutoff',\
#           property_type='Minimum difference in eigenvalues',prtvol=1,plot_figures=0)
myparam=basic_utils.set_convergence_parameter_scalar(ecuts,27.211*maxvalues,eig_ev_convtol,parameter_type='Kinetic energy cutoff',\
           property_type='Minimum difference in eigenvalues',prtvol=10,plot_figures=1)
#myparam=basic_utils.set_convergence_parameter_scalar(ecuts,maxvalues,eig_ev_convtol/27.211,parameter_type='Kinetic energy cutoff',\
#           property_type='Maximum difference in eigenvalues',prtvol=1,plot_figures=0)
#myparam=basic_utils.set_convergence_parameter_scalar(ecuts,meanvalues,eig_ev_convtol/27.211,parameter_type='Kinetic energy cutoff',\
#           property_type='Mean absolute difference in eigenvalues',prtvol=1,plot_figures=0)
#etsfIO.set_do_not_read_any()
#sys.exit()

#!/usr/bin/python

import os, string
import numpy as N
import netCDF4
import basic_utils
from imp import find_module

#############
# Variables #
#############

#Container for any group of variables
class VariableContainer:pass

#Constants
constants = VariableContainer()
constants.DO_NOT_READ = N.int(0)         #Do not read the variable or dimension or attribute
constants.DO_READ = N.int(1)             #Just read the variable or dimension or attribute in an numpy.array format
constants.DO_READ_CHECK = N.int(2)       #Read the variable or dimension or attribute and checks whether it is the same for all files
constants.DO_READ_CHECK_ERROR = N.int(3) #Read the variable or dimension or attribute and checks whether it is the same for all files. If it
                                         #is not, exit the program (error)
constants.DO_READ_STR = N.int(4)             #Just read the variable or dimension or attribute in a string format
constants.DO_READ_STR_CHECK = N.int(5)       #Read the variable or dimension or attribute and checks whether it is the same for all files
constants.DO_READ_STR_CHECK_ERROR = N.int(6) #Read the variable or dimension or attribute and checks whether it is the same for all files. If it
                                         #is not, exit the program (error)
constants.DO_READ_INT = N.int(7)             #Just read the variable or dimension or attribute in a numpy.int format
constants.DO_READ_INT_CHECK = N.int(8)       #Read the variable or dimension or attribute and checks whether it is the same for all files
constants.DO_READ_INT_CHECK_ERROR = N.int(9) #Read the variable or dimension or attribute and checks whether it is the same for all files. If it
                                         #is not, exit the program (error)
constants.DO_READ_FLOAT = N.int(10)             #Just read the variable or dimension or attribute in a numpy.float format
constants.DO_READ_FLOAT_CHECK = N.int(11)       #Read the variable or dimension or attribute and checks whether it is the same for all files
constants.DO_READ_FLOAT_CHECK_ERROR = N.int(12) #Read the variable or dimension or attribute and checks whether it is the same for all files. If it
                                         #is not, exit the program (error)
constants.CONSTANT = N.int(0)
constants.VARIABLE = N.int(1)

constants.ERROR = N.int(0)
constants.NOERROR = N.int(1)

constants.ACTIVATED = N.int(1)
constants.DEACTIVATED = N.int(0)

constants.ATOMIC_UNITS = str('atomic units')

constants.ALIGN_OPTION_MEAN_OCCEIG = N.int(1)

constants.HA_TO_EV = 27.211

#Flags for the attributes (0 : do not read, 1 : read, 2 : read and check that it is the same for all files)
attributes_flags = dict()
dimensions_flags = dict()
variables_flags = dict()
variables_attributes_flags = dict()
# Get all the flags in the dictionnaries and set them to DO_NOT_READ
flagsfile = str(find_module('etsfIO')[1][:-9]+'flags.list')
try:
    flagsreader = open(flagsfile,'r')
except IOError:
    error_message = 'The flags.list file has not been found'
    basic_utils.error_exit(error_message)
flagslines = flagsreader.readlines()
basic_utils.clean(flagslines)
flagsreader.close()
attrread = False
dimread = False
varread = False
for iline in range(len(flagslines)):
    line = flagslines[iline]
    if len(line.split()) == 0:continue
    elif line.split()[0] == 'ATTRBEGIN':attrread = True
    elif line.split()[0] == 'ATTREND':attrread = False
    elif line.split()[0] == 'DIMBEGIN':dimread = True
    elif line.split()[0] == 'DIMEND':dimread = False
    elif line.split()[0] == 'VARBEGIN':varread = True
    elif line.split()[0] == 'VAREND':varread = False
    elif attrread:
        attributes_flags[line] = constants.DO_NOT_READ
    elif dimread:
        dimensions_flags[line] = constants.DO_NOT_READ
    elif varread:
        variables_flags[line.split()[0]] = constants.DO_NOT_READ
        if len(line.split()) > 1:
            variables_attributes_flags[line.split()[0]] = dict()
            for ivarattr in range(1,len(line.split())):
                variables_attributes_flags[line.split()[0]][line.split()[ivarattr]] = constants.DO_NOT_READ
    else:
        error_message = 'The reading of the flags is corrupted'
        basic_utils.error_exit(error_message)

#attributes_flags['file_format'] = constants.DO_READ_STR
#attributes_flags['file_format_version'] = constants.DO_READ_FLOAT
#attributes_flags['Conventions'] = constants.DO_READ_STR
#attributes_flags['title'] = constants.DO_READ_STR
#attributes_flags['history'] = constants.DO_READ_STR
#
##Flags for the dimensions (0 : do not read, 1 : read, 2 : read and check that it is the same for all files)
#
#dimensions_flags['number_of_atoms'] = constants.DO_READ_INT_CHECK_ERROR
#dimensions_flags['number_of_atom_species'] = constants.DO_READ_INT_CHECK_ERROR
#dimensions_flags['number_of_cartesian_directions'] = constants.DO_READ_INT
#dimensions_flags['number_of_symmetry_operations'] = constants.DO_READ_INT
#dimensions_flags['number_of_vectors'] = constants.DO_READ_INT
#dimensions_flags['max_number_of_states'] = constants.DO_READ_INT
#dimensions_flags['number_of_kpoints'] = constants.DO_READ_INT_CHECK
#dimensions_flags['number_of_spins'] = constants.DO_READ_INT_CHECK
#dimensions_flags['number_of_spinor_components'] = constants.DO_READ_INT_CHECK
#dimensions_flags['max_number_of_coefficients'] = constants.DO_READ_INT
#dimensions_flags['real_or_complex_coefficients'] = constants.DO_READ_INT
#dimensions_flags['number_of_reduced_dimensions'] = constants.DO_READ_INT
#dimensions_flags['number_of_components'] = constants.DO_READ_INT
#
#dimensions_flags['character_string_length'] = constants.DO_NOT_READ
#dimensions_flags['number_of_wavelet_resolutions'] = constants.DO_NOT_READ
#dimensions_flags['max_number_of_basis_grid_points'] = constants.DO_NOT_READ
#dimensions_flags['complex'] = constants.DO_NOT_READ
#dimensions_flags['max_number_of_angular_momenta'] = constants.DO_NOT_READ
#dimensions_flags['max_number_of_projectors'] = constants.DO_NOT_READ
#dimensions_flags['number_of_coefficients_dielectric_function'] = constants.DO_NOT_READ
#dimensions_flags['number_of_frequencies_dielectric_function'] = constants.DO_NOT_READ
#dimensions_flags['number_of_grid_points_vector1'] = constants.DO_NOT_READ
#dimensions_flags['number_of_grid_points_vector2'] = constants.DO_NOT_READ
#dimensions_flags['number_of_grid_points_vector3'] = constants.DO_NOT_READ
#dimensions_flags['number_of_localization_regions'] = constants.DO_NOT_READ
#dimensions_flags['number_of_qpoints_dielectric_function'] = constants.DO_NOT_READ
#dimensions_flags['number_of_qpoints_gamma_limit'] = constants.DO_NOT_READ
#dimensions_flags['symbol_length'] = constants.DO_NOT_READ
#dimensions_flags['npsp'] = constants.DO_NOT_READ
#dimensions_flags['codvsnlen'] = constants.DO_NOT_READ
#dimensions_flags['rhoijdim1'] = constants.DO_NOT_READ
#dimensions_flags['psptitlen'] = constants.DO_NOT_READ
#
##Flags for the variables (0 : do not read, 1 : read, 2 : read and check that it is the same for all files)
#
#variables_flags['space_group'] = constants.DO_READ_INT
#variables_flags['primitive_vectors'] = constants.DO_READ_CHECK
#variables_flags['reduced_symmetry_matrices'] = constants.DO_READ
#variables_flags['reduced_symmetry_translations'] = constants.DO_READ
#variables_flags['atom_species'] = constants.DO_READ_CHECK_ERROR
#variables_flags['reduced_atom_positions'] = constants.DO_READ_CHECK
#variables_flags['atomic_numbers'] = constants.DO_READ_CHECK_ERROR
#variables_flags['atom_species_names'] = constants.DO_READ_STR
#variables_flags['chemical_symbols'] = constants.DO_NOT_READ
#variables_flags['reduced_coordinates_of_kpoints'] = constants.DO_READ_CHECK
#variables_flags['kpoint_weights'] = constants.DO_READ
#variables_flags['number_of_states'] = constants.DO_READ
#variables_flags['fermi_energy'] = constants.DO_READ_FLOAT
#variables_flags['smearing_scheme'] = constants.DO_READ_STR_CHECK
#variables_flags['smearing_width'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['tsmear'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['ecutdg'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['usepaw'] = constants.DO_READ_INT_CHECK_ERROR
#variables_flags['etot'] = constants.DO_READ_FLOAT
#variables_flags['kinetic_energy_cutoff'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['ecut_eff'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['ixc'] = constants.DO_READ_INT_CHECK
#variables_flags['ecutsm'] = constants.DO_READ_FLOAT_CHECK
#variables_flags['exchange_functional'] = constants.DO_READ_STR
#variables_flags['correlation_functional'] = constants.DO_READ_STR
#variables_flags['occopt'] = constants.DO_READ_INT
#
#variables_flags['valence_charges'] = constants.DO_NOT_READ
#variables_flags['number_of_wavelets'] = constants.DO_NOT_READ
#variables_flags['number_of_coefficients_per_grid_point'] = constants.DO_NOT_READ
#variables_flags['coordinates_of_basis_grid_points'] = constants.DO_NOT_READ
#variables_flags['pseudopotential_types'] = constants.DO_NOT_READ
#variables_flags['number_of_electrons'] = constants.DO_NOT_READ
#variables_flags['eigenvalues'] = constants.DO_NOT_READ
#variables_flags['occupations'] = constants.DO_NOT_READ
#variables_flags['basis_set'] = constants.DO_NOT_READ
#variables_flags['number_of_coefficients'] = constants.DO_NOT_READ
#variables_flags['reduced_coordinates_of_plane_waves'] = constants.DO_NOT_READ
#variables_flags['date'] = constants.DO_NOT_READ
#variables_flags['codvsn'] = constants.DO_NOT_READ
#variables_flags['headform'] = constants.DO_NOT_READ
#variables_flags['fform'] = constants.DO_NOT_READ
#variables_flags['intxc'] = constants.DO_NOT_READ
#variables_flags['pertcase'] = constants.DO_NOT_READ
#variables_flags['residm'] = constants.DO_NOT_READ
#variables_flags['stmbias'] = constants.DO_NOT_READ
#variables_flags['tphysel'] = constants.DO_NOT_READ
#variables_flags['pspcod'] = constants.DO_NOT_READ
#variables_flags['pspdat'] = constants.DO_NOT_READ
#variables_flags['pspso'] = constants.DO_NOT_READ
#variables_flags['pspxc'] = constants.DO_NOT_READ
#variables_flags['qptn'] = constants.DO_NOT_READ
#variables_flags['so_psp'] = constants.DO_NOT_READ
#variables_flags['symafm'] = constants.DO_NOT_READ
#variables_flags['title'] = constants.DO_NOT_READ
#variables_flags['znuclpsp'] = constants.DO_NOT_READ
#variables_flags['lmn_size'] = constants.DO_NOT_READ
#variables_flags['rhoij'] = constants.DO_NOT_READ
#variables_flags['usewvl'] = constants.DO_NOT_READ
#variables_flags['coefficients_of_wavefunctions'] = constants.DO_NOT_READ
#
##Flags for the attributes of some variables
#fermi_energy_attributes_flags = dict()
#kinetic_energy_cutoff_attributes_flags = dict()
#reduced_symmetry_matrices_attributes_flags = dict()
#reduced_symmetry_translations_attributes_flags = dict()
#number_of_states_attributes_flags = dict()
#eigenvalues_attributes_flags = dict()
#number_of_coefficients_attributes_flags = dict()
#reduced_coordinates_of_plane_waves_attributes_flags = dict()
#
#variables_attributes_flags['fermi_energy'] = fermi_energy_attributes_flags
#variables_attributes_flags['kinetic_energy_cutoff'] = kinetic_energy_cutoff_attributes_flags
#variables_attributes_flags['reduced_symmetry_matrices'] = reduced_symmetry_matrices_attributes_flags
#variables_attributes_flags['reduced_symmetry_translations'] = reduced_symmetry_translations_attributes_flags
#variables_attributes_flags['number_of_states'] = number_of_states_attributes_flags
#variables_attributes_flags['eigenvalues'] = eigenvalues_attributes_flags
#variables_attributes_flags['number_of_coefficients'] = number_of_coefficients_attributes_flags
#variables_attributes_flags['reduced_coordinates_of_plane_waves'] = reduced_coordinates_of_plane_waves_attributes_flags
#
#fermi_energy_attributes_flags['units'] = constants.DO_READ_STR
#fermi_energy_attributes_flags['scale_to_atomic_units'] = constants.DO_NOT_READ
#
#kinetic_energy_cutoff_attributes_flags['units'] = constants.DO_READ_STR
#kinetic_energy_cutoff_attributes_flags['scale_to_atomic_units'] = constants.DO_NOT_READ
#
#reduced_symmetry_matrices_attributes_flags['symmorphic'] = constants.DO_READ_STR
#
#reduced_symmetry_translations_attributes_flags['symmorphic'] = constants.DO_READ_STR
#
#number_of_states_attributes_flags['k_dependent'] = constants.DO_READ_STR
#
#eigenvalues_attributes_flags['units'] = constants.DO_READ_STR
#eigenvalues_attributes_flags['scale_to_atomic_units'] = constants.DO_NOT_READ
#
#number_of_coefficients_attributes_flags['k_dependent'] = constants.DO_READ_STR
#
#reduced_coordinates_of_plane_waves_attributes_flags['k_dependent'] = constants.DO_READ_STR
#reduced_coordinates_of_plane_waves_attributes_flags['used_time_reversal_at_gamma'] = constants.DO_NOT_READ

#List containing all the wavefunction containers
wavefunctions_list                       = list()

#####################
# Class definitions #
#####################

#WaveFunctionData object containing all the data of a wavefunction
class WaveFunctionData(object):
    def __init__(self,filename=None):
        self.variables_attributes = dict()
        if not filename == None:
            check_file_etsf_wavefunction(filename)
            ncdata = netCDF4.Dataset(filename,'r')
            for attrname in ncdata.ncattrs():
                if attrname in attributes_flags:
                    if attributes_flags[attrname]:
                        attr = getattr(ncdata,attrname)
                        attrtype = get_type(attr)
                        put_attr_in_wfkdata(self,attrname,attr,attrtype)
                else:
                    warning_message = 'Attribute "%s" from the ETSF wavefunction file is an unknown attribute' %attrname
                    basic_utils.print_warning(warning_message)
            for dimname,dimobj in ncdata.dimensions.iteritems():
                if dimname in dimensions_flags:
                    if dimensions_flags[dimname]:
                        self.put(dimname,N.int(len(dimobj)))
                else:
                    warning_message = 'Dimension "%s" from the ETSF wavefunction file is an unknown dimension' %dimname
                    basic_utils.print_warning(warning_message)
            for varname in ncdata.variables:
                varobj = ncdata.variables[varname]
                if varname in variables_flags:
                    if N.ndim(varobj) == 0:
                        vartype = get_type(varobj.getValue())
                        varndim = N.int(0)
                        varshape = ()
                    else:
                        vartype = get_type(varobj[:])
                        varndim = N.ndim(varobj[:])
                        varshape = N.shape(varobj[:])
                    if variables_flags[varname]:
                        put_var_in_wfkdata(self,varname,varobj,varndim,varshape,vartype)
                else:
                    warning_message = 'Variable "%s" from the ETSF wavefunction file is an unknown variable' %varname
                    basic_utils.print_warning(warning_message)
            ncdata.close()
    def set_eig_occ(self,eigenvalues,occupations):
        self.eigenvalues = eigenvalues
        self.occupations = occupations
        self.has_eig_occ = True
    def set_wfk_coeff(self,number_of_coefficients,reduced_coordinates_of_plane_waves,coefficients_of_wavefunctions):
        self.number_of_coefficients = number_of_coefficients
        self.reduced_coordinates_of_plane_waves = reduced_coordinates_of_plane_waves
        self.coefficients_of_wavefunctions = coefficients_of_wavefunctions
        self.has_wfk_coeff = True
    def has_eig_occ(self):
        return self.has_eig_occ
    def has_wfk_coeff(self):
        return self.has_wfk_coeff
    def has_key(self,keyname):
        return keyname in self.__dict__
    def get_from_key(self,keyname):
        return self.__dict__[keyname]
    def put(self,keyname,variable):
        self.__dict__[keyname] = variable
    def set_attribute(self,keyvar,keyattr,attribute):
        if keyvar in self.variables_attributes:
            self.variables_attributes[keyvar][keyattr] = attribute
        else:
            self.variables_attributes[keyvar] = dict()
            self.variables_attributes[keyvar][keyattr] = attribute
    def get_units(self):
        for ivar in self.variables_attributes:
            if 'units' in self.variables_attributes[ivar]:
                if self.variables_attributes[ivar]['units'] != constants.ATOMIC_UNITS:
                    error_message = 'This wavefunction file contains variables in units differing from "%s"' %constants.ATOMIC_UNITS
                    basic_utils.error_exit(error_message)
        return constants.ATOMIC_UNITS
    def get_eig_align(self,align_option=constants.ALIGN_OPTION_MEAN_OCCEIG,number_of_bands=1,bands=N.array([1],N.int)):
        if 'eigenvalues' in self.__dict__:
            if align_option == constants.ALIGN_OPTION_MEAN_OCCEIG:
                if 'alignment_mean_occeig' in self.__dict__:
                    return self.alignment_mean_occeig
                elif 'occupations' in self.__dict__:
                    occtimeseig = N.multiply(self.occupations,self.eigenvalues)
                    wtk = N.repeat(N.repeat(N.reshape(self.kpoint_weights,[1,len(self.kpoint_weights),1]),\
                                   N.shape(occtimeseig)[2],axis=2),N.shape(occtimeseig)[0],axis=0)
                    wtktimesocctimeseig = N.multiply(wtk,occtimeseig)
                    occtimeswtk = N.multiply(wtk,self.occupations)
                    self.alignment_mean_occeig = N.sum(wtktimesocctimeseig)/N.sum(occtimeswtk)
                    return self.alignment_mean_occeig
                else:
                    error_message = 'The wavefunction data that has been read does not contain occupations, it is not possible to get the alignment'
                    basic_utils.error_exit(error_message)
            else:
                error_message = 'Alignment option is invalid'
                basic_utils.error_exit(error_message)
        else:
            error_message = 'The wavefunction data that has been read does not contain eigenvalues, it is not possible to get the alignment'
            basic_utils.error_exit(error_message)
    def align(self,align_option=constants.ALIGN_OPTION_MEAN_OCCEIG,number_of_bands=1,bands=N.array([1],N.int)):
        myalign = self.get_eig_align(align_option,number_of_bands,bands)
        if 'current_alignment' in self.__dict__ and self.current_alignment == align_option and 'aligned_eigenvalues' in self.__dict__:
            return self.aligned_eigenvalues
        else:
            self.current_alignment = align_option
            self.aligned_eigenvalues = self.eigenvalues - myalign
            return self.aligned_eigenvalues
    def realspace_wfk_from_pw_coeff(self):
        if 'realspace_wfk_comp' in self.__dict__:
            return self.realspace_wfk_comp
        elif 'coefficients_of_wavefunctions' in self.__dict__:
            self.realspace_wfk_comp = N.fft.ifftn(self.coefficients_of_wavefunctions)
            return self.realspace_wfk_comp
        else:
            error_message = 'Wavefunction does not contain plane wave coefficients,\
                             the real space wavefunction cannnot be computed'
            basic_utils.error_exit(error_message)

###########
# Methods #
###########
def check_file_etsf_wavefunction(filename):
    if not (os.path.isfile(filename)):
        error_message = 'File "%s" does not exists' %filename
        basic_utils.error_exit(error_message)
    else:
        try:
            test = netCDF4.Dataset(filename)
            attributes_dict = test.__dict__
            if not ('file_format' in attributes_dict):
                error_message = 'File "%s" does not contain the "file_format" attribute' %filename
                test.close()
                basic_utils.error_exit(error_message)
            elif not (attributes_dict['file_format'] == 'ETSF Nanoquanta'):
                error_message = 'Attribute "file_format" of file "%s" is %s \
                               \nwhile it should be "ETSF Nanoquanta"' %(filename,attributes_dict['file_format'])
                test.close()
                basic_utils.error_exit(error_message)
        except RuntimeError:
            error_message = 'File "%s" cannot be opened as a netCDF file' %filename
            basic_utils.error_exit(error_message)

def get_type(variable):
    if N.isscalar(variable):
        typ = type(variable)
        if typ is N.bool or typ is bool:
            return N.bool
        elif typ is N.int or typ is N.int8 or typ is N.int16 or typ is N.int32 or typ is N.int64:
            return N.int
        elif typ is N.uint or typ is N.uint8 or typ is N.uint16 or typ is N.uint32 or typ is N.uint64:
            return N.uint
        elif typ is N.float or typ is N.float32 or typ is N.float64:
            return N.float
        elif typ is N.complex or typ is N.complex64 or typ is N.complex128:
            return N.complex
        elif typ is str or typ is unicode:
            return str
        else:
            error_message = 'Type "%s" not found in get_type(variable)' %type(variable)
            basic_utils.error_exit(error_message)
    else:
        if N.issubdtype(N.bool,variable.dtype):
            return N.bool
        elif N.issubdtype(variable.dtype,N.int):
            return N.int
        elif N.issubdtype(variable.dtype,N.uint):
            return N.uint
        elif N.issubdtype(N.float,variable.dtype):
            return N.float
        elif N.issubdtype(variable.dtype,N.complex):
            return N.complex
        elif N.issubdtype(variable.dtype,str):
            return str
        else:
            error_message = 'Type "%s" not found in get_type(variable)' %variable.dtype
            basic_utils.error_exit(error_message)

def put_attr_in_wfkdata(wfkdata,attrname,attr,attrtype,varname=None):
    if varname is None:
        if attrtype is N.bool:
            wfkdata.put(attrname,N.bool(attr))
        elif attrtype is N.int:
            wfkdata.put(attrname,N.int(attr))
        elif attrtype is N.uint:
            wfkdata.put(attrname,N.uint(attr))
        elif attrtype is N.float:
            wfkdata.put(attrname,N.float(attr))
        elif attrtype is N.complex:
            wfkdata.put(attrname,N.complex(attr))
        elif attrtype is str:
            wfkdata.put(attrname,str(attr))
        else:
            error_message = 'Global attribute "%s" of type "%s" cannot be put in a wavefunction data' %(attrname,attrtype)
            basic_utils.error_exit(error_message)
    else:
        if attrtype is N.bool:
            newattr = N.bool(attr)
        elif attrtype is N.int:
            newattr = N.int(attr)
        elif attrtype is N.uint:
            newattr = N.uint(attr)
        elif attrtype is N.float:
            newattr = N.float(attr)
        elif attrtype is N.complex:
            newattr = N.complex(attr)
        elif attrtype is str:
            newattr = str(attr)
        else:
            error_message = 'Attribute "%s" of variable "%s" of type "%s" cannot be put in a wavefunction data' %(attrname,varname,attrtype)
            basic_utils.error_exit(error_message)
        wfkdata.set_attribute(varname,attrname,newattr)

def put_var_in_wfkdata(wfkdata,varname,varobj,varndim,varshape,vartype):
    if varndim == 0:
        if vartype == N.bool:
            wfkdata.put(varname,N.bool(varobj.getValue()))
        elif vartype == N.int:
            wfkdata.put(varname,N.int(varobj.getValue()))
        elif vartype == N.uint:
            wfkdata.put(varname,N.uint(varobj.getValue()))
        elif vartype == N.float:
            wfkdata.put(varname,N.float(varobj.getValue()))
        elif vartype == N.complex:
            wfkdata.put(varname,N.complex(varobj.getValue()))
        elif vartype == str:
            wfkdata.put(varname,str(varobj.tostring()))
        else:
            error_message = 'Variable "%s" with type "%s" cannot be put in a wavefunction data' %(varname,vartype)
            basic_utils.error_exit(error_message)
    else:
        if vartype == N.bool or vartype == N.int or vartype == N.uint or vartype == N.float or vartype == N.complex:
            newarray = N.reshape(N.array(varobj,vartype),varshape)
#            newarray = N.array(varobj,vartype)
            wfkdata.put(varname,newarray)
        elif vartype == str:
            newstring = str(varobj[:].tostring())
            wfkdata.put(varname,newstring)
        else:
            error_message = 'Variable "%s" with type "%s" and shape "%s" cannot be put in a wavefunction data' %(varname,vartype,varshape)
            basic_utils.error_exit(error_message)
    if variables_flags[varname] and varname in variables_attributes_flags:
        for attrname in varobj.ncattrs():
            if attrname in variables_attributes_flags[varname]:
                if variables_attributes_flags[varname][attrname]:
                    attr = varobj.getncattr(attrname)
                    attrtype = get_type(attr)
                    put_attr_in_wfkdata(wfkdata,attrname,attr,attrtype,varname)

def get_units(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        wfklist[0].get_units()
    else:
        ref_units = wfklist[0].get_units()
        for iwfk in range(1,len(wfklist)):
            wfk_units = wfklist[iwfk].get_units()
            if wfk_units != ref_units:
                error_message = 'The wavefunction file number %s contains variables defined in other units than "atomic units"'
                basic_utils.error_exit(error_message)
    return 'atomic units'

def check_system(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.NOERROR
    else:
        ref_number_of_atoms = wfklist[0].number_of_atoms
        ref_number_of_atom_species = wfklist[0].number_of_atom_species
        ref_atom_species = wfklist[0].atom_species
        ref_atomic_numbers = wfklist[0].atomic_numbers
        ref_number_of_spins = wfklist[0].number_of_spins
        ref_number_of_spinor_components = wfklist[0].number_of_spinor_components
        ref_number_of_components = wfklist[0].number_of_components
        for iwfk in range(1,len(wfklist)):
            if (wfklist[iwfk].number_of_atoms != ref_number_of_atoms) or \
                  (wfklist[iwfk].number_of_atom_species != ref_number_of_atom_species) or \
                  (wfklist[iwfk].number_of_spins != ref_number_of_spins) or \
                  (wfklist[iwfk].number_of_spinor_components != ref_number_of_spinor_components) or \
                  (wfklist[iwfk].number_of_components != ref_number_of_components) or \
              not (N.array_equal(wfklist[iwfk].atom_species,ref_atom_species)) or \
              not (N.array_equal(wfklist[iwfk].atomic_numbers,ref_atomic_numbers)):
                error_message = 'The list contains wavefunctions of different systems : the atoms or atom species differ'
                basic_utils.error_exit(error_message)
        return constants.NOERROR

def check_primitive_vectors(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.CONSTANT
    else:
        ref_primitive_vectors = wfklist[0].primitive_vectors
        for iwfk in range(1,len(wfklist)):
            if not N.array_equal(ref_primitive_vectors,wfklist[iwfk].primitive_vectors):
                return constants.VARIABLE
        return constants.CONSTANT

def check_red_coord(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.CONSTANT
    else:
        ref_reduced_atom_positions = wfklist[0].reduced_atom_positions
        for iwfk in range(1,len(wfklist)):
            if not N.array_equal(ref_reduced_atom_positions,wfklist[iwfk].reduced_atom_positions):
                return constants.VARIABLE
        return constants.CONSTANT

def check_cart_coord(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.CONSTANT
    else:
        ref_cartesian_atom_positions = N.dot(wfklist[0].reduced_atom_positions,wfklist[0].primitive_vectors)
        for iwfk in range(1,len(wfklist)):
            cartesian_atom_positions = N.dot(wfklist[iwfk].reduced_atom_positions,wfklist[iwfk].primitive_vectors)
            if not N.array_equal(ref_cartesian_atom_positions,cartesian_atom_positions):
                return constants.VARIABLE
        return constants.CONSTANT

def check_kinetic_energy_cutoff(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.CONSTANT
    else:
        ref_kinetic_energy_cutoff = wfklist[0].kinetic_energy_cutoff
        for iwfk in range(1,len(wfklist)):
            if (wfklist[iwfk].kinetic_energy_cutoff != ref_kinetic_energy_cutoff):
                return constants.VARIABLE
        return constants.CONSTANT

def check_kpoints(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return constants.CONSTANT
    else:
        ref_number_of_kpoints = wfklist[0].number_of_kpoints
        ref_reduced_coordinates_of_kpoints =  wfklist[0].reduced_coordinates_of_kpoints
        for iwfk in range(1,len(wfklist)):
            if (wfklist[iwfk].number_of_kpoints != ref_number_of_kpoints) or \
               not N.array_equal(wfklist[iwfk].reduced_coordinates_of_kpoints,ref_reduced_coordinates_of_kpoints):
                return constants.VARIABLE
        return constants.CONSTANT

def get_indices_kinetic_energy(wfklist):
    if len(wfklist) == 0:
        error_message = 'The list does not contain any wavefunction data'
        basic_utils.error_exit(error_message)
    elif len(wfklist) == 1:
        return N.array([0],N.int)
    else:
        ecuts = N.zeros(len(wfklist),N.float)
        for iwfk in range(len(wfklist)):
            ecuts[iwfk] = wfklist[iwfk].kinetic_energy_cutoff
        indices = N.argsort(ecuts)
        return indices,ecuts[indices]

def have_keys(wfklist,keyname):
    for wfk in wfklist:
        if not wfk.has_key(keyname):
            return False
    return True

def get_values_from_key(wfklist,keyname):
    values = N.zeros(len(wfklist),N.float)
    for iwfk in range(len(wfklist)):
        values[iwfk] = wfklist[iwfk].get_from_key(keyname)
    return values

def eigenvalues_analysis(wfklist,sort_indices=None,elk_compare=constants.DEACTIVATED,bands=None,elk_bands=None):
    if sort_indices is None:
        comment_message = 'Sorting indices have been set to [0,1,2,...,length(wavefunction_list)] in eigenvalues_analysis'
        basic_utils.print_comment(comment_message)
        sort_indices = range(len(wfklist))
    if bands == None:
        bands = N.arange(N.min(wfklist[0].number_of_states))
        comment_message = 'All the bands are included in the eigenvalues analysis'
        basic_utils.print_comment(comment_message)
    conveig = wfklist[sort_indices[-1]].align()[:,:,bands]
    minvalues = N.zeros(len(wfklist),N.float)
    maxvalues = N.zeros(len(wfklist),N.float)
    meanvalues = N.zeros(len(wfklist),N.float)
    stdvalues = N.zeros(len(wfklist),N.float)
    varvalues = N.zeros(len(wfklist),N.float)
    for i_iwfk in range(len(sort_indices)):
        iwfk = sort_indices[i_iwfk]
        eig_iwfk = wfklist[iwfk].align()[:,:,bands]
        minvalues[i_iwfk] = N.min(eig_iwfk-conveig)
        maxvalues[i_iwfk] = N.max(eig_iwfk-conveig)
        meanvalues[i_iwfk] = N.mean(N.abs(eig_iwfk-conveig))
        stdvalues[i_iwfk] = N.std(eig_iwfk-conveig)
        varvalues[i_iwfk] = N.var(eig_iwfk-conveig)
    return minvalues,maxvalues,meanvalues,stdvalues,varvalues

def set_do_not_read_any():    # Sets all the flags whether to read or not an attribute, dimension or variable to DO_NOT_READ
    for attr in attributes_flags:
        attributes_flags[attr] = constants.DO_NOT_READ
    for dims in dimensions_flags:
        dimensions_flags[dims] = constants.DO_NOT_READ
    for var in variables_flags:
        variables_flags[var] = constants.DO_NOT_READ
    for var in variables_attributes_flags:
        for varattr in variables_attributes_flags[var]:
            variables_attributes_flags[var][varattr] = constants.DO_NOT_READ

def set_read_attributes():
    for attr in attributes_flags:
        attributes_flags[attr] = constants.DO_READ

def set_read_dimensions():
    for dims in dimensions_flags:
        dimensions_flags[dims] = constants.DO_READ

def set_read_variables():
    for var in variables_flags:
        variables_flags[var] = constants.DO_READ
    for var in variables_attributes_flags:
        for varattr in variables_attributes_flags[var]:
            variables_attributes_flags[var][varattr] = constants.DO_READ

def set_read_all():
    set_read_attributes()
    set_read_dimensions()
    set_read_variables()

def set_read_basic():
    set_read_attributes()
    set_read_dimensions()
    var_to_read = ['atom_species','atomic_numbers','kinetic_energy_cutoff','reduced_coordinates_of_kpoints',\
                   'primitive_vectors','reduced_atom_positions','kpoint_weights','etot']
    for varname in var_to_read:
        variables_flags[varname] = constants.DO_READ
        if varname in variables_attributes_flags:
            for varattr in variables_attributes_flags[varname]:
                variables_attributes_flags[varname][varattr] = constants.DO_READ

def set_read_eigenvalues():
    var_to_read = ['eigenvalues','occupations']
    for varname in var_to_read:
        variables_flags[varname] = constants.DO_READ
        if varname in variables_attributes_flags:
            for varattr in variables_attributes_flags[varname]:
                variables_attributes_flags[varname][varattr] = constants.DO_READ

def set_read_coefficients():
    var_to_read = ['number_of_coefficients','reduced_coordinates_of_plane_waves','coefficients_of_wavefunctions']
    for varname in var_to_read:
        variables_flags[varname] = constants.DO_READ
        if varname in variables_attributes_flags:
            for varattr in variables_attributes_flags[varname]:
                variables_attributes_flags[varname][varattr] = constants.DO_READ

#
# Default configuration parameters obtained from Abinit's Autoconf
#

# ---------------------------------------------------------------------------- #

#
# System information
#

# Architecture
abi_cpu_model   = xeon
abi_cpu_64bits  = yes
abi_cpu_spec    = intel_xeon
abi_sys_spec    = linux-x86_64

# ---------------------------------------------------------------------------- #

#
# Libraries and linking
#

AR             = ar
ARFLAGS        =      rc
ARFLAGS_64BITS = 
ARFLAGS_DEBUG  = 
ARFLAGS_OPTIM  = 
ARFLAGS_EXTRA  = 
ARFLAGS_CMD    = rc
RANLIB         = ranlib

# ---------------------------------------------------------------------------- #

#
# Preprocessing
#

CPP               = gcc -E
CPPFLAGS          = 
CPPFLAGS_DEBUG    = 
CPPFLAGS_EXTRA    = 
CPPFLAGS_HINTS    = 
CPPFLAGS_OPTIM    = 

FPP               = 
FPPFLAGS          =     
FPPFLAGS_DEBUG    = 
FPPFLAGS_EXTRA    = 
FPPFLAGS_HINTS    = 
FPPFLAGS_OPTIM    = 

INCLUDES          = 

TRUE_CPP          = cpp -P -std=c99

# ---------------------------------------------------------------------------- #

#
# C compilation
#

cc_info_string    = gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-16)

abi_cc_vendor     = gnu
abi_cc_version    = 4.8

CC                = gcc
CFLAGS            =  -g -O2 -mtune=native -march=native  
CFLAGS_64BITS     = 
CFLAGS_DEBUG      = -g
CFLAGS_EXTRA      = 
CFLAGS_HINTS      = 
CFLAGS_OPTIM      = -O2 -mtune=native -march=native

CC_LDFLAGS_64BITS = 
CC_LDFLAGS_DEBUG  = 
CC_LDFLAGS_EXTRA  = 
CC_LDFLAGS_HINTS  = 
CC_LDFLAGS_OPTIM  = 
CC_LDFLAGS        =     

CC_LIBS_DEBUG     = 
CC_LIBS_EXTRA     = 
CC_LIBS_HINTS     = 
CC_LIBS_OPTIM     = 
CC_LIBS           =     

# ---------------------------------------------------------------------------- #

#
# C++ compilation
#

cxx_info_string    = g++ (GCC) 4.8.5 20150623 (Red Hat 4.8.5-16)

abi_cxx_vendor     = gnu
abi_cxx_version    = 4.8

CXX                = g++
CXXFLAGS           =  -g -O2 -mtune=native -march=native  
CXXFLAGS_64BITS    = 
CXXFLAGS_DEBUG     = -g
CXXFLAGS_EXTRA     = 
CXXFLAGS_HINTS     = 
CXXFLAGS_OPTIM     = -O2 -mtune=native -march=native

CXX_LDFLAGS_64BITS = 
CXX_LDFLAGS_DEBUG  = 
CXX_LDFLAGS_EXTRA  = 
CXX_LDFLAGS_HINTS  = 
CXX_LDFLAGS_OPTIM  = 
CXX_LDFLAGS        =     

CXX_LIBS_DEBUG     = 
CXX_LIBS_EXTRA     = 
CXX_LIBS_HINTS     = 
CXX_LIBS_OPTIM     = 
CXX_LIBS           =     

# ---------------------------------------------------------------------------- #

#
# Fortran compilation
#

fc_info_string    = 

fc_vendor         = gnu
fc_version        = 4.8

FC                = gfortran
FCFLAGS           =  -g -ffree-line-length-none  
FCFLAGS_64BITS    = 
FCFLAGS_DEBUG     = -g
FCFLAGS_EXTRA     = 
FCFLAGS_HINTS     = -ffree-line-length-none
FCFLAGS_OPTIM     = -O2 -mtune=native -march=native

FCFLAGS_FIXEDFORM = -ffixed-form
FCFLAGS_FREEFORM  = -ffree-form

FC_LDFLAGS_64BITS = 
FC_LDFLAGS_DEBUG  = 
FC_LDFLAGS_EXTRA  = 
FC_LDFLAGS_HINTS  = 
FC_LDFLAGS_OPTIM  = 
FC_LDFLAGS        =     

FC_LIBS_DEBUG     = 
FC_LIBS_EXTRA     = 
FC_LIBS_HINTS     = 
FC_LIBS_OPTIM     = 
FC_LIBS           =       -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L. -L/gpfs/home/chuang3/codes/libxc-2.0.2/install/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../.. -lgfortran -lm -lquadmath

MODEXT            = mod

# ---------------------------------------------------------------------------- #

#
# Exports
#

# Algorithmic use flags
lib_algo_incs      = 
lib_algo_libs      = 

                    ########################################

# FFT use flags
lib_fft_incs       = 
lib_fft_libs       = 

                    ########################################

# Linear algebra use flags
lib_linalg_fcflags = 
lib_linalg_ldflags = 
lib_linalg_incs    = 
lib_linalg_libs    = -llapack  -lblas  

                    ########################################

# Mathematical use flags
lib_math_incs      = 
lib_math_libs      = 

                    ########################################

# Timer use flags
lib_timer_incs     = 
lib_timer_libs     =  -lrt

                    ########################################

# ETSF_IO use flags
lib_etsf_io_incs   = 
lib_etsf_io_libs   = 

# FOX use flags
lib_fox_incs       = 
lib_fox_libs       = 

# Libpspio use flags
lib_libpspio_incs  = 
lib_libpspio_libs  = 

# NetCDF use flags
lib_netcdf_incs    = 
lib_netcdf_libs    = 

# YAML use flags
lib_yaml_incs    = 
lib_yaml_libs    = 

                    ########################################

# LibXC use flags
lib_libxc_incs     = -I$(fallbacks_instdir)/include
lib_libxc_libs     = -L$(fallbacks_instdir)/lib -lxc

                    ########################################

# Fortran flags
FCFLAGS_ATOMPAW_EXT   =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_BIGDFT_EXT    =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_ETSF_IO_EXT   =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_FOX_EXT       =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_LIBPSPIO_EXT  =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_LIBXC_EXT     =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_LINALG_EXT    =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_NETCDF_EXT    =  -g -ffree-line-length-none   -O2 -mtune=native -march=native
FCFLAGS_WANNIER90_EXT =  -g -ffree-line-length-none   -O2 -mtune=native -march=native

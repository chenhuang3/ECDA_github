# To be installed:
noinst_LIBRARIES = libflib.a
module_DATA = $(modules)

# Source code.
if HAVE_FC_FLUSH
flush_files = flush.f90
flush_object= flush.o
else
flush_files = flush_fake.f90
flush_object= flush_fake.o
endif
if USE_MPI
mpi_source =
mpi_include =
mpi_inc=
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
mpi_inc = -I$(builddir) 
endif
if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
else
getarg_files = get_command_argument_fake.f90
endif
libflib_a_SOURCES = utils.c utilsadd.c \
	yaml_strings.f90 dictionaries_base.f90 dictionaries.f90 \
	callbacks.f90 \
	yaml_output.f90 mem_profiling.f90 dynamic_memory.f90 \
	razero.f90 getadd.f90 time.f90 f_lib_highlev.f90 \
	yaml_parser.c yaml_parse.f90 \
	$(flush_files) $(mpi_source) $(getarg_files)

DYN_INCLUDES =  deallocate-inc.f90 allocate-inc.f90 \
	deallocate-c-inc.f90 allocate-c-inc.f90 \
	f_malloc-inc.f90 f_malloc-extra-inc.f90 \
	deallocate-profile-inc.f90 allocate-profile-inc.f90 \
	malloc_templates-inc.f90 \
	f_malloc-base-inc.f90 f_malloc-null-inc.f90 \
	f_malloc-total-inc.f90 f_malloc-bound-inc.f90 \
	f_malloc-bounds-inc.f90 f_malloc-simple-inc.f90

YAML_INCLUDES = yaml_map-inc.f90 yaml_map-arr-inc.f90 yaml_map-mat-inc.f90
STRINGS_INCLUDES = yaml_toa-inc.f90 yaml_toa-arr-inc.f90
ERR_INCLUDES = get_err-inc.f90
DICT_INCLUDES = error_handling.f90 dict_add-inc.f90
GETADD_INCLUDES = getadd-c-inc.f90

EXTRA_DIST = $(DYN_INCLUDES) $(YAML_INCLUDES) $(STRINGS_INCLUDES) $(ERR_INCLUDES) $(DICT_INCLUDES) $(GETADD_INCLUDES)

# Module files.
if CAPITALIZE
modules = YAML_STRINGS.@MODULE_EXT@ \
	DICTIONARIES.@MODULE_EXT@ \
	DICTIONARIES_BASE.@MODULE_EXT@ \
	YAML_OUTPUT.@MODULE_EXT@ \
	YAML_PARSE.@MODULE_EXT@ \
	METADATA_INTERFACES.@MODULE_EXT@ \
	DYNAMIC_MEMORY.@MODULE_EXT@ \
	RANDOMDATA.@MODULE_EXT@ \
	MEMORY_PROFILING.@MODULE_EXT@ \
	EXCEPTION_CALLBACKS.@MODULE_EXT@
else
modules = yaml_strings.@MODULE_EXT@ \
	dictionaries.@MODULE_EXT@ \
	dictionaries_base.@MODULE_EXT@ \
	yaml_output.@MODULE_EXT@ \
	yaml_parse.@MODULE_EXT@ \
	metadata_interfaces.@MODULE_EXT@ \
	dynamic_memory.@MODULE_EXT@ \
	randomdata.@MODULE_EXT@ \
	memory_profiling.@MODULE_EXT@ \
	exception_callbacks.@MODULE_EXT@
endif

CLEANFILES = $(modules) *.@MODULE_EXT@ $(mpi_include)

AM_FCFLAGS = -I$(srcdir) $(mpi_inc) @LIBABINIT_INCLUDE@ @MPI_INCLUDE@
AM_CFLAGS = @LIB_YAML_CFLAGS@

# Dependencies
$(modules): $(libflib_a_OBJECTS)
yaml_strings.o: $(STRINGS_INCLUDES)
dictionaries.o: yaml_strings.o dictionaries_base.o $(DICT_INCLUDES) callbacks.o
#mem_profiling.o : error_handling.o
yaml_output.o : dictionaries.o $(YAML_INCLUDES)
yaml_parse.o : dictionaries.o
#error_handling.o: dictionaries.o callbacks.o yaml_output.o $(ERR_INCLUDES)
dynamic_memory.o: getadd.o yaml_strings.o yaml_output.o mem_profiling.o $(DYN_INCLUDES)
mem_profiling.o: yaml_output.o $(mpi_include)
f_lib_highlev.o: dictionaries.o yaml_output.o dynamic_memory.o yaml_parse.o
getadd.o: yaml_strings.o $(GETADD_INCLUDES)
$(flush_object): yaml_output.o dictionaries.o

time.o: $(mpi_include) #this is temporary, mpif.h should only be in wrappers

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_COMM_NULL=2,MPI_SUCCESS=0,MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1,MPI_LAND=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_REQUEST_NULL=1,MPI_STATUS_IGNORE=1" >> mpif.h &&\
	echo "integer, parameter :: mpi_tag_ub=1,mpi_address_kind=8,mpi_info_null=0" >> mpif.h &&\
	echo "integer, parameter :: mpi_mode_noprecede=0" >> mpif.h &&\
	echo "real(kind=8), external :: mpi_wtime" >> mpif.h


# Copy modules in a common include directory.
all: all-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done

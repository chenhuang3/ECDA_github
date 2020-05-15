#include <string.h>

#include <config.h>

/* Configure arguments. */
static char ARGS[] = " '--prefix=/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports' '--bindir=/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/bin' '--libdir=/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/lib' '--program-suffix=-abinit' '--disable-binaries' '--disable-bindings' '--enable-libbigdft' '--without-archives' '--with-moduledir=/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/include' '--enable-internal-libyaml' '--disable-shared' '--disable-internal-libxc' '--with-libxc-incs=-I/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/include' '--with-libxc-libs=-L/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/lib -lxc' 'CC=gcc' 'CFLAGS=-g -O2' 'LDFLAGS= ' 'LIBS=' 'CPPFLAGS= ' 'CXX=g++' 'CXXFLAGS=' 'FC=gfortran' 'FCFLAGS=-ffree-form -g -ffree-line-length-none   -O2 -mtune=native -march=native -g -O2 ' 'F77=gfortran' 'FFLAGS=-ffixed-form -g -ffree-line-length-none   -O2 -mtune=native -march=native -g -O2 ' 'CPP=gcc -E'";

/* Compilers. */
static char CC[]  = "gcc";
static char FC[]  = "gfortran";
static char CXX[] = "g++";

/* Compiler flags. */
static char CFLAGS[]   = "-g -O2";
static char FCFLAGS[]  = "-ffree-form -g -ffree-line-length-none   -O2 -mtune=native -march=native -g -O2 ";
static char CXXFLAGS[] = "";
static char CPPFLAGS[] = " ";

/* Linker stuff. */
static char LD[]        = "@LD@";
static char LDFLAGS[]   = "  -L$(top_builddir)/libABINIT/src ";
static char LIBS[]      = "-lrt ";
static char LIBS_DEPS[] = "-labinit -L/gpfs/home/chuang3/XCEP_codes/abinit-7.10.4_xcep_v2pbc/fallbacks/exports/lib -lxc     -llapack -lblas  $(top_builddir)/yaml-0.1.4/src/.libs/libyaml.a   ";

#define TRUNCATED "(...)"


void FC_FUNC_(bigdft_config_get_user_args, BIGDFT_CONFIG_GET_USER_ARGS)
     (char *args, unsigned int *ln, int args_ln)
{
  memset(args, ' ', sizeof(char) * *ln);
  if (sizeof(ARGS) <= *ln)
    memcpy(args, ARGS, sizeof(ARGS)-1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(args, TRUNCATED, sizeof(TRUNCATED)-1);
}

void FC_FUNC_(bigdft_config_get_compilers, BIGDFT_CONFIG_GET_COMPILERS)
     (char *cc, char *fc, char *cxx, unsigned int *ln, int cc_ln, int fc_ln, int cxx_ln)
{
  memset(cc, ' ', sizeof(char) * *ln);
  if (sizeof(CC) <= *ln)
    memcpy(cc, CC, sizeof(CC)-1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(cc, TRUNCATED, sizeof(TRUNCATED)-1);

  memset(fc, ' ', sizeof(char) * *ln);
  if (sizeof(FC) <= *ln)
    memcpy(fc, FC, sizeof(FC) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(fc, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(cxx, ' ', sizeof(char) * *ln);
  if (sizeof(CXX) <= *ln)
    memcpy(cxx, CXX, sizeof(CXX) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(cxx, TRUNCATED, sizeof(TRUNCATED) - 1);
}

void FC_FUNC_(bigdft_config_get_compiler_flags, BIGDFT_CONFIG_GET_COMPILER_FLAGS)
     (char *cflags, char *fcflags, char *cxxflags, char *cppflags, unsigned int *ln,
      int c_ln, int fc_ln, int cxx_ln, int cpp_ln)
{
  memset(cflags, ' ', sizeof(char) * *ln);
  if (sizeof(CFLAGS) <= *ln)
    memcpy(cflags, CFLAGS, sizeof(CFLAGS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(cflags, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(fcflags, ' ', sizeof(char) * *ln);
  if (sizeof(FCFLAGS) <= *ln)
    memcpy(fcflags, FCFLAGS, sizeof(FCFLAGS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(fcflags, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(cxxflags, ' ', sizeof(char) * *ln);
  if (sizeof(CXXFLAGS) <= *ln)
    memcpy(cxxflags, CXXFLAGS, sizeof(CXXFLAGS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(cxxflags, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(cppflags, ' ', sizeof(char) * *ln);
  if (sizeof(CPPFLAGS) <= *ln)
    memcpy(cppflags, CPPFLAGS, sizeof(CPPFLAGS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(cppflags, TRUNCATED, sizeof(TRUNCATED) - 1);
}

void FC_FUNC_(bigdft_config_get_linker, BIGDFT_CONFIG_GET_LINKER)
     (char *ld, char *ldflags, char *libs, char *libs_deps, unsigned int *ln,
      int ld_ln, int ldflags_ln, int libs_ln, int libs_deps_ln)
{
  memset(ld, ' ', sizeof(char) * *ln);
  if (sizeof(LD) <= *ln)
    memcpy(ld, LD, sizeof(LD) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(ld, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(ldflags, ' ', sizeof(char) * *ln);
  if (sizeof(LDFLAGS) <= *ln)
    memcpy(ldflags, LDFLAGS, sizeof(LDFLAGS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(ldflags, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(libs, ' ', sizeof(char) * *ln);
  if (sizeof(LIBS) <= *ln)
    memcpy(libs, LIBS, sizeof(LIBS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(libs, TRUNCATED, sizeof(TRUNCATED) - 1);

  memset(libs_deps, ' ', sizeof(char) * *ln);
  if (sizeof(LIBS_DEPS) <= *ln)
    memcpy(libs_deps, LIBS_DEPS, sizeof(LIBS_DEPS) - 1);
  else if (*ln > sizeof(TRUNCATED))
    memcpy(libs_deps, TRUNCATED, sizeof(TRUNCATED) - 1);
}

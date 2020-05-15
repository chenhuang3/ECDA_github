# -*- Autoconf -*-
#
# Copyright (C) 2005-2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

# Generated by make-macros-output on 2015/03/12 05:26:24 +0000

#
# Output macros for the "configure" script
#

#
# IMPORTANT NOTE
#
# This file has been automatically generated by the make-macros-output
# script. If you try to edit it, your changes will systematically be
# overwritten.
#



# ABI_OUTPUT()
# ------------
#
# Outputs configuration for ABINIT.
#
AC_DEFUN([ABI_OUTPUT],[
  dnl Config files
  AC_CONFIG_FILES([
    config.dump
    config.pc
    config.sh
    bindings/config.sh
    config/wrappers/wrap-fc
    fallbacks/config.mk
    special/Makefile
    special/special.env
    src/incs/Makefile
    src/mods/Makefile
    src/10_dumpinfo/m_build_info.F90
    Makefile
    src/Makefile
    src/01_gsl_ext/Makefile
    src/01_linalg_ext/Makefile
    src/01_macroavnew_ext/Makefile
    src/02_clib/Makefile
    src/10_defs/Makefile
    src/10_dumpinfo/Makefile
    src/11_memory_mpi/Makefile
    src/11_qespresso_ext/Makefile
    src/12_hide_mpi/Makefile
    src/14_hidewrite/Makefile
    src/15_gpu_toolbox/Makefile
    src/16_hideleave/Makefile
    src/18_timing/Makefile
    src/21_psiesta_noabirule/Makefile
    src/27_toolbox_oop/Makefile
    src/28_numeric_noabirule/Makefile
    src/32_util/Makefile
    src/41_geometry/Makefile
    src/41_xc_lowlevel/Makefile
    src/42_libpaw/Makefile
    src/42_nlstrain/Makefile
    src/42_parser/Makefile
    src/43_ptgroups/Makefile
    src/43_wvl_wrappers/Makefile
    src/44_abitypes_defs/Makefile
    src/45_geomoptim/Makefile
    src/47_xml/Makefile
    src/49_gw_toolbox_oop/Makefile
    src/51_manage_mpi/Makefile
    src/52_fft_mpi_noabirule/Makefile
    src/52_manage_cuda/Makefile
    src/53_ffts/Makefile
    src/53_spacepar/Makefile
    src/54_abiutil/Makefile
    src/56_io_mpi/Makefile
    src/56_mixing/Makefile
    src/56_recipspace/Makefile
    src/56_xc/Makefile
    src/57_iopsp_parser/Makefile
    src/57_iovars/Makefile
    src/61_ionetcdf/Makefile
    src/62_cg_noabirule/Makefile
    src/62_ctqmc/Makefile
    src/62_iowfdenpot/Makefile
    src/62_occeig/Makefile
    src/62_poisson/Makefile
    src/62_wvl_wfs/Makefile
    src/63_bader/Makefile
    src/64_atompaw/Makefile
    src/65_lotf_base/Makefile
    src/65_nonlocal/Makefile
    src/65_psp/Makefile
    src/66_fock/Makefile
    src/66_paw/Makefile
    src/66_wfs/Makefile
    src/67_common/Makefile
    src/68_dmft/Makefile
    src/68_lotf/Makefile
    src/68_recursion/Makefile
    src/68_rsprc/Makefile
    src/69_wfdesc/Makefile
    src/70_gw/Makefile
    src/71_bse/Makefile
    src/72_response/Makefile
    src/77_ddb/Makefile
    src/77_suscep/Makefile
    src/79_seqpar_mpi/Makefile
    src/83_cut3d/Makefile
    src/94_scfcv/Makefile
    src/95_drive/Makefile
    src/98_main/Makefile
    src/libs/Makefile
 ])

  dnl Commands
  AC_CONFIG_COMMANDS([dump-optim],[/bin/sh ${abinit_srcdir}/config/scripts/make-optim-dumper])
  AC_CONFIG_COMMANDS([script-perms],[chmod u+x config/wrappers/wrap-fc])
  AC_CONFIG_COMMANDS([long-lines],[/bin/sh ${abinit_srcdir}/config/scripts/shrink-src-files ${abinit_srcdir} ${abinit_builddir}])

  dnl Output everything
  AC_OUTPUT
]) # ABI_OUTPUT
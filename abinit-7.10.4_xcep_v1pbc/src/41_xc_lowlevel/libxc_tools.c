/*
 * Copyright (C) 2009-2014 ABINIT group (MT)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

/* ===============================================================
 * Set of functions interfacing the LibXC library.
 * (see http://www.tddft.org/programs/octopus/wiki/index.php/Libxc)
 * LibXC contains a set of FORTRAN APIs but all the functions
 * of the package are not present in it.
 * The following interfaces fill this gap.
 * ===============================================================
 */

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_DFT_LIBXC
#if !defined HAVE_DFT_LIBXC_GET_NUMBER

#define CC_FORTRAN_INT int

#include "string_f.h"

void FC_FUNC_(xc_f90_functional_get_name, XC_F90_FUNCTIONAL_GET_NAME)
     (CC_FORTRAN_INT *func_number, STR_F_TYPE func_string STR_ARG1)
{
  char *name;

  name = xc_functional_get_name(*func_number);
  if ( name == NULL )
  {
    name = (char *) malloc(256);
    sprintf(name, "unknown\0");
  }

  TO_F_STR1(name, func_string);
  free(name);
}

CC_FORTRAN_INT FC_FUNC_(xc_f90_functional_get_number, XC_F90_FUNCTIONAL_GET_NUMBER)
     (STR_F_TYPE func_string STR_ARG1)
{
  char *name;
  int ret;

  TO_C_STR1(func_string, name);

  ret = xc_functional_get_number(name);
  free(name);

  return (CC_FORTRAN_INT) ret;
}

#endif
#endif

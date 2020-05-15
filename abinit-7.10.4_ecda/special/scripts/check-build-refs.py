#!/usr/bin/env python
"check build refs"
#
# Copyright (C) 2011-2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

from time import gmtime,strftime

import commands
import os
import re
import sys

# ---------------------------------------------------------------------------- #
def abinit_test_generator():
  def test_func(abenv):
     "check build refs"
     try:
       return main(abenv.home_dir)
     except Exception:
       import sys
       raise sys.exc_info()[1] # Reraise current exception (py2.4 compliant)
  return {"test_func" : test_func}
#
# Main program
#

def main(home_dir=""):
  from os.path import join as pj
                                                                           
  # Check if we are in the top of the ABINIT source tree
  my_name = os.path.basename(__file__) + ".main"
  if ( not os.path.exists( pj(home_dir,"configure.ac") ) or
       not os.path.exists( pj(home_dir,"src/98_main/abinit.F90")) ):
    print "%s: You must be in the top of an ABINIT source tree." % my_name
    print "%s: Aborting now." % my_name
    sys.exit(1)

  # Init
  nerr = 0
  bex_diffs = list()
  bex_missing = list()

  bex_dir = pj(home_dir,"doc/build/config-examples")
  ref_dir = pj(home_dir,"special/buildsys/Refs")

  # Check files
  ref_list = os.listdir(ref_dir)
  ref_list.sort()
  for ref_file in ref_list:
    if os.path.exists("%s/%s" % (bex_dir,ref_file)):
      (ret,tmp) = commands.getstatusoutput("diff -q %s/%s %s/%s" % \
        (ref_dir,ref_file,bex_dir,ref_file))
      if ( ret != 0 ):
        bex_diffs.append(ref_file)
    else:
      bex_missing.append(ref_file)

  nerr = len(bex_diffs) + len(bex_missing)

  # Report any mismatch
  if nerr > 0:
    sys.stderr.write("%s: reporting wrongly generated build examples\n\n" % (os.path.basename(sys.argv[0])))
    sys.stderr.write("X: D=Difference detected / M=Missing File\n\n")
    sys.stderr.write("%s  %-64s\n" % ("X","File"))
    sys.stderr.write("%s  %s\n" % ("-","-" * 64))

    for bex in bex_diffs:
      sys.stderr.write("%s  %-64s\n" % ("D",bex))
    for bex in bex_missing:
      sys.stderr.write("%s  %-64s\n" % ("M",bex))

    sys.stderr.write("\n")

  return nerr

if __name__ == "__main__":
  if len(sys.argv) == 1: 
    home_dir = "."
  else:
    home_dir = sys.argv[1] 
                               
  exit_status = main(home_dir)
  sys.exit(exit_status)

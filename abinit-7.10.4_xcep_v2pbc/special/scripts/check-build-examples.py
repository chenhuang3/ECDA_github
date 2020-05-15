#!/usr/bin/env python
"check test farm build examples"
#
# Copyright (C) 2010-2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

# FIXME: detect duplicate definitions

from ConfigParser import ConfigParser
from time import gmtime,strftime

import commands
import os
import re
import sys

class MyConfigParser(ConfigParser):

  def optionxform(self,option):
    return str(option)

# ---------------------------------------------------------------------------- #

#
# Functions
#

env_ignore = list()
opt_ignore = ["fcflags_opt_*","status"]

def is_ignored(keyword):
  for opt in env_ignore + opt_ignore:
    if ( "*" in opt ):
      if ( re.match(opt,keyword) ):
        return True
    elif ( opt == keyword ):
        return True
  return False

def key_is_ok(mode,key):

  # Init keys to ignore
  cnf_ignore = dict()
  cnf_ignore["mpi"] = ("status","CC","CXX","FC")
  cnf_ignore["raw"] = ("status")
  cnf_ignore["serial"] = ("status","with_mpi_prefix")

  if ( key in cnf_ignore[mode] ):
    return False
  else:
    return True

# ---------------------------------------------------------------------------- #
def abinit_test_generator():
  def test_func(abenv):
     "check test farm build examples"
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
  re_env = re.compile("^[A-Z][0-9A-Z_]*")
  re_opt = re.compile("^[a-z][0-9a-z_]*")

  # Extract environment variables from config file
  cnf_env = MyConfigParser()
  cnf_env.read( pj(home_dir,"config/specs/environment.conf") )
  env_config = list()
  for env in cnf_env.sections():
    if cnf_env.get(env,"reset") == "no":
      if not is_ignored(env):
        env_config.append(env)
  env_config.sort()

  # Extract options from config file
  cnf_opt = MyConfigParser()
  cnf_opt.read( pj(home_dir,"config/specs/options.conf") )
  opt_config = list()
  opt_removed = list()
  for opt in cnf_opt.sections():
    tmp_sta = cnf_opt.get(opt,"status")
    if ( tmp_sta == "removed" or tmp_sta == "dropped" ):
      opt_removed.append(opt)
    elif "renamed" in tmp_sta:
      opt_removed.append(tmp_sta.split()[1])
      if not is_ignored(opt):
        opt_config.append(opt)
    else:
      if not is_ignored(opt):
        opt_config.append(opt)
  opt_config.sort()
  opt_removed.sort()

  # Extract information from build example config file
  cnf_bex = MyConfigParser()
  cnf_bex.read( pj(home_dir, "config/specs/testfarm.conf") )
  env_examples = list()
  opt_examples = list()
  env_dict = dict()
  opt_dict = dict()
  for bot in cnf_bex.sections():
    for var in cnf_bex.options(bot):
      if not is_ignored(var):
        if re_env.match(var):
          if not var in env_examples:
            env_examples.append(var)
          if not var in env_dict:
            env_dict[var] = list()
          env_dict[var].append(bot)
        elif re_opt.match(var):
          if not var in opt_examples:
            opt_examples.append(var)
          if not var in opt_dict:
            opt_dict[var] = list()
          opt_dict[var].append(bot)
  env_examples.sort()
  opt_examples.sort()

  # Compare environment and options
  denv_examples = [env for env in env_examples if not env in env_config]
  dopt_examples = [opt for opt in opt_examples if not opt in opt_config + opt_removed]
  dopt_removed = [opt for opt in opt_examples if opt in opt_removed]

  nerr = len(denv_examples) + len(dopt_examples) + len(dopt_removed)

  # Compare build examples and generated files
  bex_data = dict()
  dbex_files = list()
  for acf in os.listdir( pj(home_dir, "doc/build/config-examples") ):
    if ( re.match("bb_",acf) ):
      acf_section = re.sub("\.ac","",acf)
      if ( cnf_bex.has_section(acf_section) ):
        acf_data = file( pj(home_dir, "doc/build/config-examples/"+acf),"r").readlines()
        acf_dict = dict()
        for line in acf_data:
          line = re.sub("#.*","",line)
          line = line.strip()
          if ( len(line) > 0 ):
            idx = line.find("=")
            key = line[:idx]
            val = re.sub("\"","",line[idx+1:])
            acf_dict[key] = val
        bex_data[acf_section] = acf_dict
      else:
        dbex_files.append(acf_section)

  dbex_sections = [bot for bot in cnf_bex.sections() \
    if ( re.match("bb_",bot) and not bot in bex_data.keys() )]

  dbex_keys = dict()
  dbex_vals = dict()
  for bot in cnf_bex.sections():
    if re.match("bb_",bot):
      # Check for the presence of MPI options
      cnf_mode = "raw"
      if "with_mpi_prefix" in cnf_bex.options(bot):
        cnf_mode = "mpi"
      if "enable_mpi" in cnf_bex.options(bot):
        if cnf_bex.get(bot,"enable_mpi") == "no":
          cnf_mode = "serial"

      for var in cnf_bex.options(bot):
        if ( (not is_ignored(var)) and key_is_ok(cnf_mode,var) ):
          if ( var not in bex_data[bot].keys() ):
            if ( bot not in dbex_keys ):
              dbex_keys[bot] = list()
            dbex_keys[bot].append(var)
          elif ( cnf_bex.get(bot,var) != bex_data[bot][var] ):
            if ( bot not in dbex_vals ):
              dbex_vals[bot] = dict()
            dbex_vals[bot][var] = "%s != %s" % \
              (cnf_bex.get(bot,var),bex_data[bot][var])

  nerr += len(dbex_sections) + len(dbex_keys.keys()) + len(dbex_vals.keys())

  # Report any mismatch
  if nerr > 0:
    sys.stderr.write("%s: reporting use of wrong options\n\n" % (os.path.basename(sys.argv[0])))
    sys.stderr.write("X: R=Removed / U=Undefined / F=Missing File / K=Missing Key / V=Value Mismatch\n")
    sys.stderr.write("   >=Continued from previous line\n\n")
    sys.stderr.write("%s  %-24s  %-48s\n" % ("X","Wrong option","Bot"))
    sys.stderr.write("%s  %s  %s\n" % ("-","-" * 24,"-" * 48))

    for env in denv_examples:
      for bot in env_dict[env]:
        sys.stderr.write("%s  %-24s  %-48s\n" % ("U",env,bot))
    for opt in dopt_examples:
      for bot in opt_dict[opt]:
        sys.stderr.write("%s  %-24s  %-48s\n" % ("U",opt,bot))
    for opt in dopt_removed:
      for bot in opt_dict[opt]:
        sys.stderr.write("%s  %-24s  %-48s\n" % ("R",opt,bot))

    for bot in dbex_sections:
      sys.stderr.write("%s  %-24s  %-48s\n" % ("F","---",bot))
    for bot in dbex_keys:
      for key in dbex_keys[bot]:
        sys.stderr.write("%s  %-24s  %-48s\n" % ("K",key,bot))
    for bot in dbex_vals:
      for key in dbex_vals[bot].keys():
        val = "  " + dbex_vals[bot][key]
        sys.stderr.write("%s  %-24s  %-48s\n" % ("V",key,bot))
        sys.stderr.write("%s  %-24s  %-48s\n" % (">",val,""))

    sys.stderr.write("\n")

  return nerr

if __name__ == "__main__":
  if len(sys.argv) == 1: 
    home_dir = "."
  else:
    home_dir = sys.argv[1] 

  exit_status = main(home_dir)
  sys.exit(exit_status)

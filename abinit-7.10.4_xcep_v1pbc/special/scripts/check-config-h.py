#!/usr/bin/env python

import re
import os
import sys

# Init
re_srcfile = re.compile("\.(F|F90)$")
re_config  = re.compile(\
  "#if defined HAVE_CONFIG_H\n#include .config\.h.\n#endif\n",re.MULTILINE)

def main(top):
  retval     = 0
  for root,dirs,files in os.walk(top):
    # Sort dirs
    dirs.sort()

    # Check line lengths in Fortran source files
    for item in files:
      if re_srcfile.search(item):
        src_data = file("%s/%s" % (root,item),"r").read()
        src_count = len(re.findall(re_config,src_data))

        if src_count == 0:
          sys.stderr.write("%s/%s: missing include of config.h\n" % (root,item))
          retval = 1
        elif src_count > 1:
          sys.stdout.write("%s/%s: %d includes of config.h\n" % (root,item,src_count))

  return retval

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    top = "src"
  else:
    top = sys.argv[1] 

  exit_status = main(top)
  sys.exit(exit_status)

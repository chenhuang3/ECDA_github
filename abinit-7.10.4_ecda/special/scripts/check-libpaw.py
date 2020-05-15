#!/usr/bin/env python

import tempfile
from subprocess import Popen
import commands
import string
import glob,os
import re
import sys
from shutil import rmtree

# ---------------------------------------------------------------------------- #
def abitest(abenv, *args, **kwargs):
  return main(home_dir=abenv.home_dir)

def main(home_dir=""):
  # create tarball
  sys.stdout.write("Creating tarball...\n")

  cmd = "cd bindings/libpaw && make && cp *tar.gz /tmp"
  ou = tempfile.TemporaryFile()
  er = tempfile.TemporaryFile()

  process = Popen( cmd, shell=True, stdout=ou, stderr=er )
  process.wait()
  process.communicate()
  rc=process.returncode
  if rc != 0:
     ou.seek(0)
     er.seek(0)
     sys.stdout.write("%s\n" % ou.read())
     sys.stderr.write("%s\n" % er.read())
     ou.close()
     er.close()
     retval = 1
     return retval

  ou.close()
  er.close()
  sys.stdout.write(" done...\n")
	
  # test tarball
  sys.stdout.write("\nTesting tarball...\n")

  cmd = "cd /tmp;libpaw=`ls libpaw*.tar.gz`; tar xzf $libpaw; cd libpaw; make"
  ou = tempfile.TemporaryFile()
  er = tempfile.TemporaryFile()

  process = Popen( cmd, shell=True, stdout=ou, stderr=er )
  process.wait()
  process.communicate()
  rc=process.returncode
  if rc != 0:
     ou.seek(0)
     er.seek(0)
     sys.stdout.write("%s\n" % ou.read())
     sys.stderr.write("%s\n" % er.read())
     ou.close()
     er.close()
     retval = 1
     return retval

  ou.seek(0)
  er.seek(0)
  sys.stdout.write("%s\n" % ou.read())
# sys.stderr.write("%s\n" % er.read())
  ou.close()
  er.close()

  # cleaning
  retval = 0
  try:
    rmtree("/tmp/libpaw")
    f = glob.glob("/tmp/libpaw*.tar.gz")
    os.remove(f[0])
  except:
    sys.stderr.write("cleaning error")
    retval = 1
  
  return retval

if __name__ == "__main__":
  if len(sys.argv) == 1: 
    home_dir = "."
  else:
    home_dir = sys.argv[1] 

  exit_status = main(home_dir)
  sys.exit(exit_status)

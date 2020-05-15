#!/usr/bin/env python

import sys
import string
import os
import os.path
import glob
import commands

gnu_warnings = { # ( warning_string, warno, src_excluded )
    3  : ( 'Unused variable', ['12_hide_mpi','65_psp','68_dmft'] ),
    4  : ( 'Unused dummy argument',  ['65_psp'] ),
    5  : ( 'Nonstandard type declaration',  ['interfaces','28_numeric_noabirule','01_macroavnew_ext','01_linalg_ext'] ),
    6  : ( 'Same actual argument associated with INTENT', []),  
    7  : ( 'CHARACTER expression will be truncated in assignment',  ["57_iopsp_parser",] ),
    8  : ( 'Limit of 39 continuations exceeded',  [] ),
    9  : ( 'DOUBLE COMPLEX at (1) does not conform to the Fortran 95 standard',  ['interfaces','01_linalg_ext'] ),
    10 : ( 'at (1) defined but not used', [] ),
    11 : ( 'Character length of actual argument shorter than of dummy argument', [] ),
    #12 : ( 'may be used uninitialized',  [] ), FIXME Disabled cause it sigfaults
    13 : ( 'Obsolescent', [] ),
    14 : ( 'Type specified for intrinsic function', [] ),
    15 : ( 'Nonconforming tab character', [] ),
}

def abinit_suite_generator():

  def make_callable(wno):
    def test_func(abenv):
       try:
         return main(wno, home_dir=abenv.home_dir)
       except Exception:
         import sys
         raise sys.exc_info()[1] # Reraise current exception (py2.4 compliant)
    test_func.__doc__ = gnu_warnings[wno][0]
    return test_func

  warnos = gnu_warnings.keys()
  warnos.sort()
  for wno in warnos:
    yield {"test_func" : make_callable(wno)}

def usage():
    print "\n Usage: warningschk test_number \n "

def main(warno, home_dir=""):
  from os.path import join as pj
  debug = 0

  if not home_dir:
    cwd_dir = commands.getoutput('pwd')
    if os.path.isabs(sys.argv[0]):
        home_dir = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "../.."))
        inp_dir = os.path.join(home_dir, "tests/abirules/Input")
    else:
        inp_dir = os.path.join("..", "Input")
        home_dir = os.path.join(cwd_dir,"../../..")

  else:
    inp_dir = pj(home_dir, "tests", "abirules", "Input")

  #print "home_dir", home_dir
  warno = int(warno)
  Warning      = gnu_warnings[warno][0]
  src_excluded = gnu_warnings[warno][1]

  # read variable from file : warnings and src_excluded
  #try:
  #    test_number = pj(inp_dir, "warnings_"+str(warno)+".in")
  #except IndexError:
  #    usage()
  #    sys.exit(2)
  #except:
  #    print "unknown error : ", sys.exc_info()[0]
  #    raise
  #src_excluded = []
  #Warning = ""
  #try:
  #    f=open(test_number, 'r')
  #    exec f
  #    f.close()
  #except IOError:
  #    print "%s: no such file" % test_number
  #    sys.exit(4)
  #if Warning == "":
  #    print "Pattern not defined..."
  #    sys.exit(3)

  # header
  print "**********************************************************************"
  print "Warning pattern : '"+Warning+"'"
  print "**********************************************************************"
  #if len(src_excluded) > 0: print src_excluded

  makelog = pj(home_dir, "make.log")
  try:
    logfile = open(makelog)
  except IOError:
    raise

  words = []
  Buffer = []
  linec = 0
  warning_count = 0
  start = False
  for line in logfile:
      linec = linec + 1
      #print 'linec : %d' % linec
      if linec > 5 : Buffer.pop(0)
      Buffer.append(line)
      if start == False :
          # Examine the make.log file, starting with the section where the directory 10_defs was treated.
          if line.find("Making all in 10_defs") == -1 :
              continue
          else:
              #print linec
              start = True
      if line.find(Warning) != -1 :
          if debug:
              print "[DEBUG] Buffer[0] : " + string.strip(Buffer[0])          # source.F90:line.pos:
              print "[DEBUG] Buffer[2] : " + string.strip(Buffer[2])          # instruction
              print "[DEBUG] Buffer[4] : " + string.strip(Buffer[4])          # Warning: msg
          if True:
              if debug: print "[DEBUG] len of Buffer[0] : " + str(len(string.strip(Buffer[0])))
              if len(string.strip(Buffer[0])) != 0:
                  source = Buffer[0].split(":")[0]
                  if source.find('Included at'): source = source.split(" ")[-1]
                  sourceline = Buffer[0].split(":")[1]
                  try:
                      sourceline = sourceline.split(".")[0]
                  except IndexError:
                      pass
                  pattern = pj(home_dir, "src") + "/*/"+source
                  #pattern = '../../../src/*/'+source
                  #print pattern
                  path = glob.glob(pattern)
                  assert len(path) < 2
                  try:
                      source_dir = path[0].split('/')
                      if debug: print "[DEBUG] source_dir :" + source_dir[-2]
                      if src_excluded.index(source_dir[-2]) :
                          pass
                  except IndexError:
                      pass
                  except ValueError:
                      warning_count += 1
                      try:
                          print source + ' : var = ' + Buffer[4].split("'")[1] +' ['+source_dir[-2]+']'
                      except IndexError:
                          print source + ' : line = ' + sourceline +' ['+source_dir[-2]+']'
              else:
                  print " ***** Can't determine source but warning exists..."
          else:
              source = Buffer[4].split(":")[0]
              sourceline = Buffer[4].split(":")[1]
              pattern = pj(home_dir, "src") + "/*/"+source
              #pattern = '../../../src/*/'+source
              path = glob.glob(pattern)
              source_dir = path[0].split('/')
              if debug: print "[DEBUG] source_dir :" + source_dir[-2]
              try:
                  if src_excluded.index(source_dir[-2]) :
                      warning_count += 1
                      print string.strip(Buffer[4]) +' ['+source_dir[-2]+']'
              except ValueError:
                  pass

  logfile.close()

  # footer
  print "**********************************************************************"
  print "Warning count = " + str(warning_count)
  print "**********************************************************************"
  return warning_count


# ---------------------------------------------------------------------------
if __name__ == "__main__":

  warno = sys.argv[1]
  try:
    home_dir = os.path.abspath(sys.argv[2])
  except IndexError:
    home_dir = ""

  exit_status = main(warno, home_dir=home_dir) 
  sys.exit(exit_status)

#!/usr/bin/env python
"Check documentation and input variables"

import sys
import os
import os.path
import glob
import re

def usage():
    print "\n Usage: docchk \n "

def abinit_test_generator():
  def test_func(abenv):
     "Check documentation and input variables"
     top = abenv.apath_of("src")
     return main(abenv.home_dir)
  return {"test_func" : test_func}

def main(home_dir, verbose=False):

  home_dir = os.path.abspath(home_dir)

  # construct list of input keywords that appear in chkvars.F90
  chkvarsf90 = os.path.join(home_dir, "src/57_iovars/chkvars.F90")
  if (os.path.isfile(chkvarsf90)):
      varfile = open(chkvarsf90)
  else:
      print " \n File ", chkvarsf90," not found! "
      sys.exit(2)

  in_block = False
  words = []
  for line in varfile:
      if line.find("admitted variable names") > 0:
          in_block = True
      if line.find("Extra token") > 0:
          in_block = False
      if in_block == True and line.find("list_var") > 0:
          line_words=(line.split("'")[1]).split()
          for i in range(len(line_words)):
              words.append(line_words[i])

  if not words:
      print "Found empty list of words in %s " % chkvarsf90
      print "Perhaps someone changed the format of the file?"
      print "Please modify the code in " + __file__
      sys.exit(2)

  print " ============================================================= "
  print " ABINIT Input variables: Check in documentation                "
  print " ============================================================= "
  varhtml = glob.glob(os.path.join(home_dir, "doc/input_variables/var*html"))
  varkeyhr = glob.glob(os.path.join(home_dir, "doc/input_variables/keyhr.html"))
  ret_code = 0
  for iwords in range(len(words)):
      deffiles = []
      for ivarhtml in range(len(varhtml)):
          varhtmldata = open(varhtml[ivarhtml]).read()
          if words[iwords] in varhtmldata:
              deffiles.append(varhtml[ivarhtml])
      if len(deffiles) > 0:
          if verbose: print "SUCCESS: ",words[iwords]," appears in ",len(deffiles)," var*html files "
      else:
          print "FAIL: ",words[iwords]," does not appear in any var*html files "
          ret_code += 1
      deffiles = []
      for ivarkeyhr in range(len(varkeyhr)):
          varkeyhrdata = open(varkeyhr[ivarkeyhr]).read()
          if words[iwords] in varkeyhrdata:
              deffiles.append(varkeyhr[ivarkeyhr])
      if len(deffiles) > 0:
          if verbose: print "SUCCESS: ",words[iwords]," appears in ",len(deffiles)," keyhr.html file as well"
      else:
          print "FAIL: ",words[iwords]," does not appear in the central keyhr.html file "
          ret_code += 1

  print " ============================================================= "
  print " ABINIT Input variables: Check in test suite                   "
  print " ============================================================= "
  for iwords in range(len(words)):
      autotest = False
      for root, dirs, files in os.walk(os.path.join(home_dir, 'tests')):
          if root.find("Input")>0:
              for ifiles in range(len(files)):
                  testfilename = os.path.join(root,files[ifiles])
                  testfileinput = open(testfilename).read()
                  if words[iwords] in testfileinput:
                      autotest = True
                      break
          if autotest == True:
              break
      if autotest == True:
          if verbose: print "SUCCESS: ",words[iwords]," appears in automatic test suite "
      else: 
          print "FAIL: ",words[iwords]," does not appear in automatic test suite "
          ret_code += 1
                  
  # construct list of key words appearing in anaddb input
  invars9f90 = os.path.join(home_dir, "src/77_ddb/m_anaddb_dataset.F90")
  if (os.path.isfile(invars9f90)):
      varfile = open(invars9f90)
  else:
      print " \n File ", invars9f90," not found! "
      sys.exit(2)

  # Scan the source and search for the calls to intagm. Parse the arguments
  # and extract the name of the variable. The prototype of intagm is:
  #    call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
  re_call = re.compile(r'\s*call\s+intagm\((.+)\)\w*', re.I)

  words = []
  for line in varfile:
      m = re_call.match(line)
      if m: 
        tokens = m.group(1).split(",")
        assert len(tokens) == 9
        words.append(tokens[-3].replace("'","").replace('"',""))

  if not words:
      print "Found empty list of words in file %s" % invars9f90
      print "Perhaps someone changed the format of the file?"
      print "Please modify the code in " + __file__
      sys.exit(2)
  #print words

  print " ============================================================= "
  print " ANADDB Input variables: Check in documentation                "
  print " ============================================================= "
  varhtml = os.path.join(home_dir, "doc/users/anaddb_help.html")
  for iwords in range(len(words)):
      varhtmldata = open(varhtml).read()
      if words[iwords] in varhtmldata:
          if verbose: print "SUCCESS: ",words[iwords]," appears in ",varhtml
      else:
          print "FAIL: ",words[iwords]," does not appear ",varhtml
          ret_code += 1

  print " ============================================================= "
  print " ANADDB Input variables: Check in test suite                   "
  print " ============================================================= "
  for iwords in range(len(words)):
      autotest = False
      for root, dirs, files in os.walk(os.path.join(home_dir, 'tests')):
          if root.find("Input")>0:
              for ifiles in range(len(files)):
                  testfilename = os.path.join(root,files[ifiles])
                  testfileinput = open(testfilename).read()
                  if words[iwords] in testfileinput:
                      autotest = True
                      break
          if autotest == True:
              break
      if autotest == True:
          if verbose: print "SUCCESS: ",words[iwords]," appears in automatic test suite "
      else: 
          print "FAIL: ",words[iwords]," does not appear in automatic test suite "
          ret_code += 1

  return ret_code

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    home_dir = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "../.."))
  else:
    home_dir = sys.argv[1] 

  exit_status = main(home_dir, verbose=False)
  sys.exit(exit_status)

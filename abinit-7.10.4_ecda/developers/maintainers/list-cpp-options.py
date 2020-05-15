#!/usr/bin/env python

import os
import re
import sys

fortran = re.compile("\.([Ff]|[Ff]90)$")
cppline = re.compile("^#")
cppkeys = ("define .*","include.*","ifdef","ifndef","elif","^if ","else","endif","defined","undef","!","&&","\|\|","\(","\)")

def main(top):
  cppopts = dict()
  for root,dirs,files in os.walk(top):
    for src in files:
      if ( fortran.search(src) ):
        code = file(os.path.join(root,src),"r").readlines()

        for line in code:
          if ( cppline.match(line) ):
            line = re.sub("^#","",line).strip()
            for kw in cppkeys:
              line = re.sub(kw,"",line)

            line = line.split()
            for item in line:
              if ( item in cppopts ):
                cppopts[item] += 1
              else:
                cppopts[item] = 1

  print "Option                             Occurences"
  print "--------------------------------   ----------"
  names = cppopts.keys()
  names.sort()
  for opt in names:
    print "%-32s   %10d" % (opt,cppopts[opt])
  print ""

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    top = "src"
  else:
    top = sys.argv[1] 

  exit_status = main(top)
  sys.exit(exit_status)

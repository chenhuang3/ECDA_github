#!/usr/bin/env python

import os
import json
import yaml
import re
import sys
from htmlparser import MyHTMLParser
from variables import *
import ast

show_errors = False

def from_html_to_dict(text,vars):
    if(text.startswith("<a") and text.endswith("</a>")):
       parser = MyHTMLParser()
       parser.feed(text)
       return parser.value
    elif text.replace("[[","").replace("]]","") in vars.keys():
       return text
    else:
       return int(text)

def get_dimensions(text,vars):
    alltext = text.split(",")
    dims = []
    for dim in alltext:
       if dim.strip():
          dims.append(from_html_to_dict(dim.strip(),vars))
    return dims

def get_number(text):
    try:
        return int(text)
    except ValueError:
        return float(text.replace("d","E"))

def get_default(text,vars):
    text = text.strip()
    if(text.endswith(".")):
      text = text[:-1]
    if text.replace("[[","").replace("]]","") in vars.keys():
       return text
   
    units = None
    
    if "Ha" in text or "hartree" in text or "Hartree" in text:
        units = "Ha"
        text = text.replace("Ha","")
    elif "Bohr" in text:
        units = "Bohr"
        text = text.replace("Bohr","")
    elif "eV" in text:
        units = "eV"
        text = text.replace("eV","")
    elif "a.u" in text:
        units = "a.u."
        text = text.replace("a.u","")


    defs = []
    if "*" in text:
       tab = text.split("*")
       curval = get_default(tab[0],vars)
       if(isinstance(curval,int)):
           nelem = int(curval)
           elem = get_default(tab[1],vars)
    
           for i in range(nelem):
               defs.append(elem) 
       else:
           return text

    elif "," in text:
       alltext = text.split(",")
       for val in alltext:
          if val.strip():
              defs.append(get_default(val,vars))
    
    elif " " in text:
        alltext = text.split(" ")
        for val in alltext:
           if val.strip():
               defs.append(get_default(val,vars))
    else:
        return get_number(text)
    
    if(len(defs) == 1):
        defs=defs[0]
    elif(len(defs)==0):
        defs=None
    
    if(units is not None):
        return ValueWithUnit(value=defs,units=units)
        #return {"value":defs,"unit":units}
    else:
        return defs

file='ABINIT_variables.yml'

f = open(file,'r')
variables=yaml.load(f);
f.close();

file='log.txt'
log = open(file,'w')

file='err.txt'
err = open(file,'w')

CheckItems=[
['definition','Mnemonics: '],
['category','Characteristic: '],
['vartype','Variable type: '],
['default','Default is ']]

log.write("Remove useless text ...\n")
for i in sorted(variables.keys()):     
    for CI in CheckItems:
        key=CI[0]
        prefix=CI[1]
        leng=len(prefix)
        if key not in variables[i].keys():
            log.write(i+" : \""+key+"\" field added\n")
            variables[i][key] = 'None'
        elif variables[i][key][:leng]==prefix:
            log.write(i+" : \""+prefix+"\" removed\n")
            variables[i][key]=variables[i][key][leng:]  

listnames = []
catvars = dict()

for i in sorted(variables.keys()):
    listnames.append(i)
    try:
        catvars[variables[i]["section"]].append(i)
    except KeyError:
        catvars[variables[i]["section"]] = [i]
        
for i in sorted(variables.keys()):
    for key,value in variables[i].items():
        pattern = '<a href="([a-z]+).html#([a-zA-Z0-9_]+)">([a-zA-Z0-9_]+)</a>'
        for (cat,var,var2) in re.findall(pattern,variables[i][key]):
            if(var == var2):
                s = '<a href="'+str(cat)+'.html#'+str(var)+'">'+str(var)+'</a>'
                variables[i][key] = variables[i][key].replace(s,'[['+str(var)+']]')
        
        
alltypes=[
'real parameter',
'real parameters',
'real array',
'real arrays',
'real number',
'real numbers',
'real',
'reals',
'integer parameter',
'integer parameters',
'integer array',
'integer arrays',
'integer',
'integers',
'character string',
'character file name']

new_variables_list = []
new_variables = dict()

log.write("Compute new var type and dimensions\n")
for i in sorted(variables.keys()):
    oldtype = variables[i]['vartype']

    # Trying to find the dimensions between parenthesis
    idim_start = oldtype.find('(')
    idim_stop = oldtype.rfind(')')

    variables[i]["errordimensions"] = ""
 
    if idim_start != -1 and idim_stop != -1:
       try:
           dimensions = get_dimensions(oldtype[idim_start+1:idim_stop],variables)
       except ValueError,e:
           if show_errors:
           	print e," --- ",oldtype[idim_start+1:idim_stop]
           dimensions = None
           error = oldtype+" [error reading dimensions]"
           err.write(i+" : "+error+"\n")
           variables[i]["errordimensions"] += error
           
       if(oldtype[idim_stop+1:].strip()):
           text2 = oldtype[idim_stop+1:].strip()
           variables[i]['text'] = variables[i]['text']+" <br>From type : "+str(text2)
           error = "[ text after dimensions moved to text ! ]"
           err.write(i+" : "+oldtype+" [text after dimensions moved to text]"+"\n")
           variables[i]["errordimensions"] += error
    elif "array" in oldtype:
        dimensions = None
        error = oldtype+" [no dimension]"
        err.write(i+" : "+error+"\n")
        variables[i]["errordimensions"] += error
    else:
        dimensions = None


    # Look for arrays described in plain text
    foundtype = False
    if re.match('(real|integer) array of [1-9] elements',oldtype.strip()) :
        m=re.search('(?<=array of )[1-9](?= elements)',oldtype.strip())
        if m :
    	    foundtype = True
            nelem = int(m.group(0))
            dimensions = [nelem]
            log.write(i+" : array of X elements converted"+"\n")

    variables[i]['errortype'] = ""

    # Remove parenthesis, try to find type
    onlytype = oldtype[0:idim_start].strip()
    foundtype = True
    if not foundtype and not any(onlytype in tp for tp in alltypes):
        onlytype2 = onlytype.replace(i,"").strip()
        if not foundtype and not any(onlytype2 in tp for tp in alltypes):
            error = oldtype+" [unrecognized type]"
            err.write(i+" : "+error+"\n")
            variables[i]['errortype'] += error
            foundtype = False
    
    if foundtype:
    	if oldtype.strip().startswith('integer'): 
    	    variables[i]['vartype'] = 'integer'
    	    log.write(i+" : integer converted"+"\n")
    	elif oldtype.strip().startswith('real'):
    	    variables[i]['vartype'] = 'real'
    	    log.write(i+" : real converted"+"\n")
    	elif oldtype.strip().startswith('character') or \
    	       oldtype.strip().startswith('string'):
    	    variables[i]['vartype'] = 'string'
    	    log.write(i+" : string converted"+"\n")
    	else:
    	    variables[i]['vartype'] = oldtype
            error = oldtype+" [unknown type]"
            err.write(i+" : "+error+"\n")
            variables[i]['errortype'] += error

    #if dimensions is not None:
    variables[i]['dimensions'] = dimensions

    # Try to parse default value
    oldval = variables[i]['default']
    oldval = oldval.strip()
    oldval = oldval.replace(". ","")

    variables[i]['errordefault'] = ""

    try:
       defaultval = get_default(oldval,variables)
    except ValueError as e:
        if show_errors:
		print e," ::: ",oldval
        idim_start = oldval.find('(')
        idim_stop = oldval.find(')')
        
        if idim_start != -1 and idim_stop != -1:
           try:
               defaultval = get_default(oldval[idim_start+1:idim_stop],variables)
           except ValueError,e:
               if show_errors:
               	print e," === ",oldval[idim_start+1:idim_stop]
               defaultval = None
               error = oldval+" [error reading default]"
               err.write(i+" : "+error+"\n")
               variables[i]['errordefault'] += error
        else:
            defaultval = None
            error = oldval+" [error reading default]"
            err.write(i+" : "+error+"\n")
            variables[i]['errordefault'] += error

    variables[i]['default'] = defaultval
    
    # For the moment !
    variables[i]['range'] = None

    variables[i]['varname'] = i
  
    try:
      new_variables_list.append(Variable.from_array(variables[i]))
      new_variables[i] = Variable.from_array(variables[i])
    except:
      print "Error with variable : ",variables[i]
      sys.exit(-1)
    
log.close()
err.close()
    
file_yml = 'ABINIT_structvariables.yml'
f_yaml = open(file_yml,'w')
output = yaml.dump(new_variables_list, indent = 4, default_flow_style=False).replace('!!python/object:','!!')
f_yaml.write(output)
f_yaml.close()

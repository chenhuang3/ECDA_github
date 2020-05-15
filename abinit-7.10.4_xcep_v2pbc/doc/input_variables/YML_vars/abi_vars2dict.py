#!/usr/bin/env python

import os
import json
import yaml

def getVar(variables,varname,section,curVar,curVal):
    if(curVar != None or curVal != None):
        variables[varname][curVar]=curVal

maindir = '..'
sections=[ x[:-5] for x in os.listdir(maindir) if x[:3]=='var' and x[-5:]=='.html' ]

#sections = ["varbas"]

variables={}

curVar=None

for section in sections:
    file=maindir+'/'+section+'.html'
    rf=open(file)
    lines=rf.readlines()
    rf.close()
    numlines=len(lines)

    i=-1
    FormGroup=False
    Group=[]
    for line in lines:
        i=i+1
        if 'id="title"' in line:
            varname=line.split('"')[3]
            variables[varname]={}
            definition=None
            category=""
            vartype=""
            default="Default is "
            text=""
            TEXTini=False
            TEXTfin=False
            curVal=''
            curVar = None
            for j in range(i,numlines):
                tempVal = ''
                found = False
                for key in ['definition','category','vartype','text','default']:
                    if 'id="'+key+'"' in lines[j]:
                        curVar = key
                        tempVal = lines[j].split('id="'+key+'"')[1][1:].strip()
                        found = True
                        
                if (not found and curVar != None ):
                    tempVal = lines[j].strip()
                    
                if 'Go to the top' in lines[j]:
                    break
#                if 'id="definition"' in lines[j]:
#                    curVar = 'definition'
#                    tempVal = lines[j].split('id="definition"')[1][1:].strip()
#                elif 'id="category"' in lines[j]:
#                    curVar = 'category'
#                    tempVal = lines[j].
                        
                if tempVal[-7:]=='</font>':
                    if curVal != '':
                      curVal = curVal + ' ' + tempVal[:-7]
                    else:
                      curVal = tempVal[:-7]
                    getVar(variables,varname,section,curVar,curVal)
                    curVal = ''
                    curVar = None
                else:
                    if curVal != '':
                      curVal = curVal + ' ' + tempVal
                    else:
                      curVal = tempVal
                 
            variables[varname]['section']=section


f=open('ABINIT_variables.yml','w')
#f=open('ABINIT_variables.json','w')
#json.dump(variables, f,sort_keys = True, indent = 4)
yaml.dump(variables, f, indent=4, default_flow_style=False)
f.close()

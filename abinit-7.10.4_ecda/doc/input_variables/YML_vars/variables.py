#!/usr/bin/env python

import yaml

class Variable(yaml.YAMLObject):

    vartype = '' # String containing the type
    category = None # String containing the category
    definition = None # String containing the mnemonics
    dimensions = None # Array containing either int, formula or another variable
    defaultval = None # Either constant number, formula or another variable
    text = None # Description (str)
    varname = None # Name of the variable (str)
    errortype = None
    errordimensions = None
    errordefault = None
    section = None
    range = None

    yaml_tag = u'!variable'
    
    def attrs(self):
        return ['vartype','category','definition','dimensions','defaultval','text',
                'varname','errortype','errordimensions','errordefault','section']

    def __init__(self, vartype=None, category=None,
                definition=None,dimensions=None,default=None,
                text=None, varname=None,section=None,range=None,errordefault=None,
                errortype=None,errordimensions=None):

        self.vartype = vartype
        self.category = category
        self.definition = definition
        self.dimensions = dimensions
        self.defaultval = default
        self.text = text
        self.varname = varname
        self.section = section
        self.errortype = errortype
        self.errordimensions = errordimensions
        self.errordefault = errordefault
        self.range = range

    @classmethod
    def from_array(cls,array):
        return Variable(vartype=array["vartype"],category=array["category"],
                        definition=array["definition"],dimensions=array["dimensions"],
                        default=array["default"],text=array["text"],varname=array["varname"],
                        section=array["section"],range=array["range"],errordefault=array["errordefault"],errortype=array["errortype"],
                        errordimensions=array["errordimensions"])

    def __str__(self):
        return "Variable (default = "+str(self.defaultval)+")"


class ValueWithUnit(yaml.YAMLObject):
    
    value = None
    units = None
    yaml_tag = u'!valuewithunit'
    
    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units
        
    def __str__(self):
        return str(self.value)+" "+str(self.units)  
    
    def __repr__(self):
        return str(self)
    
def valuewithunit_representer(dumper, data):
    
    return dumper.represent_mapping('!valuewithunit',data.__dict__)
    
class Range(yaml.YAMLObject):
    
    start = None
    stop = None
    
    yaml_tag = u'!range'
    
    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop
        
    def isin(self, value):
        isin = True
        if(start is not None):
            isin = isin and (start <= value)
        if(stop is not None):
            isin = isin and (stop > value)
        return str(self)
    
    def __repr__(self):
        if(self.start is not None and self.stop is not None):
            return "["+str(self.start)+" .. "+str(self.stop)+"]"
        if(self.start is not None):
            return "["+self.start+"; ->"
        if(self.stop is not None):
            return "<-;"+self.stop+"]"
        else:
            return None

class ValueWithConditions(yaml.YAMLObject):
    
    yaml_tag = u'!valuewithconditions'
    
    def __repr__(self):
        return str(self.__dict__)
        
    def __str__(self):
        return str(self.defaultval)+" -- "+str(self.__dict__)
======================================
==    Abinit vars in YML  (2014)    ==
======================================
===    Y. Gillet, F. Abreu Araujo  ===
===          X. Gonze              ===
======================================

This folder contains different scripts and applets concerning the structured format of Abinit input variables.

Here under you can find a description of the procedure:

1) Execute abi_vars2dict.py to generate ABINIT_variables.yml

  This step reads the html files and create a YML file with the same data

2) Execute abi_dict2struct.py to generate ABINIT_structvariables.yml

  This step tries to "structure" the data from ABINIT_variables.yml in different fields and write ABINIT_structvariables.yml.
  Two useful files are also created:
    - err.txt : contains the errors encountered during the parse
    - log.txt : contains the different actions taken by the parser

3) The file 'abinit_vars.yml' is the "modified" version of 'ABINIT_structvariables.yml', so that we don't overwrite all the modifications when playing with python scripts.

Note : the specifications of the YML file are availables in 'specs.txt'.

-----------------------------
-  Abivars Java Previewer   -
-----------------------------

This Java Application can be run with the command:
  java -jar Abivars.jar

It permits an easy graphical view of the different variables
  and links between them

This application opens the file "abinit_vars.yml" by default. 
You can open any other file from the "File -> Open" menu.



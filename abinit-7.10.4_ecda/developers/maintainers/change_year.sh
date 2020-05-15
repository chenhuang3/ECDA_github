# Copyright (C) 1998-2014 ABINIT group (XG)
# 
# The purpose of this script is to change the copyright year
# in nearly all files in the ABINIT package. 
# First you should reapply the present script without changing it, because it has been seen that some people bring
# routines with erroneous date during a few months after the beginning of a new one ... !
#
# Then one should update the present script to put the current year (at present, valid for going from 2013 to 2014 !) !
#
# Then should be called from the top directory  (here a list of generic filenames, ordered on the basis of the alphanumeric string)
# developers/maintainers/change_year.sh *.ac */*.ac */*/*.am */*/*.bat */*/*/*.bat */*/*.c */*/*/*.c */*/*.cf */*/*.cnf */*/*.com */*/*.conf */*/*.cu */*/*.csh 
# developers/maintainers/change_year.sh */*/*.dep */*/*.dir */*.env */*/*.f90 */*/*.F90 */*/*/*.F90 *.in */*.in */*/*.in *.h */*/*.h */*/*.help */*/*/*.help */*/*.html */*/*/*.log */*/*.m */*/*/*.m */*/make* 
# developers/maintainers/change_year.sh */*/*.mk */*/*.m4 */*/*/*.m4 */*/Makefile */*/*/*.out */*/*.pl */*/*/*.pl */*/*.py */*/*/*.py */README */*/README */*/*.sav */*.sh */*/*.src */*/*/*.stdout */*/*.tex */*/*.txt */*/*_ 
#
# In the previous list, files without an extension are not treated (except the Makefile and README files), 
# and */*/*.sh are not treated (except tests/*/*.sh), because of conflict with the present script file extension !!
# Also config/scripts/abilint cannot be treated automatically ...
#
# So, also issue, one after the other (cut and paste the following):
# developers/maintainers/change_year.sh autom4te.cache/tr* config/scripts/a* config/scripts/clean* config/scripts/u* config/wrappers/wrap-fc
# developers/maintainers/change_year.sh developers/bzr_helpers/abinit-forge-branch developers/bzr_helpers/bzr-make-patch 
# developers/maintainers/change_year.sh developers/maintainers/change2.sh developers/maintainers/change.sh developers/various/change_perl.sh developers/various/fixed_to_free tests/cpu/Refs/changeref 
# developers/maintainers/change_year.sh developers/various/*.sh developers/various/fixed_to_free doc/config/scripts/make* doc/manpages/abinit.1 fallbacks/config/scripts/make* INSTALL 
# developers/maintainers/change_year.sh README tests/config/scripts/make-makefiles-tests tests/cpu/Refs/changeref scripts/configure/upgrade-build-config packages/debian/copyright 
# 
# Moreover, one should complement the present script with a search 
# grep 'past_year ABINIT' * */* */*/* */*/*/* */*/*/*/*
# and treat by hand the remaining files ...
#
#XG 100118 Remarked still many other problems with copyrights, by using the following command (replace 2013 by the present year !):
# grep -i opyright * */* */*/* */*/*/* */*/*/*/* | grep -v 2013 | grep -v '!! COPYRIGHT' | grep -v 'Oldenburg' | grep -v 'Stefan Goedecker' | grep -v 'doc/rel' | grep -v 'Remove' | grep -v 'tests/' | grep -v 'EXC group' | grep -v 'PWSCF group' | grep -v 'Makefile' | grep -v 'abinit.d' | grep -v 'fallbacks' | grep -v 'doc/features/features' | grep -v 'doc/install_notes/install' | grep -v 'COPYING' | grep -v 'gui'

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.yr*  
 sed -e 's&Copyright (c)&Copyright (C)&' $file > tmp.yrup
 sed -e 's&(C) 1987-2013 ABINIT&(C) 1987-2014 ABINIT&' tmp.yrup > tmp.yr87
 sed -e 's&(C) 1991-2013 ABINIT&(C) 1991-2014 ABINIT&' tmp.yr87 > tmp.yr91
 sed -e 's&(C) 1992-2013 ABINIT&(C) 1992-2014 ABINIT&' tmp.yr91 > tmp.yr92
 sed -e 's&(C) 1993-2013 ABINIT&(C) 1993-2014 ABINIT&' tmp.yr92 > tmp.yr93
 sed -e 's&(C) 1996-2013 ABINIT&(C) 1996-2014 ABINIT&' tmp.yr93 > tmp.yr96
 sed -e 's&(C) 1997-2013 ABINIT&(C) 1997-2014 ABINIT&' tmp.yr96 > tmp.yr97
 sed -e 's&(C) 1998-2013 ABINIT&(C) 1998-2014 ABINIT&' tmp.yr97 > tmp.yr98
 sed -e 's&(C) 1999-2013 ABINIT&(C) 1999-2014 ABINIT&' tmp.yr98 > tmp.yr99
 sed -e 's&(C) 2000-2013 ABINIT&(C) 2000-2014 ABINIT&' tmp.yr99 > tmp.yr00
 sed -e 's&(C) 2001-2013 ABINIT&(C) 2001-2014 ABINIT&' tmp.yr00 > tmp.yr01
 sed -e 's&(C) 2002-2013 ABINIT&(C) 2002-2014 ABINIT&' tmp.yr01 > tmp.yr02
 sed -e 's&(C) 2003-2013 ABINIT&(C) 2003-2014 ABINIT&' tmp.yr02 > tmp.yr03
 sed -e 's&(C) 2004-2013 ABINIT&(C) 2004-2014 ABINIT&' tmp.yr03 > tmp.yr04
 sed -e 's&(C) 2005-2013 ABINIT&(C) 2005-2014 ABINIT&' tmp.yr04 > tmp.yr05
 sed -e 's&(C) 2006-2013 ABINIT&(C) 2006-2014 ABINIT&' tmp.yr05 > tmp.yr06
 sed -e 's&(C) 2007-2013 ABINIT&(C) 2007-2014 ABINIT&' tmp.yr06 > tmp.yr07
 sed -e 's&(C) 2008-2013 ABINIT&(C) 2008-2014 ABINIT&' tmp.yr07 > tmp.yr08
 sed -e 's&(C) 2009-2013 ABINIT&(C) 2009-2014 ABINIT&' tmp.yr08 > tmp.yr09
 sed -e 's&(C) 2010-2013 ABINIT&(C) 2010-2014 ABINIT&' tmp.yr09 > tmp.yr10
 sed -e 's&(C) 2011-2013 ABINIT&(C) 2011-2014 ABINIT&' tmp.yr10 > tmp.yr11
 sed -e 's&(C) 2012-2013 ABINIT&(C) 2012-2014 ABINIT&' tmp.yr11 > tmp.yr12
#The next lines are both needed, as some developers decide to use one, and some the other ...
 sed -e 's&(C) 2013-2013 ABINIT&(C) 2013-2014 ABINIT&' tmp.yr12 > tmp.yr13
 sed -e 's&(C) 2013 ABINIT&(C) 2013-2014 ABINIT&' tmp.yr13 > tmp.yr
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.yr $file
 echo "file $file written "
done
rm -f tmp.yr*  
#chmod 755 */*/*.sh */*/*.py */*/*.pl */*/*.com config/*/make* developers/*/make*  */config/scripts/* */*.sh
chmod 755 */*/*.sh */*/*.pl config/*/make* */config/scripts/* */*.sh
chmod 755 config/scripts/* developers/bzr_helpers/* developers/*/* tests/cpu/Refs/changeref

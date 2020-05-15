#!/usr/bin/env python
"""Script for checking the ABINIT automatic tests."""
from __future__ import print_function, division #, unicode_literals

import sys
import os

from pprint import pprint
from optparse import OptionParser

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

#try:
#    import tests
#except ImportError:
# Add the directory [...]/abinit/tests to $PYTHONPATH
pack_dir, tail = os.path.split(os.path.abspath(__file__))
pack_dir, tail = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)

import tests
from tests import abitests
abenv = tests.abenv

__version__ = "0.1"
__author__ = "Matteo Giantomassi"


def check_authors(suite):
    def first_second_name(string):
        idx = string.rfind(".")
        if idx == -1:
            first, second = "", string
        else:
            first, second = string[:idx], string[idx+1:]

        return first.strip(), second.strip()

    second_names = []
    for test in suite:
        for string in test.authors:
            f, s = first_second_name(string)
            if not f and s and s != "Unknown":
                print("author(s) first name is missing in file %s, string = %s " %(test.full_id, s))
            second_names.append(s)
        #print(test.id, first_second_name(test.authors[0])[1])
    second_names = set(second_names)

    return second_names


def main():
    usage = "usage: %prog [suite_name] [options] [-h|--help] for help)"
    version = "%prog "+ str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="verbose mode")

    options, args = parser.parse_args()

    # Get the full database.
    # TODO should use with_disabled=True
    full_database = abitests.build_database(with_disabled=False)

    retcode = 0

    print("Test on stale or lost reference files...  ", end="")
    err = full_database.find_stale_or_lost_refs()

    if err:
        retcode += len(err)
        print("FAILED")
        sys.stderr.write(err)
    else:
        print("OK")

    print("Test on stale or lost inputs...  ", end="")
    err = full_database.find_stale_or_lost_inputs()

    if err:
        retcode += len(err)
        print("FAILED")
        sys.stderr.write(err)
    else:
        print("OK")

    print("Test on unknown keywords ...  ", end="")
    unknowns = full_database.find_unknown_keywords()

    if unknowns:
        retcode += len(unknowns)
        print("FAILED")
        print("The following keys are not documented:\n\t%s" % unknowns)
        #pprint(unknowns)
    else:
        print("OK")

    print("Test whether important TEST_INFO options are present", end="")
    full_database.check_testinfo_options()

    # Check authors.
    #print("Test wheter authors are defined", end="")
    #second_names = set()
    #for suite_name, suite in full_database.items():
    #  second_names = second_names.union( check_authors(suite) )
    #pprint(second_names)

    return retcode


if __name__ == "__main__":
    sys.exit(main())

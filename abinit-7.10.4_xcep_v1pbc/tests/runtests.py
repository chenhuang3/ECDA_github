#!/usr/bin/env python
"""This script executes the ABINIT suite of automatic tests."""
from __future__ import print_function, division #, unicode_literals

import sys
import os
import platform

from os.path import abspath as absp, basename
from warnings import warn

from optparse import OptionParser
from socket import gethostname

#try:
#    import tests
#except ImportError:
# Add the directory [...]/abinit/tests to $PYTHONPATH
# TODO: Use Absolute imports and rename tests --> abitests to 
# avoid possible conflicts with the packages in PYTHONPATH
# monty installs the subpackage paths and this breaks the import below
pack_dir, x = os.path.split(absp(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)

import logging
logger = logging.getLogger(__name__)

# TODO change name!
import tests
#print(dir(tests))
#print(tests.__path__)
abenv = tests.abenv
abitests = tests.abitests

from tests.pymods.devtools import number_of_cpus
from tests.pymods.testsuite import find_top_build_tree

JobRunner = tests.pymods.jobrunner.JobRunner
OMPEnvironment = tests.pymods.jobrunner.OMPEnvironment
TimeBomb = tests.pymods.jobrunner.TimeBomb
BuildEnvironment = tests.pymods.testsuite.BuildEnvironment

__version__ = "0.4"
__author__ = "Matteo Giantomassi"

_my_name = basename(__file__) + "-" + __version__


def str_examples():
    examples = """
Usage example (assuming the script is executed within a build tree):
\n
    runtests.py                      ==> Run the entire test suite with one python thread.
    runtests.py -j2                  ==> Run the entire test suite with two python threads.
    runtests.py v1 v2 -k abinit      ==> Run only the tests in v1,v2 containing abinit as keyword.
    runtests.py v3[:4] v4[45:] v5[3] ==> Run the tests in v3 from t1.in to t3.in, all the
                                         tests in v4 starting from t45.in, and test t3 in v5
                                         (Note Python indexing)
    runtests.py v3- v4- -k anaddb    ==> Run the anaddb tests, except those in v3 and v4
    runtests.py paral -n4 -c mpi.cfg ==> Run the paral tests with 4 MPI processes. The MPI
                                         environment is read from from file mpi.cfg
                                         (default runner is 'mpirun -n' if -c is not used)
    runtests.py -k GW perl-          ==> Run the tests containing the keyword 'GW', exclude those with
                                         the keyword 'perl'.
    runtests.py v1 -a Wall-          ==> Run the tests in v1 but exclude those contributed by the author 'Wall'
    runtests  v1 -t0 --force-mpirun  ==> Disable the timeout tool, use mpirun also for sequential runs.


Debugging mode:

    runtests.py v3[30] --gdb         ==> Run test v3[30] under the control of the GNU debugger gdb
    runtests.py paral[1] -n 2 --gdb  ==> Run test paral[1] with 2 MPI nodes under gdb 
                                         (will open 2 xterminal sessions, it may not work!)
    runtests.py v3[30] -V "memcheck" ==> Run test v3[30] with valgrind memcheck
    runtests.py v3[30] --pedantic    ==> Mark test as failed if stderr is not empty
    """
    return examples


def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg:
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def vararg_callback(option, opt_str, value, parser):
    """Callback for an option with variable arguments"""
    assert value is None
    value = []

    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2: break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg): break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)


def make_abinit(num_threads):
    """
    Find the top-level directory of the build tree and issue `make -j num_threads`.

    Returns:
        Exit status of the subprocess.
    """
    top = find_top_build_tree(".", with_abinit=False)
    cmd = "cd %s && make -j%d" % (top, num_threads)

    return os.system(cmd)


def main():

    usage = "usage: %prog [suite_args] [options]. Use [-h|--help] for help."
    version = "%prog " + str(__version__)

    class MyOptionParser(OptionParser):
        def print_help(self):
            OptionParser.print_help(self)
            print("\n" + str_examples())

    parser = MyOptionParser(usage=usage, version=version)

    parser.add_option("-c", "--cfg_file", dest="cfg_fname", type="string",
                      help="Read options from cfg FILE.", metavar="FILE")

    parser.add_option("--force-mpirun", default=False, action="store_true",
                      help="Force execution via mpiruner even for sequential jobs, i.e. np==1, defaults to False")

    parser.add_option("--use-mpiexec", default=False, action="store_true",
                      help="Replace mpirun with mpiexec (ignored if -c is provided)")

    parser.add_option("-n", "--num-mpi-processors", dest="mpi_nprocs", type="int", default=1,
                      help="<aximum number of MPI processes.")

    parser.add_option("-i", "--input-vars", dest="input_vars", type="string", default="",
                      help="String with the  variables (and values) that should be present in the input file " +
                      "Format: 'name1 value1, name2 value2, name3' " +
                      "If value is not given, a wild card is assumed" +
                      "Example: -i 'optdriver 3, getden' will execute only those tests where the " +
                      "input file contains optdriver with value 3, and the variable getden (irrespectively of its value)."
                      )

    parser.add_option("-j", "--jobs", dest="py_nthreads", type="int", default=1,
                      help="Number of python threads.")

    parser.add_option("-r", "--regenerate", dest="regenerate", default=False, action="store_true",
                      help="Regenerate the test suite database" +
                           "(use this option after any change of the input files of the test suite or any change of the python scripts.")

    parser.add_option("-k", "--keywords", dest="keys", default=[], action="callback", callback=vararg_callback,
                      help="Run the tests containing these keywords.")

    parser.add_option("-a", "--authors", dest="authors", default=[], action="callback", callback=vararg_callback,
                      help="Run the tests contributed by these developers.")

    parser.add_option("-t", "--timeout", dest="timeout_time", type="int", default=900,
                      help="Timeout value for Fortran executables (in seconds). -t 0 disables the timeout")

    parser.add_option("-b", "--build-tree", dest="build_dir_path", default="",
                      help="Path to the top level directory of the build tree.")

    parser.add_option("-d", "--dry-run", default=False, action="store_true",
                      help="Print list of tests and exit")

    parser.add_option("--gdb", action="store_true", 
                      help="Run the test(s) under the control of the GNU gdb debugger. " + 
                           "Support both sequential and MPI executions. In the case of MPI runs, " + 
                           "the script will open multiple instances of xterm (it may not work depending of your architecture).")

    parser.add_option("--nag", action="store_true", help="Option for developers")

    parser.add_option("--perf", default="", help="Use perf command to profile the test")

    parser.add_option("-s", "--show-info", dest="show_info", default=False, action="store_true",
                      help="Show information on the test suite (keywords, authors...) and exit")

    parser.add_option("-l", "--list-tests-info", dest="list_info", default=False, action="store_true",
                      help="List the tests in test suite (echo description section in ListOfFile files) and exit")

    parser.add_option("-m", "--make", dest="make", type="int", default=0,
                      help="Find the abinit build tree, and compile to code with 'make -j#MAKE' before running the tests")

    parser.add_option("-w", "--workdir", dest="workdir", type="string", default="",
                      help="Directory where the test suite results will be produced")

    parser.add_option("-o", "--omp_num-threads", dest="omp_nthreads", type="int", default=0,
                      help="Number of OMP threads to use (set the value of the env variable OMP_NUM_THREADS.\n" +
                           "Not compatible with -c. Use the cfg file to specify the OpenMP runtime variables.\n")

    parser.add_option("-p", "--patch", dest="patch", type="str", default="",
                      help="Patch the reference files of the tests with the status specified by -p. Diff tool can be specified via $PATCHER.\n" +  
                           "Examples: -p failed will patch the failed tests. Status can be concatenated by '+' e.g. -p failed+passed")

    parser.add_option("-e", "--edit", dest="edit", type="str", default="",
                      help="Edit the input files of the tests with the specified status. Use $EDITOR as editor.\n" +  
                           "Examples: -i failed will edit the input files of the the failed tests. Status can be concatenated by '+' e.g. failed+passed")

    parser.add_option("--stderr", type="str", default="",
                      help="Edit the stderr files of the tests with the specified status. Use $EDITOR as editor.\n" +
                           "Examples: --stderr failed will edit the error files of the the failed tests. Status can be concatenated by '+' e.g. failed+passed")

    parser.add_option("-v", "--verbose", dest="verbose", action="count", default=0, # -vv --> verbose=2
                      help='Verbose, can be supplied multiple times to increase verbosity')

    parser.add_option("-V", "--valgrind_cmdline", type="str", default="", 
                      help="Run test(s) under the control of valgrind.\n" +
                           "Examples: runtests.py -V memcheck or runtests.py -V 'memcheck -v' to pass options to valgrind")

    parser.add_option("--Vmem", action="store_true", 
                      help="Shortcut to run test(s) under the control of valgrind memcheck:\n"+
                           "Use --leak-check=full --show-reachable=yes --track-origins=yes")

    parser.add_option("--pedantic", action="store_true", help="Mark test(s) as failed if stderr is not empty.") 

    parser.add_option("--erase-files", dest="erase_files", type="int", default=2,
                      help=("0 => Keep all files.\n" +
                            "1 => Remove files but only if the test passed or succeeded.\n"+
                            "2 => Remove files even if the test failed.\n" +
                            "default=2\n") )

    parser.add_option("--make-html-diff", dest="make_html_diff", type="int", default=0,
                      help=("0 => Do not produce diff files in HTML format\n" +
                            "1 => Produce HTML diff but only if the test failed\n" +
                            "2 => Produce HTML diff independently of the final status of the test.\n" +
                            "default=0\n") )

    parser.add_option("--sub-timeout", dest="sub_timeout", type="int", default=30,
                      help="Timeout (s) for small subprocesses (fldiff.pl, python functions)")

    parser.add_option("--with-pickle", type="int",  default=1, 
                      help="Save test database in pickle format (default: True).")

    parser.add_option('--loglevel', default="ERROR", type="str",
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Parse command line.
    options, suite_args = parser.parse_args()
                                              
    if options.show_info:
        abitests.show_info()
        return 0

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    ncpus_detected = max(1, number_of_cpus())
    system, node, release, version, machine, processor = platform.uname()

    mpi_nprocs = options.mpi_nprocs
    omp_nthreads = options.omp_nthreads
    py_nthreads = options.py_nthreads

    if options.verbose:
        print("Running on %s -- system %s -- ncpus %s -- Python %s -- %s" % (
              gethostname(), system, ncpus_detected, platform.python_version(), _my_name))
        print("mpi_nprocs %s, omp_nthreads %s, py_nthreads %s" % (mpi_nprocs, omp_nthreads, py_nthreads))

    # Compile the code before running the tests.
    if options.make:
        retcode = make_abinit(options.make)
        if retcode: return retcode

    # Initialize info on the build. User's option has the precedence.
    build_dir_path = os.path.curdir
    if options.build_dir_path:
        build_dir_path = absp(options.build_dir_path)

    build_env = BuildEnvironment(build_dir_path)

    timeout_time = options.timeout_time

    if timeout_time > 0 and build_env.has_bin("timeout"):
        # Run executables under the control of timeout.
        timeout_path = build_env.path_of_bin("timeout")
        #timeout_signal = ""
        timebomb = TimeBomb(timeout_time, exec_path=timeout_path)
    else:
        warn("Cannot find timeout executable at: %s" % build_env.path_of_bin("timeout"))
        timebomb = TimeBomb(timeout_time)

    # ------------------------------------------------
    # Initialize the jobrunner for the (MPI|seq) mode
    # ------------------------------------------------
    if options.cfg_fname:
        # read the [mpi] and the [openmp] sections from the external cfg file.
        assert omp_nthreads == 0
        cfg_fname = options.cfg_fname
        logger.info("Initalizing JobRunner from cnf file: %s" % cfg_fname)
        runner = JobRunner.fromfile(cfg_fname, timebomb=timebomb)

    else:
        if mpi_nprocs == 1 and not options.force_mpirun:
            logger.info("Initalizing JobRunner for sequential runs.")
            runner = JobRunner.sequential(timebomb=timebomb)
        else:
            logger.info("Initalizing JobRunner assuming generic_mpi. [-c option not provided]")
            runner = JobRunner.generic_mpi(use_mpiexec=options.use_mpiexec, timebomb=timebomb)

        if omp_nthreads > 0:
            omp_env = OMPEnvironment(OMP_NUM_THREADS = omp_nthreads)
            runner.set_ompenv(omp_env)

    # Valgrind support.
    if options.valgrind_cmdline:
        runner.set_valgrind_cmdline(options.valgrind_cmdline)

    # Valgrind shortcuts.
    if options.Vmem:
        runner.set_valgrind_cmdline("memcheck --leak-check=full --show-reachable=yes --track-origins=yes")

    if runner.has_valgrind: 
        cmd = "valgrind --tool=%s " % runner.valgrind_cmdline
        print("Will invoke valgrind with:\n %s" % cmd)

    # Debugging with GNU gdb
    if options.gdb:
        runner.set_debugger("gdb")

    # Profiling with perf.
    if options.perf:
        runner.set_perf_command(options.perf)

    # Select tests according to the input variables.
    # Note that the parser is very primitive and it does not
    # have acces to the default values used by the codes.
    ivars = None
    if options.input_vars:
        ivars = {}
        string = options.input_vars
        if "," in string:
            tokens = string.split(",")
        else:
            tokens = [string,]
        for tok in tokens:
            keyval = tok.split()
            if len(keyval) == 1:
                ivars[keyval[0]] = None
            elif len(keyval) == 2:
                k, v = keyval[0], int(keyval[1])
                ivars[k] = v
            else:
                raise ValueError("Don't know how to interpret string: %s" % tok)

    try:
        test_suite = abitests.select_tests(suite_args, regenerate=options.regenerate,
                                           keys=options.keys, authors=options.authors,
                                           ivars=ivars, with_pickle=options.with_pickle)

    except Exception as exc:
        raise
        show_examples_and_exit(str(exc))

    if not len(test_suite):
        print("test_suite is empty. No test fulfills the requirements specified by the user!")

    workdir = options.workdir
    if not workdir:
        workdir = "Test_suite"

    # Create workdir.
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    else:
        print("%s directory already exists! Old files will be removed" % workdir)

    # Run the final test suites.
    if omp_nthreads == 0:
        msg = ("Running ntests = %s, MPI_nprocs = %s, py_nthreads = %s..."
               % (test_suite.full_length, mpi_nprocs, py_nthreads))
    else:
        msg = ("Running ntests = %s, MPI_nprocs = %s, OMP_nthreads %s, py_nthreads = %s..."
               % (test_suite.full_length, mpi_nprocs, omp_nthreads, py_nthreads))

    print(msg)

    if options.list_info:
        with open("ListOfTests.html", "w") as fh:
            fh.write(test_suite.make_listoftests(width=160, html=True))

        with open("ListOfTests.txt", "w") as fh:
            fh.write(test_suite.make_listoftests(width=100, html=False))

        sys.exit(1)

    if options.dry_run:
        print("Dry-run mode, print list of tests and exit(0)")
        for i, test in enumerate(test_suite):
            print("%d) %s" % (i, test))
        sys.exit(0)

    # If np > 1, use dynamic runmode.
    runmode = "static"
    if mpi_nprocs > 1: runmode = "dynamic"

    results = test_suite.run_tests(build_env, workdir, runner, 
                                   nprocs=mpi_nprocs, 
                                   nthreads=py_nthreads,
                                   runmode=runmode,
                                   erase_files=options.erase_files,
                                   make_html_diff=options.make_html_diff,
                                   sub_timeout=options.sub_timeout,
                                   pedantic=options.pedantic)

    def parse_stats(stats):
        if "+" in stats:
            stats = stats.split("+")
        elif stats == "all":
            stats = ["failed", "passed", "succeeded"]
        else:
            # String
            stats = [stats]

        # TODO Use BaseTest class attribute
        _possible_status = ["failed", "passed", "succeeded", "skipped", "disabled",]
        for s in stats:
            if s not in _possible_status:
                raise ValueError("%s is not a valid status" % s)

        return stats

    # Edit input files.
    if options.edit:
        for status in parse_stats(options.edit):
            print("Editing input files of tests with status %s" % status)
            results.edit_inputs(status=status)

    # Edit error files.
    if options.stderr:
        for status in parse_stats(options.stderr):
            print("Opening stderror files of tests with status %s" % status)
            results.inspect_stderrs(status=status)

    # Patch reference files.
    if options.patch:
        for status in parse_stats(options.patch):
            print("Patching tests with status %s" % status)
            results.patch_refs(status=status)

    if options.nag:
        for test in results.failed_tests:
            traces = test.get_backtraces()
            for trace in traces:
                print(trace)
                trace.edit_source()

    logger.info("Execution completed.")

    print("Test suite results in HTML format are available in %s" % (
        os.path.join(workdir,"suite_report.html")))

    try:
        return results.nfailed
    except AttributeError:
        return 99


if __name__ == "__main__":
    sys.exit(main())

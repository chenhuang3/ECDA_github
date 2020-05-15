"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_MPI",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: This suite contains tests executed with different numbers of MPI processes.
is_multi_parallel = True

subsuites = [
#"tgswl", 
"dfpt",
"gspw",
#"mbt",
#"moldyn",
#"string",
]

#: List of input files
inp_files = [
#"gswvl_01.in",
#"gswvl_02.in",
#"tdfpt_01.in",
#"tdfpt_02.in",
"tdfpt_03.in",
"tdfpt_04.in",
#"tdfpt_03PAW.in,"
#"tdfpt_04PAW.in",
"tgspw_01.in",
#"tgspw_02.in",
#"tgspw_03.in",
#"tgspw_04.in",
#"tgspw_05.in",
#"tmbt_1.in",
#"tmbt_2.in",
#"tmbt_3.in",
#"tmbt_4.in",
#"tmbt_5.in",
#"tmoldyn_01.in",
#"tmoldyn_02.in",
#"tmoldyn_03.in",
#"tmoldyn_04.in",
#"tmoldyn_05.in",
#"tmoldyn_06.in",
#"tmoldyn_07.in",
#"tstring_01.in",
#"tstring_02.in",
#"tstring_03.in",
#"tstring_04.in",
]

sections={
'varbas': {
'title': "BASIC variables",
'header': "ABINIT, basic input variables:",
'text': """This document lists and provides the description
of the name (keywords) of the "basic" input
variables to be used in the main input file of the abinit code.
""",
},
'vardev': {
'title': "DEVELOPMENT variables",
'header': "ABINIT, developper input variables:",
'text': """This document lists and provides the description
of the name (keywords) of the input
variables "for developpers" to be used in the
main input file of the abinit code.
""",
},
'varff': {
'title': "FINITE_FIELD variables",
'header': "ABINIT, applied finite field calculation variables:",
'text': """This document lists and provides the description
of the name (keywords) of the input
variables to be used in the main input file of the abinit code,
to define finite fields and associated polarisation (including Berry phase calculations).
""",
},
'varfil': {
'title': "FILES variables",
'header': "ABINIT, file handling input variables:",
'text': """This document lists and provides the description
of the name (keywords) of "file handling" input
variables to be used in the main input file of the abinit code.
""",
},
'vargeo': {
'title': "GEOMETRY input variables",
'header': "ABINIT, geometry builder + symmetry related input variables:",
'text': """This document lists and provides the description
of the name (keywords) of all geometry builder + symmetry related input
variables to be used in the main input file of the abinit code.
""",
},
'vargs': {
'title': "GROUND-STATE variables",
'header': "ABINIT, ground-state calculation variables:",
'text': """This document lists and provides the description
of the name (keywords) of miscellaneous (ground state-related) input
variables to be used in the main input file of the abinit code.
""",
},
'vargw': {
'title': "GW variables",
'header': "ABINIT, GW input variables:",
'text': """This document lists and provides the description
of the name (keywords) of the GW input
variables to be used in the main input file of the abinit code.

<br>The new user is advised to read first the
  <a href="../users/new_user_guide.html">new user's guide</a>,
  then the <a href="../users/abinit_help.html">abinit help file</a>,
  before reading the present file. It will be easier to discover the
  present file with the help of the two lessons
  <a href="../tutorial/lesson_gw1.html">lesson_gw1</a> and
  <a href="../tutorial/lesson_gw2.html">lesson_gw2</a> of the tutorial.
  Other input variables directly related to the GW computation are :
  <a href="varfil.html#getkss">getkss</a>&nbsp;&nbsp;
  <a href="varfil.html#getscr">getscr</a>&nbsp;&nbsp;
  <a href="vargs.html#optdriver">optdriver</a>.
  In case of parallel execution, the input variable
  <a href="varpar.html#gwpara">gwpara</a>
  is used to choose between the different parallelization levels presently implemented
  (k-points or bands).
<br>
  Note also the <B>Mrgscr</B> utility
  which allows one to merge dielectric matrices calculated at different q-points
  into a unique screening file that can be used for subsequent GW calculations.
""",
},
'varint': {
'title': "INTERNAL variable",
'header': "ABINIT internal variables:",
'text': """This document describes several important internal variables, to which
the user has no direct access through a keyword,
but that are derived from the input variables at the time
of their processing and used internally. Their
value is fixed for a specific dataset.
They are present in the
dtset array, in addition to the input variables that can be
directly addressed by the user.
""",
},
'varpar': {
'title': "PARALLELISATION variables",
'header': "ABINIT parallelisation input variables:",
'text': """This document lists and provides the description
of the name (keywords) of parallelisation input
variables to be used in the main input file of the abinit code.
""",
},
'varpaw': {
'title': "PAW variables",
'header': "ABINIT, PAW input variables:",
'text': """This document lists and provides the description
of the name (keywords) of the input
variables for runs based on the Projector Augmented Waves methodology,
to be used in the
main input file of the abinit code.
""",
},
'varrf': {
'title': "RESPONSE FUNCTION variables",
'header': "ABINIT, response function input variables:",
'text': """>This document lists and provides the description
of the name (keywords) of the response function input
variables to be used in the main input file of the abinit code.
""",
},
'varrlx': {
'title': "STRUCTURE OPTIMIZATION variables",
'header': "ABINIT, structural optimization input variables:",
'text': """This document lists and provides the description
of the name (keywords) of structural optimization input
variables to be used in the main input file of the abinit code.
""",
},
'varw90': {
'title': "Wannier90 interfac",
'header': "ABINIT Wannier90 interface input variables:",
'text': """This document lists and provides the description
of the name (keywords) of parallelisation input
variables to be used in the main input file of the abinit code.
""",
}
}

AC_INIT([bagel],[2.0.0],[paboyle@ph.ed.ac.uk])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

AC_MSG_NOTICE([Hello, world])

AC_LANG([C++])
AC_PROG_CXX([g++])
AC_PROG_RANLIB

dnl allow setting of alignment
AC_ARG_ENABLE(isize,
   AC_HELP_STRING([--enable-isize=ALIGNMENT],
	 [ Force the integer size used for pointers - necessary as cross compilers used 4byte is default]),
   [ ac_ISIZE=${enable_isize} ],
   [ ac_ISIZE=4 ]
)

AC_ARG_ENABLE(itype,
   AC_HELP_STRING([--enable-itype=type-name],
	 [ Force the integer type used to contain pointers - necessary as cross compilers used 4byte is default]),
   [ ac_ITYPE=${enable_itype} ],
   [ ac_ITYPE="uint32_t" ]
)

AC_ARG_ENABLE(issize,
   AC_HELP_STRING([--enable-issize=ALIGNMENT],
	 [ Force the integer size used for non-pointers]),
   [ ac_ISSIZE=${enable_issize} ],
   [ ac_ISSIZE=4 ]
)


AC_ARG_ENABLE(istype,
   AC_HELP_STRING([--enable-istype=type-name],
	 [ Force the integer type used for non-pointers ]),
   [ ac_ISTYPE=${enable_istype} ],
   [ ac_ISTYPE="uint32_t" ]
)

AC_ARG_ENABLE(ifmt,
   AC_HELP_STRING([--enable-ifmt=fmt],
	 [forces a llx printf format to be used for bagel ints in noarch]),
   [ ac_IFMT=${enable_ifmt} ],
   [ ac_IFMT="x" ]
)


AC_ARG_ENABLE(isfmt,
   AC_HELP_STRING([--enable-isfmt=fmt],
	 [forces a llx printf format to be used for bagel ints in noarch]),
   [ ac_ISFMT=${enable_isfmt} ],
   [ ac_ISFMT="x" ]
)

AC_CONFIG_HEADERS([lib/bagel_config_internal.h])
AC_SUBST(ac_ISIZE)
AC_SUBST(ac_ITYPE)
AC_SUBST(ac_IFMT)
AC_SUBST(ac_ISSIZE)
AC_SUBST(ac_ISTYPE)
AC_SUBST(ac_ISFMT)

AC_MSG_NOTICE([ac_IFMT is ${ac_IFMT} ])
AC_MSG_NOTICE([ac_ISIZE is ${ac_ISIZE} ])
AC_MSG_NOTICE([ac_ITYPE is ${ac_ITYPE} ])

dnl Initialise automake

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(lib/bagel_int.h)
AC_CONFIG_FILES(bagel-config)
AC_OUTPUT

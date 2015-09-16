AC_INIT([bfm],[2.0.0],[paboyle@ph.ed.ac.uk])
AC_CONFIG_AUX_DIR(config)

dnl boilerplace
AC_LANG([C++])
AC_PROG_CXX([g++])
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, [ar])

dnl init automake
AM_INIT_AUTOMAKE

dnl set up assembler for automake
AM_PROG_AS

dnl user options go to wfm_options_internal.h
AC_CONFIG_HEADERS(bfm/bfm_options.h)

dnl user must specify where the bagel installation is 
AC_ARG_WITH(bagel,
   AC_HELP_STRING([--with_bagel=BAGELDIR],
                  [ Use BAGEL installed in BAGELDIR ] ),
   [ bagel_bindir=${with_bagel}/bin
     bagel_CXX=`${bagel_bindir}/bagel-config --cxx`
     bagel_CXXFLAGS=`${bagel_bindir}/bagel-config --cxxflags`
     bagel_LDFLAGS=`${bagel_bindir}/bagel-config --ldflags`
     bagel_LIBS=`${bagel_bindir}/bagel-config --libs`
     bagel_ITYPE=`${bagel_bindir}/bagel-config --itype`
     bagel_ISTYPE=`${bagel_bindir}/bagel-config --istype`
     bagel_CXXFLAGS="${bagel_CXXFLAGS} -g"
     bagel_LDFLAGS="${bagel_LDFLAGS} -g"
     AC_SUBST(bagel_CXX)
     AC_SUBST(bagel_CXXFLAGS)
     AC_SUBST(bagel_LDFLAGS)
     AC_SUBST(bagel_LIBS)
     AC_SUBST(bagel_ITYPE)
     AC_SUBST(bagel_ISTYPE)
    ],
)


AC_MSG_NOTICE( [Bagel CXX:      ${bagel_CXX} ] )
AC_MSG_NOTICE( [Bagel CXXFLAGS: ${bagel_CXXFLAGS} ] )
AC_MSG_NOTICE( [Bagel LDFLAGS:  ${bagel_LDFLAGS} ] )
AC_MSG_NOTICE( [Bagel LIBS:     ${bagel_LIBS} ] )

dnl set the target cpu
AC_ARG_ENABLE(target-cpu,
   AC_HELP_STRING([--enable-target-cpu=CPU],
       [ Create code for target (bgl|noarch|ppc440) default is bgl ]),
       [ ac_TARGET_CPU=${enable_target_cpu} ],	
       [ ac_TARGET_CPU=yes ]
)
bagel_MUNGED_TARGET=${ac_TARGET_CPU};
AC_SUBST(bagel_MUNGED_TARGET)

dnl set the default comms type
AC_ARG_ENABLE(comms,
   AC_HELP_STRING([--enable-comms=COMMSTYPE],
       [ Default comms type (QMP|IROIRO|FAKE) default is FAKE ]),
       [ ac_COMMS=${enable_comms} ],	
       [ ac_COMMS=FAKE ]
)

dnl override dslash to SPI comms
AC_ARG_ENABLE(spidslash,
   AC_HELP_STRING([--enable-spidslash=yes/no],
       [ Use bgq spi calls in dslash ]),
       [ ac_SPIDSLASH=${enable_spidslash} ],	
       [ ac_SPIDSLASH=no ]
)
case ${ac_SPIDSLASH} in
     no)
     echo BFM OPTIONS: BGQ spi dslash DISABLED
     ;;
     yes)
     AC_DEFINE([BFM_BGQ_SPI_DSLASH],[1],[BFM_BGQ_SPI_DSLASH enabled] )
     echo  BFM OPTIONS: BGQ spi dslash is ENABLED
     ;;   
     *)
     echo  BFM OPTIONS: BGQ spi dslash disabled
     ;;
esac


dnl
AC_ARG_ENABLE(chroma-regression,
   AC_HELP_STRING([--enable-chroma-regression],
                  [ Check solutions against chroma. Requires chroma linking and availability ] ),
   [ ac_chroma_regression=${enable_chroma_regression} ],
   [ ac_chroma_regression=no ]
)
case ${ac_chroma_regression} in
     no)
     echo  BFM OPTIONS: Chroma regression DISABLED
     ;;
     yes)
     AC_DEFINE([BFM_CHROMA_BINDING],[1],[BFM_CHROMA_BINDING enabled] )
     echo  BFM OPTIONS: Chroma regression ENABLED
     ;;   
     *)
     echo  BFM OPTIONS: Chroma regression DISABLED
     ;;
esac

dnl set the thread model
AC_ARG_ENABLE(thread-model,
   AC_HELP_STRING([--enable-thread-model=single|openmp|spi],
       [ Use threading model default is single thread ]),
       [ ac_THREAD_MODEL=${enable_thread_model} ],	
       [ ac_THREAD_MODEL=single ]
)
case ${ac_THREAD_MODEL} in
     single)
     echo BFM OPTIONS: SINGLE threaded model
       AC_DEFINE([THREAD_MODEL_SINGLE],[1],[single threaded] )
     ;;   
     openmp)
       echo BFM OPTIONS: OpenMP thread model
       AC_DEFINE([THREAD_MODEL_OPENMP],[1],[openmp threads] )
     ;;
     spi)
       echo BFM OPTIONS: BGQ spi thread model
       AC_DEFINE([THREAD_MODEL_SPI],[1],[bgq spi] )
     ;;
     *)
     AC_MSG_ERROR([${ac_THREAD_MODEL} unsupported --enable-thread-model option]);
     ;;
esac


dnl user must specify where the bagel installation is 
AC_ARG_ENABLE(qdp,
   AC_HELP_STRING([--enable-qdp],
                  [ Compile the QDP binding for BAGEL BFM ] ),
   [ ac_QDP=${enable_qdp} ],
   [ ac_QDP=yes ]
)
case ${ac_QDP} in
     no)
     echo BFM OPTIONS: QDP binding is DISABLED
     CHROMA_LIBS=
     ;;   
     yes)
       echo BFM OPTIONS: QDP binding is ENABLED
       AC_DEFINE([BFM_QDP_BINDING],[1],[BFM_QDP_BINDING enabled] )
       CHROMA_LIBS=`chroma-config --libs`
       echo Chroma Libraries $CHROMA_LIBS
     ;;
     *)
     AC_MSG_ERROR([BAGEL ${ac_QDP} unsupported --enable-qdp option]);
     ;;
esac

dnl Enable qpx code paths through xlc
AC_ARG_ENABLE(xlcqpx,
   AC_HELP_STRING([--enable-xlcqpx],
                  [ Compile the QPX xlc codes ] ),
   [ ac_XLCQPX=${enable_xlcqpx} ],
   [ ac_XLCQPX=no ]
)
case ${ac_XLCQPX} in
     yes)
     AC_DEFINE([USE_XLC_OPTIMISED_CODE],[1],[XLC OPTIMISED CODE PATHS ENABLED])
     echo "Using XLC QPX code paths"
     ;;   
     no)
     echo "NOT using XLC QPX code paths"
     ;;
     *)  
     ;;
esac

AM_CONDITIONAL(BUILD_XLCQPX,[ test  "X${ac_XLCQPX}X" = "XyesX"])

dnl Set compile flags
case ${ac_COMMS} in
     QMP|qmp)
       AC_DEFINE([BFM_QMP],[1],[ NOARCH target] )
       echo QMP comms
     ;;
     IROIRO|iroiro)
       AC_DEFINE([BFM_IROIRO],[1],[ NOARCH target] )
       echo IROIRO comms
     ;;
     spi|SPI)
       AC_DEFINE([BFM_QMP],[1],[ NOARCH target] )
       AC_DEFINE([BFM_BGQ_SPI_DSLASH],[1],[ NOARCH target] )
       echo spi comms
     ;;
     FAKE|fake)
       AC_DEFINE([BFM_FAKE],[1],[ NOARCH target] )
       echo FAKE comms
     ;;
      *)
     AC_MSG_ERROR([BAGEL Target comms ${ac_COMMS} unsupported ]);
     ;;
esac

AM_CONDITIONAL(BUILD_OPENMP,[ test  "X${ac_THREAD_MODEL}X" = "XopenmpX"])
AM_CONDITIONAL(BUILD_SPI,[ test  "X${ac_THREAD_MODEL}X" = "XspiX"])

AM_CONDITIONAL(BUILD_CHROMA_CONFIG,[ test  "X${ac_chroma_regression}X" = "XyesX"])
AM_CONDITIONAL(BUILD_QDP,[ test  "X${ac_QDP}X" = "XyesX"])

AM_CONDITIONAL(BUILD_COMM_QMP,[ test  "X${ac_COMMS}X" == "XQMPX" ])
AM_CONDITIONAL(BUILD_COMM_IROIRO,[ test  "X${ac_COMMS}X" == "XIROIROX" ])
AM_CONDITIONAL(BUILD_COMM_SPI,[ test  "X${ac_COMMS}X" == "XspiX" ])
AM_CONDITIONAL(BUILD_COMM_FAKE,[ test  "X${ac_COMMS}X" == "XFAKEX" ])

AM_CONDITIONAL(BUILD_BFM_BGL,[ test "X${ac_TARGET_CPU}X" == "XbglX" ])
AM_CONDITIONAL(BUILD_BFM_BGQ,[ test "X${ac_TARGET_CPU}X" == "XbgqX" ])
AM_CONDITIONAL(BUILD_BFM_KNC,[ test "X${ac_TARGET_CPU}X" == "XkncX" ])
AM_CONDITIONAL(BUILD_BFM_NOARCH,[ test "X${ac_TARGET_CPU}X" == "XnoarchX" ])
AM_CONDITIONAL(BUILD_BFM_NOARCHBGQ,[ test "X${ac_TARGET_CPU}X" == "XnoarchbgqX" ])


bagel_NSIMD=`${bagel_bindir}/bagel-config --nsimd ${bagel_MUNGED_TARGET}`
AC_DEFINE_UNQUOTED([BAGEL_NSIMD],${bagel_NSIMD},[Width of Complex SIMD arithmetic])
AC_SUBST(bagel_NSIMD)
AC_SUBST(CHROMA_LIBS)

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(bfm/Makefile)
AC_CONFIG_FILES(bfm/bagel_int.h)
AC_CONFIG_FILES(code_gen/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(tests/wilson_nochroma/Makefile)
AC_CONFIG_FILES(tests/wilson_vec/Makefile)
AC_CONFIG_FILES(tests/wilson_cg/Makefile)
AC_CONFIG_FILES(tests/wilson/Makefile)
AC_CONFIG_FILES(tests/clover/Makefile)
AC_CONFIG_FILES(tests/clover_threaded/Makefile)
AC_CONFIG_FILES(tests/clover_cg/Makefile)
AC_CONFIG_FILES(tests/wilson_pthread_cg/Makefile)
AC_CONFIG_FILES(tests/wilson_singlenode/Makefile)
AC_CONFIG_FILES(tests/wilson_firmware_singlenode/Makefile)
AC_CONFIG_FILES(tests/wilson_firmware_singlenode/prepare.sh)
AC_CONFIG_FILES(tests/wilson_pthread_singlenode/Makefile)
AC_CONFIG_FILES(tests/wilson_cg_threaded/Makefile)
AC_CONFIG_FILES(tests/wilsontm_cg_threaded/Makefile)
AC_CONFIG_FILES(tests/dwf_singlenode/Makefile)
AC_CONFIG_FILES(tests/dwf_firmware_singlenode/Makefile)
AC_CONFIG_FILES(tests/dwf_firmware_singlenode/prepare.sh)
AC_CONFIG_FILES(tests/dwf_firmware_cg/Makefile)
AC_CONFIG_FILES(tests/dwf_firmware_cg/prepare.sh)
AC_CONFIG_FILES(tests/dwf_nochroma/Makefile)
AC_CONFIG_FILES(tests/dwf_cg_nochroma/Makefile)
AC_CONFIG_FILES(tests/dwf_mooee/Makefile)
AC_CONFIG_FILES(tests/dwf_mooeeinv/Makefile)
AC_CONFIG_FILES(tests/dwf_prec/Makefile)
AC_CONFIG_FILES(tests/dwf_cg/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_dperp/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_unprec/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_unprec_spi/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_cg/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_cg_threaded/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_cg_threaded_spi/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_prop/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_comms/Makefile)
AC_CONFIG_FILES(tests/dwf_rb5d_multishift/Makefile)
AC_CONFIG_FILES(tests/dwf_4d_operator/Makefile)
AC_CONFIG_FILES(tests/dwf_4d_mres/Makefile)
AC_CONFIG_FILES(tests/dwf_4d_bulk/Makefile)
AC_CONFIG_FILES(tests/dwf_cg_fpga/Makefile)
AC_CONFIG_FILES(tests/linalg/Makefile)
AC_CONFIG_FILES(tests/majority_vote/Makefile)
AC_CONFIG_FILES(tests/eguchi_kawamoto/Makefile)
AC_CONFIG_FILES(tests/eguchi_kawamoto_cg/Makefile)
AC_CONFIG_FILES(tests/general_5d_init_test/Makefile)
AC_CONFIG_FILES(tests/general_5d_munprec/Makefile)
AC_CONFIG_FILES(tests/general_5d_mprec/Makefile)
AC_CONFIG_FILES(tests/general_5d_cg/Makefile)
AC_CONFIG_FILES(tests/general_5d_cg_prec/Makefile)
AC_CONFIG_FILES(tests/general_5d_mres/Makefile)
AC_CONFIG_FILES(tests/general_5d_forceterm/Makefile)
AC_CONFIG_FILES(tests/general_5d_deltaL/Makefile)
AC_CONFIG_FILES(tests/general_5d_deltaL_nersc/Makefile)
AC_CONFIG_FILES(tests/inexact_deflation/Makefile)
AC_CONFIG_FILES(tests/inexact_deflation_qdp/Makefile)
AC_CONFIG_FILES(tests/ibm_qcd_qpx_diag/Makefile)
AC_CONFIG_FILES(tests/benchmark/Makefile)
AC_CONFIG_FILES(tests/dwf_deflated_cg/Makefile)
AC_CONFIG_FILES(tests/BfmMultiGrid/Makefile)
AC_CONFIG_FILES(tests/BfmMultiGridDev/Makefile)

AC_OUTPUT

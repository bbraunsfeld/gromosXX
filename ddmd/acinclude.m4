dnl Local function, as GROMOSXX depends
dnl on the STL. -> Does the C++ compiler
dnl support the STL to the degree necessary?
dnl 
AC_DEFUN([AC_CV_CXX_VERSION_OK],
  [AC_CACHE_CHECK(whether the compiler supports the STL,
   ac_cv_cxx_version_ok,
     [AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([#include <vector>],[
	std::vector<int> v; v.push_back(0);return 0;],
         ac_cv_cxx_version_ok=yes, ac_cv_cxx_version_ok=no)
      AC_LANG_RESTORE
  ])
])

AC_DEFUN([AC_PROG_CXX_MWERKS],
[AC_CACHE_CHECK(whether we are using Metrowerks Codewarrior C++, MWERKS_CXX,
[cat > conftest.c <<EOF
#ifdef __MWERKS__
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  MWERKS_CXX=yes
  compiler=mwerks
else
  MWERKS_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_MSVCPP],
[AC_CACHE_CHECK(whether we are using Visual C++, MSVCPP_CXX,
[cat > conftest.c <<EOF
#if defined(_MSC_VER) && !defined(__MWERKS__) && !defined(__ICL)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  MSVCPP_CXX=yes
  compiler=msvcpp
else
  MSVCPP_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_SGICC],
[AC_CACHE_CHECK(whether we are using SGI MIPSpro C++, SGI_CXX,
[cat > conftest.c <<EOF
# if defined(__sgi) && !defined(__GNUC__)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  SGI_CXX=yes
  compiler=sgicc
else
  SGI_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_SUNCC],
[AC_CACHE_CHECK(whether we are using Sun C++, SUN_CXX,
[cat > conftest.c <<EOF
# if defined(__SUNPRO_CC) 
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  SUN_CXX=yes
  compiler=suncc
else
  SUN_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_INTELCC],
[AC_CACHE_CHECK(whether we are using Intel C++, INTEL_CXX,
[cat > conftest.c <<EOF
# if defined(__ICC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  INTEL_CXX=yes
  compiler=intelcc
else
  INTEL_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_KAICC],
[AC_CACHE_CHECK(whether we are using KAI C++, KAI_CXX,
[cat > conftest.c <<EOF
# if defined(__KCC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  KAI_CXX=yes
  compiler=kaicc
else
  KAI_CXX=no
fi])])

AC_DEFUN([AC_BZ_SET_COMPILER],
  [cxxwith=`echo $1 | sed -e 's/ /@/'`
   case "$cxxwith" in
     *:*@*)                 # Full initialization syntax
       CXX=`echo "$cxxwith" | sed  -n -e 's/.*:\(.*\)@.*/\1/p'`
       CXXFLAGS=`echo "$cxxwith" | sed  -n -e 's/.*:.*@\(.*\)/\1/p'`
     ;;
     *:*)                   # Simple initialization syntax
       CXX=`echo "$cxxwith" | sed  -n -e 's/.*:\(.*\)/\1/p'`
       CXXFLAGS="$CXXFLAGS $3"
     ;;
     *)                     # Default values
       CXX=$2
       CXXFLAGS="$CXXFLAGS $3"
     ;;
   esac])

dnl Determine a Fortran 77 compiler to use.  If `F77' is not already set
dnl in the environment, check for `g77', `f77' and `f2c', in that order.
dnl Set the output variable `F77' to the name of the compiler found.
dnl 
dnl If using `g77' (the GNU Fortran 77 compiler), then `AC_PROG_F77'
dnl will set the shell variable `G77' to `yes', and empty otherwise.  If
dnl the output variable `FFLAGS' was not already set in the environment,
dnl then set it to `-g -02' for `g77' (or `-O2' where `g77' does not
dnl accept `-g').  Otherwise, set `FFLAGS' to `-g' for all other Fortran
dnl 77 compilers.
dnl 
dnl AC_PROG_F77()
AC_DEFUN([AC_MTL_PROG_F77],
[AC_BEFORE([$0], [AC_PROG_CPP])dnl
if test -z "$F77"; then
  AC_CHECK_PROGS(F77, g77 f77 f2c)
    test -z "$F77" && AC_MSG_WARN([no acceptable Fortran 77 compiler found in \$PATH])
fi

AC_PROG_F77_WORKS
AC_PROG_F77_GNU

if test $ac_cv_prog_g77 = yes; then
  G77=yes
dnl Check whether -g works, even if FFLAGS is set, in case the package
dnl plays around with FFLAGS (such as to build both debugging and
dnl normal versions of a library), tasteless as that idea is.
  ac_test_FFLAGS="${FFLAGS+set}"
  ac_save_FFLAGS="$FFLAGS"
  FFLAGS=
  AC_PROG_F77_G
  if test "$ac_test_FFLAGS" = set; then
    FFLAGS="$ac_save_FFLAGS"
  elif test $ac_cv_prog_f77_g = yes; then
    FFLAGS="-g -O2"
  else
    FFLAGS="-O2"
  fi
else
  G77=
  test "${FFLAGS+set}" = set || FFLAGS="-g"
fi
])

dnl @synopsis AC_CXX_LIB_BLITZ([optional-string "required"])
dnl
dnl Check whether Blitz++ is installed.
dnl Blitz++ is available at http://oonumerics.org/blitz.
dnl
dnl   Set the path for Blitz++  with the option
dnl      --with-blitz[=DIR]
dnl   Blitz headers should be under DIR/include
dnl   Blitz library should be under DIR/lib
dnl   Then try to compile and run a simple program with a Blitz Array
dnl   Optional argument `required' triggers an error if Blitz++ not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_BLITZ],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Blitz++ template library
available at http://oonumerics.org/blitz
When installed give the directory of installation with the option
  --with-blitz@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_BLITZ],[

AC_ARG_WITH(blitz,
AS_HELP_STRING([--with-blitz@<:@=DIR@:>@],[Set the path for Blitz++]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_BLITZ
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
	fi
	LIBS="-lblitz"

	AC_CACHE_CHECK([whether Blitz++ is installed],ac_cxx_lib_blitz,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <blitz/array.h>
]],[[
blitz::Array<int,1> x(10);
x = blitz::tensor::i;
	]])],[ac_cxx_lib_blitz=yes],[ac_cxx_lib_blitz=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_blitz" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="-I$withval/include"
			BLITZ_LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
dnl			AC_SUBST(BLITZ_CPPFLAGS)
			AC_SUBST(BLITZ_LDFLAGS)
		fi
		BLITZ_LIB="-lblitz"
		AC_SUBST(BLITZ_LIB)
	else
		if test "$1" = required ; then
			AC_MSG_ERROR_BLITZ
		fi
	fi

fi

])

dnl @synopsis AC_CXX_LIB_GSL([optional-string "required"])
dnl
dnl Check whether Gnu Scientific Library (GSL) is installed.
dnl GSL is available from
dnl www.gnu.org
dnl
dnl   Set the path for GSL  with the option
dnl      --with-gsl[=DIR]
dnl   GSL headers should be under DIR/include
dnl   GSL library should be under DIR/lib
dnl   Then try to compile and run a simple program with a gsl random number
dnl   Optional argument `required' triggers an error if GSL not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_GSL],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Gnu Scientific Library (GSL)
When installed give the directory of installation with the option
  --with-gsl@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_GSL],[

AC_ARG_WITH(gsl,
AS_HELP_STRING([--with-gsl@<:@=DIR@:>@],[Set the path for GSL]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_GSL
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib"
	fi
	LIBS="-lgsl -lgslcblas"

	AC_CACHE_CHECK([whether Gnu Scientific Library is installed],ac_cxx_lib_gsl,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <gsl/gsl_matrix.h>
]],[[
gsl_matrix * mat = gsl_matrix_alloc(3,3);
gsl_matrix_set_zero(mat);
gsl_matrix_free(mat);
	]])],[ac_cxx_lib_gsl=yes],[ac_cxx_lib_gsl=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_gsl" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="$CPPFLAGS -I$withval/include"
			GSL_LDFLAGS="-L$withval/lib"
			AC_SUBST(GSL_LDFLAGS)
		fi
		GSL_LIBS="-lgsl -lgslcblas"
		AC_SUBST(GSL_LIBS)

   		AC_DEFINE_UNQUOTED([HAVE_GSL],[],[Gnu Scientific Library])

	else
		if test "$1" = required ; then
			AC_MSG_ERROR_GSL
		fi
	fi

fi

])


dnl check for GSL
AC_DEFUN([AM_PATH_GSL],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
  min_gsl_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
  no_gsl=""
  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
  else
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
       gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
       gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
       gsl_micro_version=0
    fi

    if test "x$enable_gsltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSL_CFLAGS"
      LIBS="$LIBS $GSL_LIBS"

      rm -f conf.gsltest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gsltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gsl_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gsl_version");
     exit(1);
   }

   if (($gsl_major_version > major) ||
      (($gsl_major_version == major) && ($gsl_minor_version > minor)) ||
      (($gsl_major_version == major) && ($gsl_minor_version == minor) && ($gsl_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      printf("\n*** 'gsl-config --version' returned %d.%d.%d, but the minimum version\n", $gsl_major_version, $gsl_minor_version, $gsl_micro_version);
      printf("*** of GSL required is %d.%d.%d. If gsl-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If gsl-config was wrong, set the environment variable GSL_CONFIG\n");
      printf("*** to point to the correct copy of gsl-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      exit(1);
    }
}

],, no_gsl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gsl" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GSL_CONFIG" = "no" ; then
       echo "*** The gsl-config script installed by GSL could not be found"
       echo "*** If GSL was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GSL_CONFIG environment variable to the"
       echo "*** full path to gsl-config."
     else
       if test -f conf.gsltest ; then
        :
       else
          echo "*** Could not run GSL test program, checking why..."
          CFLAGS="$CFLAGS $GSL_CFLAGS"
          LIBS="$LIBS $GSL_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GSL or finding the wrong"
          echo "*** version of GSL. If it is not finding GSL, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GSL was incorrectly installed"
          echo "*** or that you have moved GSL since it was installed. In the latter case, you"
          echo "*** may want to edit the gsl-config script: $GSL_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     ifelse([$3], , :, [$3])
  fi
  GSL_CXXFLAGS="$GSL_CFLAGS"
  AC_SUBST(GSL_CXXFLAGS)
  AC_SUBST(GSL_LIBS)
dnl just add it...
  CXXFLAGS="$CXXFLAGS $GSL_CXXFLAGS"
  rm -f conf.gsltest
])


AC_INIT(configure.ac)

AC_CANONICAL_HOST
echo Host: $host_os
case $host_os in
	*mingw* ) CPPFLAGS="$CPPFLAGS -DWINDOWS_PATHS";;
esac

#This definition is used when autoconf is invoked
#It does not need to be correct when running the created configure
m4_define(M4LIBCMATRIX,../libcmatrix)

if test "x$CPPFLAGS" = x; then
	CPPFLAGS="-I."
fi
m4_include(M4LIBCMATRIX/configure-compiler.ac)

AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(int)

if test "$ac_cv_sizeof_size_t" != "$ac_cv_sizeof_int"; then
AC_MSG_NOTICE([size_t and unsigned int differ in size.])
CPPFLAGS="$CPPFLAGS -DNEED_SEPARATE_COUNT_T"
else
AC_MSG_NOTICE([size_t and unsigned int are same size.])
fi 

# Moved forward - autoconf weirdness
# Add note to prevent libcmatrix being added to LIBS
AC_CHECK_LIB(cmatrix,main,AC_MSG_NOTICE([libcmatrix library found]),[AC_MSG_ERROR([libcmatrix library not found
Check LDFLAGS is defined correctly])])

m4_include(M4LIBCMATRIX/configure-libraries.ac)

AC_CHECK_HEADERS(MetaPropagation.h,,AC_MSG_FAILURE([libcmatrix header files (e.g. MetaPropagation.h) not found
Check CPPFLAGS is defined to include the relevant directory]))

AC_CHECK_HEADERS(boost/spirit.hpp,,AC_MSG_FAILURE([boost::spirit header files not found
Check full Boost or Spirit is installed and CPPFLAGS is defined to include the directory containing boost]))

AC_CHECK_HEADERS(CrystalSystem.h,USEPERIODIC="yes",USEPERIODIC="no")
if test $USEPERIODIC = no; then
	CXXFLAGS="$CXXFLAGS -DNOPERIODIC";
	AC_MSG_WARN([Periodic spin system support not enabled (missing CrystalSystem.h)])
fi

AC_CHECK_HEADERS(sys/resource.h,CXXFLAGS="$CXXFLAGS -DHAVE_SYS_RESOURCE")

dnl This check is currently unreliable. Assuming C++11 behaviour
dnl AC_CHECK_DECL(std::isnan,CXXFLAGS="$CXXFLAGS -DISNAN_IN_STD",,[#include <cmath>])

AC_ARG_ENABLE(static,[  --enable-static  use static linkage],LDFLAGS="$LDFLAGS -static -static-libgcc -static-libstdc++",)
AC_ARG_ENABLE(CXX11,[  --disable-CXX11  disable C++11 features],,CPPFLAGS="$CPPFLAGS -std=c++11 -DNMRSIM_USE_CXX11")

dnl AC_CHECK_HEADERS(dlfcn.h,USEMODULE="yes",USEMODULE="no")
dnl #AC_ARG_WITH(module,[  --with-module  enable loadable modules],USEMODULE="yes",USEMODULE="no")
dnl if test $USEMODULE = yes; then
dnl 	CXXFLAGS="$CXXFLAGS -DENABLE_DL"
dnl	AC_CHECK_LIB(dl,dlopen,LDFLAGS="$LDFLAGS -ldl",AC_MSG_WARN([libdl not found
dnl LDFLAGS may not be set correctly]))
dnl else
dnl 	AC_MSG_WARN([dlopen header file (dlfcn.h) not found 
dnl Loadable modules will not be supported
dnl or check CPPFLAGS is defined to include the relevant directory])
dnl fi

AC_CHECK_FUNC(strtok_r,CXXFLAGS="$CXXFLAGS -DHAVE_STRTOK_R")
AC_CHECK_FUNC(setenv,,CXXFLAGS="$CXXFLAGS -DNOSETENV")

dnl AC_ARG_WITH(debugging,
dnl 	AS_HELP_STRING([--with-debugging],[enable debugging (default=no)]),,with_debugging=no)
dnl if test x$with_debugging = yes; then
dnl 	DEBUGFLAGS="-g"
dnl #	POSTBUILD=
dnl 	AC_MSG_NOTICE([Debugging: yes])
dnl else
dnl 	DEBUGFLAGS="-DNDEBUG"
dnl #	POSTBUILD="strip \$@@EXEEXT@"
dnl 	AC_MSG_NOTICE([Debugging: no])
dnl fi
dnl AC_SUBST(DEBUGFLAGS)
#AC_SUBST(POSTBUILD)

if test x$FOUNDMINUIT = xyes; then
	CXXFLAGS="$CXXFLAGS -DUSE_MINUIT"
fi
if test x$NEEDG2C = xyes; then
	LIBS="$LIBS -lg2c"
fi

#AC_LIBLTDL_CONVENIENCE
#AC_SUBST(LTDLINCL)
#AC_SUBST(LIBLTDL)
#AC_LIBTOOL_DLOPEN
#AC_LIBTOOL_WIN32_DLL
#AC_PROG_LIBTOOL
#AC_SUBST(LIBTOOL_DEPS)
#AC_CONFIG_SUBDIRS(libltdl)

AC_SUBST(OPTFLAGS)
AC_SUBST(LDFLAGS)
AC_OUTPUT(Makefile extras/Makefile)

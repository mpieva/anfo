#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([anfo], [svn], [udo_stenzel@eva.mpg.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/util.cc])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PATH_PROG([PROTOC], [protoc])
AC_ARG_VAR([PROTOC], [location of protobuf compiler])

# Checks for libraries.
AC_SEARCH_LIBS([exp], [m])
AC_SEARCH_LIBS([log1p], [m])
AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthread=1])
AC_CHECK_LIB([Judy], [Judy1Set])
AC_CHECK_LIB([z], [deflate])
AC_CHECK_LIB([bz2], [BZ2_bzCompress])
AC_CHECK_LIB([popt], [poptGetContext])
AC_CHECK_LIB([protobuf], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdint.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_VOLATILE
AC_HEADER_STDBOOL
AC_SYS_LARGEFILE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

AC_COMPUTE_INT( SIZEOF_INT_P, sizeof( int* ) )
AC_SUBST( SIZEOF_INT_P )

# Checks for library functions.
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([memset munmap pow strerror])

AC_PATH_PROG(GUILE_CONFIG,guile-config)
if test "$GUILE_CONFIG" ; then
	GUILE_FLAGS
fi

AM_CONDITIONAL([PTHREADS],[test "$have_pthread"])
AM_CONDITIONAL([GUILE],[test "$GUILE_CONFIG"])

AC_CONFIG_FILES([
		 Makefile
		 man/Makefile
		 src/Makefile
		 ])
AC_OUTPUT

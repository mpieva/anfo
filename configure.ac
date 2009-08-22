#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([anfo], esyscmd([svnversion -n]), [udo_stenzel@eva.mpg.de])
AC_CONFIG_SRCDIR([util.cc])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PATH_PROG([PROTOC], [protoc])
AC_PATH_PROG([DOXYGEN], [doxygen])
AC_SUBST(PROTOC)
AC_SUBST(DOXYGEN)

# Checks for libraries.
AC_SEARCH_LIBS([exp], [m])
AC_SEARCH_LIBS([log1p], [m])
AC_SEARCH_LIBS([pthread_create], [pthread])
AC_CHECK_LIB([Judy], [Judy1Set])
AC_CHECK_LIB([z], [deflate])
AC_CHECK_LIB([bz2], [BZ2_bzCompress])
AC_CHECK_LIB([popt], [poptGetContext])
AC_CHECK_LIB([protobuf], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdint.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_C_VOLATILE

AC_COMPUTE_INT( SIZEOF_INT_P, sizeof( int* ) )
AC_SUBST( SIZEOF_INT_P )

# Checks for library functions.
AC_REPLACE_FNMATCH
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_CHECK_FUNCS([memset munmap pow strerror])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

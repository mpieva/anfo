#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(anfo, 0.7.x, udo_stenzel@eva.mpg.de)
AC_CONFIG_SRCDIR([util.cc])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_ARG_VAR([PROTOC], [path to google protobuf compiler])
AC_PATH_PROG([PROTOC], [protoc])

# Checks for libraries.
AC_CHECK_LIB([Judy], [Judy1Set])
AC_CHECK_LIB([bz2], [BZ2_bzCompress])
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([popt], [poptGetContext])
# FIXME: Replace `main' with a function in `-lprotobuf':
AC_CHECK_LIB([protobuf], [main])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([z], [deflate])

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

# Checks for library functions.
AC_FUNC_ALLOCA
AC_REPLACE_FNMATCH
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_CHECK_FUNCS([memset munmap pow strerror])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

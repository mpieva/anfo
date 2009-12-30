#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT([anfo],[esyscmd([svnversion -n])],[udo_stenzel@eva.mpg.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/util.cc])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
LT_INIT
AC_PROG_CXX
AC_PROG_INSTALL

AC_ARG_VAR([PROTOC], [location of protobuf compiler])
AC_PATH_PROG([PROTOC], [protoc], AC_MSG_ERROR([protoc not found.  Need to install protobuf-compiler?]))

# Checks for libraries.
AC_LANG([C++])
AC_SEARCH_LIBS([exp], [m])
AC_SEARCH_LIBS([log1p], [m])

AC_CHECK_LIB([Judy], [Judy1Set])
AC_CHECK_LIB([popt], [poptGetContext])
AC_CHECK_LIB([z], [deflate])
AC_CHECK_LIB([bz2], [BZ2_bzCompress])
AC_CHECK_LIB([protobuf], [main])
AC_CHECK_LIB([elk], [Elk_Init])

AC_CHECK_HEADER([Judy.h], [], [AC_MSG_ERROR([Judy.h not found.  Need to install libjudy-dev?])])
AC_CHECK_HEADERS([popt.h], [], [AC_MSG_ERROR([popt.h not found.  Need to install libpopt-dev?])])
AC_CHECK_HEADERS([google/protobuf/message.h], [], [AC_MSG_ERROR([message.h not found.  Need to install libprotobuf-dev?])])
AC_CHECK_HEADERS([zlib.h], [], [AC_MSG_WARN([zlib.h not found.  Want to install zlib1g-dev?])])
AC_CHECK_HEADERS([bzlib.h], [], [AC_MSG_WARN([bzlib.h not found.  Want to install libbz2-dev?])])

PKG_CHECK_MODULES([protobuf], [protobuf >= 2.1.0])

AC_ARG_ENABLE( [elk], 
      AS_HELP_STRING([--disable-elk], [do not build Elk bindings]),
      [], [enable_elk=check] ) 

AS_IF( [test "x$enable_elk" = check],
       [AC_CHECK_HEADERS( [elk/scheme.h], [],
	  	          [AC_MSG_WARN( [elk/scheme.h not found.  Want to install libelk0-dev?] )] )],
       [test "x$enable_elk" != xno],
       [AC_CHECK_HEADERS( [elk/scheme.h], [],
			  [AC_MSG_ERROR( [elk/scheme.h not found, but Elk was requested.] )] )],
       [] )

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

# Checks for library functions.
AC_CHECK_FUNCS([memset munmap pow strerror])

AC_CONFIG_FILES([
		 Makefile
		 man/Makefile
		 src/Makefile
		 ])
AC_OUTPUT

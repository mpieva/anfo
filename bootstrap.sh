#! /bin/sh
autoheader -Wall
aclocal -Wall
libtoolize --automake
automake -a -Wall
autoconf -Wall

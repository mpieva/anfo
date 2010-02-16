#! /bin/sh
autoheader -Wall
aclocal -Wall
libtoolize -c --automake
automake -ac -Wall
autoconf -Wall

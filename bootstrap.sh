#! /bin/sh
echo 'autoheader...' && autoheader -Wall && \
echo 'aclocal...' && aclocal -Wall && \
echo 'libtoolize' && libtoolize -c --automake && \
echo 'automake...' && automake -ac -Wall && \
echo 'autoconf...' && autoconf -Wall

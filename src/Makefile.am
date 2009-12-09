bin_PROGRAMS =  fa2dna  dnaindex  anfo-sge  anfo-tool
libexec_PROGRAMS =  index-test  martin-table  file-info
EXTRA_PROGRAMS =  anfo
lib_LTLIBRARIES = libanfo.la
dist_pkgdata_DATA = config.proto output.proto 
dist_noinst_DATA = anfo.scm.in
pkginclude_HEADERS = align.h anfo_common.h compress_stream.h conffile.h \
		     ducttape.h index.h logdom.h misc_streams.h output_streams.h \
		     sequence.h stream.h trim.h util.h 

BUILT_SOURCES = config.pb.h output.pb.h config.pb.cc output.pb.cc
CLEANFILES = config.pb.h output.pb.h config.pb.cc output.pb.cc anfo.scm

LDADD = libanfo.la

libanfo_la_LDFLAGS = -version-info 0:0:0

libanfo_la_SOURCES = config.pb.h output.pb.h \
		     align.cc anfo_common.cc config.pb.cc ducttape.cc \
		     index.cc misc_streams.cc output.pb.cc output_streams.cc \
		     sequence.cc stream.cc util.cc \
		     judy++.h fastlz.h fastlz.c

anfo_SOURCES = anfo.cc queue.h 
anfo_sge_SOURCES = anfo-sge.cc 
anfo_tool_SOURCES = anfo-tool.cc
dnaindex_SOURCES = dnaindex.cc 
fa2dna_SOURCES = fa2dna.cc 
file_info_SOURCES = file-info.cc
index_test_SOURCES = index-test.cc
martin_table_SOURCES = martin-table.cc 


if PTHREADS
  bin_PROGRAMS += anfo
endif

if ELK
  libanfo_la_SOURCES += anfo-elk.cc
  schemedir = $(dir $(ELK))../share/elk
  scheme_DATA = anfo.scm
  LDADD += -lelk
endif

anfo.scm: anfo.scm.in
	sed -e "s#__LIBDIR__#$(libdir)#g" > anfo.scm < anfo.scm.in

config.pb.cc: config.proto
	$(PROTOC) --cpp_out=. config.proto

output.pb.cc: output.proto
	$(PROTOC) --cpp_out=. output.proto

config.pb.h: config.pb.cc
	@:

output.pb.h: output.pb.cc
	@:


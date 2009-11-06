bin_PROGRAMS =  fa2dna  dnaindex  anfo-sge  anfo-tool  file-info
noinst_PROGRAMS =  index-test  martin-table
EXTRA_PROGRAMS =  anfo  anfo-shell
lib_LTLIBRARIES = libanfo.la
dist_pkgdata_DATA = config.proto output.proto 
pkginclude_HEADERS = align.h anfo_common.h compress_stream.h conffile.h \
		     ducttape.h index.h logdom.h output_streams.h \
		     sequence.h stream.h trim.h util.h 

if GUILE
  bin_PROGRAMS += anfo-shell
endif

if PTHREADS
  bin_PROGRAMS += anfo
endif

BUILT_SOURCES = config.pb.h output.pb.h config.pb.cc output.pb.cc
CLEANFILES = config.pb.h output.pb.h config.pb.cc output.pb.cc

AM_CXXFLAGS = $(GUILE_CFLAGS)
LDADD = libanfo.la

libanfo_la_LDFLAGS = -version-info 0:0:0

libanfo_la_SOURCES = config.pb.h output.pb.h judy++.h \
		     align.cc anfo_common.cc config.pb.cc ducttape.cc \
		     index.cc output.pb.cc output_streams.cc sequence.cc \
		     stream.cc util.cc 

anfo_SOURCES = anfo.cc queue.h 
anfo_sge_SOURCES = anfo-sge.cc 
anfo_tool_SOURCES = anfo-tool.cc
dnaindex_SOURCES = dnaindex.cc 
fa2dna_SOURCES = fa2dna.cc 
file_info_SOURCES = file-info.cc
index_test_SOURCES = index-test.cc
martin_table_SOURCES = martin-table.cc 

anfo_shell_LDADD = $(GUILE_LDFLAGS) $(LDADD)
anfo_shell_SOURCES = anfo-shell.cc


config.pb.cc: config.proto
	$(PROTOC) --cpp_out=. config.proto

output.pb.cc: output.proto
	$(PROTOC) --cpp_out=. output.proto

config.pb.h: config.pb.cc
	@:

output.pb.h: output.pb.cc
	@:


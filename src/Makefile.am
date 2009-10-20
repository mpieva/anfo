bin_PROGRAMS =  fa2dna  dnaindex  anfo-sge  anfo-tool  file-info
noinst_PROGRAMS =  index-test  
EXTRA_PROGRAMS =  anfo  anfo-shell

if GUILE
bin_PROGRAMS += anfo-shell
endif
if PTHREADS
bin_PROGRAMS += anfo
endif

BUILT_SOURCES = config.pb.cc output.pb.cc
CLEANFILES = config.pb.h output.pb.h config.pb.cc output.pb.cc

AM_CXXFLAGS = $(GUILE_CFLAGS)

fa2dna_SOURCES = \
    align.h conffile.h config.pb.h \
    index.h sequence.h util.h \
    \
    align.cc \
    anfo_common.cc config.pb.cc ducttape.cc fa2dna.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

dnaindex_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc \
    anfo_common.cc config.pb.cc dnaindex.cc ducttape.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

anfo_sge_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc anfo-sge.cc \
    anfo_common.cc config.pb.cc ducttape.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

anfo_SOURCES = \
    config.proto output.proto \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc anfo.cc \
    anfo_common.cc config.pb.cc ducttape.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

anfo_shell_LDADD = $(GUILE_LDFLAGS)
anfo_shell_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc anfo-shell.cc \
    anfo_common.cc config.pb.cc ducttape.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

anfo_tool_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc anfo-tool.cc \
    anfo_common.cc config.pb.cc ducttape.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

file_info_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc \
    anfo_common.cc config.pb.cc ducttape.cc file-info.cc \
    index.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

index_test_SOURCES = \
    align.h anfo_common.h compress_stream.h conffile.h config.pb.h \
    ducttape.h index.h judy++.h logdom.h output.pb.h output_streams.h \
    queue.h sequence.h stream.h trim.h util.h \
    \
    align.cc \
    anfo_common.cc config.pb.cc ducttape.cc \
    index.cc index-test.cc output.pb.cc output_streams.cc \
    sequence.cc stream.cc util.cc 

config.pb.cc: config.proto
	$(PROTOC) --cpp_out=. config.proto

output.pb.cc: output.proto
	$(PROTOC) --cpp_out=. output.proto

config.pb.h: config.pb.cc
	@:

output.pb.h: output.pb.cc
	@:


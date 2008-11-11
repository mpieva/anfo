svnversion := $(shell svnversion)
version ?= $(svnversion)

PROTOC ?= protoc
LDLIBS += -lprotobuf -lpopt
# CXXFLAGS += -O3 -Wall -MMD
CXXFLAGS = -ggdb -Wall -MMD -DVERSION='"$(version)"'

TARGETS := fa2dna dnaindex index_test 
# hg18.dna chr21.dna chr21_10.idx hg18_10.idx

all: $(TARGETS)

dnaindex: util.o Index.o metaindex.pb.o
fa2dna: util.o metaindex.pb.o
index_test: util.o Index.o

index_test: util.h Index.h
fa2dna.o: util.h metaindex.pb.h
dnaindex.o: Index.h util.h metaindex.pb.h
util.o: util.h
Index.o: Index.h

hg18.dna:
	./fa2dna -g hg18 -d "Homo Sapiens genome, revision 18" -c index.txt -o $@ -v \
	    	/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa

chr21.dna: 
	./fa2dna -g chr21 -d "Homo Sapiens chromosome 21" -c index.txt -o $@ -v \
		/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

%_10.idx: %.dna 
	./dnaindex $< $@

%: %.o
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@
 
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.pb.cc %.pb.h: %.proto
	$(PROTOC) --cpp_out=$(*D) $<

-include *.d

.PRECIOUS: %.o %.cc %.h
.PHONY: clean all
.SUFFIXES:

clean:
	-rm $(TARGETS)
	-rm *.o *.pb.cc *.pb.h


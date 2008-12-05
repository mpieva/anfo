svnversion := $(shell svnversion)
version ?= $(svnversion)

PROTOC ?= protoc
LDLIBS += -lprotobuf -lpopt -lJudy
# CXXFLAGS += -O3 -Wall -MMD
CXXFLAGS = -ggdb -Wall -MMD -DVERSION='"$(version)"'

TARGETS := fa2dna dnaindex index_test 
DATABASES := hg18.dna chr21.dna chr21_10.idx hg18_10.idx

all: $(TARGETS)
dbs: $(DATABASES)

OBJECTS := util.o index.o metaindex.pb.o conffile.o

fa2dna: $(OBJECTS)
dnaindex: $(OBJECTS)
index_test: $(OBJECTS)

hg18.dna: fa2dna
	./fa2dna -g hg18 -d "Homo Sapiens genome, revision 18" -c index.txt -o $@ -v \
	    	/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa

chr21.dna: fa2dna
	./fa2dna -g chr21 -d "Homo Sapiens chromosome 21" -c index.txt -o $@ -v \
		/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

hg18_10.idx: hg18.dna dnaindex
	./dnaindex -g hg18 -o $@ -c index.txt -v -s 10

chr21_10.idx: chr21.dna dnaindex
	./dnaindex -g chr21 -o $@ -c index.txt -v -s 10

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
	-rm *.o *.d *.pb.cc *.pb.h


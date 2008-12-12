svnversion := $(shell svnversion)
version ?= $(svnversion)

PROTOC ?= protoc
LDLIBS += -lprotobuf -lpopt -lJudy
CXXFLAGS += -O3 -Wall -MMD -DVERSION='"$(version)"'
# CXXFLAGS = -ggdb -Wall -MMD -DVERSION='"$(version)"'

TARGETS := fa2dna dnaindex index_test 
DATABASES := ../data/hg18.dna ../data/chr21.dna ../data/chr21_10.idx ../data/hg18_10.idx

all: $(TARGETS)
dbs: $(DATABASES)

OBJECTS := util.o index.o metaindex.pb.o conffile.o

fa2dna: $(OBJECTS)
dnaindex: $(OBJECTS)
index_test: $(OBJECTS)

../data/hg18.dna: fa2dna
	./fa2dna -g hg18 -d "Homo Sapiens genome, revision 18" -c index.txt -o ../data/$@ -v \
	    	/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa

../data/chr21.dna: fa2dna
	./fa2dna -g chr21 -d "Homo Sapiens chromosome 21" -c index.txt -o ../data/$@ -v \
		/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

../data/hg18_10.idx: ../data/hg18.dna dnaindex
	./dnaindex -g hg18 -o ../data/$@ -c index.txt -v -s 10

../data/chr21_10.idx: ../data/chr21.dna dnaindex
	./dnaindex -g chr21 -o ../data/$@ -c index.txt -v -s 10

%: %.o 
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@
 
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.pb.cc %.pb.h: %.proto
	$(PROTOC) --cpp_out=$(*D) $<

-include *.d

.PRECIOUS: %.o %.cc %.h
.PHONY: clean all doc
.SUFFIXES:

doc:
	doxygen

clean:
	-rm $(TARGETS)
	-rm *.o *.d *.pb.cc *.pb.h


svnversion := $(shell svnversion)
version ?= $(svnversion)

PROTOC ?= protoc
LDLIBS += -lprotobuf -lpopt -lJudy
# LDFLAGS += -p
CXXFLAGS += -Wall -MMD -DVERSION='"$(version)"'
CXXFLAGS += -ggdb 
# CXXFLAGS += -O3

TARGETS := fa2dna dnaindex index_test file-info anfo-standalone
# DATABASES := ../data/chr21.dna ../data/chr21_10.idx ../data/hg18.dna ../data/hg18_10.idx
DATABASES := ../data/chr21.dna ../data/chr21_10.idx 

all: tags $(TARGETS)
dbs: $(DATABASES)

OBJECTS := config.pb.o output.pb.o util.o index.o conffile.o sequence.o

fa2dna: $(OBJECTS)
dnaindex: $(OBJECTS)
file-info: $(OBJECTS)
index_test: $(OBJECTS)
anfo-standalone: $(OBJECTS)

../data/hg18.dna: fa2dna
	./fa2dna -g hg18 -d "Homo Sapiens genome, revision 18" -O ${@D} -v \
	    	/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa

../data/chr21.dna: fa2dna
	./fa2dna -g chr21 -d "Homo Sapiens chromosome 21" -O ${@D} -v \
		/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

../data/hg18_10.idx: ../data/hg18.dna dnaindex
	./dnaindex -O ${@D} -G ${<D} -g ${<F} -v -s 10

../data/chr21_10.idx: ../data/chr21.dna dnaindex
	./dnaindex -O ${@D} -G ${<D} -g ${<F} -v -s 10

%: %.o 
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@
 
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.pb.cc %.pb.h: %.proto
	$(PROTOC) --cpp_out=$(*D) $<

-include *.d

.PRECIOUS: %.o %.pb.cc %.pb.h
.PHONY: clean all doc
.SUFFIXES:

doc:
	doxygen

tags: *.cc *.h
	ctags -R --exclude=*.pb.cc --exclude=*.pb.h 

clean:
	-rm $(TARGETS)
	-rm *.o *.d *.pb.cc *.pb.h


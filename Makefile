svnversion := $(shell svnversion)
version ?= $(svnversion)

PROTOC ?= protoc
LDLIBS += -lprotobuf -lpopt -lJudy
# LDFLAGS += -p
CXXFLAGS += -O3 -Wall -MMD -DVERSION='"$(version)"'
# CXXFLAGS += -O3 -Wall -MMD -DVERSION='"$(version)"'
# CXXFLAGS = -ggdb -Wall -MMD -DVERSION='"$(version)"'

TARGETS := fa2dna dnaindex index_test
# TARGETS := fa2dna dnaindex index_test anfo-standalone
# DATABASES := ../data/chr21.dna ../data/chr21_10.idx ../data/hg18.dna ../data/hg18_10.idx
DATABASES := ../data/chr21.dna ../data/chr21_10.idx ../data/hg18.dna

all: tags $(TARGETS)
dbs: $(DATABASES)

OBJECTS := config.pb.o output.pb.o util.o index.o conffile.o sequence.o

fa2dna: $(OBJECTS)
dnaindex: $(OBJECTS)
index_test: $(OBJECTS)
anfo-standalone: $(OBJECTS)

../data/hg18.dna: fa2dna
	./fa2dna -g hg18 -d "Homo Sapiens genome, revision 18" -c $(@:.dna=.cfg) -o $@ -v \
	    	/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa

../data/chr21.dna: fa2dna
	./fa2dna -g chr21 -d "Homo Sapiens chromosome 21" -o $@ -v \
		/mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

../data/hg18_10.idx: ../data/hg18.dna dnaindex
	./dnaindex -o $@ -g $(<:.dna=.cfg) -c $(@:.idx=.cfg) -v -s 10

../data/chr21_10.idx: ../data/chr21.dna dnaindex
	./dnaindex -o $@ -g $(<:.dna=.cfg) -c $(@:.idx=.cfg) -v -s 10

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


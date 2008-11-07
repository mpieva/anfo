CXXFLAGS += -O3 -Wall
# CXXFLAGS = -ggdb

TARGETS := fa2dna dnaindex index_test hg18.dna chr21.dna chr21_10.idx hg18_10.idx

all: $(TARGETS)

dnaindex: util.o Index.o
fa2dna: util.o
index_test: util.o Index.o

index_test: util.h Index.h
fa2dna.o: util.h
dnaindex.o: Index.h util.h
util.o: util.h
Index.o: Index.h

hg18.dna:
	cat /mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr*.fa | ./fa2dna $@ 

chr21.dna: 
	./fa2dna $@ < /mnt/sequencedb/ucsc/goldenPath/hg18/bigZips/chr21.fa

%_10.idx: %.dna 
	./dnaindex $< $@

%: %.o
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@
 
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.SUFFIXES:
.PHONY: clean all

clean:
	-rm $(TARGETS)
	-rm *.o

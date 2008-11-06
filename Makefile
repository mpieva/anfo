CXXFLAGS += -O3 -Wall
# CXXFLAGS = -ggdb

TARGETS := fa2dna dnaindex index_test

all: $(TARGETS)

dnaindex: util.o Index.o
fa2dna: util.o
index_test: util.o Index.o

index_test: util.h Index.h
fa2dna.o: util.h
dnaindex.o: Index.h util.h
util.o: util.h
Index.o: Index.h

%: %.o
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@
 
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.SUFFIXES:
.PHONY: clean all

clean:
	-rm $(TARGETS)
	-rm *.o

CXXFLAGS += -O3 -Wall

TARGETS := fa2dna dnaindex

all: $(TARGETS)

fa2dna: fa2dna.cc util.cc 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

dnaindex: dnaindex.cc util.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@


.SUFFIXES:
.PHONY: clean all

clean:
	-rm $(TARGETS)
	-rm *.o

all: make_ancestral_genome

Config.hs: ../config.proto
	hprotoc -I .. -v config.proto

make_ancestral_genome: make_ancestral_genome.hs Config.hs
	ghc -Wall -O2 --make $< -XFlexibleInstances -XMultiParamTypeClasses -XDeriveDataTypeable 

clean:
	-rm make_ancestral_genome Config.hs
	-rm *.hi *.o 
	-rm -r Config


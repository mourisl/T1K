CXX = g++
CXXFLAGS= -O3 -g #-pg #-Wall #-O3
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = 

all: fastq-extractor bam-extractor genotyper

genotyper: Genotyper.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

bam-extractor: BamExtractor.o
	if [ ! -f ./samtools-0.1.19/libbam.a ] ; \
	        then \
		                cd samtools-0.1.19 ; make ;\
	fi ;
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS) -lbam

fastq-extractor: FastqExtractor.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)


Genotyper.o: Genotyper.cpp Genotyper.hpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h KmerCode.hpp
BamExtractor.o: BamExtractor.cpp alignments.hpp defs.h SeqSet.hpp
FastqExtractor.o: FastqExtractor.cpp ReadFiles.hpp defs.h SeqSet.hpp BarcodeCorrector.hpp SimpleVector.hpp
#Alignment.o: Alignment.cpp Alignment.h SimpleVector.h defs.h StatsTests.h KmerTree.h ReadSet.h KmerIndex.h poa.h

clean:
	rm -f *.o *.gch genotyper bam-extractor fastq-extractor

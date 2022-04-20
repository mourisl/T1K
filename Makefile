CXX = g++
CXXFLAGS= -O3 -g #-pg #-Wall #-O3
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19 -I./htslib-1.15.1/ -L./htslib-1.15.1/
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = 

# Change to 1 if need to parse CRAM file
htslib=0
BAMFLAGS = 
ifeq ($(htslib),1)
	CXXFLAGS += -DHTSLIB
	BAMFLAGS += -lhts
else
	BAMFLAGS += -lbam
endif

#asan=1
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LDFLAGS+=-fsanitize=address -ldl
endif


all: fastq-extractor bam-extractor genotyper analyzer

genotyper: Genotyper.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

analyzer: Analyzer.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

bam-extractor: BamExtractor.o
	if [ $(htslib) -eq 1 ] ; then \
		if [ ! -f ./htslib-1.15.1/libhts.a ] ; \
			then \
				cd htslib-1.15.1 ; \
				cd htscodecs ; autoreconf -i ; ./configure ; make ; \
				cd .. \
				autoreconf -i ; ./configure ; make ; \
		fi \
	else \
		if [ ! -f ./samtools-0.1.19/libbam.a ] ; \
	        then \
		                cd samtools-0.1.19 ; make ;\
		fi  \
	fi ; 
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS) $(BAMFLAGS)

fastq-extractor: FastqExtractor.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)


Genotyper.o: Genotyper.cpp Genotyper.hpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h KmerCode.hpp KmerCount.hpp
Analyzer.o: Analyzer.cpp Genotyper.hpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h KmerCode.hpp KmerCount.hpp VariantCaller.hpp BarcodeSummary.hpp
BamExtractor.o: BamExtractor.cpp alignments.hpp defs.h SeqSet.hpp
FastqExtractor.o: FastqExtractor.cpp ReadFiles.hpp defs.h SeqSet.hpp BarcodeCorrector.hpp SimpleVector.hpp
#Alignment.o: Alignment.cpp Alignment.h SimpleVector.h defs.h StatsTests.h KmerTree.h ReadSet.h KmerIndex.h poa.h

clean:
	rm -f *.o *.gch genotyper analyzer bam-extractor fastq-extractor

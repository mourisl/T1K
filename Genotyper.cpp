#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>

#include <vector>
#include <pthread.h>

#include "Genotyper.hpp"

char usage[] = "./genotyper [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the reference genome sequence\n"
		"\t[Read file]\n"
		"\t-u STRING: path to single-end read file\n"
		"\t-1 STRING -2 STRING: path to paired-end files\n" 
		"Optional:\n"
		"\t-a STRING: path to the abundance file\n"
		"\t-t INT: number of threads (default: 1)\n"
		"\t-o STRING: output prefix (defult: t1k)\n"
		"\t-n INT: maximal number of alleles per read (default: 2000)\n"
		"\t-s FLOAT: filter alignments with alignment similarity less than specified value (defalut: 0.8)\n"
		"\t--alleleWhitelist STRING: only consider read aligned to the listed allele sereies. (default: not used)"
		"\t--barcode STRING: path to the barcode file\n"
		"\t--frac FLOAT: filter if abundance is less than the frac of dominant allele (default: 0.15)\n"
		"\t--cov FLOAT: filter genes with average coverage less than the specified value (default: 1.0)\n"
		"\t--crossGeneRate FLOAT: the effect from other gene's expression (default: 0.04)\n"
		"\t--relaxIntronAlign: allow one more mismatch in intronic alignment (default: false)\n"
		"\t--alleleDigitUnits INT: the number of units in genotyping result (default: automatic)\n"
		"\t--alleleDelimiter CHR: the delimiter character for digit unit (default: automatic)\n"
    "\t--outputReadAssignment: output the allele assignment for each read to prefix_assign.tsv file (default: not used)\n"
		;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[4] = {'A', 'C', 'G', 'T'} ;

static const char *short_options = "f:a:u:1:2:o:t:n:s:b:" ;
static struct option long_options[] = {
	{ "frac", required_argument, 0, 10000 },
	{ "cov", required_argument, 0, 10001 },
	{ "crossGeneRate", required_argument, 0, 10002 },
	{ "barcode", required_argument, 0, 10003 },
	{ "relaxIntronAlign", no_argument, 0, 10004 },
	{ "alleleDigitUnits", required_argument, 0, 10005 },  
	{ "alleleDelimiter", required_argument, 0, 10006 },  
	{ "alleleWhitelist", required_argument, 0, 10007 },  
	{ "outputReadAssignment", no_argument, 0, 10008 },  
	{(char *)0, 0, 0, 0}
} ;

struct _genotypeRead
{
	char *id ;
	char *seq ;
	char *qual ;

	int strand ;
	bool hasN ;

	int barcode ;
	int umi ;

	int mate ; // 0-first mate, 1-second mate
	int idx ; 
	
	int info ;
	bool fragmentAssigned ;

	bool operator<(const struct _genotypeRead &b) const	
	{
		return strcmp(seq, b.seq) < 0 ;
	}
} ;

struct _assignReadsThreadArg
{
	int tid ;
	int threadCnt ;

	SeqSet *pRefSet ;
	Genotyper *pGenotyper ;	
	std::vector<struct _genotypeRead> *pReads ;
	std::vector< std::vector<struct _overlap> *> *pReadAssignments ;
} ;

struct _readAssignmentToFragmentAssignmentThreadArg
{
	int tid; 
	int threadCnt ;

	int start ;
	int end ;

	bool hasMate ;

	SeqSet *pRefSet ;
	Genotyper *pGenotyper ;
	std::vector<struct _genotypeRead> *pReads1 ;
	std::vector<struct _genotypeRead> *pReads2 ;	
	std::vector< std::vector<struct _overlap> *> *pReadAssignments ;
} ;

char buffer[10241] = "" ;

void PrintLog( const char *fmt, ... )
{
	va_list args ;
	va_start( args, fmt ) ;
	vsprintf( buffer, fmt, args ) ;

	time_t mytime = time(NULL) ;
	struct tm *localT = localtime( &mytime ) ;
	char stime[500] ;
	strftime( stime, sizeof( stime ), "%c", localT ) ;
	fprintf( stderr, "[%s] %s\n", stime, buffer ) ;
}

void *AssignReads_Thread( void *pArg )
{
	struct _assignReadsThreadArg &arg = *( (struct _assignReadsThreadArg *)pArg ) ;
	int start, end ;
	int i, j, k ;
	int totalReadCnt = arg.pReads->size() ;
	start = totalReadCnt / arg.threadCnt * arg.tid ;
	end = totalReadCnt / arg.threadCnt * ( arg.tid + 1 ) ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = totalReadCnt ;
	struct _overlap assign ;
	
	std::vector<struct _genotypeRead> &reads = *(arg.pReads) ;
	std::vector< std::vector<struct _overlap> *> &readAssignments = *(arg.pReadAssignments);
	SeqSet &refSet = *(arg.pRefSet) ;
	std::vector<struct _overlap> *assignments = NULL ;
	int weight = 0 ;
	for ( i = start ; i < end ; )
	{
			for (j = i + 1 ; j < end ; ++j)
				if (strcmp(reads[j].seq, reads[i].seq) != 0)
					break ;
			assignments = new std::vector<struct _overlap> ;
			refSet.AssignRead(reads[i].seq, -1, j - i, *assignments) ;
			//if (arg.tid == 0 && i % 10000 == 0)
			//	printf("%d\n", i * arg.threadCnt) ;
			for (k = i ; k < j ; ++k)
				readAssignments[k] = assignments ;

			i = j ;
	}
	pthread_exit( NULL ) ;
}

void *ReadAssignmentToFragmentAssignment_Thread(void *pArg)
{
	struct _readAssignmentToFragmentAssignmentThreadArg &arg = *((struct _readAssignmentToFragmentAssignmentThreadArg *)pArg) ;
	int i, j ;
	int totalReadCnt = arg.pReads1->size() ;
	//if ( arg.tid == arg.threadCnt - 1 )
	//	end = totalReadCnt ;
	struct _overlap assign ;
	
	std::vector<struct _genotypeRead> &reads1 = *(arg.pReads1) ;
	std::vector<struct _genotypeRead> &reads2 = *(arg.pReads2) ;
	std::vector< std::vector<struct _overlap> *> &readAssignments = *(arg.pReadAssignments);
	SeqSet &refSet = *(arg.pRefSet) ;
	std::vector<struct _fragmentOverlap > fragmentAssignment ;
	
	//for (i = start ; i < end ; ++i)
	for (i = arg.tid + arg.start ; i < arg.end ; i += arg.threadCnt)
	{
		bool hasN = reads1[i].hasN ;
		if (arg.hasMate && reads2[i].hasN) 
			hasN = true ;

		if (!arg.hasMate)	
			refSet.ReadAssignmentToFragmentAssignment( readAssignments[ reads1[i].info ], NULL, reads1[i].barcode, hasN, fragmentAssignment) ;
		else
			refSet.ReadAssignmentToFragmentAssignment( readAssignments[reads1[i].info], readAssignments[reads2[i].info],
					reads1[i].barcode, hasN, fragmentAssignment) ;
		arg.pGenotyper->SetReadAssignments(i, fragmentAssignment ) ;	
		if (fragmentAssignment.size() > 0)
			reads1[i].fragmentAssigned = true ;
	}
  pthread_exit( NULL ) ;
}

int main(int argc, char *argv[])
{
	int i, j, k ;
	int c, option_index = 0 ;
	
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}
	
	char outputPrefix[1024] = "t1k" ;
	
	Genotyper genotyper(11) ;
	ReadFiles reads ;
	ReadFiles mateReads ;
	ReadFiles barcodeFile ;
	char *refFile = NULL ;
	bool hasMate = false ;
	bool hasBarcode = false ;
	bool hasUmi = false ;
	std::vector<struct _genotypeRead> reads1 ;
	std::vector<struct _genotypeRead> reads2 ;
	int threadCnt = 1 ;
	int maxAssignCnt = 2000 ;
	FILE *fpAbundance = NULL ;
	FILE *fpAlleleWhitelist = NULL ;
	FILE *fpOutput ;
	double filterFrac = 0.15 ;
	double filterCov = 1.0 ;
	double crossGeneRate = 0.04 ;
	double filterAlignmentSimilarity = 0.8 ;
	bool keepMissingBarcode = false ;
	bool relaxIntronAlign = false ;
	int alleleDigitUnits = -1 ;
	char alleleDelimiter = '\0' ;
  bool outputReadAssignment = false ;
	std::map<std::string, int> barcodeStrToInt ;
	std::vector<std::string> barcodeIntToStr ;

	while (1)	
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;

		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			//seqSet.InputRefFa( optarg ) ;
			//genotyper.InitRefSet( optarg ) ;
			refFile = strdup(optarg) ;
		}
		else if ( c == 'a' )
		{
			fpAbundance = fopen( optarg, "r" ) ; 
		}
		else if ( c == 'u' )
		{
			reads.AddReadFile( optarg, false ) ;
		}
		else if ( c == '1' )
		{
			reads.AddReadFile( optarg, true ) ;
		}
		else if ( c == '2' )
		{
			mateReads.AddReadFile( optarg, true ) ;
			hasMate = true ;
		}
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;
		}
		else if ( c == 't' )
		{
			threadCnt = atoi( optarg ) ;
		}
		else if ( c == 'n' )
		{
			maxAssignCnt = atoi(optarg) ;
		}
		else if ( c == 's' )
		{
			filterAlignmentSimilarity = atof(optarg) ;
		}
		else if ( c == 10000 ) // --frac
		{
			filterFrac = atof(optarg) ;
		}
		else if ( c == 10001 ) // --cov
		{
			filterCov = atof(optarg) ;
		}
		else if ( c == 10002 ) // --crossGeneRate
		{
			crossGeneRate = atof(optarg) ;
		}
		else if ( c == 10003 )
		{
			barcodeFile.AddReadFile(optarg, false) ;
			hasBarcode = true ;
		}
		else if ( c == 10004 )
		{
			relaxIntronAlign = true ;
		}
		else if ( c == 10005 )
		{
			alleleDigitUnits = atoi(optarg) ;
		}
		else if ( c == 10006 )
		{
			alleleDelimiter = optarg[0] ;
		}
		else if ( c == 10007 )
		{
			fpAlleleWhitelist = fopen(optarg, "r") ;
		}
    else if ( c == 10008 ) 
    {
      outputReadAssignment = true ;
    }
		else
		{
			fprintf( stderr, "%s", usage ) ;
			return EXIT_FAILURE ;
		}
	}
	if (refFile == NULL)
	{
		fprintf( stderr, "Need to use -f to specify the reference sequences.\n" );
		return EXIT_FAILURE;
	}

	genotyper.SetFilterFrac(filterFrac) ;
	genotyper.SetFilterCov(filterCov) ;
	genotyper.SetCrossGeneRate(crossGeneRate) ;
	genotyper.SetAlleleNameStructure(alleleDigitUnits, alleleDelimiter) ;
	genotyper.InitRefSet(refFile) ;
	if (fpAlleleWhitelist != NULL)
	{
		genotyper.SetAlleleWhitelist(fpAlleleWhitelist) ;
		fclose(fpAlleleWhitelist) ;
	}
	
	SeqSet &refSet = genotyper.refSet ;
	if ( refSet.Size() == 0 )
	{
		fprintf( stderr, "Need to use -f to specify the reference sequences.\n" );
		return EXIT_FAILURE;
	}
	refSet.SetRefSeqSimilarity(filterAlignmentSimilarity) ;
	refSet.SetRelaxIntronAlign(relaxIntronAlign) ;
	if (threadCnt > 1)
		refSet.InitPthread() ;
	int alleleCnt = refSet.Size() ;

	//else
	//{
	//	fprintf( stderr, "Need to use -a to specify the abundance estimation.\n" );
	//	return EXIT_FAILURE;
	//}

	// Read in the sequencing data	
	i = 0;
	int maxReadLength = 0 ;
	while (reads.Next())
	{
		struct _genotypeRead nr ;
		struct _genotypeRead mateR ;
		int barcode = -1 ;
		int umi = -1 ;
		
		if ( hasBarcode )
		{
			barcodeFile.Next() ;
			
			if ( !strcmp( barcodeFile.seq, "missing_barcode" ) && !keepMissingBarcode )
			{
				if ( hasMate )
					mateReads.Next() ;
				continue ;
			}

			std::string s( barcodeFile.seq ) ;
			if ( barcodeStrToInt.find( s ) != barcodeStrToInt.end() )
				barcode = barcodeStrToInt[s] ;
			else
			{
				barcode = barcodeIntToStr.size() ;
				barcodeStrToInt[s] = barcode ;
				barcodeIntToStr.push_back( s ) ;
			}
		}

		nr.seq = strdup(reads.seq) ;
		nr.id = strdup(reads.id) ;
		if (reads.qual != NULL)
			nr.qual = strdup(reads.qual);
		else
			nr.qual = NULL;
		nr.barcode = barcode ;
		nr.umi = umi ;
		nr.idx = i ;
		nr.mate = 0 ;
		nr.fragmentAssigned = false ;
		nr.hasN = false ;
		for (j = 0 ; nr.seq[j] ; ++j)
			if (nr.seq[j] == 'N')
			{
				nr.hasN = true ;
				break ;
			}
		if (strlen(nr.seq) > maxReadLength)
			maxReadLength = strlen(nr.seq) ;
		reads1.push_back(nr);		

		if (hasMate)
		{
			mateReads.Next() ;
			mateR.seq = strdup(mateReads.seq);
			mateR.id = strdup(mateReads.id);
			mateR.barcode = barcode;
			mateR.umi = umi;
			if (mateReads.qual != NULL)
				mateR.qual = strdup( mateReads.qual );
			else
				mateR.qual = NULL;
			mateR.idx = i ;	
			mateR.mate = 1 ;
			mateR.fragmentAssigned = false ;
			mateR.hasN = false ;
			for (j = 0 ; mateR.seq[j] ; ++j)
				if (mateR.seq[j] == 'N')
				{
					mateR.hasN = true ;
					break ;
				}
			if (strlen(mateR.seq) > maxReadLength)
				maxReadLength = strlen(mateR.seq) ;
			reads2.push_back(mateR);			
		}
		++i;
	}
	genotyper.SetReadLength(maxReadLength) ;

	int readCnt = reads1.size() ;
	//std::vector<struct _fragmentOverlap> alleleAssignments ;
	genotyper.InitReadAssignments(readCnt, maxAssignCnt) ;	
	PrintLog("Found %d read fragments. Start read assignment.", readCnt) ;
	
	// Put all the ends into one big array to reuse alignment information.
	std::vector<struct _genotypeRead> allReads ; 
	allReads = reads1 ;
	allReads.insert(allReads.end(), reads2.begin(), reads2.end()) ;
	std::sort(allReads.begin(), allReads.end()) ;
	std::vector< std::vector<struct _overlap> *> readAssignments ;
	
	int allReadCnt = allReads.size() ; 
	int alignedFragmentCnt = 0 ;
	//std::vector< std::vector<struct _fragmentOverlap> > fragmentAssignments ;
	readAssignments.resize(allReadCnt) ;
	//fragmentAssignments.resize(readCnt) ;

	if (threadCnt <= 1)
	{
		std::vector<struct _overlap> *assignments = NULL ;
		for (i = 0 ; i < allReadCnt ; )
		{
			for (j = i + 1 ; j < allReadCnt ; ++j)
				if (strcmp(allReads[j].seq, allReads[i].seq) != 0)
					break ;
			assignments = new std::vector<struct _overlap> ;
			refSet.AssignRead(allReads[i].seq, -1, j - i, *assignments) ;
			//if (arg.tid == 0 && i % 10000 == 0)
			//	printf("%d\n", i * arg.threadCnt) ;
			for (k = i ; k < j ; ++k)
				readAssignments[k] = assignments ;

			i = j ;
		}
	}
	else
	{
		pthread_t *threads = new pthread_t[ threadCnt ] ;
		struct _assignReadsThreadArg *args = new struct _assignReadsThreadArg[threadCnt] ;
		pthread_attr_t attr ;

		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;

		for ( i = 0 ; i < threadCnt ; ++i )
		{
			args[i].tid = i ;
			args[i].threadCnt = threadCnt ;
			args[i].pRefSet = &refSet ;
			args[i].pGenotyper = &genotyper ;
			args[i].pReads = &allReads ;
			args[i].pReadAssignments = &readAssignments ;
			pthread_create( &threads[i], &attr, AssignReads_Thread, (void *)( args + i ) ) ;
		}

		for ( i = 0 ; i < threadCnt ; ++i )
			pthread_join( threads[i], NULL ) ;

		pthread_attr_destroy(&attr);
		delete[] threads ;
		delete[] args ;
	}
	PrintLog("Finish read end assignments.") ;
	
	// Matching up read ends to form fragment assignment
	for (i = 0 ; i < allReadCnt ; ++i)
	{
		if (allReads[i].mate == 0)		
		{
			reads1[ allReads[i].idx ].info = i ; // map to the read assignment part
		}
		else if (allReads[i].mate == 1)
		{
			reads2[ allReads[i].idx ].info = i ;
		}
	}
	
	const int coalesceSize = 500000 ;
  fpOutput = NULL ;
  if (outputReadAssignment)
  {
    sprintf(buffer, "%s_assign.tsv", outputPrefix) ;
    fpOutput = fopen(buffer, "w") ;
  }

	if (threadCnt == 1)
	{
		int coalesceStart = 0 ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			std::vector<struct _fragmentOverlap> fragmentAssignment ;
			bool hasN = reads1[i].hasN ;
			if (hasMate && reads2[i].hasN) 
				hasN = true ;
			
			if (!hasMate)	
				refSet.ReadAssignmentToFragmentAssignment( readAssignments[ reads1[i].info ], NULL, reads1[i].barcode, hasN, fragmentAssignment) ;
			else
				refSet.ReadAssignmentToFragmentAssignment( readAssignments[reads1[i].info], readAssignments[reads2[i].info],
						reads1[i].barcode, hasN, fragmentAssignment) ;
#ifdef DEBUG
			std::vector<struct _fragmentOverlap> &assignments = fragmentAssignment ;
			for (int j = 0 ; j < assignments.size() ; ++j)
				printf("%s\t%d\t%s\t%d\t%lf. %d %d. %d %d. %d %d\n", reads1[i].id, assignments[j].seqIdx, refSet.GetSeqName(assignments[j].seqIdx),
						assignments[j].matchCnt, assignments[j].similarity, assignments[j].overlap1.matchCnt, assignments[j].overlap2.matchCnt, 
						assignments[j].overlap1.seqStart, assignments[j].overlap1.seqEnd,
						assignments[j].overlap2.seqStart, assignments[j].overlap2.seqEnd);
#endif
			genotyper.SetReadAssignments(i, fragmentAssignment ) ;	
      if (outputReadAssignment)
      {
        std::vector<struct _readAssignment> assignments = genotyper.GetReadAssignments(i) ;
        for (int j = 0 ; j < assignments.size() ; ++j)
          fprintf(fpOutput, "%s\t%s\t%d\t%d\n", 
              reads1[i].id, refSet.GetSeqName(assignments[j].alleleIdx),
              assignments[j].start, assignments[j].end) ;
      }

			if (fragmentAssignment.size() > 0)
				reads1[i].fragmentAssigned = true ;

			if (i > 0 && i % coalesceSize == 0)
			{
				alignedFragmentCnt += genotyper.CoalesceReadAssignments(coalesceStart, i) ;
				coalesceStart = i + 1 ;
			}
		}
		alignedFragmentCnt += genotyper.CoalesceReadAssignments(coalesceStart, readCnt - 1) ;
	}
	else
	{
		pthread_t *threads = new pthread_t[ threadCnt ] ;
		struct _readAssignmentToFragmentAssignmentThreadArg *args = new struct _readAssignmentToFragmentAssignmentThreadArg[threadCnt] ;
		pthread_attr_t attr ;

		pthread_attr_init( &attr ) ;
		pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
		int start, end ;
		for (start = 0 ; start < readCnt ; start += coalesceSize)
		{
			end = start + coalesceSize - 1 <= (readCnt - 1) ? start + coalesceSize - 1 : (readCnt - 1) ;
			for ( i = 0 ; i < threadCnt ; ++i )
			{
				args[i].tid = i ;
				args[i].threadCnt = threadCnt ;
				args[i].pRefSet = &refSet ;
				args[i].pGenotyper = &genotyper ;
				args[i].pReads1 = &reads1 ;
				args[i].pReads2 = &reads2 ;
				args[i].start = start ;
				args[i].end = end + 1 ;  // the thread part is [begin, end)
				args[i].pReadAssignments = &readAssignments ;
				args[i].hasMate = hasMate ;
				pthread_create( &threads[i], &attr, ReadAssignmentToFragmentAssignment_Thread, (void *)( args + i ) ) ;
			}

			for ( i = 0 ; i < threadCnt ; ++i )
				pthread_join( threads[i], NULL ) ;
			
      if (outputReadAssignment)
      {
        for (i = start ; i < readCnt && i < start + coalesceSize ; ++i)
        {
          std::vector<struct _readAssignment> assignments = genotyper.GetReadAssignments(i) ;
          for (int j = 0 ; j < assignments.size() ; ++j)
            fprintf(fpOutput, "%s\t%s\t%d\t%d\n", 
                reads1[i].id, refSet.GetSeqName(assignments[j].alleleIdx),
                assignments[j].start, assignments[j].end) ;
        }
      }
			alignedFragmentCnt += genotyper.CoalesceReadAssignments(start, end) ;
		}
		pthread_attr_destroy(&attr);
		delete[] threads ;
		delete[] args ;
	}
	// Release the memory for read end assignment
	for (i = 0 ; i < allReadCnt ;)
	{
		for (j = i + 1 ; j < allReadCnt ; ++j)
			if (readAssignments[j] != readAssignments[i])
				break ;
		delete readAssignments[i] ;
		i = j ;
	}
  if (outputReadAssignment)
    fclose(fpOutput) ;

	genotyper.FinalizeReadAssignments() ;
	PrintLog( "Finish read fragment assignments. %d read fragments can be assigned (average %.2lf alleles/read).", 
			alignedFragmentCnt, genotyper.GetAverageReadAssignmentCnt()) ;

	// Get some global abundance information, 
	// need for allele selection.
	if (fpAbundance)
		genotyper.InitAlleleAbundance(fpAbundance) ;			
	else
	{
		int emIterCnt = genotyper.QuantifyAlleleEquivalentClass() ;
		PrintLog( "Finish allele quantification in %d EM iterations.", emIterCnt) ;
	}
	genotyper.RemoveLowLikelihoodAlleleInEquivalentClass() ;

	// Main function to do genotyping
	genotyper.SelectAllelesForGenes() ;
	
	// Output the genotype information
	int geneCnt = genotyper.GetGeneCnt() ;
	char *bufferAllele[3] ;
	bufferAllele[0] = new char[20 * refSet.Size() + 40] ;
	bufferAllele[1] = new char[20 * refSet.Size()+ 40] ;
	bufferAllele[2] = new char[20 * refSet.Size()+ 40] ;
	sprintf(buffer, "%s_genotype.tsv", outputPrefix) ;
	fpOutput = fopen(buffer, "w") ;	
	for (i = 0 ; i < geneCnt ; ++i)
	{
		int calledAlleleCnt = genotyper.GetAlleleDescription(i, bufferAllele[0], bufferAllele[1],
        bufferAllele[2]) ;
		fprintf(fpOutput, "%s\t%d", genotyper.GetGeneName(i), calledAlleleCnt) ;
		for (j = 0 ; j < 3 ; ++j)
		{
			fprintf(fpOutput, "\t%s", bufferAllele[j]);
		}
		fprintf(fpOutput, "\n") ;	
	}
	delete[] bufferAllele[0] ;
	delete[] bufferAllele[1] ;
  delete[] bufferAllele[2] ;
	fclose(fpOutput) ;
	
	// Output the representative alleles
	sprintf(buffer, "%s_allele.tsv", outputPrefix) ;
	genotyper.OutputRepresentativeAlleles(buffer) ;

	// Output the read that got assigned
	if (hasMate)
		sprintf(buffer, "%s_aligned_1.fa", outputPrefix) ;
	else
		sprintf(buffer, "%s_aligned.fa", outputPrefix) ;
	fpOutput = fopen(buffer, "w") ;
	for (i = 0 ; i < readCnt ; ++i)	
	{
		if (reads1[i].fragmentAssigned)
		{
			fprintf(fpOutput, ">%s\n%s\n", reads1[i].id, reads1[i].seq) ;
		}
	}
	fclose(fpOutput) ;

	if (hasMate)
	{
		sprintf(buffer, "%s_aligned_2.fa", outputPrefix) ;
		fpOutput = fopen(buffer, "w") ;
		for (i = 0 ; i < readCnt ; ++i)	
		{
			if (reads1[i].fragmentAssigned)
			{
				fprintf(fpOutput, ">%s\n%s\n", reads2[i].id, reads2[i].seq) ;
			}
		}
		fclose(fpOutput) ;
	}

  if (hasBarcode)
	{
		sprintf(buffer, "%s_aligned_bc.fa", outputPrefix) ;
		fpOutput = fopen(buffer, "w") ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			if (reads1[i].fragmentAssigned)
				fprintf(fpOutput, ">%s\n%s\n", reads1[i].id, barcodeIntToStr[reads1[i].barcode].c_str());
		}
	}
	
	free(refFile) ;
	for ( i = 0 ; i < readCnt ; ++i )
	{
		free( reads1[i].id ) ;
		free( reads1[i].seq ) ;
		if ( reads1[i].qual != NULL )
			free( reads1[i].qual ) ;
	
		if (hasMate)
		{
			free( reads2[i].id ) ;
			free( reads2[i].seq ) ;
			if ( reads2[i].qual != NULL )
				free( reads2[i].qual ) ;
		}
	}
	PrintLog("Genotyping finishes.") ;
	return 0;
}


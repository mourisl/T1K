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
		;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[4] = {'A', 'C', 'G', 'T'} ;

static const char *short_options = "f:a:u:1:2:o:t:" ;
static struct option long_options[] = {
	{(char *)0, 0, 0, 0}
} ;

struct _genotypeRead
{
	char *id ;
	char *read ;
	char *qual ;

	int strand ;

	int barcode ;
	int umi ;
} ;

struct _assignReadsThreadArg
{
	int tid ;
	int threadCnt ;

	SeqSet *pRefSet ;
	Genotyper *pGenotyper ;	
	std::vector<struct _genotypeRead> *pReads1 ;
	std::vector<struct _genotypeRead> *pReads2 ;	
	std::vector< std::vector<struct _fragmentOverlap> > *pReadAssignments ;
	
	bool hasMate ;
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
	int i ;
	int totalReadCnt = arg.pReads1->size() ;
	start = totalReadCnt / arg.threadCnt * arg.tid ;
	end = totalReadCnt / arg.threadCnt * ( arg.tid + 1 ) ;
	if ( arg.tid == arg.threadCnt - 1 )
		end = totalReadCnt ;
	struct _overlap assign ;
	
	std::vector<struct _genotypeRead> &reads1 = *(arg.pReads1) ;
	std::vector<struct _genotypeRead> &reads2 = *(arg.pReads2) ;
	std::vector< std::vector<struct _fragmentOverlap> > &readAssignments = *(arg.pReadAssignments);
	SeqSet &refSet = *(arg.pRefSet) ;
	for ( i = start ; i < end ; ++i )
	{
			std::vector<struct _fragmentOverlap> assignments ;
			if (arg.hasMate)
				refSet.AssignRead(reads1[i].read, reads2[i].read, reads1[i].barcode, assignments) ;
			else
				refSet.AssignRead(reads1[i].read, NULL, reads1[i].barcode, assignments) ;
			arg.pGenotyper->SetReadAssignments(i, assignments) ;
			readAssignments[i] = assignments ;
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
	
	char outputPrefix[1024] = "kir" ;
	
	Genotyper genotyper(17) ;
	ReadFiles reads ;
	ReadFiles mateReads ;
	bool hasMate = false ;
	std::vector<struct _genotypeRead> reads1 ;
	std::vector<struct _genotypeRead> reads2 ;
	int threadCnt = 1 ;
	FILE *fpAbundance = NULL ;

	while (1)	
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;

		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			//seqSet.InputRefFa( optarg ) ;
			genotyper.InitRefSet( optarg ) ;
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
		else
		{
			fprintf( stderr, "%s", usage ) ;
			return EXIT_FAILURE ;
		}
	}
	SeqSet &refSet = genotyper.refSet ;

	if ( refSet.Size() == 0 )
	{
		fprintf( stderr, "Need to use -f to specify the receptor genome sequence.\n" );
		return EXIT_FAILURE;
	}

	int alleleCnt = refSet.Size() ;

	//else
	//{
	//	fprintf( stderr, "Need to use -a to specify the abundance estimation.\n" );
	//	return EXIT_FAILURE;
	//}

	// Read in the sequencing data	
	i = 0;
	while (reads.Next())
	{
		struct _genotypeRead nr ;
		struct _genotypeRead mateR ;
		int barcode = -1 ;
		int umi = -1 ;
		
		nr.read = strdup(reads.seq) ;
		nr.id = strdup(reads.id) ;
		if (reads.qual != NULL)
			nr.qual = strdup(reads.qual);
		else
			nr.qual = NULL;
		nr.barcode = barcode ;
		nr.umi = umi ;

		reads1.push_back(nr);		

		if (hasMate)
		{
			mateReads.Next() ;
			mateR.read = strdup(mateReads.seq);
			mateR.id = strdup(mateReads.id);
			mateR.barcode = barcode;
			mateR.umi = umi;
			if (mateReads.qual != NULL)
				mateR.qual = strdup( mateReads.qual );
			else
				mateR.qual = NULL;
			
			reads2.push_back(mateR);			
		}
		++i;
	}

	int readCnt = reads1.size() ;
	//std::vector<struct _fragmentOverlap> alleleAssignments ;
	genotyper.InitReadAssignments(readCnt) ;	
	PrintLog("Found %d read fragments. Start read assignment.", readCnt) ;
	
	int alignedReadCnt = 0 ;
	std::vector< std::vector<struct _fragmentOverlap> > readAssignments ;
	readAssignments.resize(readCnt) ;

	if (threadCnt <= 1)
	{
		for (i = 0 ; i < readCnt ; ++i)
		{
			std::vector<struct _fragmentOverlap> assignments ;
			if (hasMate)
				refSet.AssignRead(reads1[i].read, reads2[i].read, reads1[i].barcode, assignments) ;
			else
				refSet.AssignRead(reads1[i].read, NULL, reads1[i].barcode, assignments) ;
			
			readAssignments[i] = assignments ;

			//	++alignedReadCnt ;
			genotyper.SetReadAssignments(i, assignments) ;
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
			args[i].pReads1 = &reads1 ;
			args[i].pReads2 = &reads2 ;
			args[i].pReadAssignments = &readAssignments ;
			args[i].hasMate = hasMate ;
			pthread_create( &threads[i], &attr, AssignReads_Thread, (void *)( args + i ) ) ;
		}

		for ( i = 0 ; i < threadCnt ; ++i )
			pthread_join( threads[i], NULL ) ;


		delete[] threads ;
		delete[] args ;
	}
	alignedReadCnt = genotyper.FinalizeReadAssignments() ;
	PrintLog("Finish read assignments. %d read fragments can be assigned.", alignedReadCnt) ;

#ifdef DEBUG
	for (i = 0 ; i < readCnt ; ++i)
	{
		std::vector<struct _fragmentOverlap> &assignments = readAssignments[i] ;
		for (int j = 0 ; j < assignments.size() ; ++j)
			printf("%s\t%d\t%s\t%d\t%lf. %d %d\n", reads1[i].id, assignments[j].seqIdx, refSet.GetSeqName(assignments[j].seqIdx),
					assignments[j].matchCnt, assignments[j].similarity, assignments[j].overlap1.matchCnt, assignments[j].overlap2.matchCnt);
	}
#endif

	// Get some global abundance information, 
	// need for allele selection.
	if (fpAbundance)
		genotyper.InitAlleleAbundance(fpAbundance) ;			
	else
		genotyper.QuantifyAlleleEquivalentClass() ;
	genotyper.RemoveLowLikelihoodAlleleInEquivalentClass() ;

	// Main function to do genotyping
	genotyper.SelectAllelesForGenes() ;
		
	// Find the representative alleles
	// Output the genotype information
	int geneCnt = genotyper.GetGeneCnt() ;
	char *bufferAllele[2] ;
	bufferAllele[0] = new char[20 * refSet.Size() + 40] ;
	bufferAllele[1] = new char[20 * refSet.Size()+ 40] ;
	for (i = 0 ; i < geneCnt ; ++i)
	{
		int calledAlleleCnt = genotyper.GetAlleleDescription(i, bufferAllele[0], bufferAllele[1]) ;
		printf("%s\t%d", genotyper.GetGeneName(i), calledAlleleCnt) ;
		for (j = 0 ; j < calledAlleleCnt ; ++j)
		{
			printf("\t%s", bufferAllele[j]);
		}
		printf("\n") ;	
	}
	delete[] bufferAllele[0] ;
	delete[] bufferAllele[1] ;


	for ( i = 0 ; i < readCnt ; ++i )
	{
		free( reads1[i].id ) ;
		free( reads1[i].read ) ;
		if ( reads1[i].qual != NULL )
			free( reads1[i].qual ) ;
	
		if (hasMate)
		{
			free( reads2[i].id ) ;
			free( reads2[i].read ) ;
			if ( reads2[i].qual != NULL )
				free( reads2[i].qual ) ;
		}
	}
	PrintLog("Genotyping finishes.") ;
	return 0;
}


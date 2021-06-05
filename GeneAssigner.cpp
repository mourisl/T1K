#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>

#include <vector>

char usage[] = "./gene-assigner [OPTIONS]:\n"
		"Required:\n"
		"\t-f STRING: fasta file containing the reference genome sequence\n"
		"\t[Read file]\n"
		"\t-u STRING: path to single-end read file\n"
		"\t-1 STRING -2 STRING: path to paired-end files\n" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[4] = {'A', 'C', 'G', 'T'} ;

static const char *short_options = "f:u:1:2:o:t:" ;

int main(int argc, char *argv[])
{
	int i, j, k ;
	int c, option_index = 0 ;
	
	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
	}
	
	SeqSet refSet(9) ;
	char outputPrefix[1024] = "kir" ;

	ReadFiles reads ;
	ReadFiles mateReads ;
	bool hasMate = false ;
	std::vector<struct _Read> read1 ;
	std::vector<struct _Read> read2 ;

	while (1)	
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;

		if ( c == -1 )
			break ;

		if ( c == 'f' )
		{
			//seqSet.InputRefFa( optarg ) ;
			refSet.InputRefFa( optarg ) ;
		}
		else if ( c == 'u' )
		{
			reads.AddReadFile( optarg, false ) ;
		}
		else if ( c == '1' )
		{
			reads.AddReadFile( optarg, true ) ;
			hasMate = true ;
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

	if ( refSet.Size() == 0 )
	{
		fprintf( stderr, "Need to use -f to specify the receptor genome sequence.\n" );
		return EXIT_FAILURE;
	}

	i = 0;
	while (reads.Next())
	{
		struct _Read nr ;
		struct _Read mateR ;
		nr.read = strdup(reads.seq) ;
		nr.id = strdup(reads.id) ;
		if (reads.qual != NULL)
			nr.qual = strdup(reads.qual);
		else
			nr.qual = NULL;
		reads1.push_back(nr);
		
		if (hasMate)
		{
			mateR.read = strdup(mateReads.seq);
			mateR.id = strdup(mateReads.id);
			mateR.barcode = barcode;
			mateR.umi = umi;
			if (mateReads.qual != NULL)
				mateR.qual = strdup( mateReads.qual );
			else
				mateR.qual = NULL;
			
			read2.push_back(mateR);			
		}
		++i;
	}

	

	return 0;


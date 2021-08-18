// The call handles kmer count

#ifndef _LSONG_KMERCOUNT_HEADER
#define _LSONG_KMERCOUNT_HEADER

#include <map>
#include <algorithm>

#include "KmerCode.hpp"

#define KCOUNT_HASH_MAX 103

class KmerCount
{
private:
	std::map<uint64_t, int> *count ;
	int kmerLength ;
	KmerCode kmerCode ;
	int maxReadLen ;

	int *c ;

	int GetHash( uint64_t k )
	{
		return k % KCOUNT_HASH_MAX ;
	}
public:
	KmerCount( int k ): kmerCode( k ) 
	{ 
		kmerLength = k ; 
		maxReadLen = -1 ;
		c = NULL ;
		count = new std::map<uint64_t, int>[KCOUNT_HASH_MAX] ;
	}
	
	KmerCount(): kmerCode(31)
	{ 
		kmerLength = 31 ; 
		maxReadLen = -1 ;
		c = NULL ;
		count = new std::map<uint64_t, int>[KCOUNT_HASH_MAX] ;
	}

	~KmerCount() 
	{
		if ( c != NULL )
			delete[] c ;
		if ( count != NULL )
			delete[] count ;
	}

	int AddCount( char *read )
	{
		int i ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return 0 ;

		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;

		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( kmerCode.IsValid() )
			{
				uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
				++count[ GetHash(kcode) ][ kcode ] ;
				/*if ( count[ GetHash( kcode ) ][ kcode ] >= 500 )
				{
					printf( "%s\n", read + i - kmerLength + 1 ) ;
				}*/
			}
		}
		
		if ( len > maxReadLen )
			maxReadLen = len ;
		return 1 ;
	}

	void AddCountFromFile( char *file )
	{
		FILE *fp = fopen( file, "r" ) ;
		char buffer[100] ;
		int i ;

		while ( fscanf( fp, "%s", buffer ) != EOF )
		{
			int c = atoi( &buffer[ 1 ] ) ;
			fscanf( fp, "%s", buffer ) ;
			if ( c <= 1 )
				continue ;

			kmerCode.Restart() ;
			for ( i = 0 ; buffer[i] ; ++i )
				kmerCode.Append( buffer[i] ) ;
			uint64_t kcode = kmerCode.GetCode() ;
			count[ GetHash( kcode ) ][ kcode ] = c ;
		}
		fclose( fp ) ;
	}

	void Output( char *file )
	{
		int i, j ;
		FILE *fp = fopen( file, "r" ) ;
		char *buffer = new char[kmerLength + 1] ;
		for ( i = 0 ; i < KCOUNT_HASH_MAX ; ++i )	
		{
			for ( std::map<uint64_t, int>::iterator it = count[i].begin() ; it != count[i].end() ; ++it )
			{
				if ( it->second <= 1 )
					continue ;
				
				for ( j = 0 ; j < kmerLength ; ++j )
				{
					buffer[j] = nucToNum[ ( it->first >> ( 2 * j ) ) & 3 ] ;  
				}
				buffer[j] = '\0' ;
				printf( ">%d\n%s\n", it->second, buffer ) ;
			}
		}
		delete[] buffer ;
	}

	void SetBuffer( int sz )
	{
		maxReadLen = sz ;
		if ( c == NULL )
			 c = new int[ sz ] ;
	}
	
	int GetCount( char *kmer )
	{
		int i ;
		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength ; ++i )
			kmerCode.Append( kmer[i] ) ;
		if ( kmerCode.IsValid() )
		{
			uint64_t kcode = kmerCode.GetCanonicalKmerCode() ;
			int key = GetHash( kcode ) ;
			if ( count[ key ].find( kcode ) == count[key].end() )
				return 0 ;
			else
				return count[ GetHash( kcode ) ][ kcode ] ;
		}
		else
			return 0 ;
	}

	void Release()
	{
		if ( c != NULL )
			delete[] c ;
		if ( count != NULL )
			delete[] count ;

		c = NULL ;
		count = NULL ;
	}

	// Jaccard index
	double GetCountSimilarity(const KmerCount &b)
	{
		int i ;
		int countA = 0 ;
		int countB = 0 ;
		int sharedCount = 0 ;
		int tmp ;
		for (i = 0 ; i < KCOUNT_HASH_MAX ; ++i)
		{
			for (std::map<uint64_t, int>::iterator it = count[i].begin() ; it != count[i].end() ; ++it)
			{
				countA += it->second ;
				if (b.count[i].find(it->first) != b.count[i].end())
				{
					tmp = b.count[i][it->first] ;
					sharedCount += (countA < tmp ? countA : tmp) ;
				}
			}

			for (std::map<uint64_t, int>::iterator it = b.count[i].begin() ; it != b.count[i].end() ; ++it)
			{
				countB += it->second ;
			}
		}

		return (double)sharedCount / (countA + countB - sharedCount) ;
	}
} ;

#endif

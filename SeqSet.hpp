// The data structure holds the set of sequences
#ifndef _MOURISL_SEQSET_HEADER
#define _MOURISL_SEQSET_HEADER

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <string>

#include "SimpleVector.hpp"
#include "KmerIndex.hpp"
#include "ReadFiles.hpp"
#include "AlignAlgo.hpp"

struct _seqWrapper
{
	char *name ;
	char *consensus ; // This should be handled by malloc/free.
	int consensusLen ;
	int effectiveLen ; // the length that should be used for abundance estimation
	SimpleVector<struct _posWeight> posWeight ;
	bool isRef ; // this is from reference.

	int minLeftExtAnchor, minRightExtAnchor ; // only overlap with size larger than this can be counted as valid extension.

	struct _triple info[3] ; // For storing extra information. for ref, info[0,1] contains the coordinate for CDR1,2 and info[2].a for CDR3
					// In extending seqs with mate pair information, these are used to store rough V, J, C read coordinate.	
	int barcode ; // transformed barcode. -1: no barcode

	bool operator<( const struct _seqWrapper &b )	const
	{
		int i ;
		int weightA = 0 ;
		int weightB = 0 ;
		for ( i = 0 ; i < consensusLen ; ++i )
			weightA += posWeight[i].Sum() ;
		for ( i = 0 ; i < b.consensusLen ; ++i )
			weightB += b.posWeight[i].Sum() ;

		if ( weightA != weightB )
			return weightA > weightB ;
		else
			return consensusLen > b.consensusLen ;
	}
} ;

struct _hit
{
	struct _indexInfo indexHit ;
	int readOffset ;
	int strand ; // -1: different strand, 1: same strand. When strand==-1, the readOffset is the offset in the rcSeq.
	
	int repeats ; // record how many times this hit with other index part.

	bool operator<( const struct _hit &b ) const
	{
		if ( strand != b.strand )
			return strand < b.strand ;
		else if ( indexHit.idx != b.indexHit.idx )
			return indexHit.idx < b.indexHit.idx ;
		else if ( readOffset != b.readOffset )
			return  readOffset < b.readOffset ;
		else if ( indexHit.offset != b.indexHit.offset )
			return indexHit.offset < b.indexHit.offset ;
		
		return false ;
	}
} ;

struct _overlap
{
	int seqIdx ;
	int readStart, readEnd ; // When strand ==-1, the start,end is the offset in the rcSeq.
	int seqStart, seqEnd ;
	int strand ;
	
	int matchCnt ; // The number of matched bases, count TWICE.
	double similarity ;

	SimpleVector<struct _pair> *hitCoords ;
	SimpleVector<int> *info ; // store extra informations 
	int infoFromHits ; // some information obtained GetOverlapFromHits

	bool operator<( const struct _overlap &b ) const
	{
		// The overlap with more matched bases should come first.
		//if ( matchCnt > b.matchCnt + 2 || matchCnt < b.matchCnt - 2 )
		if ( matchCnt != b.matchCnt )
			return matchCnt > b.matchCnt ;
		else if ( similarity != b.similarity )
			return similarity > b.similarity ; 
		else if ( readEnd - readStart != b.readEnd - b.readStart )
			return readEnd - readStart > b.readEnd - b.readStart ;
		else if ( seqIdx != b.seqIdx )
			return seqIdx < b.seqIdx ;
		else if ( strand !=  b.strand )
			return strand < b.strand ;
		else if ( readStart != b.readStart )
			return readStart < b.readStart ;
		else if ( readEnd != b.readEnd )
			return readEnd < b.readEnd ;
		else if ( seqStart != b.seqStart )
			return seqStart < b.seqStart ;
		else 
			return seqEnd < b.seqEnd ; 

		return false ;
	}

	double UpdateSimilarity( int rlen, int slen, int mcnt )
	{
		double origLen = matchCnt / similarity ;
		similarity = ( matchCnt + mcnt ) / ( origLen + rlen + slen ) ;
		matchCnt += mcnt ;
	}
} ;

struct _fragmentOverlap 
{
	int seqIdx ;
	int seqStart, seqEnd ;

	int matchCnt ;
	double similarity ;

	bool hasMatePair ;

	struct _overlap overlap1 ; // for read 1
	struct _overlap overlap2 ; // for read 2
	
	bool operator<( const struct _fragmentOverlap &b ) const
	{
		if ( matchCnt != b.matchCnt )
			return matchCnt > b.matchCnt ;
		else if ( similarity != b.similarity )
			return similarity > b.similarity ; 
		return overlap1 < b.overlap1 ;
	}
} ;

// This order works better against reference set, because it seems works better for the 5' end site
struct _sortOverlapOnRef
{
	bool operator() (const struct _overlap &a, const struct _overlap &b) const
	{
		// The overlap with more matched bases should come first.
		//if ( a.matchCnt > b.matchCnt + 2 || a.matchCnt < b.matchCnt - 2 )
		if ( a.matchCnt != b.matchCnt )
			return a.matchCnt > b.matchCnt ;
		else if ( a.similarity != b.similarity )
			return a.similarity > b.similarity ; 
		else if ( a.readEnd - a.readStart != b.readEnd - b.readStart )
			return a.readEnd - a.readStart > b.readEnd - b.readStart ;
		else if ( a.strand !=  b.strand )
			return a.strand < b.strand ;
		else if ( a.seqStart != b.seqStart )
			return a.seqStart < b.seqStart ;
		else if ( a.seqEnd != b.seqEnd )
			return a.seqEnd < b.seqEnd ; 
		else if ( a.readStart != b.readStart )
			return a.readStart < b.readStart ;
		else if ( a.readEnd != b.readEnd )
			return a.readEnd < b.readEnd ;
		else 
			return a.seqIdx < b.seqIdx ;

		return false ;
	}
} ;

struct _assignRead
{
	char *id ;
	char *read ;
	int barcode ;
	int umi ;

	int info ; 
	struct _overlap overlap ;
} ;

class SeqSet
{
private:
	std::vector<struct _seqWrapper> seqs ;
	KmerIndex seqIndex ;
	int kmerLength ;
	int radius ;
	int hitLenRequired ;
	int gapN ;
	bool isLongSeqSet ; // Whether this seq set is built from long reads. Long reads may require more drastic filtration.

	// Some threshold
	double novelSeqSimilarity ;
	double refSeqSimilarity ;

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		if ( p1.b != p2.b )
			return p1.b < p2.b ;
		else
			return p1.a < p2.a ;
	}
	
	static bool CompSortPairAInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.a < p2.a ;
	}

	static bool CompSortOverlapsOnReadCoord( const struct _overlap &a, const struct _overlap &b )
	{
		return a.readStart < b.readStart ; 
	}

	static bool CompSortAssignedReadById( const struct _assignRead &a, const struct _assignRead &b )
	{
		return strcmp( a.id, b.id ) < 0 ;
	}
		
	static bool CompSortOverlapByCoord( const struct _overlap &a, const struct _overlap &b )	
	{
		if ( a.seqIdx != b.seqIdx )
			return a.seqIdx < b.seqIdx ;
		else if ( a.readStart != b.readStart )
			return a.readStart < b.readStart ;
		else
			return a.readEnd < b.readEnd ;
	}

	static bool CompSortHitCoordDiff( const struct _triple &a, const struct _triple &b )
	{
		if ( a.c != b.c )
			return a.c < b.c ;
		else if ( a.b != b.b )
			return a.b < b.b ;
		else
			return a.a < b.a ;
	}

	bool IsReverseComplement( char *a, char *b )
	{
		int i, j ;
		int len = strlen( a ) ;
		if ( len != strlen( b) ) 
			return false ;
		for ( i = 0, j = len - 1 ; i < len ; ++i, --j )
			if ( a[i] == 'N' && b[j] == 'N' )
				continue ;
			else if ( a[i] != 'N' && b[j] != 'N' )
			{
				if ( 3 - nucToNum[ a[i] - 'A' ] != nucToNum[ b[j] - 'A' ] )
					return false ;
			}
			else
				return false ;
		return true ;
	}

	void Reverse( char *r, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			r[i] = seq[len - 1 - i] ;  
		r[i] = '\0' ;
	}
	
	bool IsPosWeightCompatible( const struct _posWeight &a, const struct _posWeight &b )
	{
		int sumA = a.Sum() ;
		int sumB = b.Sum() ;
	
		if ( sumA == 0 || sumB == 0 
			|| ( sumA < 3 * a.count[0] && sumB < 3 * b.count[0] ) 
			|| ( sumA < 3 * a.count[1] && sumB < 3 * b.count[1] )
			|| ( sumA < 3 * a.count[2] && sumB < 3 * b.count[2] )
			|| ( sumA < 3 * a.count[3] && sumB < 3 * b.count[3] ) )
			return true ;
		return false ;	
	}

	bool IsOverlapIntersect( const struct _overlap &a, const struct _overlap &b )
	{
		if ( a.seqIdx == b.seqIdx && 
			( ( a.seqStart <= b.seqStart && a.seqEnd >= b.seqStart ) 
			|| ( b.seqStart <= a.seqStart && b.seqEnd >= a.seqStart ) ) )
			return true ;
		return false ;
	}

	// Return the first index whose hits.a is smaller or equal to valA
	int BinarySearch_LIS( int top[], int size, int valA, SimpleVector<struct _pair> &hits )
	{
		int l = 0, r = size - 1 ;
		int m ;
		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( valA == hits[ top[m] ].a )
			{
				return m ;
			}
			else if ( valA < hits[ top[m] ].a )
			{
				r = m - 1 ;
			}
			else
			{
				l = m + 1 ;
			}
		}
		return l - 1 ;
	}

	// The O(nlogn) method for solving LIS problem, suppose there are n hits.
	// Return the LIS, the LIS's length is returned by the function
	int LongestIncreasingSubsequence( SimpleVector<struct _pair> &hits, SimpleVector<struct _pair> &LIS ) 
	{
		// Only use the first hit of each qhit
		// Bias towards left

		int i, j, k ;
		int ret = 0 ;
		int size = hits.Size() ;

		int *record = new int[size] ; // The index of the selected hits
		int *top = new int[size] ; // record the index of the hits with smallest valB of the corresponding LIS length. kind of the top element.
		int *link = new int[size] ; // used to retrieve the LIS

		int rcnt = 1 ;
		record[0] = 0 ;
		for ( i = 1 ; i < size ; ++i )
		{
			//if ( hits[i].b == hits[i - 1].b )
			//	continue ;
			record[rcnt] = i ;
			++rcnt ;
		}
		top[0] = 0 ;
		link[0] = -1 ;
		ret = 1 ;
		for ( i = 1 ; i < rcnt ; ++i )
		{
			int tag = 0 ;
			if ( hits[ top[ ret - 1 ] ].a <= hits[ record[i] ].a )
				tag = ret - 1 ;
			else
				tag = BinarySearch_LIS( top, ret, hits[ record[i] ].a, hits ) ;			
			
			if ( tag == -1 )
			{
				top[0] = record[i] ;
				link[ record[i] ] = -1 ;
			}
			else if ( hits[ record[i] ].a > hits[ top[tag] ].a )
			{
				if ( tag == ret - 1 )
				{
					top[ret] = record[i] ;
					++ret ;
					link[ record[i] ] = top[tag] ;
				}
				else if ( hits[ record[i] ].a < hits[ top[tag + 1] ].a )
				{
					top[ tag + 1 ] = record[i] ;
					link[ record[i] ] = top[tag] ;
				}
			}
		}


		k = top[ret - 1] ;
		for ( i = ret - 1 ; i >= 0 ; --i )
		{
			LIS.PushBack( hits[k] ) ;
			k = link[k] ;	
		}
		LIS.Reverse() ;
		//for ( i = 0 ; i < ret ; ++i )
		//	LIS.PushBack( hits[ top[i] ] ) ;
		
		// Remove elements with same b.
		if ( ret > 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < ret ; ++i )
			{
				if ( LIS[i].b == LIS[k - 1].b )
					continue ;
				LIS[k] = LIS[i] ;
				++k ;
			}
			ret = k ;
		}

		delete []top ;
		delete []record ;
		delete []link ;

		return ret ;
	}

	void GetAlignStats( char *align, bool update, int &matchCnt, int &mismatchCnt, int &indelCnt)
	{
		int k ;
		if ( !update )
		{
			matchCnt = mismatchCnt = indelCnt = 0 ;
		}

		for ( k = 0 ; align[k] != -1 ; ++k )
		{
			if ( align[k] == EDIT_MATCH )
				++matchCnt ;
			else if ( align[k] == EDIT_MISMATCH )
				++mismatchCnt ;
			else 
				++indelCnt ;
		}
	}


	bool IsOverlapLowComplex( char *r, struct _overlap &o )
	{
		int cnt[4] = {0, 0, 0, 0} ;
		int i ;
		for ( i = o.readStart ; i <= o.readEnd ; ++i )
		{
			if ( r[i] == 'N' )
				continue ;
			++cnt[ nucToNum[ r[i] - 'A' ] ] ;
		}
		int len = o.readEnd - o.readStart + 1 ;
		int lowCnt = 0 ; 
		int lowTotalCnt = 0 ;
		for ( i = 0 ; i < 4 ; ++i )
		{
			if ( cnt[i] <= 2 )
			{
				++lowCnt ;
				lowTotalCnt += cnt[i] ;
			}
		}
		if ( lowTotalCnt * 7 >= o.readEnd - o.readStart + 1 )
			return false ;

		if ( lowCnt >= 2 )
			return true ;
		return false ;
	}

	// o is from the mate1 set. Test whether o's mate pair could be missing due
	// to too short reference sequence.
	bool TruncatedMatePairOverlap(struct _overlap &o, struct _overlap &compMate1, struct _overlap &compMate2)
	{
		if (o.seqIdx == -1 || compMate1.seqIdx == -1 || compMate2.seqIdx == -1)
			return false ;

		if (o.strand == 1) // mate is on the firght
		{
			if (seqs[o.seqIdx].consensusLen - 1 < o.seqEnd + compMate2.seqEnd - compMate1.seqEnd)
				return true ;
		}
		else if (o.strand == -1)
		{
			if (o.seqStart - (compMate1.seqStart - compMate2.seqStart) < 0)
				return true ;
		}

		return false ;	
	}

public:
	SeqSet( int kl ) 
	{
		kmerLength = kl ;
		radius = 10 ;
		hitLenRequired = 31 ;
		isLongSeqSet = false ;

		novelSeqSimilarity = 0.9 ;
		refSeqSimilarity = 0.8 ; 
	}

	~SeqSet() 
	{
		int size ;
		int i ;
		size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].consensus != NULL )
				free( seqs[i].consensus ) ;
			if ( seqs[i].name != NULL )
				free( seqs[i].name ) ;	
		}
	}

	int Size()
	{
		return seqs.size() ;
	}

	int SetRadius( int r )  
	{
		return radius = r ;
	}

	int SetHitLenRequired( int l )
	{
		return hitLenRequired = l ;
	}

	char *GetSeqName( int seqIdx ) 
	{
		return seqs[ seqIdx ].name ;
	}

	char *GetSeqConsensus( int seqIdx )
	{
		return seqs[ seqIdx ].consensus ;
	}

	int GetSeqConsensusLen( int seqIdx )
	{
		return seqs[seqIdx].consensusLen ;
	}

	int GetSeqEffectiveLen( int seqIdx )
	{
		return seqs[seqIdx].effectiveLen ;
	}
	
	// Input some baseline sequence to match against.
	void InputRefFa( char *filename ) 
	{
		int i, j, k ;
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		KmerCode kmerCode( kmerLength ) ;
		while ( fa.Next() )
		{
			// Insert the kmers 
			struct _seqWrapper ns ;
			ns.name = strdup( fa.id ) ;
			ns.isRef = true ;

			int id = seqs.size() ;
			seqs.push_back( ns ) ;

			struct _seqWrapper &sw = seqs[id] ;
			int seqLen = strlen( fa.seq ) ;
			sw.consensus = strdup( fa.seq ) ;	
			sw.consensusLen = strlen( sw.consensus );
			sw.effectiveLen = sw.consensusLen ; //effectiveLen ; 	
			sw.barcode = -1 ;
			seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, sw.consensusLen, id ) ;
		}
	}

	// Compute the length of hit from the read, take the overlaps of kmer into account 
	int GetTotalHitLengthOnRead( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;
		//for ( i = 0 ; i < hitSize ; ++i )
		//	printf( "%d %d\n", i, hits[i].readOffset) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
			{
				if ( hits[j].readOffset > hits[j - 1].readOffset + kmerLength - 1 )	
					break ;
			}

			ret += hits[j - 1].readOffset - hits[i].readOffset + kmerLength ;

			i = j ;
		}
		return ret ;
	}
	
	int GetTotalHitLengthOnSeq( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;

		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].indexHit.offset > hits[j - 1].indexHit.offset + kmerLength - 1 )
					break ;
			ret += hits[j - 1].indexHit.offset - hits[i].indexHit.offset + kmerLength ;
			i = j ;
		}
		return ret ;
	}

	int GetHitsFromRead( char *read, char *rcRead, int strand, int barcode, bool allowTotalSkip, SimpleVector<struct _hit> &hits, SimpleVector<bool> *puse )
	{
		int i, j ;
		int len = strlen( read ) ;

		KmerCode kmerCode( kmerLength ) ;
		KmerCode prevKmerCode( kmerLength ) ;

		// Locate the hits from the same-strand case.
		//int skipLimit = 3 ;
		int skipLimit = kmerLength / 2 ; 

		int skipCnt = 0 ;
		int downSample = 1 ;
		if ( len > 200 && isLongSeqSet )
		{
			downSample = 1 + len / 200 ;
			//skipLimit /= downSample ;
			//if ( skipLimit < 2 )
			//	skipLimit = 2 ; 
		}

		if ( strand != -1 )
		{
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( read[i] ) ;

			for ( ; i < len ; ++i )
			{
				kmerCode.Append( read[i] ) ;
				if ( downSample > 1 && ( i - kmerLength + 1 ) % downSample != 0 )
					continue ;

				if ( i == kmerLength - 1 || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;
					if ( size >= 100 && puse == NULL /*&& barcode == -1*/ && i != kmerLength - 1 && i != len - 1 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}

					if ( size >= 100 && allowTotalSkip )
						continue ;

					skipCnt = 0 ;
					int repeats = size ;
					if ( puse != NULL )
					{
						repeats = 0 ;
						for ( j = 0 ; j < size ; ++j )
						{
							if ( !puse->Get( indexHit[j].idx ) ) 
								continue ;
							++repeats ;
						}
					}

					if ( barcode != -1 )
						repeats = 1 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						nh.indexHit = indexHit[j] ;
						nh.readOffset = i - kmerLength + 1 ;
						nh.strand = 1 ;
						nh.repeats = repeats ;
						if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
							continue ;
						if ( barcode != -1 && seqs[ indexHit[j].idx ].barcode != barcode )
							continue ;
						hits.PushBack( nh ) ;
					}
				}

				prevKmerCode = kmerCode ;
			}
		}
		// Locate the hits from the opposite-strand case.
		ReverseComplement( rcRead, read, len ) ;		

		if ( strand != 1 )
		{
			kmerCode.Restart() ;
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( rcRead[i] ) ;

			skipCnt = 0 ; 
			for ( ; i < len ; ++i )
			{
				kmerCode.Append( rcRead[i] ) ;
				if ( downSample > 1 && ( i - kmerLength + 1 ) % downSample != 0 )
					continue ;
				if ( i == kmerLength - 1 || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;

					if ( size >= 100 && puse == NULL /*&& barcode == -1*/ && i != kmerLength - 1 && i != len - 1 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}
					if ( size >= 100 && allowTotalSkip )
						continue ;

					skipCnt = 0 ;
					
					int repeats = size ;
					if ( puse != NULL )
					{
						repeats = 0 ;
						for ( j = 0 ; j < size ; ++j )
						{
							if ( !puse->Get( indexHit[j].idx ) ) 
								continue ;
							++repeats ;
						}
					}

					if ( barcode != -1 )
						repeats = 1 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						nh.indexHit = indexHit[j] ;
						nh.readOffset = i - kmerLength + 1 ;
						nh.strand = -1 ;
						nh.repeats = repeats ;
						if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
							continue ;
						if ( barcode != -1 && seqs[ indexHit[j].idx ].barcode != barcode )
							continue ;

						hits.PushBack( nh ) ;
						/*if ( seqs[indexHit[j].idx].name == NULL)
						{
							printf( "%d %d\n", indexHit[j].idx, indexHit[j].offset ) ;
						}*/
						assert( seqs[indexHit[j].idx].name != NULL ) ;
					}
				}

				prevKmerCode = kmerCode ;
			}
		}
		return hits.Size() ;
	}

	// Use the hits to extract overlaps from SeqSet
	int GetOverlapsFromHits( SimpleVector<struct _hit> &hits, int hitLenRequired, int filter, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _triple> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;
		int maxReadOffset = -1 ;

		for (i = 0 ; i < hitSize ; ++i)
			if (hits[i].readOffset > maxReadOffset)
				maxReadOffset = hits[i].readOffset ;
		int *readOffsetUsed = new int[maxReadOffset + 1] ;
		
		// Compute the minHitRequired. 
		// NOTE: each strand should have its own minHitRequired, it could be that on one strand,
		//    each hit is matched to too many places and the skip hits mechanism is triggered.
		int novelMinHitRequired[2] = {3, 3} ;
		int refMinHitRequired[2] = {3, 3} ;
		bool removeOnlyRepeats[2] = {false, false} ; // Remove the hits on a seq that are all repeats hit.
		int possibleOverlapCnt[2] = {0, 0} ;
		if ( filter == 1 )
		{
			int longestHits[2] = {0, 0} ;
			for ( i = 0 ; i < hitSize ; ++i )
			{
				int isPlusStrand = ( 1 + hits[i].strand ) / 2 ; 
				for ( j = i + 1 ; j < hitSize ; ++j )
					if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
						break ;
				if ( !seqs[ hits[i].indexHit.idx].isRef )
				{
					if ( j - i > novelMinHitRequired[ isPlusStrand ] )
						++possibleOverlapCnt[ isPlusStrand ] ;
					if ( j - i > longestHits[ isPlusStrand ] )
						longestHits[ isPlusStrand] = j - i  ;
				}
				
				if ( !removeOnlyRepeats[isPlusStrand] )
				{
					int cnt = 0 ;
					for ( k = i ; k < j ; ++k )
						if ( hits[k].repeats <= 10000 )
							++cnt ;
					if ( cnt >= novelMinHitRequired[ isPlusStrand ] )
					{
						removeOnlyRepeats[ isPlusStrand ] = true ;
					}
				}

				i = j ;
			}
			// filter based on the repeatability of overlaps.
			for ( i = 0 ; i <= 1 ; ++i )
			{
				if ( possibleOverlapCnt[i] > 100000 )
					novelMinHitRequired[i] = longestHits[i] * 0.75 ;
				else if ( possibleOverlapCnt[i] > 10000 )
					novelMinHitRequired[i] = longestHits[i] / 2 ;
				else if ( possibleOverlapCnt[i] > 1000 )
					novelMinHitRequired[i] = longestHits[i] / 3 ;
				else if ( possibleOverlapCnt[i] > 100 )
					novelMinHitRequired[i] = longestHits[i] / 4 ;
			}
		}

		//if ( novelMinHitRequired > 3 )
		//	printf( "novelMinHitRequired=%d\n", novelMinHitRequired ) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
					break ;

			int minHitRequired = novelMinHitRequired[ ( 1 + hits[i].strand ) / 2 ] ;
			if ( seqs[ hits[i].indexHit.idx].isRef )
				minHitRequired = refMinHitRequired[ ( 1 + hits[i].strand ) / 2 ];
			
			//[i,j) holds the hits onto the same seq on the same strand.	
			if ( j - i < minHitRequired )
			{
				i = j ;
				continue ;
			}
			
			if ( removeOnlyRepeats[( 1 + hits[i].strand ) / 2] )
			{
				bool hasUnique = false ;
				for ( k = i ; k < j ; ++k )
				{
					if ( hits[k].repeats <= 10000 )
					{
						hasUnique = true ;
						break ;
					}
				}
				if ( !hasUnique )
				{
					i = j ;
					continue ;
				}
			}

			hitCoordDiff.Clear() ;
			for ( k = i ; k < j ; ++k )
			{
				struct _triple nh ;
				nh.a = hits[k].readOffset ;
				nh.b = hits[k].indexHit.offset ;
				nh.c = hits[k].readOffset - hits[k].indexHit.offset ;
				hitCoordDiff.PushBack( nh ) ;
			}
			std::sort( hitCoordDiff.BeginAddress(), hitCoordDiff.EndAddress(), CompSortHitCoordDiff ) ;

			// Pick the best concordant hits.
			int s, e ;
			int adjustRadius = radius ;
			if ( !seqs[ hits[i].indexHit.idx ].isRef )
				adjustRadius = 0 ;
			
			int currCoordDiffCnt = 0 ;
			int currCoordDiff = 0 ;
			int dominantCoordDiff = 0 ;
			int dominantCoordDiffCnt = 0 ;

			for ( s = 0 ; s < j - i ; )
			{
				int diffSum = 0 ;
				currCoordDiff = hitCoordDiff[s].c ;
				currCoordDiffCnt = 1 ;
				dominantCoordDiffCnt = 0 ;
				readOffsetUsed[hitCoordDiff[s].a] = -1 ; 
				for ( e = s + 1 ; e < j - i ; ++e )
				{
					int diff = hitCoordDiff[e].c - hitCoordDiff[e - 1].c ;
					if ( diff < 0 )
						diff = -diff ;
					
					if ( diff > adjustRadius ) 
						break ;
					
					if (diff == 0)
						++currCoordDiffCnt ;
					else
					{
						if (currCoordDiffCnt > dominantCoordDiffCnt)
						{
							dominantCoordDiff = currCoordDiff ;
							dominantCoordDiffCnt = currCoordDiffCnt ;
						}
						currCoordDiff = hitCoordDiff[e].c ;
						currCoordDiffCnt = 1 ;
					}
					//printf("%d %d\n", currCoordDiffCnt, hitCoordDiff[e].c) ;
					
					readOffsetUsed[hitCoordDiff[e].a] = -1 ; 
					diffSum += diff ; 
				}
				if (currCoordDiffCnt > dominantCoordDiffCnt)
				{
					dominantCoordDiff = currCoordDiff ;
					dominantCoordDiffCnt = dominantCoordDiffCnt ;
				}

				//printf( "%d %d: %d %d\n", i, j, s, e ) ;
				if ( e - s < minHitRequired 
					|| ( e - s ) * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}


				if ( removeOnlyRepeats[( 1 + hits[i].strand ) / 2] )
				{
					bool hasUnique = false ;
					for ( k = s ; k < e ; ++k )
					{
						if ( hits[k].repeats <= 10000 )
						{
							hasUnique = true ;
							break ;
						}
					}
					if ( !hasUnique )
					{
						s = e ;
						continue ;
					}
				}

				// [s, e) holds the candidate in the array of hitCoordDiff 
				concordantHitCoord.Clear() ;
				
				for ( k = s ; k < e ; ++k )
				{
					struct _pair nh ;
					nh.a = hitCoordDiff[k].a ;
					nh.b = hitCoordDiff[k].b ;
					concordantHitCoord.PushBack( nh ) ;
				}
				
				if ( adjustRadius > 0 )
				{
					for (k = 0 ; k < e - s ; ++k)
						if (readOffsetUsed[ concordantHitCoord[k].a ] == -1
								|| readOffsetUsed[concordantHitCoord[k].a] > ABS(concordantHitCoord[k].a - concordantHitCoord[k].b - dominantCoordDiff))
						{
							readOffsetUsed[concordantHitCoord[k].a] = ABS(concordantHitCoord[k].a - concordantHitCoord[k].b - dominantCoordDiff ) ;
						}
					int l = 0 ;
					for (k = 0 ; k < e - s ; ++k)
					{
						if (ABS(concordantHitCoord[k].a - concordantHitCoord[k].b - dominantCoordDiff) == readOffsetUsed[concordantHitCoord[k].a] )
						{
							concordantHitCoord[l] = concordantHitCoord[k] ;
							++l ;
						}
					}
					concordantHitCoord.Resize( l ) ;
					std::sort( concordantHitCoord.BeginAddress(), concordantHitCoord.EndAddress(), CompSortPairBInc ) ;
				}
				//for ( k = 0 ; k < e - s ; ++k )	
				//	printf( "%d (%d-%d): %d %s %d %d\n", i, s, e, hits[i].indexHit.idx, seqs[ hits[i].indexHit.idx ].name, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;

				// Compute the longest increasing subsequence.
				//printf( "lis for %d (%d %d; %d %d). strand=%d (%d)\n", e - s, i, j, s, e, hits[i].strand, seqs.size() ) ;
				hitCoordLIS.Clear() ;
				int lisSize = LongestIncreasingSubsequence( concordantHitCoord, hitCoordLIS ) ; 
				if ( lisSize * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				// Rebuild the hits.
				int lisStart = 0 ;
				int lisEnd = lisSize - 1 ;
				// Ignore long insert gaps.
				if ( isLongSeqSet )
				{
					int maxGap = 2 * hitLenRequired + 3 * kmerLength ;
					if ( filter == 0 )//&& possibleOverlapCnt[( 1 + hits[i].strand ) / 2] > 1000 )
						maxGap *= 4 ;
					if ( maxGap < 200 )
						maxGap = 200 ;
					int max = -1 ;
					for ( k = 0 ; k < lisSize ; )
					{
						int l ;
						for ( l = k + 1 ; l < lisSize ; ++l )
						{
							if ( hitCoordLIS[l].a - hitCoordLIS[l - 1].a > maxGap )
								break ;
						}
						if ( l - k > max )
						{
							max = l - k ;
							lisStart = k ;
							lisEnd = l - 1 ;
						}

						k = l ;	
					}
				}

				finalHits.Clear() ;
				for ( k = lisStart ; k <= lisEnd ; ++k )
				{
					struct _hit nh = hits[i];
					nh.readOffset = hitCoordLIS[k].a ;
					nh.indexHit.offset = hitCoordLIS[k].b ;
					//if (seqs.size() == 1 )
					//	printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.strand ) ;
					finalHits.PushBack( nh ) ;
				}
				lisSize = lisEnd - lisStart + 1 ;

				int hitLen = GetTotalHitLengthOnRead ( finalHits ) ;
				if ( hitLen < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				else if ( GetTotalHitLengthOnSeq( finalHits ) < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				
				struct _overlap no ;
				no.seqIdx = hits[i].indexHit.idx ;
				no.readStart = finalHits[0].readOffset ;
				no.readEnd = finalHits[ lisSize - 1 ].readOffset + kmerLength - 1 ;
				no.strand = finalHits[0].strand ;
				no.seqStart = finalHits[0].indexHit.offset ;
				no.seqEnd = finalHits[ lisSize - 1 ].indexHit.offset + kmerLength - 1 ;
				no.matchCnt = 2 * hitLen ;
				no.similarity = 0 ;

				if ( !seqs[ no.seqIdx ].isRef && hitLen * 2 < no.seqEnd - no.seqStart + 1 )
				{
					s = e ; 
					continue ;
				}
				no.hitCoords = new SimpleVector<struct _pair> ;
				no.hitCoords->Reserve( lisSize ) ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _pair nh ;
					nh.a = finalHits[k].readOffset ;
					nh.b = finalHits[k].indexHit.offset ;
					no.hitCoords->PushBack( nh ) ;
				}
				overlaps.push_back( no ) ;

				s = e ;
			} // iterate through concordant hits.
			i = j ;
		}
		delete []readOffsetUsed ;
		return overlaps.size() ;
	}

	void SortHits( SimpleVector<struct _hit> &hits, bool alreadyReadOrder )
	{
		int i, k ;
		if ( hits.Size() > 2 * seqs.size() && alreadyReadOrder ) 
		{
			// Bucket sort.
			int hitCnt = hits.Size() ;
			int seqCnt = seqs.size() ;
			SimpleVector<struct _hit> *buckets[2] ;
			buckets[0] = new SimpleVector<struct _hit>[seqCnt] ;
			buckets[1] = new SimpleVector<struct _hit>[seqCnt] ;

			for ( i = 0 ; i < hitCnt ; ++i )
			{
				int tag = hits[i].strand == 1 ? 1 : 0 ;
				buckets[tag][ hits[i].indexHit.idx ].PushBack( hits[i] ) ;
			}
			
			hits.Clear() ;
			for ( k = 0 ; k <= 1 ; ++k )
			{
				for ( i = 0 ; i < seqCnt ; ++i )
				{
					hits.PushBack( buckets[k][i] ) ;
				}
			}

			delete[] buckets[0] ;
			delete[] buckets[1] ;
		}
		else
			std::sort( hits.BeginAddress(), hits.EndAddress() ) ;
	}

	// Obtain the overlaps, each overlap further contains the hits induce the overlap. 
	// Return: the number of overlaps.
	int GetOverlapsFromRead( char *read, int strand, int barcode, std::vector<struct _overlap> &overlaps)
	{
		int i, j, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return -1 ;
		
		int overlapCnt = 0 ;
		SimpleVector<struct _hit> hits ;
		char *rcRead =  new char[len + 1] ;
		
		GetHitsFromRead( read, rcRead, strand, barcode, false, hits, NULL ) ;
		SortHits( hits, true ) ;
		// Find the overlaps.
		//if ( seqs.size() == 1 )
		//	for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//		printf( "- %d %s %d %d\n", it->readOffset, seqs[ it->indexHit.idx ].name, it->indexHit.offset, it->strand ) ;
		//if ( seqs.size() == 1 )
		//	for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//		printf( "- %d %s %d %d\n", it->readOffset, seqs[ it->indexHit.idx ].name, it->indexHit.offset, it->strand ) ;
		overlapCnt = GetOverlapsFromHits( hits, hitLenRequired, 0, overlaps ) ;
		std::sort( overlaps.begin(), overlaps.end() ) ;
			
		//for ( i = 0 ; i < overlapCnt ; ++i )
		//	printf( "%d: %d %s %d. %d %d %d %d. %d\n", i, overlaps[i].seqIdx,seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].matchCnt ) ; 
		if ( overlapCnt > 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < overlapCnt ; ++i )
			{
				if ( overlaps[i].strand != overlaps[0].strand )
				{
					delete overlaps[i].hitCoords ;
					overlaps[i].hitCoords = NULL ;
					continue ;		
				}
				if ( i != k )
					overlaps[k] = overlaps[i] ;
				++k ;
			}

			overlaps.resize( k ) ;
			overlapCnt = k ;
		}
		/*for ( i = 1 ; i < overlapCnt ; ++i )
		{
			if ( overlaps[i].strand != overlaps[i - 1].strand )
			{
				overlaps.clear() ;
				for ( i = 0 ; i < overlapCnt ; ++i )
					delete overlaps[i].hitCoords ;
				return 0 ;
			}
		}*/
		ReverseComplement( rcRead, read, len ) ;		
		
		int firstRef = -1 ;
		int bestNovelOverlap = -1 ;
		SimpleVector<struct _pair> readOverlapRepresentatives ; // The non-subset best overlaps
		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *r ;
			if ( overlaps[i].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			overlaps[i].infoFromHits = i ;

			SimpleVector<struct _pair> &hitCoords = *overlaps[i].hitCoords ; 	
			int hitCnt = hitCoords.Size() ;
			int matchCnt = 0, mismatchCnt = 0, indelCnt = 0  ;
			double similarity = 1 ;
			
			//if ( overlaps[i].seqIdx == 4 )
			//	printf( "%d %d %d. %d %d\n", seqs[ overlaps[i].seqIdx ].isRef, overlapCnt, bestNovelOverlap,
			//			overlaps[i].matchCnt, overlaps[ bestNovelOverlap ].matchCnt ) ;
			// Some fast pre-filters
			if ( seqs[ overlaps[i].seqIdx ].isRef )
			{
				if ( firstRef == -1 )
					firstRef = i ;
				/*else
				{
					if ( overlaps[i].matchCnt < 0.9 * overlaps[ firstRef ].matchCnt )
					{
						overlaps[i].similarity = 0 ;  // No need to look into this.
						continue ;
					}
				}*/
			}
			
			matchCnt += 2 * kmerLength ;
			char *align = new char[ overlaps[i].readEnd - overlaps[i].readStart + 1 + 
				overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 1] ;
			for ( j = 1 ; j < hitCnt ; ++j )
			{
				if ( hitCoords[j - 1].b - hitCoords[j - 1].a == hitCoords[j].b - hitCoords[j].a )
				{
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ; 
						
						if ( seqs[ overlaps[i].seqIdx ].isRef  )
						{
							//printf( "Use ref %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;
						}
						else
						{
						//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							//printf( "Use novel %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;

							/*if ( seqs.size() == 1 )
							{
								AlignAlgo::VisualizeAlignment( 
										seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
										hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
										r + hitCoords[j - 1].a + kmerLength, 
										hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
										align ) ;

							}*/
						}

						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;
						

						if ( ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef ) && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}
					}
				}
				else
				{
					if ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef )
					{
						similarity = 0 ;
						break ;
					}

					//printf( "%d %d=>%d %d\n", hitCoords[j - 1].a, hitCoords[j - 1].b, hitCoords[j].a, hitCoords[j].b ) ;
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 < hitCoords[j].b )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ; //+ kmerLength ;
						// Make the two kmer hit match on coordinate.
						indelCnt += ( hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) + 
							( hitCoords[j].a + kmerLength - hitCoords[j - 1].a )  ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 < hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += kmerLength + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ( hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) +
							( hitCoords[j].b + kmerLength - hitCoords[j - 1].b ) ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a &&
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += ( hitCoords[j].a - hitCoords[j - 1].a ) + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * MIN( hitCoords[j].a - hitCoords[j - 1].a, hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ABS( ( hitCoords[j].a - hitCoords[j].b ) - 
							( hitCoords[j - 1].a - hitCoords[j - 1].b ) ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ;
						 
						if ( seqs[ overlaps[i].seqIdx ].isRef )
						{
							//printf( "Use ref2 %d %d.\n", hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
							//	hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) ) ;
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						else
						{
							//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						
						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;
						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

					}
				}
			} // for j
			delete[] align ;
			
			//printf( "%d: %d %d %d %lf\n", overlaps[i].seqIdx, matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, similarity ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( similarity == 1 )
				overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
			else
				overlaps[i].similarity = 0 ;
			
			if ( IsOverlapLowComplex( r, overlaps[i]) )
				overlaps[i].similarity = 0 ;
			
		  //printf( "%d %s: %d %d %d %lf. %d %d\n", overlaps[i].seqIdx, seqs[overlaps[i].seqIdx].name, matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, overlaps[i].similarity, hitCnt, overlaps[i].matchCnt ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity > 0 )
			{
				if ( bestNovelOverlap == -1 || overlaps[i] < overlaps[ bestNovelOverlap ] ) // the less than means has higher priority
				{
					bestNovelOverlap = i ;
				}
			}
			/*if ( overlaps[i].similarity > 0 )
			{
				printf( "%d: %d %d %d %d %d %lf\n", matchCnt, overlaps[i].seqIdx, overlaps[i].readStart, overlaps[i].readEnd, 
							overlaps[i].seqStart, overlaps[i].seqEnd, similarity ) ;
			}
			assert( overlaps[i].similarity <= 1 ) ;*/

			if ( overlaps[i].similarity > 0 )
			{
				int size = readOverlapRepresentatives.Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					int k = readOverlapRepresentatives[j].a ;
					if ( overlaps[i].readStart >= overlaps[k].readStart 
						&& overlaps[i].readEnd <= overlaps[k].readEnd )
						break ;		
				}

				if ( j >= size )
				{
					struct _pair np ;
					np.a = i ;
					np.b = 0 ;
					readOverlapRepresentatives.PushBack( np ) ;
				}
			}
		} // for i
		delete[] rcRead ;

		// Release the memory for hitCoords.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			overlaps[i].hitCoords->Release() ;
			delete overlaps[i].hitCoords ;
			overlaps[i].hitCoords = NULL ;
		}

		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < refSeqSimilarity )
				continue ;
			else if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < novelSeqSimilarity )
				continue ;

			//printf( "%d %s: %d-%d %d-%d: %d %d %lf\n", overlaps[i].seqIdx, seqs[overlaps[i].seqIdx].name, 
			//	overlaps[i].readStart, overlaps[i].readEnd, 
			//	overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].matchCnt, overlaps[i].strand, overlaps[i].similarity ) ;
			overlaps[k] = overlaps[i] ;
			++k ;
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		
		//printf( "return: %d\n", overlapCnt) ;
		return overlapCnt ;
	}

	// Test whether the read share a kmer hit on the seqs.
	bool HasHitInSet( char *read, char *rcRead )
	{
		int i, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return false ;

		SimpleVector<struct _hit> hits ;	
		GetHitsFromRead( read, rcRead, 0, -1, false, hits, NULL  ) ;

		int hitCnt = hits.Size() ;
		if ( hitCnt == 0 )
			return false ;

		// Bucket sort.
		int seqCnt = seqs.size() ;
		SimpleVector<struct _hit> *buckets[2] ;
		buckets[0] = new SimpleVector<struct _hit>[seqCnt] ;
		buckets[1] = new SimpleVector<struct _hit>[seqCnt] ;

		for ( i = 0 ; i < hitCnt ; ++i )
		{
			int tag = hits[i].strand == 1 ? 1 : 0 ;
			buckets[tag][ hits[i].indexHit.idx ].PushBack( hits[i] ) ;
		}
		
		// Find the best bucket.
		int max = -1 ;
		int maxTag = -1 ;
		int maxSeqIdx = -1 ;
		for ( k = 0 ; k <= 1 ; ++k )
		{
			for ( i = 0 ; i < seqCnt ; ++i )
			{
				int size = buckets[k][i].Size() ;
				if ( size > 0 && size > max )
				{
					maxTag = k ;
					maxSeqIdx = i ;
					max = size ;
				}
			}
		}
		
		std::vector<struct _overlap> overlaps ;
		GetOverlapsFromHits( buckets[maxTag][maxSeqIdx], hitLenRequired, 0, overlaps ) ;
		delete[] buckets[0] ;
		delete[] buckets[1] ;
		//printf( "%d %d\n", hitCnt, overlaps.size() ) ;	
		for ( i = 0 ; i < overlaps.size() ; ++i )
		{
			//printf( "%s\n", seqs[ overlaps[i].seqIdx ].name ) ;
			delete overlaps[i].hitCoords ;
		}
		if ( overlaps.size() == 0 )
			return false ;
		/*printf( "%s %d %d %lf\n", seqs[ overlaps[0].seqIdx ].name, overlaps[0].readStart, overlaps[0].readEnd, overlaps[0].similarity ) ;
		for ( i = 0 ; i < overlaps[0].hitCoords->Size() ; ++i )
			printf( "%d %d\n", overlaps[0].hitCoords->Get(i).a, overlaps[0].hitCoords->Get(i).b ) ;*/
		return true ;
	}

	// Extend the overlap to include the overhang parts and filter the overlaps if the overhang does not match well.
	// return: whether this is a valid extension or not
	int ExtendOverlap( char *r, int len, struct _seqWrapper &seq,
		char *align, struct _overlap &overlap, struct _overlap &extendedOverlap )
	{
		// Check whether the overhang part is compatible with each other or not.
		// Extension to 5'-end ( left end )
		int matchCnt, mismatchCnt, indelCnt ;
		int leftOverhangSize = MIN( overlap.readStart, overlap.seqStart ) ;
		int ret = 1 ;
		int i, k ;
		int goodLeftOverhangSize = 0 ;
		AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqStart - leftOverhangSize,
				leftOverhangSize, 
				r + overlap.readStart - leftOverhangSize, leftOverhangSize, align ) ;
		GetAlignStats( align, false, matchCnt, mismatchCnt, indelCnt ) ;
		
		/*for ( i = 0 ; align[i] != -1 ; ++i )
			;
		int tmpMatchCnt = 0  ;
		for ( i = i - 1, k = 1 ; i >= 0 ; --i, ++k )
		{
			if ( align[i] == EDIT_MATCH )
			{
				++tmpMatchCnt ;
				if ( tmpMatchCnt > 0.75 * k )
					goodLeftOverhangSize = k ;
			}
			else if ( align[i] != EDIT_MISMATCH )
				break ;
		}*/
		// Extension to 3'-end ( right end )
		int rightOverhangSize = MIN( len - 1 - overlap.readEnd, seq.consensusLen - 1 - overlap.seqEnd ) ;
		int goodRightOverhangSize = 0 ;
		AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqEnd + 1, 
				rightOverhangSize,
				r + overlap.readEnd + 1, rightOverhangSize, align ) ;
		int oldIndelCnt = indelCnt ;
		GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;

		extendedOverlap.seqIdx = overlap.seqIdx ;
		extendedOverlap.readStart = overlap.readStart - leftOverhangSize ;
		extendedOverlap.readEnd = overlap.readEnd + rightOverhangSize ;
		extendedOverlap.seqStart = overlap.seqStart - leftOverhangSize ;
		extendedOverlap.seqEnd = overlap.seqEnd + rightOverhangSize ;
		extendedOverlap.strand = overlap.strand ;	
		extendedOverlap.matchCnt = 2 * matchCnt + overlap.matchCnt ;
		extendedOverlap.similarity = (double)( 2 * matchCnt + overlap.matchCnt ) / 
			( extendedOverlap.readEnd - extendedOverlap.readStart + 1 + extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ) ;	
		//printf("%d %d %d %d. %d\n", extendedOverlap.readStart, extendedOverlap.readEnd, extendedOverlap.seqStart, extendedOverlap.seqEnd,
		//		extendedOverlap.matchCnt);		
		if (extendedOverlap.similarity < refSeqSimilarity)
			ret = 0 ;
		/*if (leftOverhangSize + rightOverhangSize > 10
				&& (double)matchCnt / (leftOverhangSize + rightOverhangSize) < 0.5)
			ret = 0 ;*/
		/*else
		{
			AlignAlgo::GlobalAlignment(seq.consensus + extendedOverlap.seqStart, extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ,
					r, len, align) ;
			AlignAlgo::VisualizeAlignment(seq.consensus + extendedOverlap.seqStart, extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ,
					r, len, align) ;
		}*/
		return ret ;
	}


	void ReverseComplement( char *rcSeq, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
		{
			if ( seq[len - 1 - i] != 'N' )
				rcSeq[i] = numToNuc[ 3 - nucToNum[seq[len - 1 - i] - 'A'] ];
			else
				rcSeq[i] = 'N' ;
		}
		rcSeq[i] = '\0' ;
	}
	
	// Find the seq id this read belongs to.
	int AssignRead( char *read, char *read2, int barcode, std::vector<struct _fragmentOverlap> &assign )
	{
		assign.clear() ;
		
		int i, j ;
		std::vector<struct _overlap> overlaps ;
		std::vector<struct _overlap> overlaps2 ;
		
		int overlapCnt = GetOverlapsFromRead( read, 0, barcode, overlaps ) ;
		int overlapCnt2 = 0 ;
		if (read2 != NULL)
		{
			overlapCnt2 = GetOverlapsFromRead(read2, 0, barcode, overlaps2 ) ;
		}
		//printf( "%d %d\n", overlapCnt, mateOverlapCnt ) ;
		//printf( "%d %s\n%d %s\n", overlaps[0].strand, reads[i].seq, mateOverlaps[0].strand, reads[i + 1].seq ) ;

		if ( overlapCnt <= 0 || seqs.size() == 0
				|| (read2 && overlapCnt2 <=0 ) )
		{
			return -1 ;
		}
		
		std::sort( overlaps.begin(), overlaps.end() ) ;

		int len = strlen( read ) ;
		int len2 = 0 ;
		char *rc = new char[len + 1] ;
		char *rc2 = NULL ;

		ReverseComplement( rc, read, len ) ;

		if (read2)
		{
			len2 = strlen(read2) ;
			rc2 = new char[len + 1] ;
			ReverseComplement( rc2, read2, len2 ) ;
		}

		char *r = read ;
		if ( overlaps[0].strand == -1 )
			r = rc ;

		std::vector<struct _overlap> extendedOverlaps ;
		std::vector<struct _overlap> extendedOverlaps2 ;
		struct _overlap eOverlap ;
		
		char *align = new char[ 2 * len + 2 ] ;
		
		int extendCnt = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			//printf( "%d %d: %d-%d %d-%d %lf\n", i, overlaps[i].seqIdx, overlaps[i].readStart, overlaps[i].readEnd,
			//		overlaps[i].seqStart, overlaps[i].seqEnd, overlaps[i].similarity) ;
			if ( overlaps[i].readEnd - overlaps[i].readStart > len / 2 
					&& ExtendOverlap( r, len, seqs[ overlaps[i].seqIdx ], align, overlaps[i], eOverlap ) == 1 )
			{
				extendedOverlaps.push_back(eOverlap) ;
			}
		}

		if (read2 != NULL)
		{
			r = read2 ;
			if (overlaps2[0].strand == -1)
				r = rc2 ;
			for (i = 0 ; i < overlapCnt2 ; ++i)
			{
				if ( overlaps2[i].readEnd - overlaps2[i].readStart > len2 / 2 
						&& ExtendOverlap(r, len, seqs[overlaps2[i].seqIdx], align, overlaps2[i], eOverlap) == 1)
				{
					extendedOverlaps2.push_back(eOverlap) ;
				}
			}
		}

		// Match mate pairs
		SimpleVector<struct _pair> fragments ;
		overlapCnt = extendedOverlaps.size() ;
		overlapCnt2 = extendedOverlaps2.size() ;
		if (read2 == NULL)
		{
			for (i = 0 ; i < overlapCnt ; ++i)
			{
				struct _pair nf ;
				nf.a = i ;
				nf.b = -1 ;
				fragments.PushBack(nf) ;
			}
		}
		else
		{
			// TODO: improve the efficiency
			//printf("%d %d\n", overlapCnt, overlapCnt2);
			//printf("%d %d\n", extendedOverlaps[0].strand, extendedOverlaps2[0].strand);
			for (i = 0 ; i < overlapCnt ; ++i)
			{
				for (j = 0 ; j < overlapCnt2 ; ++j)
				{
					// compatible mate pairs
					if ( extendedOverlaps[i].strand == extendedOverlaps2[j].strand
							|| extendedOverlaps[i].seqIdx != extendedOverlaps2[j].seqIdx)
						continue ;
					if ((extendedOverlaps[i].strand == 1 && extendedOverlaps[i].seqStart < extendedOverlaps2[j].seqStart)
							|| (extendedOverlaps[i].strand == -1 && extendedOverlaps[i].seqStart > extendedOverlaps2[j].seqStart))
					{
						struct _pair nf ;
						nf.a = i ;
						nf.b = j ;
						//printf("hi: %d %d\n", i, j);
						fragments.PushBack(nf) ;
					}
				}
			}
		}

		int fragmentCnt = fragments.Size() ;
		std::map<int, int> seqIdxToOverlapIdx ;
		for (i = 0 ; i < fragmentCnt ; ++i)
		{
			struct _fragmentOverlap fragmentOverlap ;
			struct _overlap &o = extendedOverlaps[fragments[i].a] ;

			fragmentOverlap.matchCnt = o.matchCnt ;
			fragmentOverlap.similarity = o.similarity ;
			fragmentOverlap.seqIdx = o.seqIdx ;
			fragmentOverlap.seqStart = o.seqStart ;
			fragmentOverlap.seqEnd = o.seqEnd ;
			fragmentOverlap.hasMatePair = false ;
			fragmentOverlap.overlap1 = o ;
			
			if (fragments[i].b >= 0)
			{
				struct _overlap &o2 = extendedOverlaps2[fragments[i].b] ;
				fragmentOverlap.matchCnt += o2.matchCnt ;
				if (o.strand == 1)
					fragmentOverlap.seqEnd = o2.seqEnd ;
				else
					fragmentOverlap.seqStart = o2.seqStart ;

				// TODO: there might be other better ways to combine the mate pairs.
				fragmentOverlap.similarity = (double)fragmentOverlap.matchCnt / 
					(o.readEnd - o.readStart + 1 + o2.readEnd - o2.readStart + 1 + 
					 o.seqEnd - o.seqStart + 1 + o2.seqEnd - o2.seqStart + 1 ) ;
			
				fragmentOverlap.hasMatePair = true ;
				fragmentOverlap.overlap2 = o2 ;
			}	
			
			if (seqIdxToOverlapIdx.find(fragmentOverlap.seqIdx) != seqIdxToOverlapIdx.end())
			{
				// Note < here is for ranking, so smaller has higher rank
				if (fragmentOverlap < assign[seqIdxToOverlapIdx[fragmentOverlap.seqIdx]])
				{
					assign[seqIdxToOverlapIdx[fragmentOverlap.seqIdx]] = fragmentOverlap ; 
				}	
			}
			else
			{
				assign.push_back( fragmentOverlap ) ;
				seqIdxToOverlapIdx[fragmentOverlap.seqIdx] = assign.size() - 1;
			}	
		}

		std::sort(assign.begin(), assign.end()) ;
		int assignCnt = assign.size() ;
		if (assignCnt > 0 && assign[0].similarity < refSeqSimilarity)
		{
			assign.clear() ;
			assignCnt = 0 ;
		}
		for (i = 1 ; i < assignCnt ; ++i)
		{
			if (assign[i].matchCnt < assign[0].matchCnt
					|| assign[i].similarity < assign[0].similarity ) // TODO: maybe allow more difference
			{
				assign.resize(i) ;
				break ;
			}	
		}

		// Check whether there is better alignment but mate could not be aligned due to truncated reference gene (e.g. UTR).
		if (assign.size() > 0 && read2 != NULL)
		{
			struct _overlap representative ;
			representative = assign[0].overlap1 ;
			bool filter = false ;
			for (i = 0 ; i < overlapCnt && !filter ; ++i)
			{
				if (extendedOverlaps[i].matchCnt > representative.matchCnt
						|| (extendedOverlaps[i].matchCnt == representative.matchCnt
							&& extendedOverlaps[i].similarity > representative.similarity) 
						&& seqIdxToOverlapIdx.find(extendedOverlaps[i].seqIdx) == seqIdxToOverlapIdx.end())
				{
					if (TruncatedMatePairOverlap(extendedOverlaps[i], assign[0].overlap1, assign[0].overlap2))
					{
						filter = true ;
					}
				}
			}
			// Better read 2
			representative = assign[0].overlap2 ;
			for (i = 0 ; i < overlapCnt2 && !filter ; ++i)
			{
				if (extendedOverlaps2[i].matchCnt > representative.matchCnt
						|| (extendedOverlaps2[i].matchCnt == representative.matchCnt
							&& extendedOverlaps2[i].similarity > representative.similarity) 
						&& seqIdxToOverlapIdx.find(extendedOverlaps2[i].seqIdx) == seqIdxToOverlapIdx.end())
				{
					if (TruncatedMatePairOverlap(extendedOverlaps2[i], assign[0].overlap2, assign[0].overlap1))
					{
						filter = true ;
					}
				}
			}

			if (filter)
				assign.clear() ;
		}


		delete[] rc ;
		delete[] align ;
		if (read2 != NULL)
		{
			delete[] rc2 ;
		}
		return assign.size() ;
	}

	int GetSeqNameToIdxMap(std::map<std::string, int>& nameToIdx)
	{
		int seqCnt = seqs.size() ;
		int i ;
		for (i = 0 ; i < seqCnt ; ++i)
		{
			std::string s(seqs[i].name) ;
			nameToIdx[s] = i ;
		}
		return seqCnt ;
	}
}	;


#endif

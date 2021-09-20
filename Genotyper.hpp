#ifndef _MOURISL_GENOTYPER
#define _MOURISL_GENOTYPER

#include <math.h>

#include "SeqSet.hpp"

#include "defs.h"
#include "KmerCount.hpp"
#include "SimpleVector.hpp"

#define GENETYPE_KIR 0
#define GENETYPE_HLA 1

struct _alleleInfo
{
	int majorAlleleIdx ; 
	int geneIdx ; // which gene this allele belongs to

	int	alleleRank ; // -1 not select, 0-first allele, 1-second alelle
	int genotypeQuality ; // assignment quality.
	double abundance ;

	int equivalentClass ; // the class id for the alleles with the same set of read alignment
	double ecAbundance ; // the abundance for the equivalent class.
} ;

struct _readGroupInfo
{
	double count ; // number of reads this group contains.
} ;

struct _readAssignment
{
	int alleleIdx ;
	int start, end ;
	double weight ; 
	double qual ;

	bool operator<(const struct _readAssignment &b) const
	{
		return alleleIdx < b.alleleIdx ;
	}
} ;

class Genotyper
{
private:
	void ParseAlleleName(char *allele, char *gene, char *majorAllele)
	{
		int i, j ;
		int parseType = 1 ; // 1-first three; 2-colon(HLA)
		strcpy(gene, allele) ;
		strcpy(majorAllele, allele) ;
		for (i = 0 ; allele[i] ; ++i)
			if (allele[i] == ':')
				parseType = 2 ;

		if (parseType == 1)
		{
			for (i = 0 ; allele[i] ; ++i)
			{
				if (allele[i] == '*')
					break ;
			}
			gene[i] = '\0' ;
			for (j = 0 ; j <= 3 && allele[i + j] ; ++j)
				;
			majorAllele[i + j] = '\0' ;
		}
		else if (parseType == 2)
		{
			for (i = 0 ; allele[i] ; ++i)
			{
				if (allele[i] == '*')
					break ;
			}
			gene[i] = '\0' ;
			int k = 0 ;
			for (j = i ; allele[j] ; ++j)
			{
				if (allele[j] == ':')	
				{
					++k ;
					if (k >= 3)
						break ;
				}
			}
			majorAllele[j] = '\0' ;
		}
	}
	
	static bool CompSortPairByBDec( const struct _pair &p1, const struct _pair &p2 )
	{
		if (p1.b != p2.b)
			return p2.b < p1.b ;
		return p1.a < p2.a ;
	}

	static bool CompSortPairIntDoubleBDec( const struct _pairIntDouble &p1, const struct _pairIntDouble &p2 )
	{
		if (p2.b != p1.b)
			return p2.b < p1.b ;
		return p1.a < p2.a ;
	}

	static bool CompSortDoubleDec(const double &a, const double &b)
	{
		return b < a ;
	}
	
	bool IsAssignedReadTheSame(const std::vector<struct _pair> &l1, const std::vector<struct _pair> &l2)
	{
		int cnt1 = l1.size() ;
		int cnt2 = l2.size() ;
		int i ;
		if (cnt1 != cnt2) 
			return false ;
		// The read id in each vector should be sorted
		for (i = 0 ; i < cnt1 ; ++i)
		{
			if (l1[i].a != l2[i].a 
					|| readAssignments[l1[i].a][l1[i].b].qual != readAssignments[l2[i].a][l2[i].b].qual)
				return false ;
		}
		return true ;
	}

	bool IsReadAssignmentTheSame(const std::vector<struct _readAssignment> &a1, const std::vector<struct _readAssignment> &a2)
	{
		int cnt1 = a1.size() ;
		int cnt2 = a2.size() ;
		int i ;
		if (cnt1 != cnt2) 
			return false ;
		// The read id in each vector should be sorted
		for (i = 0 ; i < cnt1 ; ++i)
		{
			if (a1[i].alleleIdx != a2[i].alleleIdx
					|| a1[i].qual != a2[i].qual)
				return false ;
		}
		return true ;
	}

	bool IsReadsInAlleleIdxOptimal(const std::vector<struct _pair> &readsInAllele, int k)
	{
		if (readAssignments[readsInAllele[k].a][readsInAllele[k].b].qual == 1)
			return true ;
		return false ;
	}

	double ReadAssignmentWeight(const struct _fragmentOverlap &o)
	{
		double ret = 1 ;
		
		double similarity = o.similarity ; //o.overlap1.similarity ;
		//if (o.hasMatePair && o.overlap2.similarity < similarity)
		//	similarity = o.overlap2.similarity ;
		
		if (similarity < 0.85)
			ret = 0.01 ;
		else if (similarity < 0.9)
			ret = 0.1 ;
		else if (similarity < 0.95)
			ret = 0.5 ;
		//else if (similarity < 1)
		//	ret = 0.5 ;
		
		return ret ;
	}

	int Rand()
	{
		return randomSeed = (48271 * randomSeed) & 0x7fffffff ;
	}

	double alnorm ( double x, bool upper )

		//****************************************************************************80
		//
		//  Purpose:
		//
		//    ALNORM computes the cumulative density of the standard normal distribution.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license. 
		//
		//  Modified:
		//
		//    17 January 2008
		//
		//  Author:
		//
		//    Original FORTRAN77 version by David Hill.
		//    C++ version by John Burkardt.
		//
		//  Reference:
		//
		//    David Hill,
		//    Algorithm AS 66:
		//    The Normal Integral,
		//    Applied Statistics,
		//    Volume 22, Number 3, 1973, pages 424-427.
		//
		//  Parameters:
		//
		//    Input, double X, is one endpoint of the semi-infinite interval
		//    over which the integration takes place.
		//
		//    Input, bool UPPER, determines whether the upper or lower
		//    interval is to be integrated:
		//    .TRUE.  => integrate from X to + Infinity;
		//    .FALSE. => integrate from - Infinity to X.
		//
		//    Output, double ALNORM, the integral of the standard normal
		//    distribution over the desired interval.
		//
	{
		double a1 = 5.75885480458;
		double a2 = 2.62433121679;
		double a3 = 5.92885724438;
		double b1 = -29.8213557807;
		double b2 = 48.6959930692;
		double c1 = -0.000000038052;
		double c2 = 0.000398064794;
		double c3 = -0.151679116635;
		double c4 = 4.8385912808;
		double c5 = 0.742380924027;
		double c6 = 3.99019417011;
		double con = 1.28;
		double d1 = 1.00000615302;
		double d2 = 1.98615381364;
		double d3 = 5.29330324926;
		double d4 = -15.1508972451;
		double d5 = 30.789933034;
		double ltone = 7.0;
		double p = 0.398942280444;
		double q = 0.39990348504;
		double r = 0.398942280385;
		bool up;
		double utzero = 18.66;
		double value;
		double y;
		double z;

		up = upper;
		z = x;

		if ( z < 0.0 )
		{
			up = !up;
			z = - z;
		}

		if ( ltone < z && ( ( !up ) || utzero < z ) )
		{
			if ( up )
			{
				value = 0.0;
			}
			else
			{
				value = 1.0;
			}
			return value;
		}

		y = 0.5 * z * z;

		if ( z <= con )
		{
			value = 0.5 - z * ( p - q * y 
					/ ( y + a1 + b1 
						/ ( y + a2 + b2 
							/ ( y + a3 ))));
		}
		else
		{
			value = r * exp ( - y ) 
				/ ( z + c1 + d1 
						/ ( z + c2 + d2 
							/ ( z + c3 + d3 
								/ ( z + c4 + d4 
									/ ( z + c5 + d5 
										/ ( z + c6 ))))));
		}

		if ( !up )
		{
			value = 1.0 - value;
		}

		return value;
	}

	double EMupdate(double *ecAbundance0, double *ecAbundance1, double *ecReadCount, const std::vector<std::vector<struct _pairID> > &readGroupToAlleleEc, const SimpleVector<struct _readGroupInfo> readGroupInfo, const double *ecLength)
	{
		int ecCnt = equivalentClassToAlleles.size() ;
		int rgCnt = readGroupToAlleleEc.size() ;
		int i, j ;
		// E-step: find the expected number of reads
		memset(ecReadCount, 0, sizeof(double) * ecCnt) ;
		for (i = 0 ; i < rgCnt ; ++i)
		{
			double psum	= 0 ;
			int size = readGroupToAlleleEc[i].size() ;
			for (j = 0 ; j < size ; ++j)
			{
				int ecIdx = readGroupToAlleleEc[i][j].a ;
				double qual = readGroupToAlleleEc[i][j].b ;
				psum += ecAbundance0[ecIdx] / ecLength[ecIdx] * qual ;
			}
			if (psum == 0)	
				psum = 1 ;
			for (j = 0 ; j < size ; ++j)
			{
				int ecIdx = readGroupToAlleleEc[i][j].a ;
				double qual = readGroupToAlleleEc[i][j].b ;
				ecReadCount[ecIdx] += readGroupInfo[i].count * (ecAbundance0[ecIdx] * qual / ecLength[ecIdx] / psum)  ;
			}
		}

		// M-step: recompute the abundance
		double diffSum = 0 ;
		double normalization = 0 ;
		for (i = 0 ; i < ecCnt ; ++i)
			normalization += ecReadCount[i] ; // ecLength[i] ;

		for (i = 0 ; i < ecCnt ; ++i)
		{
			double tmp = ecReadCount[i] / normalization ;
			//printf("%d %s %d: %lf %lf %lf. %lf\n", i, refSet.GetSeqName(equivalentClassToAlleles[i][0]), equivalentClassToAlleles[i].size(), tmp, ecReadCount[i], ecLength[i], ecAbundance[i]) ;

			diffSum += ABS(tmp - ecAbundance0[i]) ;
			ecAbundance1[i] = tmp ;
		}
		return diffSum ;
	}

	// Compute the update coeffcient alpha in the SQUREEM paper
	double SQUAREMalpha(double *t0, double *t1, double *t2, int n)
	{
		int i ;
		double sqrSumR = 0 ;
		double sqrSumV = 0 ;
		for (i = 0 ; i < n ; ++i)
		{
			sqrSumR += (t1[i] - t0[i]) * (t1[i] - t0[i]) ;
			sqrSumV += (t2[i] - 2 * t1[i] + t0[i]) * (t2[i] - 2 * t1[i] + t0[i]) ;
		}
		return -sqrt(sqrSumR) / sqrt(sqrSumV) ;
	}


	int readCnt ;
	int totalReadCnt ;
	int maxAssignCnt ;
	std::vector< std::vector<struct _pair> >	readsInAllele ; // a-read idx. b-the index withint the readAssignment[a]
	std::vector< std::vector<struct _readAssignment> > readAssignments ; // Coalesce reads assigned to the same alleles 
	std::map<int, std::vector<int> > readAssignmentsFingerprintToIdx ;
	std::vector< std::vector<struct _readAssignment> > allReadAssignments ;
	std::vector< std::vector<int> > equivalentClassToAlleles ;
	std::vector< std::vector<struct _pair> >	selectedAlleles ; // a-allele name, b-which allele (0,1)
	
	// variables for allele, majorAllele and genes	
	char *geneBuffer ;
	char *majorAlleleBuffer ;	

	SimpleVector<struct _alleleInfo> alleleInfo ;
	std::map<std::string, int> majorAlleleNameToIdx ;
	std::map<std::string, int> geneNameToIdx ;
	std::vector<int> majorAlleleSize ; // number of alleles in the major allele
	std::vector<std::string> geneIdxToName ;
	std::vector<std::string> majorAlleleIdxToName ;
	int geneCnt ;
	int majorAlleleCnt ;
	int alleleCnt ;

	// variables for abundance
	SimpleVector<double> geneAbundance ;
	SimpleVector<double> majorAlleleAbundance ;
	SimpleVector<double> geneMaxMajorAlleleAbundance ;

	int64_t randomSeed ;

	int geneType ;  // not used actually

	// variables for filter
	double filterFrac ;
	double filterCov ;
	double crossGeneRate ;
	double **geneSimilarity ;

public:
	SeqSet refSet ;
	
	Genotyper(int kmerLength):refSet(kmerLength) 
	{
		geneBuffer = new char[256] ;
		majorAlleleBuffer = new char[256] ;
		alleleCnt = majorAlleleCnt = geneCnt = readCnt = 0 ;
		maxAssignCnt = 2000 ;
		randomSeed = 17 ;
		geneType = GENETYPE_KIR ;

		filterFrac = 0.15 ;
		filterCov = 1.0 ;
		crossGeneRate = 0.005 ;

		geneSimilarity = NULL ;
	}

	~Genotyper() 
	{
		delete[] geneBuffer ;
		delete[] majorAlleleBuffer ;

		int i ;
		for (i = 0 ; i < geneCnt ; ++i)
			delete[] geneSimilarity[i] ;
		delete[] geneSimilarity ;
	}

	void SetGeneType(int g)
	{
		geneType = g ;
	}

	void SetFilterFrac(double f)
	{
		filterFrac = f ;
	}

	void SetFilterCov(double c)
	{
		filterCov = c ;
	}

	void SetCrossGeneRate(double r)
	{
		crossGeneRate = r ;
	}

	void InitRefSet(char *filename)
	{
		int i, j ;

		//refSet.InputRefFa(filename) ;
		std::map<std::string, int> usedSeq ;
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		while ( fa.Next() )
		{
			std::string seq( fa.seq );
			if (usedSeq.find(seq) != usedSeq.end())
			{
				refSet.UpdateSeqWeight( usedSeq[seq], 1) ;
			}
			else 
			{
				usedSeq[seq] = refSet.InputRefSeq(fa.id, fa.seq, 1);
			}
		}

		alleleCnt = refSet.Size() ;
		alleleInfo.ExpandTo(alleleCnt) ;
		
		for ( i = 0 ; i < alleleCnt ; ++i )	
		{
			char *allele = refSet.GetSeqName(i) ;
			ParseAlleleName( allele, geneBuffer, majorAlleleBuffer ) ;

			std::string sGene(geneBuffer) ;
			std::string sMajorAllele(majorAlleleBuffer) ;
			if (geneNameToIdx.find(sGene) == geneNameToIdx.end())
			{
				geneNameToIdx[sGene] = geneCnt ;
				geneIdxToName.push_back(sGene) ;
				++geneCnt ;
			}
			if (majorAlleleNameToIdx.find(sMajorAllele) == majorAlleleNameToIdx.end())
			{
				majorAlleleNameToIdx[sMajorAllele] = majorAlleleCnt ;
				majorAlleleIdxToName.push_back(sMajorAllele) ;
				majorAlleleSize.push_back(0);
				++majorAlleleCnt ;
			}

			
			alleleInfo[i].abundance = 0 ;
			alleleInfo[i].geneIdx = geneNameToIdx[sGene] ;
			alleleInfo[i].majorAlleleIdx = majorAlleleNameToIdx[sMajorAllele] ;
			alleleInfo[i].alleleRank = -1 ;
			alleleInfo[i].abundance = 0 ;
			alleleInfo[i].genotypeQuality = -1 ;
			majorAlleleSize[ alleleInfo[i].majorAlleleIdx ] += refSet.GetSeqWeight(i) ;
		}
		
		geneSimilarity = new double*[geneCnt] ;
		KmerCount *kmerProfiles = new KmerCount[geneCnt];
		for (i = 0 ; i < geneCnt ; ++i)
		{
			geneSimilarity[i] = new double[geneCnt] ;
			// Select the lexcographically smallest
			int minTag = -1 ;
			for (j = 0 ; j < alleleCnt ; ++j)
			{
				if (alleleInfo[j].geneIdx != i)
					continue ;
				if (minTag == -1 || strcmp(refSet.GetSeqConsensus(j), refSet.GetSeqConsensus(minTag) ) < 0)
					minTag = j ;
			}

			kmerProfiles[i].AddCount(refSet.GetSeqConsensus(minTag)) ;
		}

		for (i = 0 ; i < geneCnt ; ++i)
		{
			geneSimilarity[i][i] = 1.0 ;
			for (j = i + 1 ; j < geneCnt ; ++j)
			{
				geneSimilarity[i][j] = geneSimilarity[j][i] = kmerProfiles[i].GetCountSimilarity( kmerProfiles[j] ) ;
				//printf("%d %d %lf\n", i, j, geneSimilarity[i][j]) ;
			}
		}
		delete[] kmerProfiles ;
	}

	void InitReadAssignments(int totalReadCnt, int maxAssignCnt)
	{
		maxAssignCnt = maxAssignCnt ;
		readCnt = 0 ;

		allReadAssignments.resize(totalReadCnt) ;
		readAssignments.clear() ;
		readAssignmentsFingerprintToIdx.clear() ;
		readsInAllele.resize(alleleCnt) ;

		int i ;
		for (i = 0 ; i < totalReadCnt ; ++i)
			allReadAssignments[i].clear() ;
		for (i = 0 ; i < alleleCnt ; ++i)
			readsInAllele[i].clear() ;

		this->totalReadCnt = totalReadCnt ;
	}

	void SetReadAssignments(int readId, const std::vector<struct _fragmentOverlap> &assignment)
	{
		int i ;
		int assignmentCnt = assignment.size() ;
		allReadAssignments[readId].clear() ;
		if (maxAssignCnt > 0 && assignmentCnt > maxAssignCnt)
			return ;

		/*for (i = 1 ; i < assignmentCnt ; ++i)
			if ( alleleInfo[assignment[i].seqIdx].geneIdx != alleleInfo[assignment[i - 1].seqIdx].geneIdx)
				return ;*/
		double weightFactor = 1.0 ;
		/*for (i = 1 ; i < assignmentCnt ; ++i)
			if (assignment[i].overlap1.matchCnt != assignment[i - 1].overlap1.matchCnt ||
					assignment[i].overlap2.matchCnt != assignment[i - 1].overlap2.matchCnt)
			{
				weightFactor = 0.1 ;
				return;
				break ;
			}*/
		for (i = 0 ; i < assignmentCnt ; ++i)
		{
			if (refSet.IsFragmentSpanSeparator(assignment[i]))
				return;
		}

		for (i = 0; i < assignmentCnt; ++i)
		{
			struct _readAssignment na ;
			na.alleleIdx = assignment[i].seqIdx ;
			na.start = assignment[i].seqStart ;
			na.end = assignment[i].seqEnd ;
			na.weight = ReadAssignmentWeight(assignment[i]);
			na.qual = assignment[i].qual ;
			allReadAssignments[readId].push_back(na) ;
		}
	}
		
	// Coalesce the [begin,end] all reads to the read assignment
	int CoalesceReadAssignments(int begin, int end)
	{
		int i, j, k ;
		int ret = 0 ;
		for (i = begin ; i <= end && i < totalReadCnt ; ++i)
		{
			const int FINGERPRINT_MAX = 20000003 ;
			int size = allReadAssignments[i].size() ;
			if (size == 0)
				continue ;
			++ret ;
			//std::sort(allReadAssignments[i].begin(), allReadAssignments[i].end()) ;
			int fingerprint = 0 ;
			for (j = 0 ; j < size ; ++j)
			{
				k = allReadAssignments[i][j].alleleIdx ;
				fingerprint = (fingerprint * (int64_t)alleleCnt + k) % FINGERPRINT_MAX ;
			}
			
			int addTo = -1 ;
			if (readAssignmentsFingerprintToIdx.find(fingerprint) == readAssignmentsFingerprintToIdx.end())
				addTo = -1 ;
			else
			{
				std::vector<int> assignmentsIdx = readAssignmentsFingerprintToIdx[fingerprint] ;
				int idxSize = assignmentsIdx.size() ;
				addTo = -1 ;
				for (j = 0 ; j < idxSize ; ++j)
				{
					if (IsReadAssignmentTheSame(allReadAssignments[i], 
								readAssignments[ assignmentsIdx[j] ] ) )
					{
						addTo = assignmentsIdx[j] ;
						break ;
					}
				}
			}

			if (addTo == -1)
			{
				readAssignments.push_back( allReadAssignments[i] ) ;
				readAssignmentsFingerprintToIdx[fingerprint].push_back(readCnt) ;
				++readCnt ;				
			}
			else
			{
				for (j = 0 ; j < size ; ++j)
				{
					if (allReadAssignments[i][j].qual == 1)
					{
						if (allReadAssignments[i][j].start < readAssignments[addTo][j].start)
							readAssignments[addTo][j].start = allReadAssignments[i][j].start ;
						if (allReadAssignments[i][j].end < readAssignments[addTo][j].end)
							readAssignments[addTo][j].end = allReadAssignments[i][j].start ;
					}
					readAssignments[addTo][j].weight += allReadAssignments[i][j].weight ;
					// The read assignment the same test makes sure they have the same quality
				}
			}
		}
		// Release the memory space for allReadAssignments
		for (i = begin ; i <= end && i < totalReadCnt ; ++i)
		{
			std::vector<struct _readAssignment>().swap(allReadAssignments[i]) ;
		}
		return ret ;
	}
	
	// Build the read in allele list
	// Ret: the number of read got assigned
	int FinalizeReadAssignments()
	{
		int i, j ;
		int ret = 0 ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			int assignmentCnt = readAssignments[i].size() ;
			std::sort(readAssignments[i].begin(), readAssignments[i].end()) ;
			if (assignmentCnt > 0)
				++ret ;
			for (j = 0; j < assignmentCnt; ++j)
			{
				struct _pair np ;
				np.a = i ;
				np.b = j ;
				readsInAllele[readAssignments[i][j].alleleIdx].push_back(np) ;
			}
		}
		
		BuildAlleleEquivalentClass() ;
		return ret ;
	}

	double GetAverageReadAssignmentCnt()
	{
		int i ;
		double sum = 0 ;
		double cnt = 0 ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			if (readAssignments[i].size() > 0)
			{
				sum += readAssignments[i].size() ;
				++cnt ;
			}
		}
		return sum / cnt ;
	}

	void SetAlleleAbundance(double *ecReadCount, double *ecLength)
	{
		int i, j, k ;
		int ecCnt = equivalentClassToAlleles.size();
		if (ecReadCount != NULL)
		{
			for (i = 0 ; i < alleleCnt ; ++i)
				alleleInfo[i].abundance = alleleInfo[i].ecAbundance = 0 ;

			for (i = 0 ; i < ecCnt ; ++i)
			{
				int size = equivalentClassToAlleles[i].size() ;
				double abund = 0 ; 
				//k = equivalentClassToAlleles[i][0] ;
				//printf("%d %d %s %lf %d %d\n", i, k, refSet.GetSeqName(k), emEcReadCount[j][i], refSet.GetSeqConsensusLen(k),readsInAllele[k].size()) ;
				abund += ecReadCount[i] ;
				//printf("%lf\n", abund) ;
				abund = abund / ecLength[i] * 1000.0 ; // FPK
				for (j = 0 ; j < size ; ++j)
				{
					k = equivalentClassToAlleles[i][j] ;
					//alleleInfo[k].abundance = ecAbundance[i] / size * effectiveReadCnt ;
					//alleleInfo[k].ecAbundance = ecAbundance[i] * effectiveReadCnt ;
					alleleInfo[k].abundance = abund / size ;
					alleleInfo[k].ecAbundance = abund ;
					//printf("%d %d %s %lf %d %d\n", i, k, refSet.GetSeqName(k), (double)abund, refSet.GetSeqConsensusLen(k),readsInAllele[k].size()) ;
				}
				//printf("%lf %lf\n", ecAbundance[i], ecAbundance[i] * effectiveReadCnt) ;
			}
		}


		// Set major allele and gene abundances
		// Init other useful abundance data
		geneAbundance.ExpandTo(geneCnt) ;
		geneAbundance.SetZero(0, geneCnt) ;
		majorAlleleAbundance.ExpandTo(majorAlleleCnt) ;
		majorAlleleAbundance.SetZero(0, majorAlleleCnt) ;
		geneMaxMajorAlleleAbundance.ExpandTo(geneCnt) ;
		geneMaxMajorAlleleAbundance.SetZero(0, geneCnt) ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			majorAlleleAbundance[ alleleInfo[i].majorAlleleIdx ] += alleleInfo[i].abundance ;
			geneAbundance[ alleleInfo[i].geneIdx ] += alleleInfo[i].abundance ;

			/*if (alleleInfo[i].abundance > geneMaxMajorAlleleAbundance[ alleleInfo[i].geneIdx ])
				{
				geneMaxMajorAlleleAbundance[ alleleInfo[i].geneIdx ] = alleleInfo[i].abundance ;
				}*/
		}

		for (i = 0 ; i < alleleCnt ; ++i)
		{
			double abund = majorAlleleAbundance[ alleleInfo[i].majorAlleleIdx ] ;
			if (abund > geneMaxMajorAlleleAbundance[ alleleInfo[i].geneIdx ])
				geneMaxMajorAlleleAbundance[ alleleInfo[i].geneIdx ] = abund ;
		}
	}

	void InitAlleleAbundance(FILE *fp)
	{
		int i, j ;	
		char buffer[256] ;
		double abundance ;
		double count ;
		int tmp ;

		std::map<std::string, int> refNameToIdx ;
		refSet.GetSeqNameToIdxMap(refNameToIdx) ;

		fscanf(fp, "%s %s %s %s %s", buffer, buffer, buffer, buffer, buffer) ; // header

		while (fscanf( fp, "%s %d %d %lf %lf", buffer, &tmp, &tmp, &count, &abundance) != EOF)
		{
			std::string s(buffer)	;
			struct _pairIntDouble np ;
			np.a = refNameToIdx[buffer] ;
			np.b = count ;//abundance ;
			alleleInfo[np.a].abundance = np.b ;
		}
		fclose(fp) ;
		
		int ecCnt = equivalentClassToAlleles.size() ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			int size = equivalentClassToAlleles[i].size() ;
			double total = 0 ;
			for (j = 0 ; j < size ; ++j)
				total += alleleInfo[equivalentClassToAlleles[i][j]].abundance ;
			for (j = 0 ; j < size ; ++j)
				alleleInfo[equivalentClassToAlleles[i][j]].ecAbundance = total ;
		}
		
		SetAlleleAbundance(NULL, NULL) ;
	}
	
	int GetGeneAlleleTypes(int geneIdx)
	{
		if ( selectedAlleles[geneIdx].size() == 0 )
			return  0 ;
		else
		{
			//return selectedAlleles[geneIdx].back().b + 1;
			int size = selectedAlleles[geneIdx].size() ;
			int i, ret = 0 ;
			for (i = 0 ; i < size ; ++i)
			{
				if (selectedAlleles[geneIdx][i].b > ret)
					ret = selectedAlleles[geneIdx][i].b ;
			}
			return ret + 1 ;
		}
	}

	// Return: the number of allele equivalent class
	int BuildAlleleEquivalentClass()
	{
		SimpleVector<struct _pair> alleleFingerprint ; // a-alleleId, b-finger print from the read it covered
		int i, j ;
		const int FINGERPRINT_MAX = 1000003 ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			struct _pair np ;
			np.a = i ;
			np.b = -1 ;
			int size = readsInAllele[i].size() ;
			alleleInfo[i].equivalentClass = -1 ;
			if (readsInAllele[i].size() > 0)
			{
				np.b = 0 ;
				for ( j = 0 ; j < size ; ++j)
				{
					np.b = ((uint32_t)np.b * readCnt + readsInAllele[i][j].a) % FINGERPRINT_MAX ;
				}
			}	
			alleleFingerprint.PushBack(np) ;
		}

		std::sort(alleleFingerprint.BeginAddress(), alleleFingerprint.EndAddress(), CompSortPairByBDec ) ;
		
		int ecCnt = 0 ;
		equivalentClassToAlleles.clear() ;

		if (alleleCnt == 0 || alleleFingerprint[0].b == -1)
			return 0 ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			bool newEc = true ;
			if (alleleFingerprint[i].b == -1)
				break ;

			for (j = i - 1 ; j >= 0 ; --j)
			{
				if (alleleFingerprint[i].b != alleleFingerprint[j].b)
					break ;
				if (IsAssignedReadTheSame(readsInAllele[ alleleFingerprint[i].a ], 
							readsInAllele[ alleleFingerprint[j].a]))
						//&& refSet.GetSeqConsensusLen(alleleFingerprint[i].a) == refSet.GetSeqConsensusLen(alleleFingerprint[j].a))
				{
					newEc = false ;
					break ;
				}
			}
			int alleleIdx = alleleFingerprint[i].a ;	
			if (newEc)
			{
				equivalentClassToAlleles.push_back( std::vector<int>() ) ;	
				equivalentClassToAlleles.back().push_back( alleleIdx ) ;
				alleleInfo[ alleleIdx ].equivalentClass = ecCnt ;
				++ecCnt ;
			}
			else
			{
				int ecIdx = alleleInfo[alleleFingerprint[j].a].equivalentClass ;
				equivalentClassToAlleles[ecIdx].push_back(alleleIdx) ;
				alleleInfo[alleleIdx].equivalentClass = ecIdx ;
			}

		}
		RemoveLowMAPQAlleleInEquivalentClass() ;

		return ecCnt ;
	}
	
	// return : EM algorithm iterations
	int QuantifyAlleleEquivalentClass()
	{
		int i, j, k ;
		int t ; // iteration for EM algorithm.
		int ecCnt = equivalentClassToAlleles.size() ;	
		int ret = 0 ;	
		
		// Convert readassignment_to_allele to readassignment_to_alleleEquivalentClass
		std::vector< std::vector<struct _pairID> > readGroupToAlleleEc ;
		SimpleVector<struct _readGroupInfo> readGroupInfo ; // the read represent the read group, so we can easily obtain the assignment information.
		readGroupToAlleleEc.resize(readCnt) ;
		readGroupInfo.ExpandTo(readCnt) ;
		std::map<int, int> ecUsed ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			double count = 0 ;
			int size = readAssignments[i].size() ;
			count = readAssignments[i][0].weight ;
			for (j = 1 ; j < size ; ++j)
				if (readAssignments[i][j].weight > count)
					count = readAssignments[i][j].weight ;
			readGroupInfo[i].count = count ;
		}
		for (i = 0 ; i < readCnt ; ++i)
		{
			int size = readAssignments[i].size() ;
			ecUsed.clear() ;
			for (j = 0 ; j < size ; ++j)
			{
				int ecIdx = alleleInfo[readAssignments[i][j].alleleIdx].equivalentClass ;
				if (ecUsed.find(ecIdx) == ecUsed.end())
				{
					struct _pairID np ;
					np.a = ecIdx ;
					np.b = readAssignments[i][j].qual ;
					ecUsed[ecIdx] = readGroupToAlleleEc[i].size() ;
					readGroupToAlleleEc[i].push_back(np) ;
				}
				else
				{
					// Should not happen though, the equivalent class
					// makes sure the quality score is the same.
					int k = ecUsed[ecIdx] ;
					if (readAssignments[i][j].qual > readGroupToAlleleEc[i][k].b)
						readGroupToAlleleEc[i][k].b = readAssignments[i][j].qual ;
				}
			}
		}

		// Start the EM algorithm
		double *emResults ;
		const int maxEMIterations = 1000 ;
		double *ecAbundance0 = new double[ecCnt] ;
		double *ecAbundance1 = new double[ecCnt] ;
		double *ecAbundance2 = new double[ecCnt] ;
		double *ecAbundance3 = new double[ecCnt] ;
		double *ecReadCount = new double[ecCnt] ;
		double *ecLength = new double[ecCnt] ; // the sequence length for equivalent class.

		for (i = 0 ; i < ecCnt ; ++i)
		{
			int size = equivalentClassToAlleles[i].size() ;
			ecLength[i] = refSet.GetSeqEffectiveLen(equivalentClassToAlleles[i][0]) ;
			for (j = 1 ; j < size ; ++j)
			{
				int len = refSet.GetSeqEffectiveLen(equivalentClassToAlleles[i][j]) ;
				if (len < ecLength[i])
					ecLength[i] = len ;
			}
		}
		
		ret = 0 ;
		const int maskRound = 10 ; // Mask low abundant alleles every 10 round
		for (i = 0 ; i < ecCnt ; ++i)
		{
			//ecAbundance[i] = Rand()%7 + 1 ; //1.0 / ecCnt ;
			//ecAbundance[i] = equivalentClassToAlleles[i].size() + (Rand()%3 - 1) * 0.5 ; //1.0 / ecCnt ;
			//ecAbundance0[i] = equivalentClassToAlleles[i].size() ;
			ecAbundance0[i] = 0;
			for (j = 0 ; j < (int)equivalentClassToAlleles[i].size() ; ++j)
				ecAbundance0[i] += refSet.GetSeqWeight(equivalentClassToAlleles[i][j]);
			/*ecAbundance0[i] = 0 ;
				for (j = 0 ; j < (int)equivalentClassToAlleles[i].size() ; ++j)
				ecAbundance0[i] += majorAlleleSize[alleleInfo[equivalentClassToAlleles[i][j]].majorAlleleIdx] ;*/
		}

		for (t = 0 ; t < maxEMIterations ; ++t)
		{
			++ret ;
			double diffSum = EMupdate(ecAbundance0, ecAbundance1, ecReadCount,
					readGroupToAlleleEc, readGroupInfo, ecLength) ;
			diffSum = EMupdate(ecAbundance1, ecAbundance2, ecReadCount,
					readGroupToAlleleEc, readGroupInfo, ecLength) ;

			double alpha = SQUAREMalpha(ecAbundance0, ecAbundance1, ecAbundance2, ecCnt) ;
			//memcpy(ecAbundance0, ecAbundance1, sizeof(double) * ecCnt) ;
			for (i = 0 ; i < ecCnt ; ++i)
			{
				ecAbundance3[i] = ecAbundance0[i] 
					- 2 * alpha * (ecAbundance1[i] - ecAbundance0[i])
					+ alpha * alpha * (ecAbundance2[i] - 2 * ecAbundance1[i] + ecAbundance0[i]) ;
				if (ecAbundance3[i] < 0)
					ecAbundance3[i] = 0 ;
			}
			EMupdate(ecAbundance3, ecAbundance1, ecReadCount,
					readGroupToAlleleEc, readGroupInfo, ecLength	) ;

			diffSum = 0 ;
			for (i = 0 ; i < ecCnt ; ++i)
			{
				diffSum += ABS(ecAbundance1[i] - ecAbundance0[i]) ;
				ecAbundance0[i] = ecAbundance1[i] ;
				//printf("%d %s %d: %lf %lf. %lf\n", i, refSet.GetSeqName(equivalentClassToAlleles[i][0]), equivalentClassToAlleles[i].size(), ecReadCount[i], ecLength[i], ecAbundance0[i]) ;
			}

			if (diffSum < 1e-5 && t < maxEMIterations - 2)
				t = maxEMIterations - 2 ; // Force one more iteration
			
			if (t > 0 && t % maskRound == 0)
			{
				// Filter the low abundant ones
				SetAlleleAbundance(ecReadCount, ecLength) ;

				for (i = 0 ; i < alleleCnt ; ++i)			
				{
					if (majorAlleleAbundance[alleleInfo[i].majorAlleleIdx] 
							< filterFrac * 0.5 * geneMaxMajorAlleleAbundance[alleleInfo[i].geneIdx])
					{
						alleleInfo[i].abundance = 0 ;
						alleleInfo[i].ecAbundance = 0 ;
					}
				}

				// Reset ecAbundance0
				for (i = 0 ; i < ecCnt ; ++i)
				{
					k = equivalentClassToAlleles[i][0] ;
					ecAbundance0[i] = alleleInfo[k].ecAbundance ;
				}
			}
		}

		SetAlleleAbundance(ecReadCount, ecLength);

		delete[] ecLength ;	
		delete[] ecAbundance0 ;
		delete[] ecAbundance1 ;
		delete[] ecAbundance2 ;
		delete[] ecAbundance3 ;
		delete[] ecReadCount ;

		return ret ;
	}

	void RemoveLowMAPQAlleleInEquivalentClass()
	{
		int i, j ;
		double *alleleReadQual = new double[alleleCnt];
		memset(alleleReadQual, 0, sizeof(double) * alleleCnt) ;
		for (i = 0 ; i < readCnt; ++i)
		{
			int size = readAssignments[i].size() ;
			for (j = 0 ; j < size ; ++j)		
				alleleReadQual[readAssignments[i][j].alleleIdx] += readAssignments[i][j].qual ;
		}
		
		int ecCnt = equivalentClassToAlleles.size() ;
		std::vector<int> keptAlleles ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			keptAlleles.clear() ; 
			int size = equivalentClassToAlleles[i].size() ;
			double maxQualSum = -1 ;

			for (j = 0 ; j < size ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				if (alleleReadQual[alleleIdx] > maxQualSum)
					maxQualSum = alleleReadQual[alleleIdx] ;
			}

			for (j = 0 ; j < size ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				if (alleleReadQual[alleleIdx] == maxQualSum)
					keptAlleles.push_back(alleleIdx) ;
			}

			equivalentClassToAlleles[i] = keptAlleles ;
		}

		delete[] alleleReadQual ;
	}

	// Based on the read coverage, remove the ones that are not likely to be true
	void RemoveLowLikelihoodAlleleInEquivalentClass()
	{
		int i, j, k ;
		int ecCnt = equivalentClassToAlleles.size() ;
		std::vector<int> keptAlleles ;
		SimpleVector<int> minStarts ;
		SimpleVector<int> maxEnds ;
		SimpleVector<double> likelihoods ;
		std::map<int, int> alleleIdxToIdx ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			keptAlleles.clear() ; 
			int size = equivalentClassToAlleles[i].size() ;
			likelihoods.ExpandTo(size) ;
			likelihoods.SetZero(0, size) ;
			minStarts.ExpandTo(size) ;
			maxEnds.ExpandTo(size) ;
			alleleIdxToIdx.clear() ;
			
			for (j = 0 ; j < size ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				minStarts[j] = refSet.GetSeqConsensusLen(alleleIdx) ;
				maxEnds[j] = -1 ;
				alleleIdxToIdx[alleleIdx] = j ;
			}
			
			// Setting up the range of coverage for each allele
			int representAlleleIdx = equivalentClassToAlleles[i][0] ;
			int assignedReadCnt = readsInAllele[representAlleleIdx].size() ;
			for (j = 0 ; j < assignedReadCnt ; ++j)
			{
				int readIdx = readsInAllele[representAlleleIdx][j].a ;
				int readToAllelesCnt = readAssignments[readIdx].size() ;
				for (k = 0 ; k < readToAllelesCnt ; ++k)
				{
					struct _readAssignment &assign = readAssignments[readIdx][k] ;
					if (alleleIdxToIdx.find(assign.alleleIdx) == alleleIdxToIdx.end())
						continue ;
					int idx = alleleIdxToIdx[assign.alleleIdx] ;
					if (assign.start < minStarts[idx])
						minStarts[idx] = assign.start ;
					if (assign.end > maxEnds[idx])
						maxEnds[idx] = assign.end ;
				}
			}

			// Compute the likelihood
			double maxLikelihood = -1 ;
			int maxTag = -1 ;
			for (j = 0 ; j < size ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				int len = refSet.GetSeqConsensusLen(alleleIdx) ;
				int effectiveLen = maxEnds[j] - minStarts[j] + 1 ;
				if (effectiveLen > len)
				{
					//printf("Wrong: %d %d %d %lf\n", minStarts[j], maxEnds[j], len, alleleInfo[alleleIdx].ecAbundance) ;
					effectiveLen = len ;
				}
				double ll = pow(double(effectiveLen) / len, alleleInfo[alleleIdx].ecAbundance) ;
				if (ll > maxLikelihood)
				{
					maxLikelihood = ll ;
					maxTag = j ;
				}
				likelihoods[j] = ll ;
			}
			
			// Store the kept alleles
			double cutoff = 0.05 ;
			for (j = 0 ; j < size ; ++j)
			{
				if (likelihoods[j] / maxLikelihood >= cutoff || likelihoods[j] == maxLikelihood)
				{
					keptAlleles.push_back(equivalentClassToAlleles[i][j]) ;
					/*printf("Kept: %d %d %s: %lf %lf\n", i, equivalentClassToAlleles[i][j],
							refSet.GetSeqName(equivalentClassToAlleles[i][j]),
							likelihoods[j], maxLikelihood);*/
				}
				/*else
				{
					printf("Filtered: %d %d %s: %lf %lf\n", i, equivalentClassToAlleles[i][j],
							refSet.GetSeqName(equivalentClassToAlleles[i][j]),
							likelihoods[j], maxLikelihood);
				}*/
			}
			equivalentClassToAlleles[i] = keptAlleles ;
		}
	}

	void SelectAllelesForGenes() // main function for genotyping
	{
		int i, j, k ;

		SimpleVector<bool> readCovered ;
		readCovered.ExpandTo(readCnt) ;
		readCovered.SetZero(0, readCnt) ;

		selectedAlleles.resize(geneCnt) ;
		
		// Compute the abundance for equivalent class
		SimpleVector<struct _pairIntDouble> ecAbundanceList ;
		int ecCnt = equivalentClassToAlleles.size() ;
		
		for (i = 0 ; i < ecCnt ; ++i)
		{
			struct _pairIntDouble np ;
			np.a = i ;
			np.b = alleleInfo[equivalentClassToAlleles[i][0]].ecAbundance ;//abundance ;
			ecAbundanceList.PushBack(np) ;
		}
		std::sort(ecAbundanceList.BeginAddress(), ecAbundanceList.EndAddress(), CompSortPairIntDoubleBDec) ;

		// equivalent classes are sorted by the abundance 
		SimpleVector<int> genesToAdd ;
		SimpleVector<int> allelesToAdd ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			int ec = ecAbundanceList[i].a ;
			
			/*double maxWeight = 0 ;
			for (j = 0 ; j < ecCnt ; ++j)
			{
				const std::vector<int> &readList = readsInAllele[equivalentClassToAlleles[j][0]] ;
				int covered = 0 ;
				int readListSize = readList.size() ;
				for (k = 0 ; k < readListSize ; ++k)
				{
					if (readCovered[readList[k]])
						++covered ;
				}

				double weight = (readListSize - covered) * alleleInfo[equivalentClassToAlleles[j][0]].ecAbundance ;
				//weight = alleleInfo[equivalentClassToAlleles[j][0]].ecAbundance ;
				if (weight > maxWeight)
				{
					maxWeight = weight ;
					ec = j ;
				}
			}*/
			int size = equivalentClassToAlleles[ec].size() ;
			int alleleIdx = equivalentClassToAlleles[ec][0] ;
			
			//if (maxWeight <= 1e-6)
			//	break ;
			
			if (alleleInfo[alleleIdx].ecAbundance <= 1e-6)
				break ;

			//if (ecAbundanceList[i].b <= 1e-6)
			//	break ;

			// Check whether there is uncovered reads.
			double covered = 0;
			const std::vector<struct _pair> &readList = readsInAllele[alleleIdx] ;
			int readListSize = readList.size() ;
			double totalAssignedWeight = 0 ;
			for (j = 0 ; j < readListSize ; ++j)
			{
				if (!IsReadsInAlleleIdxOptimal(readList, j))
						continue ;

				double weight = readAssignments[readList[j].a][0].weight ;
				if (readCovered[readList[j].a])
					covered += weight ;
				totalAssignedWeight += weight ;
			}
#ifdef DEBUG
			printf("%d %s %lf %lf %lf\n", alleleIdx, refSet.GetSeqName(alleleIdx), alleleInfo[alleleIdx].ecAbundance, covered, totalAssignedWeight ) ;
#endif
			//if (covered == totalAssignedWeight) // no uncovered reads
			//	continue ;
			// Add these alleles to the gene allele
			genesToAdd.Clear() ;
			allelesToAdd.Clear() ;
			for (j = 0 ; j < size; ++j)	
			{
				alleleIdx = equivalentClassToAlleles[ec][j] ;
				int geneIdx = alleleInfo[alleleIdx].geneIdx ;
				
				int selectedAllelesCnt = selectedAlleles[geneIdx].size() ;
				/*for (k = 0 ; k < selectedAllelesCnt ; ++k)
				{
					if (alleleInfo[alleleIdx].majorAlleleIdx == 
							alleleInfo[ selectedAlleles[geneIdx][k].a ].majorAlleleIdx)
					{
						allelesToAdd.PushBack(alleleIdx) ;
						break ;
					}
				}
				if (k < selectedAllelesCnt)
					continue ;*/

				// geneMaxAllele is at allele level, ecAbundance is at equivalent class level
				if (alleleInfo[alleleIdx].ecAbundance < filterFrac * geneMaxMajorAlleleAbundance[geneIdx]
						/*&& (totalAssignedWeight - covered < geneMaxMajorAlleleAbundance[geneIdx] / alleleInfo[alleleIdx].ecAbundance / filterFrac
							|| totalAssignedWeight - covered < geneMaxMajorAlleleAbundance[geneIdx] * filterFrac)*/)				
					continue ;
				if (covered == totalAssignedWeight 
						&& (alleleInfo[alleleIdx].ecAbundance < 0.25 * geneMaxMajorAlleleAbundance[geneIdx]
						|| selectedAlleles[geneIdx].size() == 0 || alleleInfo[alleleIdx].ecAbundance < 0.5 * alleleInfo[selectedAlleles[geneIdx].back().a].ecAbundance ) ) 
					continue ;
				/*if (GetGeneAlleleTypes(geneIdx) >= 2)
				{
					// If there is already a good amount haplottypes, we need a harsher cutoff
					int selectedAlleleSize = selectedAlleles[geneIdx].size() ;
					for (k = 0 ; k < selectedAlleleSize ; ++k)
					{
						if (alleleInfo[alleleIdx].ecAbundance > 
								0.25 * majorAlleleAbundance[ alleleInfo[ selectedAlleles[geneIdx][k].a ].majorAlleleIdx ] ) 
						{
							break ;
						}
					}
					if (k >= selectedAlleleSize)
						continue ;
				}*/

				int tmp = genesToAdd.Size() ;
				for (k = 0 ; k < tmp ; ++k)
					if (genesToAdd[k] == geneIdx)
						break ;
				if (k >= tmp)
					genesToAdd.PushBack(geneIdx) ;
				allelesToAdd.PushBack(alleleIdx) ;
			}

			int allelesToAddSize = allelesToAdd.Size() ;
			int quality = 60 ;
			if (genesToAdd.Size() > 1)
			{
				// quality set to 0
				quality = 0 ;
			}
			
			if (genesToAdd.Size() > 0)
			{
				for (j = 0 ; j < readListSize; ++j)
				{
					if (readAssignments[readList[j].a][readList[j].b].qual == 1)
						readCovered[readList[j].a] = true ;
				}
			}
			std::map<int ,int> geneAlleleTypes ;
			for (j = 0 ; j < allelesToAddSize ; ++j)
			{
				alleleIdx = allelesToAdd[j] ;
				int geneIdx = alleleInfo[alleleIdx].geneIdx ;
				int majorAlleleIdx = alleleInfo[alleleIdx].majorAlleleIdx ;
				int alleleRank = -1 ;

				int selectedAlleleSize = selectedAlleles[geneIdx].size() ;
				for (k = 0 ; k < selectedAlleleSize ; ++k)
				{
					if (alleleInfo[selectedAlleles[geneIdx][k].a].majorAlleleIdx == majorAlleleIdx)
					{
						alleleRank = selectedAlleles[geneIdx][k].b ;
						break ;
					}
				}
				if (alleleRank == -1)
				{
					if (geneAlleleTypes.find(geneIdx) != geneAlleleTypes.end()) 
						alleleRank = geneAlleleTypes[geneIdx] ;
					else
					{
						alleleRank = GetGeneAlleleTypes(geneIdx) ;
						geneAlleleTypes[geneIdx] = alleleRank ;
					}
				}
				alleleInfo[alleleIdx].genotypeQuality = quality ;
				alleleInfo[alleleIdx].alleleRank = alleleRank ;
				//printf("%s %lf %d\n", refSet.GetSeqName(alleleIdx), alleleRank, alleleInfo[alleleIdx].ecAbundance ) ;
				if (alleleInfo[alleleIdx].ecAbundance < filterFrac * geneMaxMajorAlleleAbundance[geneIdx])
					alleleInfo[alleleIdx].genotypeQuality = 0 ; 

				struct _pair np ;
				np.a = alleleIdx ;
				np.b = alleleRank ;	
				
				selectedAlleles[geneIdx].push_back(np) ;
			}
		}
		
		// Go through each gene with more than 2 alleles
		int *readCoverage = new int[readCnt] ;
		
		int iter = 0 ;
		const int iterMax = 1000 ;
		int totalCoveredReadCnt = 0 ;
		
		memset(readCoverage, 0, sizeof(int) * readCnt) ;
		std::map<int, int> usedEc ;
		for (i = 0 ; i < geneCnt ; ++i)
		{
			int selectedAlleleCnt = selectedAlleles[i].size() ;
			for (j = 0 ; j < selectedAlleleCnt ; ++j)
			{
				if (selectedAlleles[i][j].b > 1)
					continue ;
				int alleleIdx = selectedAlleles[i][j].a ;
				if (usedEc.find(alleleInfo[alleleIdx].equivalentClass) != usedEc.end())
					continue ;
				usedEc[alleleInfo[alleleIdx].equivalentClass] = 1 ;

				int size = readsInAllele[alleleIdx].size() ;
				for (int r = 0 ; r < size ; ++r)
				{
					if (!IsReadsInAlleleIdxOptimal(readsInAllele[alleleIdx], r))
						continue ;

					if (readCoverage[readsInAllele[alleleIdx][r].a] == 0)
						++totalCoveredReadCnt ;
					++readCoverage[readsInAllele[alleleIdx][r].a] ;
				}
			}
		}
		for (iter = 0 ; iter < iterMax ; ++iter)
		{
			int updatedGeneCnt = 0 ;
			for (i = 0 ; i < geneCnt ; ++i)
			{
				int alleleTypeCnt = GetGeneAlleleTypes(i) ;
				if (alleleTypeCnt <= 2)
					continue ;
				std::map<int, int> coveredReads ;
				std::map<int, int> coveredReadsFromA ;
				SimpleVector<struct _pair> bestTypes ;

				int selectedAlleleCnt = selectedAlleles[i].size() ;
				double maxCover = 0 ;	
				double maxCoverAbundance = 0 ;
				int alleleJ, alleleK ; // representative allele for type J and K 

				// Remove the effects of current gene
				usedEc.clear() ;
				for (j = 0 ; j < selectedAlleleCnt ; ++j)
				{
					if (selectedAlleles[i][j].b > 1)
						continue ;
					int alleleIdx = selectedAlleles[i][j].a ;
					if (usedEc.find(alleleInfo[alleleIdx].equivalentClass) != usedEc.end())
						continue ;
					usedEc[alleleInfo[alleleIdx].equivalentClass] = 1 ;

					int size = readsInAllele[alleleIdx].size() ;
					for (int r = 0 ; r < size ; ++r)
					{
						if (!IsReadsInAlleleIdxOptimal(readsInAllele[alleleIdx], r))
							continue ;
						--readCoverage[readsInAllele[alleleIdx][r].a] ;
					}
				}

				for (j = 0 ; j < alleleTypeCnt - 1 ; ++j)
				{
					int l ;
					usedEc.clear() ;
					coveredReadsFromA.clear() ;
						
					for (l = 0 ; l < selectedAlleleCnt ; ++l)
					{
						if (selectedAlleles[i][l].b != j)
							continue ;
						int alleleIdx = selectedAlleles[i][l].a ;

						if (usedEc.find(alleleInfo[alleleIdx].equivalentClass) != usedEc.end())
							continue ;
						usedEc[alleleInfo[alleleIdx].equivalentClass] = 1 ;

						int r ;
						int size = readsInAllele[alleleIdx].size() ;
						for (r = 0 ; r < size ; ++r)
							if (readCoverage[readsInAllele[alleleIdx][r].a] == 0
									&& IsReadsInAlleleIdxOptimal(readsInAllele[alleleIdx], r))
								coveredReadsFromA[readsInAllele[alleleIdx][r].a] |= 1 ;
						alleleJ = l ;
					}
					for (k = j + 1 ; k < alleleTypeCnt ; ++k)
					{
						coveredReads = coveredReadsFromA ;
						for (l = 0 ; l < selectedAlleleCnt ; ++l)
						{
							if (selectedAlleles[i][l].b != k)
								continue ;
							int alleleIdx = selectedAlleles[i][l].a ;

							if (usedEc.find(alleleInfo[alleleIdx].equivalentClass) != usedEc.end())
								continue ;
							usedEc[alleleInfo[alleleIdx].equivalentClass] = 1 ;

							int r ;
							int size = readsInAllele[alleleIdx].size() ;
							for (r = 0 ; r < size ; ++r)
								if (readCoverage[readsInAllele[alleleIdx][r].a] == 0
										&& IsReadsInAlleleIdxOptimal(readsInAllele[alleleIdx], r))
									coveredReads[readsInAllele[alleleIdx][r].a] |= 3 ;
							alleleK = l ;
						}

						struct _pair np ;
						np.a = j ;
						np.b = k ;
						double abundanceSum = 0 ;
						double abundanceJ = 0 ;
						double abundanceK = 0 ;
						for (l = 0 ; l < selectedAlleleCnt ; ++l)
						{
							if (selectedAlleles[i][l].b == j) 
							{
								abundanceSum += alleleInfo[selectedAlleles[i][l].a].abundance ;
								abundanceJ += alleleInfo[selectedAlleles[i][l].a].abundance ;
							}		
							else if (selectedAlleles[i][l].b == k)
							{
								abundanceSum += alleleInfo[selectedAlleles[i][l].a].abundance ;
								abundanceK += alleleInfo[selectedAlleles[i][l].a].abundance ;
							}
						}	
						//coveredReadCnt += sqrt(abundanceJ) + sqrt(abundanceK) ;
						abundanceSum = abundanceJ * abundanceK ;
						
						double coveredReadCnt = 0 ; //coveredReads.size() ;
						for (std::map<int ,int>::iterator it = coveredReads.begin(); it != coveredReads.end() ; ++it)
						{
							/*if (it->second & 1)
								coveredReadCnt += sqrt(readAssignments[it->first][0].weight) * abundanceJ / (abundanceJ + abundanceK); // the read must have some assignment to be here.
							if (it->second & 3)
								coveredReadCnt += sqrt(readAssignments[it->first][0].weight) * abundanceK / (abundanceJ + abundanceK); // the read must have some assignment to be here.*/
							coveredReadCnt += readAssignments[it->first][0].weight ;
						}
#ifdef DEBUG
						printf("Further selection %s %s %lf %lf %.2lf\n", refSet.GetSeqName(selectedAlleles[i][alleleJ].a), refSet.GetSeqName(selectedAlleles[i][alleleK].a), abundanceJ, abundanceK, coveredReadCnt) ;
#endif
						if (coveredReadCnt > maxCover
								|| (coveredReadCnt == maxCover && abundanceSum > maxCoverAbundance))
						{
							maxCover = coveredReadCnt ;
							maxCoverAbundance = abundanceSum ;
							bestTypes.Clear() ;
							bestTypes.PushBack(np) ;
						}
						else if (coveredReadCnt == maxCover )
						{
							bestTypes.PushBack(np) ;
						}
					} // for k-alleleII
				} // for j-alleleI

				struct _pair bestType = bestTypes[0] ;
				// Rearrange the first two allele as the best selections
				if (bestType.a != 0 || bestType.b != 1)
				{	
					++updatedGeneCnt ;
					for (j = 0 ; j < selectedAlleleCnt ; ++j)
					{
						int newAlleleRank ;
						if (selectedAlleles[i][j].b == bestType.a)
							newAlleleRank = 0 ;
						else if (selectedAlleles[i][j].b == bestType.b)
							newAlleleRank = 1 ;
						else if (selectedAlleles[i][j].b < bestType.a)
							newAlleleRank = selectedAlleles[i][j].b + 2 ;
						else if (selectedAlleles[i][j].b < bestType.b)
							newAlleleRank = selectedAlleles[i][j].b + 1 ;
						else 
							continue ;

						//selectedAlleles[i][j] = selectedAlleles[i][j] ;
						//printf("%d %d %d\n", i, j, newAlleleRank) ;
						selectedAlleles[i][j].b = newAlleleRank ;
						alleleInfo[ selectedAlleles[i][j].a ].alleleRank = newAlleleRank ;
					} // for j-selected allele
				} // for if need update selected alleles

				// Update read coverage.
				usedEc.clear() ;
				for (j = 0 ; j < selectedAlleleCnt ; ++j)
				{
					if (selectedAlleles[i][j].b > 1)
						continue ;
					int alleleIdx = selectedAlleles[i][j].a ;
					if (usedEc.find(alleleInfo[alleleIdx].equivalentClass) != usedEc.end())
						continue ;
					usedEc[alleleInfo[alleleIdx].equivalentClass] = 1 ;

					int size = readsInAllele[alleleIdx].size() ;
					for (int r = 0 ; r < size ; ++r)
					{
						if (IsReadsInAlleleIdxOptimal(readsInAllele[alleleIdx], r))
							++readCoverage[readsInAllele[alleleIdx][r].a] ;
					}
				} // for j- selected allele 
			} // for i-geneCnt

			if (updatedGeneCnt == 0)
				break ;
		} // for iter: global iterations

		// Set the genes with too few abundance to quality 0.
		double *geneAbundances = new double[geneCnt] ;
		double totalGeneAbundance = 0 ;
		for (i = 0 ; i < geneCnt ; ++i)
		{
			int size = selectedAlleles[i].size() ;
			geneAbundances[i] = 0 ;
			for (j = 0 ; j < size ; ++j)
				geneAbundances[i] += alleleInfo[selectedAlleles[i][j].a].abundance ;
			totalGeneAbundance += geneAbundances[i] ;
		}
		/*std::sort(geneAbundances, geneAbundances + geneCnt, CompSortDoubleDec) ;
		double geneAbundanceCutoff = geneAbundances[0] / 10.0 ;
		if (geneCnt > 5)
			geneAbundanceCutoff = geneAbundances[1] / 10.0 ;

		for (i = 0 ; i < geneCnt ; ++i)
		{
			int size = selectedAlleles[i].size() ;
			double abund = 0 ;
			for (j = 0 ; j < size ; ++j)
				abund += alleleInfo[selectedAlleles[i][j].a].abundance ;
			if (abund < geneAbundanceCutoff)
			{
				for (j = 0 ; j < size ; ++j)
					alleleInfo[selectedAlleles[i][j].a].genotypeQuality = 0 ;
			}
		}*/

		// Compute the quality score statistically
		double crossAlleleRate = 0.01 ;
		for (i = 0 ; i < geneCnt ; ++i)
		{
			int type = 0 ;
			std::vector<double> alleleRankAbund ;
			int size = selectedAlleles[i].size() ;
			alleleRankAbund.clear() ;
			int rankCnt = GetGeneAlleleTypes(i) ;
			for (j = 0 ; j < rankCnt ; ++j)
				alleleRankAbund.push_back(0) ;
			
			for (j = 0 ; j < size ; ++j)
				alleleRankAbund[ selectedAlleles[i][j].b ] += alleleInfo[selectedAlleles[i][j].a].abundance ;
			
			double crossGeneNoise = 0 ;
			for (j = 0 ; j < geneCnt ; ++j)
			{
				if (i == j)
					continue ;
				crossGeneNoise += crossGeneRate * (1 + geneSimilarity[i][j]) * geneAbundances[j];
			}
			
			for (j = 0 ; j < rankCnt ; ++j)
			{
				double nullMean = (geneAbundances[i] - alleleRankAbund[j]) * crossAlleleRate + crossGeneNoise ;
				//printf("0: %d %lf %lf %lf\n", i, geneAbundances[i], totalGeneAbundance, alleleRankAbund[j]) ;
				double score = 0 ;
				if (alleleRankAbund[j])
					score = -log(alnorm(2 * (sqrt(alleleRankAbund[j]) - sqrt(nullMean)), true) /** geneCnt * 2*/)/log(double(10.0)) ;
				
				//printf("1: %d %lf %lf %lf\n", i, score, alleleRankAbund[j], nullMean) ;
				if (score > 60)
					score = 60 ;
				if (score < 0)
					score = 0 ;
				
				if (alleleRankAbund[j] < filterCov)
						score = 0 ;

				for (k = 0 ; k < size; ++k)
					if (selectedAlleles[i][k].b == j && alleleInfo[selectedAlleles[i][k].a].genotypeQuality > 0)
					{
						alleleInfo[selectedAlleles[i][k].a].genotypeQuality = (int)score;
					}
			}
		}
		 
		
		delete[] readCoverage ;
		delete[] geneAbundances ;
	}

	int GetGeneCnt()
	{
		return geneCnt ;
	}

	const char *GetGeneName(int geneIdx)
	{
		return geneIdxToName[geneIdx].c_str() ;	
	}

	int GetAlleleDescription(int geneIdx, char *allele1, char *allele2)
	{
		int i, k ;
		int type ;
		SimpleVector<int> selectedMajorAlleles ;
		SimpleVector<bool> used ;
		selectedMajorAlleles.Reserve(majorAlleleCnt) ;
		used.ExpandTo(majorAlleleCnt) ;
		int ret = 0 ;

		used.SetZero(0, majorAlleleCnt);
		int qualities[2] = {-1, -1} ;
		for (type = 0; type <= 1; ++type)	
		{
			double abundance = 0 ;
			char *buffer = allele1 ;
			if (type == 1)
				buffer = allele2 ;
			buffer[0] = '\0' ;

			int size = selectedAlleles[geneIdx].size() ;
			selectedMajorAlleles.Clear() ;

			qualities[type] = -1 ;
			if (type == 1 && qualities[0] == 0)
				used.SetZero(0, majorAlleleCnt) ;
			for (i = 0; i < size; ++i)
			{
				k = selectedAlleles[geneIdx][i].a ; 
				if (selectedAlleles[geneIdx][i].b != type)
					continue ;
				int majorAlleleIdx = alleleInfo[k].majorAlleleIdx ;
				abundance += alleleInfo[k].abundance ;
				if (!used[ majorAlleleIdx ])
				{
					qualities[type] = alleleInfo[k].genotypeQuality ;

					ret = type + 1 ;
					if (buffer[0])
					{
						sprintf(buffer + strlen(buffer), ",%s", majorAlleleIdxToName[majorAlleleIdx].c_str()) ;
					}
					else
						strcpy(buffer, majorAlleleIdxToName[majorAlleleIdx].c_str()) ;
					used[majorAlleleIdx] = 1 ;
				}	
			}
			if (qualities[type] >= 0)
				sprintf(buffer + strlen(buffer), "\t%lf\t%d", abundance, qualities[type]) ;
		}
		return ret ;
	}
} ;

#endif

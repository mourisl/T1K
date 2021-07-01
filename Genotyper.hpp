#ifndef _MOURISL_GENOTYPER
#define _MOURISL_GENOTYPER

#include <math.h>

#include "SeqSet.hpp"

#include "defs.h"
#include "SimpleVector.hpp"

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
	int representative ;
	double count ; // number of reads this group contains.
} ;

struct _readAssignment
{
	int alleleIdx ;
	int start, end ;

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
		strcpy(gene, allele) ;
		strcpy(majorAllele, allele) ;
		for (i = 0 ; allele[i] ; ++i)
		{
			if (allele[i] == '*')
				break ;
		}
		strcpy(gene, allele) ;
		gene[i] = '\0' ;
		for (j = 0 ; j <= 3 && allele[i + j] ; ++j)
			;
		majorAllele[i + j] = '\0' ;
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
	
	bool IsAssignedReadTheSame(const std::vector<int> &l1, const std::vector<int> &l2)
	{
		int cnt1 = l1.size() ;
		int cnt2 = l2.size() ;
		int i ;
		if (cnt1 != cnt2) 
			return false ;
		// The read id in each vector should be sorted
		for (i = 0 ; i < cnt1 ; ++i)
		{
			if (l1[i] != l2[i])
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
			if (a1[i].alleleIdx != a2[i].alleleIdx)
				return false ;
		}
		return true ;
	}

	int Rand()
	{
		return randomSeed = (48271 * randomSeed) & 0x7fffffff ;
	}

	int readCnt ;
	std::vector< std::vector<int> >	readsInAllele ;
	std::vector< std::vector<struct _readAssignment> > readAssignments ;
	std::vector< std::vector<int> > equivalentClassToAlleles ;
	std::vector< std::vector<struct _pair> >	selectedAlleles ; // a-allele name, b-which allele (0,1)

	// variables for allele, majorAllele and genes	
	char *geneBuffer ;
	char *majorAlleleBuffer ;	

	SimpleVector<struct _alleleInfo> alleleInfo ;
	std::map<std::string, int> majorAlleleNameToIdx ;
	std::map<std::string, int> geneNameToIdx ;
	std::vector<std::string> geneIdxToName ;
	std::vector<std::string> majorAlleleIdxToName ;
	int geneCnt ;
	int majorAlleleCnt ;
	int alleleCnt ;

	// variables for abundance
	SimpleVector<double> geneAbundance ;
	SimpleVector<double> majorAlleleAbundance ;
	SimpleVector<double> geneMaxAlleleAbundance ;

	int64_t randomSeed ;
public:
	SeqSet refSet ;
	
	Genotyper(int kmerLength):refSet(kmerLength) 
	{
		geneBuffer = new char[256] ;
		majorAlleleBuffer = new char[256] ;
		alleleCnt = majorAlleleCnt = geneCnt = readCnt = 0 ;
		randomSeed = 17 ;
	}
	~Genotyper() 
	{
		delete[] geneBuffer ;
		delete[] majorAlleleBuffer ;
	}
	
	void InitRefSet(char *filename)
	{
		int i ;

		refSet.InputRefFa(filename) ;

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
				++majorAlleleCnt ;
			}

			
			alleleInfo[i].abundance = 0 ;
			alleleInfo[i].geneIdx = geneNameToIdx[sGene] ;
			alleleInfo[i].majorAlleleIdx = majorAlleleNameToIdx[sMajorAllele] ;
			alleleInfo[i].alleleRank = -1 ;
			alleleInfo[i].abundance = 0 ;
			alleleInfo[i].genotypeQuality = -1 ;
		}
	}

	void InitReadAssignments(int readCnt)
	{
		readAssignments.resize(readCnt) ;
		readsInAllele.resize(alleleCnt) ;

		int i ;
		for (i = 0 ; i < readCnt ; ++i)
			readAssignments[i].clear() ;
		for (i = 0 ; i < alleleCnt ; ++i)
			readsInAllele[i].clear() ;

		this->readCnt = readCnt ;
	}

	void SetReadAssignments(int readId, const std::vector<struct _fragmentOverlap> &assignment)
	{
		int i ;
		int assignmentCnt = assignment.size() ;
		readAssignments[readId].clear() ;
		for (i = 0; i < assignmentCnt; ++i)
		{
			struct _readAssignment na ;
			na.alleleIdx = assignment[i].seqIdx ;
			na.start = assignment[i].seqStart ;
			na.end = assignment[i].seqEnd ;
			readAssignments[readId].push_back(na) ;
		}
	}
	
	// Build the read in allele list
	void FinalizeReadAssignments()
	{
		int i, j ;
		for (i = 0 ; i < readCnt ; ++i)
		{
			int assignmentCnt = readAssignments[i].size() ;
			std::sort(readAssignments[i].begin(), readAssignments[i].end()) ;
			for (j = 0; j < assignmentCnt; ++j)
			{
				readsInAllele[readAssignments[i][j].alleleIdx].push_back(i) ;
			}
		}

		BuildAlleleEquivalentClass() ;
	}

	void SetMajorAlleleAndGeneAbundance()
	{
		int i ;
		// Init other useful abundance data
		geneAbundance.ExpandTo(geneCnt) ;
		geneAbundance.SetZero(0, geneCnt) ;
		majorAlleleAbundance.ExpandTo(majorAlleleCnt) ;
		majorAlleleAbundance.SetZero(0, majorAlleleCnt) ;
		geneMaxAlleleAbundance.ExpandTo(geneCnt) ;
		geneMaxAlleleAbundance.SetZero(0, geneCnt) ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			majorAlleleAbundance[ alleleInfo[i].majorAlleleIdx ] += alleleInfo[i].abundance ;
			geneAbundance[ alleleInfo[i].geneIdx ] += alleleInfo[i].abundance ;

			if (alleleInfo[i].abundance > geneMaxAlleleAbundance[ alleleInfo[i].geneIdx ])
			{
				geneMaxAlleleAbundance[ alleleInfo[i].geneIdx ] = alleleInfo[i].abundance ;
			}
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
		
		SetMajorAlleleAndGeneAbundance() ;
	}
	
	int GetGeneAlleleTypes(int geneIdx)
	{
		if ( selectedAlleles[geneIdx].size() == 0 )
			return  0 ;
		else
			return selectedAlleles[geneIdx].back().b + 1;
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
					np.b = ((uint32_t)np.b * readCnt + readsInAllele[i][j]) % FINGERPRINT_MAX ;
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
		return ecCnt ;
	}

	void QuantifyAlleleEquivalentClass()
	{
		int i, j, k ;
		int t ; // iteration for EM algorithm.
		int ecCnt = equivalentClassToAlleles.size() ;	
		
		
		SimpleVector<struct _pair> readFingerprint ;
		const int FINGERPRINT_MAX = 1000003 ;
		
		// First rebuild the read groups: the reads mapping to the same set of alleles
		for (i = 0 ; i < readCnt ; ++i)
		{
			struct _pair np ;
			np.a = i ;
			np.b = -1 ;
			readFingerprint.PushBack(np) ;
		}

		for (i = 0 ; i < readCnt ; ++i)	
		{
			int size = readAssignments[i].size() ;
			if (size == 0)
				continue ;
			readFingerprint[i].b = 0 ;
			for (j = 0 ; j < size ; ++j)
			{
				k = readAssignments[i][j].alleleIdx ;
				readFingerprint[i].b = (readFingerprint[i].b * (int64_t)readCnt + k) % FINGERPRINT_MAX ;
			}
		}
	
		SimpleVector<int> readToReadGroup ;
		SimpleVector<struct _readGroupInfo> readGroupInfo ; // the read represent the read group, so we can easily obtain the assignment information.
		
		readToReadGroup.ExpandTo(readCnt) ;
		std::sort(readFingerprint.BeginAddress(), readFingerprint.EndAddress(), CompSortPairByBDec) ;

		int rgCnt = 0 ;
		int effectiveReadCnt = 0 ; // the number of reads covering alleles
		for (i = 0 ; i < readCnt ; ++i)
		{
			bool newRg = true ; // rg: read group
			if (readFingerprint[i].b == -1)
			{
				readToReadGroup[ readFingerprint[i].a ] = -1 ;
				continue ;
			}
			for (j = i - 1 ; j >= 0 ; --j)
			{
				if (readFingerprint[i].b != readFingerprint[j].b)
					break ;
				if (IsReadAssignmentTheSame(readAssignments[readFingerprint[i].a],
							readAssignments[readFingerprint[j].a]))
				{
					newRg = false ;
					break ;
				}
			}
			
			int readId = readFingerprint[i].a ;
			if (newRg)
			{
				readToReadGroup[readId] = rgCnt ;
				struct _readGroupInfo nrg ;
				nrg.representative = readId ;
				nrg.count = 1 ;
				readGroupInfo.PushBack(nrg) ;
				++rgCnt ;
			}
			else
			{
				readToReadGroup[readId] = readToReadGroup[ readFingerprint[j].a ] ;
				++readGroupInfo[ readToReadGroup[readId] ].count ;
			}
			++effectiveReadCnt ;
		}
		
		// Convert readgroup_to_allele to readgroup_to_alleleEquivalentClass
		std::vector< std::vector<int> > readGroupToAlleleEc ;
		readGroupToAlleleEc.resize(rgCnt) ;
		std::map<int, bool> ecUsed ;
		for (i = 0 ; i < rgCnt ; ++i)
		{
			int readId = readGroupInfo[i].representative ;
			int size = readAssignments[readId].size() ;
			ecUsed.clear() ;
			for (j = 0 ; j < size ; ++j)
			{
				int ecIdx = alleleInfo[readAssignments[readId][j].alleleIdx].equivalentClass ;
				if (ecUsed.find(ecIdx) == ecUsed.end())
				{
					readGroupToAlleleEc[i].push_back(ecIdx) ;
					ecUsed[ecIdx] = true ;
				}
			}
		}

		// Start the EM algorithm
		double *emResults ;
		const int maxEMIterations = 1000 ;
		const int maxRandIterations = 1 ;
		double *ecAbundance = new double[ecCnt] ;
		double *ecReadCount = new double[ecCnt] ;
		double *ecLength = new double[ecCnt] ; // the sequence length for equivalent class.
		int randIter ;

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
		
		double **emEcReadCount = new double*[maxRandIterations] ;
		for (randIter = 0 ; randIter < maxRandIterations ; ++randIter)
		{
			for (i = 0 ; i < ecCnt ; ++i)
			{
				//ecAbundance[i] = Rand()%7 + 1 ; //1.0 / ecCnt ;
				ecAbundance[i] = 1.0 / ecCnt ;
			}

			for (t = 0 ; t < maxEMIterations ; ++t)
			{
				// E-step: find the expected number of reads
				memset(ecReadCount, 0, sizeof(double) * ecCnt) ;
				for (i = 0 ; i < rgCnt ; ++i)
				{
					double psum	= 0 ;
					int size = readGroupToAlleleEc[i].size() ;
					for (j = 0 ; j < size ; ++j)
						psum += ecAbundance[readGroupToAlleleEc[i][j]] ;
					if (psum == 0)	
						psum = 1 ;
					for (j = 0 ; j < size ; ++j)
					{
						int ecIdx = readGroupToAlleleEc[i][j] ;
						ecReadCount[ecIdx] += readGroupInfo[i].count * ecAbundance[ecIdx] / psum ;
					}
				}

				// M-step: recompute the abundance
				double diffSum = 0 ;
				double normalization = 0 ;
				for (i = 0 ; i < ecCnt ; ++i)
					normalization += ecReadCount[i] / ecLength[i] ;

				for (i = 0 ; i < ecCnt ; ++i)
				{
					double tmp = ecReadCount[i] / effectiveReadCnt ;
					tmp = ecReadCount[i] / ecLength[i] / normalization ;
					//printf("%d: %lf %lf %lf. %lf\n", i, tmp, ecReadCount[i], ecLength[i], ecAbundance[i]) ;
					diffSum += ABS(tmp - ecAbundance[i]) ;
					ecAbundance[i] = tmp ;
				}
				//printf("%lf\n", diffSum) ;
				if (diffSum < 1e-3 && t < maxEMIterations - 2)
					t = maxEMIterations - 2 ; // Force one more iteration
			}

			emEcReadCount[randIter] = new double[ecCnt] ;
			memcpy(emEcReadCount[randIter], ecReadCount, sizeof(double) * ecCnt) ;
		}

		for (i = 0 ; i < alleleCnt ; ++i)
			alleleInfo[i].abundance = alleleInfo[i].ecAbundance = 0 ;

		for (i = 0 ; i < ecCnt ; ++i)
		{
			int size = equivalentClassToAlleles[i].size() ;
			double abund = 0 ; //emEcReadCount[0][i] ;
			for (j = 0 ; j < maxRandIterations ; ++j)
			{
				//if (abund < emEcReadCount[j][i])
				k = equivalentClassToAlleles[i][0] ;
				//printf("%d %d %s %lf %d %d\n", i, k, refSet.GetSeqName(k), emEcReadCount[j][i], refSet.GetSeqConsensusLen(k),readsInAllele[k].size()) ;
				abund += emEcReadCount[j][i] ;
			}
			//printf("%lf\n", abund) ;
			abund /= maxRandIterations ;
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
		SetMajorAlleleAndGeneAbundance() ;
		
		for (i = 0 ; i < maxRandIterations ; ++i)
			delete[] emEcReadCount[i] ;
		delete[] emEcReadCount ;

		delete[] ecLength ;	
		delete[] ecAbundance ;
		delete[] ecReadCount ;
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
				int readIdx = readsInAllele[representAlleleIdx][j] ;
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

	void QuantifyAlleleEquivalentClass_test()
	{
		int i, j, k ;
		int t ; // iteration for EM algorithm.
		int ecCnt = equivalentClassToAlleles.size() ;	
		
		SimpleVector<struct _pair> readFingerprint ;
		const int FINGERPRINT_MAX = 1000003 ;
		
		// First rebuild the read groups: the reads mapping to the same set of alleles
		for (i = 0 ; i < readCnt ; ++i)
		{
			struct _pair np ;
			np.a = i ;
			np.b = -1 ;
			readFingerprint.PushBack(np) ;
		}

		for (i = 0 ; i < readCnt ; ++i)	
		{
			int size = readAssignments[i].size() ;
			if (size == 0)
				continue ;
			readFingerprint[i].b = 0 ;
			for (j = 0 ; j < size ; ++j)
			{
				k = readAssignments[i][j].alleleIdx ;
				readFingerprint[i].b = (readFingerprint[i].b * (int64_t)readCnt + k) % FINGERPRINT_MAX ;
			}
		}
	
		SimpleVector<int> readToReadGroup ;
		SimpleVector<struct _readGroupInfo> readGroupInfo ; // the read represent the read group, so we can easily obtain the assignment information.
		
		readToReadGroup.ExpandTo(readCnt) ;
		std::sort(readFingerprint.BeginAddress(), readFingerprint.EndAddress(), CompSortPairByBDec) ;

		int rgCnt = 0 ;
		int effectiveReadCnt = 0 ; // the number of reads covering alleles
		for (i = 0 ; i < readCnt ; ++i)
		{
			bool newRg = true ; // rg: read group
			if (readFingerprint[i].b == -1)
			{
				readToReadGroup[ readFingerprint[i].a ] = -1 ;
				continue ;
			}
			for (j = i - 1 ; j >= 0 ; --j)
			{
				if (readFingerprint[i].b != readFingerprint[j].b)
					break ;
				if (IsReadAssignmentTheSame(readAssignments[readFingerprint[i].a],
							readAssignments[readFingerprint[j].a]))
				{
					newRg = false ;
					break ;
				}
			}
			
			int readId = readFingerprint[i].a ;
			if (newRg)
			{
				readToReadGroup[readId] = rgCnt ;
				struct _readGroupInfo nrg ;
				nrg.representative = readId ;
				nrg.count = 1 ;
				readGroupInfo.PushBack(nrg) ;
				++rgCnt ;
			}
			else
			{
				readToReadGroup[readId] = readToReadGroup[ readFingerprint[j].a ] ;
				++readGroupInfo[ readToReadGroup[readId] ].count ;
			}
			++effectiveReadCnt ;
		}
		
		// Convert readgroup_to_allele to readgroup_to_alleleEquivalentClass
		std::vector< std::vector<int> > readGroupToAlleleEc ;
		readGroupToAlleleEc.resize(rgCnt) ;
		std::map<int, bool> ecUsed ;
		for (i = 0 ; i < rgCnt ; ++i)
		{
			int readId = readGroupInfo[i].representative ;
			int size = readAssignments[readId].size() ;
			ecUsed.clear() ;
			for (j = 0 ; j < size ; ++j)
			{
				int ecIdx = alleleInfo[readAssignments[readId][j].alleleIdx].equivalentClass ;
				if (ecUsed.find(ecIdx) == ecUsed.end())
				{
					readGroupToAlleleEc[i].push_back(ecIdx) ;
					ecUsed[ecIdx] = true ;
				}
			}
		}

		// Start the EM algorithm
		const int maxEMIterations = 1000 ;
		double *ecAbundance = new double[ecCnt] ;
		double *ecReadCount = new double[ecCnt] ;
		double *ecLength = new double[ecCnt] ; // the sequence length for equivalent class.
		ecUsed.clear() ;

		double diffSum = 0 ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			ecUsed[i] = false ;
			int size = equivalentClassToAlleles[i].size() ;
			ecLength[i] = refSet.GetSeqEffectiveLen(equivalentClassToAlleles[i][0]) ;
			for (j = 1 ; j < size ; ++j)
			{
				int len = refSet.GetSeqEffectiveLen(equivalentClassToAlleles[i][j]) ;
				if (len < ecLength[i])
					ecLength[i] = len ;
			}
		}
		
		int selectionIter = 0 ;
		for (i = 0 ; i < alleleCnt ; ++i)
			alleleInfo[i].abundance = alleleInfo[i].ecAbundance = 0 ;
		
		for (selectionIter = 0 ; selectionIter < ecCnt ; ++selectionIter)
		{
			for (i = 0 ; i < ecCnt ; ++i)
				ecAbundance[i] = 1.0 / ecCnt ;
			for (t = 0 ; t < maxEMIterations ; ++t)
			{
				// E-step: find the expected number of reads
				memset(ecReadCount, 0, sizeof(double) * ecCnt) ;
				for (i = 0 ; i < rgCnt ; ++i)
				{
					double psum	= 0 ;
					int size = readGroupToAlleleEc[i].size() ;
					for (j = 0 ; j < size ; ++j)
					{
						if (ecUsed[readGroupToAlleleEc[i][j]])
							continue ;
						psum += ecAbundance[readGroupToAlleleEc[i][j]] ;
					}
					for (j = 0 ; j < size ; ++j)
					{
						int ecIdx = readGroupToAlleleEc[i][j] ;
						if (ecUsed[ecIdx])
							continue ;
						ecReadCount[ecIdx] += readGroupInfo[i].count * ecAbundance[ecIdx] / psum ;
					}
				}

				// M-step: recompute the abundance
				diffSum = 0 ;
				double normalization = 0 ;
				for (i = 0 ; i < ecCnt ; ++i)
					normalization += ecReadCount[i] / ecLength[i] ;

				for (i = 0 ; i < ecCnt ; ++i)
				{
					double tmp = ecReadCount[i] / effectiveReadCnt ;
					tmp = ecReadCount[i] / ecLength[i] / normalization ;
					//printf("%lf %lf %d. %lf\n", tmp, ecReadCount[i], effectiveReadCnt, ecAbundance[i]) ;
					diffSum += ABS(tmp - ecAbundance[i]) ;
					ecAbundance[i] = tmp ;
				}
				//printf("%lf\n", diffSum) ;
				if (diffSum < 1e-3 && t < maxEMIterations - 2)
					t = maxEMIterations - 2 ; // Force one more iteration
			}

			double maxAbund = 0 ;
			int maxEc = 0 ;
			for (i = 0 ; i < ecCnt ; ++i)		
			{
				if (ecReadCount[i] > maxAbund)
				{
					maxAbund = ecReadCount[i] ;
					maxEc = i ;
				}
			}
			if (maxAbund < 1)
				break ;
			
			// adjust the read count
			int maxEcAlleleIdx = equivalentClassToAlleles[maxEc][0] ;
			int alleleReadCnt = readsInAllele[maxEcAlleleIdx].size() ;
			std::map<int, int> visitedRg ;
			for (j = 0 ; j < alleleReadCnt; ++j)
			{
				visitedRg[readToReadGroup[readsInAllele[maxEcAlleleIdx][j]]] = 1 ;
			}

			double updatedAbund = 0 ;
			for (std::map<int, int>::iterator it = visitedRg.begin() ; it != visitedRg.end() ; ++it)
			{
				int rg = it->first ;
				double psum	= 0 ;
				int size = readGroupToAlleleEc[rg].size() ;
				for (j = 0 ; j < size ; ++j)
				{
					if (ecUsed[readGroupToAlleleEc[rg][j]])
						continue ;
					psum += ecAbundance[readGroupToAlleleEc[rg][j]] ;
				}

				double partAbund = readGroupInfo[rg].count * ecAbundance[maxEc] / psum ;
				if (ecAbundance[maxEc] / psum < 0.95)
					partAbund = readGroupInfo[rg].count * 0.95 ;

				readGroupInfo[rg].count -= partAbund ;
				updatedAbund += partAbund ;
			}
			
			//printf("%d %lf %lf %s\n", selectionIter, maxAbund, updatedAbund, refSet.GetSeqName(equivalentClassToAlleles[maxEc][0])) ;
			int ecAlleleCnt = equivalentClassToAlleles[maxEc].size() ;
			for (j = 0 ; j < ecAlleleCnt ; ++j)
			{
				k = equivalentClassToAlleles[maxEc][j] ;
				alleleInfo[k].abundance = updatedAbund / ecAlleleCnt ;
				alleleInfo[k].ecAbundance = updatedAbund ;
			}
			ecUsed[maxEc] = true ;
		}
		delete[] ecLength ;	
		delete[] ecAbundance ;
		delete[] ecReadCount ;
		SetMajorAlleleAndGeneAbundance() ;
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
			int covered = 0;
			const std::vector<int> &readList = readsInAllele[alleleIdx] ;
			int readListSize = readList.size() ;
			for (j = 0 ; j < readListSize ; ++j)
			{
				if (readCovered[readList[j]])
					++covered ;
			}
			//printf("%d %s %lf %d %d\n", alleleIdx, refSet.GetSeqName(alleleIdx), alleleInfo[alleleIdx].ecAbundance, covered, readListSize ) ;
			if (covered == readListSize) // no uncovered reads
				continue ;
			// Add these alleles to the gene allele
			genesToAdd.Clear() ;
			allelesToAdd.Clear() ;
			for (j = 0 ; j < size; ++j)	
			{
				alleleIdx = equivalentClassToAlleles[ec][j] ;
				int geneIdx = alleleInfo[alleleIdx].geneIdx ;
				
				if (alleleInfo[alleleIdx].ecAbundance < 0.1 * geneMaxAlleleAbundance[geneIdx])				
					continue ;
				//if (GetGeneAlleleTypes(geneIdx) >= 2)
				//	continue ;

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
					readCovered[readList[j]] = true ;
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
				struct _pair np ;
				np.a = alleleIdx ;
				np.b = alleleRank ;	
				
				selectedAlleles[geneIdx].push_back(np) ;
			}
		}
		
		// Go through each gene with more than 2 alleles
		for (i = 0 ; i < geneCnt ; ++i)
		{
			int alleleTypeCnt = GetGeneAlleleTypes(i) ;
			if (alleleTypeCnt <= 2)
				continue ;
			std::map<int, int> usedEc ;
			std::map<int, int> coveredReads ;
			std::map<int, int> coveredReadsFromA ;
			SimpleVector<struct _pair> bestTypes ;
			
			int selectedAlleleCnt = selectedAlleles[i].size() ;
			int maxCover = 0 ;	
			for (j = 0 ; j < alleleTypeCnt - 1 ; ++j)
			{
				int l ;
				usedEc.clear() ;
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
						coveredReadsFromA[readsInAllele[alleleIdx][r]] = 1 ;
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
							coveredReads[readsInAllele[alleleIdx][r]] = 1 ;
					}
					
					int coveredReadCnt = coveredReads.size() ;
					struct _pair np ;
					np.a = j ;
					np.b = k ;
					//printf("Further selection %d %d %d\n", j, k, coveredReadCnt) ;
					if (coveredReadCnt > maxCover )
					{
						maxCover = coveredReadCnt ;
						bestTypes.Clear() ;
						bestTypes.PushBack(np) ;
					}
					else if (coveredReadCnt == maxCover)
					{
						bestTypes.PushBack(np) ;
					}
				} // for k-alleleII
			} // for j-alleleI

			struct _pair bestType = bestTypes[0] ;
			k = 0 ;
			for (j = 0 ; j < selectedAlleleCnt ; ++j)
			{
				if (selectedAlleles[i][j].b == bestType.a
						|| selectedAlleles[i][j].b == bestType.b)
				{
					int newAlleleRank = 0;
					if (selectedAlleles[i][j].b == bestType.b)
						newAlleleRank = 1 ;
					selectedAlleles[i][k] = selectedAlleles[i][j] ;
					selectedAlleles[i][k].b = newAlleleRank ;
					alleleInfo[ selectedAlleles[i][k].a ].alleleRank = newAlleleRank ;
					++k ;
				}
			} // for j-selected allele
		} // for i-geneCnt

		// Set the genes with too few abundance to quality 0.
		double *geneAbundances = new double[geneCnt] ;
		for (i = 0 ; i < geneCnt ; ++i)
		{
			int size = selectedAlleles[i].size() ;
			geneAbundances[i] = 0 ;
			for (j = 0 ; j < size ; ++j)
				geneAbundances[i] += alleleInfo[selectedAlleles[i][j].a].abundance ;
		}
		std::sort(geneAbundances, geneAbundances + geneCnt, CompSortDoubleDec) ;
		double geneAbundanceCutoff = geneAbundances[1] / 10.0 ;

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
		}

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
			if (qualities >= 0)
				sprintf(buffer + strlen(buffer), "\t%lf\t%d", abundance, qualities[type]) ;
		}
		return ret ;
	}
} ;

#endif

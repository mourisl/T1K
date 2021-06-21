#ifndef _MOURISL_GENOTYPER
#define _MOURISL_GENOTYPER

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
	int count ; // number of reads this group contains.
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
public:
	SeqSet refSet ;
	
	Genotyper(int kmerLength):refSet(kmerLength) 
	{
		geneBuffer = new char[256] ;
		majorAlleleBuffer = new char[256] ;
		alleleCnt = majorAlleleCnt = geneCnt = readCnt = 0 ;
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

	void InitAlleleAbundance(FILE *fp)
	{
		// TODO: calculate abundance with out own method
		int i ;	
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

		const int maxEMIterations = 1000 ;
		double *ecAbundance = new double[ecCnt] ;
		double *ecReadCount = new double[ecCnt] ;
		
		double diffSum = 0 ;
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
					psum += ecAbundance[readGroupToAlleleEc[i][j]] ;
				
				for (j = 0 ; j < size ; ++j)
				{
					int ecIdx = readGroupToAlleleEc[i][j] ;
					ecReadCount[ecIdx] += readGroupInfo[i].count * ecAbundance[ecIdx] / psum ;
				}
			}

			// M-step: recompute the abundance
			diffSum = 0 ;
			for (i = 0 ; i < ecCnt ; ++i)
			{
				double tmp = ecReadCount[i] / effectiveReadCnt ;
				//printf("%lf %lf %d. %lf\n", tmp, ecReadCount[i], effectiveReadCnt, ecAbundance[i]) ;
				diffSum += ABS(tmp - ecAbundance[i]) ;
				ecAbundance[i] = tmp ;
			}
			//printf("%lf\n", diffSum) ;
			if (diffSum < 1e-3)
				break ;
		}

		for (i = 0 ; i < alleleCnt ; ++i)
			alleleInfo[i].abundance = alleleInfo[i].ecAbundance = 0 ;

		for (i = 0 ; i < ecCnt ; ++i)
		{
			int size = equivalentClassToAlleles[i].size() ;
			for (j = 0 ; j < size ; ++j)
			{
				alleleInfo[i].abundance = ecAbundance[i] / size ;
				alleleInfo[i].ecAbundance = ecAbundance[i] ;
			}
		}

		delete[] ecAbundance ;
		delete[] ecReadCount ;
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
		QuantifyAlleleEquivalentClass() ;
		
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
			int size = equivalentClassToAlleles[ec].size() ;
			int alleleIdx = equivalentClassToAlleles[ec][0] ;

			// Check whether there is uncovered reads.
			int covered = 0;
			const std::vector<int> &readList = readsInAllele[alleleIdx] ;
			int readListSize = readList.size() ;
			for (j = 0 ; j < readListSize ; ++j)
			{
				if (readCovered[readList[j]])
					++covered ;
			}
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
				if (GetGeneAlleleTypes(geneIdx) >= 2)
					continue ;

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
				int alleleRank = 0 ;
				if (geneAlleleTypes.find(geneIdx) != geneAlleleTypes.end()) 
					alleleRank = geneAlleleTypes[geneIdx] ;
				else
				{
					alleleRank = GetGeneAlleleTypes(geneIdx) ;
					geneAlleleTypes[geneIdx] = alleleRank ;
				}
				alleleInfo[alleleIdx].genotypeQuality = quality ;
				alleleInfo[alleleIdx].alleleRank = alleleRank ;
				
				struct _pair np ;
				np.a = alleleIdx ;
				np.b = alleleRank ;	
				
				selectedAlleles[geneIdx].push_back(np) ;
			}
		}
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

		for (type = 0; type <= 1; ++type)	
		{
			double abundance = 0 ;
			char *buffer = allele1 ;
			if (type == 1)
				buffer = allele2 ;
			buffer[0] = '\0' ;

			int size = selectedAlleles[geneIdx].size() ;
			selectedMajorAlleles.Clear() ;
			used.SetZero(0, majorAlleleCnt);

			int quality = -1 ;
			for (i = 0; i < size; ++i)
			{
				k = selectedAlleles[geneIdx][i].a ; 
				if (selectedAlleles[geneIdx][i].b != type)
					continue ;
				quality = alleleInfo[k].genotypeQuality ;

				ret = type + 1 ;
				abundance = alleleInfo[k].ecAbundance ;
				int majorAlleleIdx = alleleInfo[k].majorAlleleIdx ;
				if (!used[ majorAlleleIdx ] )
				{
					if (buffer[0])
					{
						sprintf(buffer + strlen(buffer), ",%s", majorAlleleIdxToName[majorAlleleIdx].c_str()) ;
					}
					else
						strcpy(buffer, majorAlleleIdxToName[majorAlleleIdx].c_str()) ;
					used[majorAlleleIdx] = 1 ;
				}	
			}
			if (quality >= 0)
				sprintf(buffer + strlen(buffer), "\t%lf\t%d", abundance, quality) ;
		}
		return ret ;
	}
} ;

#endif

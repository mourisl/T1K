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
	double ecAbundance ; // the sum of abundance from the equivalent class.
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

	int readCnt ;
	std::vector< std::vector<int> >	readsInAllele ;
	std::vector< std::vector<int> > readAssignments ;
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
	}

	void SetReadAssignments(int readId, const std::vector<struct _fragmentOverlap> &assignment)
	{
		int i ;
		int assignmentCnt = assignment.size() ;
		readAssignments[readId].clear() ;
		for (i = 0; i < assignmentCnt; ++i)
		{
			readAssignments[readId].push_back(assignment[i].seqIdx) ;
		}
	}

	void FlushReadAssignments(int readIdOffset)
	{
		int i, j ;
		int size = readAssignments.size() ;
		for (i = 0 ; i < size ; ++i)
		{
			int assignmentCnt = readAssignments[i].size() ;
			for (j = 0; j < assignmentCnt; ++j)
			{
				readsInAllele[readAssignments[i][j]].push_back(i + readIdOffset) ;
			}
		}
		if (readIdOffset + size > readCnt)
			readCnt = readIdOffset + size ;
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

	void SelectAllelesForGenes() // main function for genotyping
	{
		int i, j, k ;

		SimpleVector<bool> readCovered ;
		readCovered.ExpandTo(readCnt) ;
		readCovered.SetZero(0, readCnt) ;

		selectedAlleles.resize(geneCnt) ;
		
		SimpleVector<struct _pairIntDouble> alleleAbundanceList ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			struct _pairIntDouble np ;
			np.a = i ;
			np.b = alleleInfo[i].abundance ;//abundance ;
			alleleAbundanceList.PushBack(np) ;
		}
		std::sort(alleleAbundanceList.BeginAddress(), alleleAbundanceList.EndAddress(), CompSortPairIntDoubleBDec) ;

		// Build equivalent class: alleles covered by the same set of reads and have the same abundance.
		int ecCnt = 0 ; // equivalent class count
		std::vector< std::vector<int> > equivalentClassToAlleles ;
		for (i = 0 ; i < alleleCnt ; )
		{
			for (j = i + 1; j < alleleCnt; ++j)
				if (alleleAbundanceList[i].b != alleleAbundanceList[j].b)
					break ;
			
			int alleleIdx = alleleAbundanceList[i].a ;
			if (alleleInfo[alleleIdx].abundance <= 0)
				break ;

			alleleInfo[alleleIdx].equivalentClass = ecCnt ;
			equivalentClassToAlleles.push_back(std::vector<int>()) ;
			equivalentClassToAlleles.back().push_back(alleleIdx) ;
			++ecCnt ;
			for (k = i + 1; k < j ; ++k)
			{
				int l ;
				alleleIdx = alleleAbundanceList[k].a ;
				for (l = k - 1 ; l >= i ; --l)
				{
					int alleleIdx2 = alleleAbundanceList[l].a ;
					if (IsAssignedReadTheSame(readsInAllele[alleleIdx], readsInAllele[alleleIdx2]))
					{
						int ec = alleleInfo[alleleIdx2].equivalentClass ;
						equivalentClassToAlleles[ec].push_back(alleleIdx) ;
						alleleInfo[alleleIdx].equivalentClass = ec ;
						break ;
					}
				}
				if (l < i)
				{
					alleleInfo[alleleIdx].equivalentClass = ecCnt ;
					equivalentClassToAlleles.push_back(std::vector<int>()) ;
					equivalentClassToAlleles.back().push_back(alleleIdx) ;
					++ecCnt ;
				}
			}

			i = j ;
		}
		// Compute the abundance for equivalent class
		SimpleVector<struct _pairIntDouble> ecAbundanceList ;
		for (i = 0 ; i < ecCnt ; ++i)
		{
			int ecSize = equivalentClassToAlleles[i].size() ;
			int ecAbundance = 0 ;
			for (j = 0 ; j < ecSize ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				ecAbundance += alleleInfo[alleleIdx].abundance ;
			}
			for (j = 0 ; j < ecSize ; ++j)
			{
				int alleleIdx = equivalentClassToAlleles[i][j] ;
				alleleInfo[alleleIdx].ecAbundance = ecAbundance ; 
			}
			
			struct _pairIntDouble np ;
			np.a = i ;
			np.b = ecAbundance ;//abundance ;
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

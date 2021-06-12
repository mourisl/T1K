#ifndef _MOURISL_GENOTYPER
#define _MOURISL_GENOTYPER

#include "SeqSet.hpp"

#include "defs.h"
#include "SimpleVector.hpp"

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
	std::vector< std::vector<struct _pair> >	selectedAllele ; // a-allele name, b-which allele (0,1)

	// variables for allele, majorAllele and genes	
	char *geneBuffer ;
	char *majorAlleleBuffer ;	

	SimpleVector<struct _pairIntDouble> alleleAbundanceList ;
	SimpleVector<double> alleleAbundance ;
	std::map<std::string, int> majorAlleleNameToIdx ;
	std::map<std::string, int> geneNameToIdx ;
	SimpleVector<int> alleleToGene ;
	SimpleVector<int> alleleToMajorAllele ;
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

			alleleToGene.PushBack(geneNameToIdx[sGene]) ;
			alleleToMajorAllele.PushBack(majorAlleleNameToIdx[sMajorAllele]) ;
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
		alleleAbundanceList.Reserve(alleleCnt) ;
		alleleAbundance.ExpandTo(alleleCnt) ;
	
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
			alleleAbundanceList.PushBack(np) ;
			alleleAbundance[np.a] = np.b ;
		}
		fclose(fp) ;
	
		std::sort(alleleAbundanceList.BeginAddress(), alleleAbundanceList.EndAddress(), CompSortPairIntDoubleBDec) ;

		// Init other useful abundance data
		geneAbundance.ExpandTo(geneCnt) ;
		geneAbundance.SetZero(0, geneCnt) ;
		majorAlleleAbundance.ExpandTo(majorAlleleCnt) ;
		majorAlleleAbundance.SetZero(0, majorAlleleCnt) ;
		geneMaxAlleleAbundance.ExpandTo(geneCnt) ;
		geneMaxAlleleAbundance.SetZero(0, geneCnt) ;
		for (i = 0 ; i < alleleCnt ; ++i)
		{
			majorAlleleAbundance[ alleleToMajorAllele[i] ] += alleleAbundance[i] ;
			geneAbundance[ alleleToGene[i] ] += alleleAbundance[i] ;

			if (alleleAbundance[i] > geneMaxAlleleAbundance[ alleleToGene[i] ])
			{
				geneMaxAlleleAbundance[ alleleToGene[i] ] = alleleAbundance[i] ;
			}
		}
	}

	void SelectAllelesForGenes()
	{
		int i, j, k ;

		SimpleVector<bool> readCovered ;
		readCovered.ExpandTo(readCnt) ;
		readCovered.SetZero(0, readCnt) ;

		selectedAllele.resize(geneCnt) ;

		for (i = 0 ; i < alleleCnt ; ++i)
		{
			k = alleleAbundanceList[i].a ;
			ParseAlleleName( refSet.GetSeqName(k), geneBuffer, majorAlleleBuffer) ;
			if (alleleAbundanceList[i].b <= 0)	
				break ;
			const std::vector<int> &readList = readsInAllele[k] ;
			int size = readList.size() ;
			int covered = 0 ;
			//printf("%d %d %d %d\n", i, k, ) ;
			for (j = 0 ; j < size ; ++j)
			{
				if (!readCovered[ readList[j] ])
					++covered ;
			}
			//printf("%s %lf\n", refSet.GetSeqName(k), alleleAbundance[k]);

			int geneIdx = alleleToGene[k] ;
			if (majorAlleleAbundance[ alleleToMajorAllele[k] ] < 0.1 * geneMaxAlleleAbundance[geneIdx]) // TODO: set a good cut off
				continue ;

			int selectedCnt = selectedAllele[geneIdx].size() ;
			int max = -1 ;
			for (j = 0 ; j < selectedCnt ; ++j)
				if (selectedAllele[geneIdx][j].b > max)
					max = selectedAllele[geneIdx][j].b ;
			if ( max >= 2)
				continue ;
			if (selectedCnt > 0 && alleleAbundance[ selectedAllele[geneIdx][selectedCnt - 1].a ] > alleleAbundance[k] )
				continue ;
			if (covered > 0) 
			{
				// check whether we need to
				struct _pair np ;
				np.a = k ;
				np.b = max + 1 ;	
				selectedAllele[geneIdx].push_back(np);
				for (j = 0 ; j < size ; ++j)
				{
					if (!readCovered[ readList[j] ])
					{
						readCovered[readList[j]] = true ;
					}
				}
			}
			else
			{
				// Check whethere this is equivalent to some added alleles
				for (j = selectedCnt - 1 ; j >= 0 ; --j)
				{
					int sId = selectedAllele[geneIdx][j].a ; //selected id
					if (alleleAbundance[sId] != alleleAbundance[k])
						break ;
					if (IsAssignedReadTheSame( readsInAllele[sId], readsInAllele[k]))
					{
						struct _pair np ;
						np.a = k ;
						np.b = selectedAllele[geneIdx][j].b ;
						selectedAllele[geneIdx].push_back(np) ;
						break ;
					}
				} // for selected alleles
			} // else cover==0
		} // for i - alleleCnt
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

			int size = selectedAllele[geneIdx].size() ;
			selectedMajorAlleles.Clear() ;
			used.SetZero(0, majorAlleleCnt);
			for (i = 0; i < size; ++i)
			{
				k = selectedAllele[geneIdx][i].a ; 
				if (selectedAllele[geneIdx][i].b != type)
					continue ;
				ret = type + 1 ;
				abundance += alleleAbundance[k] ;
				int majorAlleleIdx = alleleToMajorAllele[k] ;
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
			sprintf(buffer + strlen(buffer), "\t%lf", abundance) ;
		}
		return ret ;
	}
} ;

#endif

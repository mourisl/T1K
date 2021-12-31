#ifndef _MOURISL_BARCODESUMMARY
#define _MOURISL_BARCODESUMMARY

#include <stdio.h>
#include <map>
#include <vector>

#include "SeqSet.hpp"
#include "VariantCaller.hpp"

class BarcodeSummary
{
private:
	SeqSet &refSet ;
	std::map<int, std::vector<struct _pairID> >	barcodeAlleleCount ; //a-unique count, b-fractioned count
public:
	BarcodeSummary(SeqSet &inRefSet):refSet(inRefSet)
	{
				
	}

	~BarcodeSummary() {}
	
	void AddFragment(char *read1, char *read2, int barcode, VariantCaller *pVariantCaller, std::vector<struct _fragmentOverlap> &fragmentAssignments)
	{
		int i ;
		if (barcodeAlleleCount.find(barcode) == barcodeAlleleCount.end())
		{
			std::vector<struct _pairID>	nv ;
			int alleleCnt = refSet.Size() ;
			for (i = 0 ; i < alleleCnt ; ++i)
			{
				struct _pairID np ;
				np.a = 0 ;
				np.b = 0 ;
				nv.push_back(np) ;
			}

			barcodeAlleleCount[barcode] = nv ;
		}

		std::vector<struct _fragmentOverlap> adjustedAssignments ;
		if (pVariantCaller == NULL)
			adjustedAssignments = fragmentAssignments ;
		else
			adjustedAssignments = pVariantCaller->AdjustFragmentAssignment(read1, read2, fragmentAssignments) ;

		int assignCnt = adjustedAssignments.size() ;
		for (i = 0 ; i < assignCnt ; ++i)
		{
			barcodeAlleleCount[barcode][adjustedAssignments[i].seqIdx].b += 1.0 / assignCnt ;
			if (assignCnt == 1)
			{
				++barcodeAlleleCount[barcode][adjustedAssignments[i].seqIdx].a ;
			}
		}
	}

	void Output(std::vector<std::string> &barcodeIntToStr, FILE *fp)
	{
		int i ;
		int alleleCnt = refSet.Size() ;
		// print header
		fprintf(fp, "#barcode") ;
		for ( i = 0 ; i < alleleCnt ; ++i)	
			fprintf(fp, "\t%s", refSet.GetSeqName(i)) ;
		for ( i = 0 ; i < alleleCnt ; ++i)	
			fprintf(fp, "\t%s_uniq", refSet.GetSeqName(i)) ;
		fprintf(fp, "\n") ;
		for (std::map<int, std::vector<struct _pairID> >::iterator it = barcodeAlleleCount.begin() ;
				it != barcodeAlleleCount.end() ; ++it)
		{
			fprintf(fp, "%s", barcodeIntToStr[it->first].c_str()) ;
			for (i = 0 ; i < alleleCnt ; ++i)
				fprintf(fp, "\t%lf", it->second[i].b) ;
			for (i = 0 ; i < alleleCnt ; ++i)
				fprintf(fp, "\t%d", it->second[i].a) ;
			fprintf(fp, "\n") ;
		}
	}
} ;

#endif

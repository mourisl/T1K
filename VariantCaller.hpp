#ifndef _MOURISL_VARIANT_CALLER
#define _MOURISL_VARIANT_CALLER

#include "SeqSet.hpp"
#include "Genotyper.hpp"

struct _variant
{
	int seqIdx ;
	int refStart, refEnd ;
	char ref[10] ;
	char var[10] ;
	double allSupport ;
	double varSupport ;
	double varUniqSupport ;
	
	int varGroupId ; 
	int outputGroupId ; // 0-best variants, 1- equal best variants
	int qual ;
} ;

struct _baseVariant
{
	double count[4] ;
	double uniqCount[4] ;
	double unweightedCount[4] ;
	struct _pairIntDouble alignInfo[4] ; // information for best alignment. 
	bool exon ;
	int candidateId ; // -1: not a variant candidate. 
	std::vector<int> finalVariantIds ; // the id in the final variant table

	double AllCountSum()
	{
		return count[0] + count[1] + count[2] + count[3] ;
	}
	
	double UniqCountSum()
	{
		return uniqCount[0] + uniqCount[1] + uniqCount[2] + uniqCount[3] ;
	}

	double UnweightedCountSum()
	{
		return unweightedCount[0] + unweightedCount[1] + unweightedCount[2] + unweightedCount[3] ;
	}

	double IsGoodAssignment(int matchCnt, double similarity)
	{
		int i ;
		for (i = 0 ; i < 4 ; ++i)
			if (matchCnt < alignInfo[i].a - 4) 
				return false ;
		return true ;
	}
} ;

struct _adjFragmentToBaseVariant
{
	int seqIdx ;
	int refPos ;
	char nuc[5] ; // which nucleotide this fragment support
	int weight ;

	int next ;	
} ;

struct _adjBaseVariantToFragment
{
	int fragIdx ;
	char nuc[5] ;
	int next ;	
} ;

// used to determine how to group variants
struct _adjBaseVariantToBaseVariant 
{
	int varIdx ;
	double weight ;
	bool rootCandidate ; // whether this variant is inferred from read coverage.
	int next ;
} ;

struct _enumVarResult
{
	double bestCover ; // the number of covered results
	int usedVarCnt ; // the number of introduced variant
	SimpleVector<char> bestEnumVariants ;
	SimpleVector<char> equalBestEnumVariants ;
} ;


class VariantCaller
{
private:
	SeqSet &refSet ;
	std::vector< SimpleVector<struct _baseVariant> > baseVariants ;
	std::vector<double> seqAbundance ;
	SimpleVector<struct _pair> candidateVariants ; // a: seqidx, b: refpos
	SimpleVector<int> candidateVariantGroupId ; // variant id to group id.
	SimpleVector<int> seqCopy ; // 1-homozygous, 2-heterzygous
	std::vector<struct _variant> finalVariants ; 
	void UpdateBaseVariantFromOverlap(char *read, double weight, bool filterLowQual, struct _overlap o)
	{
		if (o.seqIdx == -1)
			return ;
		int i, k ;
		int readLen = strlen(read)	;
		char *r = read ;
		if (o.strand == -1)
		{
			r = strdup(read) ;
			refSet.ReverseComplement(r, read, readLen) ;
		}
		char *align = o.align ;
		if (align == NULL)
		{
		  char *align = new char[ 3 * readLen + 2 ] ;
			AlignAlgo::GlobalAlignment( refSet.GetSeqConsensus(o.seqIdx) + o.seqStart, 
				o.seqEnd - o.seqStart + 1,
				r + o.readStart,
				o.readEnd - o.readStart + 1, align) ;
		}

		int refPos = o.seqStart ;
		int readPos = o.readStart ;

		refPos = o.seqStart ;
		readPos = o.readStart ;
		for ( k = 0 ; align[k] != -1 ; ++k )
		{
			if ( align[k] == EDIT_MATCH || align[k] == EDIT_MISMATCH)
			{
				/*if (weight > 0 
						&& (o.seqIdx == 15 || o.seqIdx == 16)&& refPos == 3490)
				{
					printf("%d %d %lf %c %s\n", o.seqIdx, align[k], weight, r[readPos], read);
				}*/
				if (filterLowQual && !baseVariants[o.seqIdx][refPos].IsGoodAssignment(o.matchCnt, o.similarity))
					continue ;
				if (r[readPos] == 'N')
					continue;

				int nucIdx = nucToNum[r[readPos] - 'A'] ;
				if (weight == 1)
					baseVariants[o.seqIdx][refPos].uniqCount[nucIdx] += weight ;
				baseVariants[o.seqIdx][refPos].count[ nucIdx ] += 1;
				baseVariants[o.seqIdx][refPos].unweightedCount[ nucIdx ] += 1 ;

				if (o.matchCnt > baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].a )
				{
					baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].a = o.matchCnt ;
					baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].b = o.similarity ;
				}
				else if (o.matchCnt == baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].a 
						&& o.similarity > baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].b)
				{
					baseVariants[o.seqIdx][refPos].alignInfo[nucIdx].b = o.similarity ;
				}
			}
			//TODO: handle indels
			//if (isValidDiff[refPos].exon && (align[k] == EDIT_DELETE) )
			//	printf("%s\n%s\n", seq.consensus + o.seqStart, r) ;
			if (align[k] != EDIT_INSERT)
				++refPos ;
			if (align[k] != EDIT_DELETE)
				++readPos ;
		}
		//char *align = new char[ 3 * readLen + 2 ] ;

		if (o.strand == -1)
			free(r);
	}

	/*int GetCandidateVariantGroup(int gid)
	{
		if (candidateVariantGroupId[gid] != gid)
			return candidateVariantGroup[gid] = GetCandidateVariantGroup(candidateVariantGroup[gid]) ;
		return gid ;
	}*/
	
	bool containCandidateVar(int start, int end, SimpleVector<int> &candidateVarAccuCount)
	{
		/*if (start == 0)
		{
			if (candidateVarAccuCount[start] > 0
					|| candidateVarAccuCount[end] != candidateVarAccuCount[start])
				return true ;
		}
		else
		{
			if (candidateVarAccuCount[end] != candidateVarAccuCount[start - 1])
				return true ;
		} 
		return false ;*/

		// candidateVarAccuCount starts from position 0 for value 0
		return candidateVarAccuCount[start] != candidateVarAccuCount[end + 1] ;
	}
	
	struct _overlap SelectOverlapFromFragmentOverlap(int k, struct _fragmentOverlap &frag)
	{
		if (k == 0)
			return frag.overlap1 ;
		else if (k == 1)
			return frag.overlap2 ;	
		return frag.overlap1 ;
	}

	/*char *SelectRead(char *r1, char *r2, struct _fragmentOverlap &frag)
	{
	}*/
	
	void ComputeCandidateVarAccuCount(int seqIdx, SimpleVector<int> &candidateVarAccuCount)
	{
		int i ;
		int len = baseVariants[seqIdx].Size() ;
		candidateVarAccuCount.ExpandTo(len + 1) ;
		candidateVarAccuCount[0] = 0 ;
		for (i = 0 ; i < len ; ++i)
		{
			if (baseVariants[seqIdx][i].candidateId != -1)
				candidateVarAccuCount[i + 1] = candidateVarAccuCount[i] + 1 ;
			else
				candidateVarAccuCount[i + 1] = candidateVarAccuCount[i] ;
		}
	}
public:
	VariantCaller(SeqSet &inRefSeq):refSet(inRefSeq) 
	{
		int i, j ;
		int seqCnt = refSet.Size() ;
		baseVariants.resize(seqCnt) ;
		for (i = 0 ; i < seqCnt ; ++i)
		{
			int len = refSet.GetSeqConsensusLen(i) ;
			baseVariants[i].ExpandTo(len) ;
			baseVariants[i].SetZero(0, len) ;
			for (j = 0 ; j < len ; ++j)
			{
				baseVariants[i][j].exon = refSet.IsPosInExon(i, j) ;
				baseVariants[i][j].candidateId = -1 ;
			}
		}
	} 
	~VariantCaller() {} 

	void SetSeqAbundance(Genotyper &genotyper) 
	{
		int seqCnt = refSet.Size() ;
		int i ;
		seqAbundance.resize(seqCnt) ;
		for (i = 0 ; i < seqCnt ; ++i)
		{
			seqAbundance[i] = genotyper.GetAlleleAbundance(i) ;
			//printf("%d %s %lf\n", i, refSet.GetSeqName(i), seqAbundance[i]) ;
		}
		std::map<int, int> geneAlleleCount ;
		for (i = 0 ; i < seqCnt ; ++i)
			geneAlleleCount[ genotyper.GetAlleleGeneIdx(i) ] += 1 ;
		seqCopy.ExpandTo(seqCnt) ;
		for (i = 0 ; i < seqCnt ; ++i)
			seqCopy[i] = geneAlleleCount[ genotyper.GetAlleleGeneIdx(i)] ;
	}
	

	// updateType: 0-weight, 1-alignInfo
	void UpdateBaseVariantFromFragmentOverlap(char *read1, char *read2, int updateType, std::vector<struct _fragmentOverlap> &fragmentAssignment)
	{
		int i ;
		int assignCnt = fragmentAssignment.size() ;
		bool filterLowQual = true ;
		double totalWeight = 0 ;
		for (i = 0 ; i < assignCnt ; ++i)
			totalWeight += seqAbundance[fragmentAssignment[i].seqIdx] ;

		for (i = 0 ; i < assignCnt ; ++i)
		{
			struct _fragmentOverlap &fragOverlap = fragmentAssignment[i] ;
			int seqIdx = fragOverlap.seqIdx ;
			double weight = seqAbundance[seqIdx] / totalWeight ;
			if (updateType == 1)
			{
				filterLowQual = false ;
				weight = 0 ;
			}
			if (fragOverlap.hasMatePair)
			{
				UpdateBaseVariantFromOverlap(read1, weight, filterLowQual, fragOverlap.overlap1) ;
				UpdateBaseVariantFromOverlap(read2, weight, filterLowQual, fragOverlap.overlap2) ;
			}
			else
			{
				if (!fragOverlap.o1FromR2)
					UpdateBaseVariantFromOverlap(read1, weight, filterLowQual, fragOverlap.overlap1) ;
				else
					UpdateBaseVariantFromOverlap(read2, weight, filterLowQual, fragOverlap.overlap1) ;
			}
		}
	}

	void FindCandidateVariants()
	{
		int i, j, k ;
		candidateVariants.Clear() ;  		
		int seqCnt = refSet.Size() ;
		const int countThreshold = 5 ;

		for (i = 0 ; i < seqCnt ; ++i)
		{
			int len = baseVariants[i].Size() ;
			const char *s = refSet.GetSeqConsensus(i) ;
			double factor = 0.5 ;
			//if (seqCopy[i] <= 1)
			//	factor = 0.25 ;
			for (j = 0 ; j < len ; ++j)
			{
				double refCount = baseVariants[i][j].count[ nucToNum[ s[j] - 'A' ]] ;
				for (k = 0 ; k < 4 ; ++k)
				{
					if (baseVariants[i][j].count[k] >= countThreshold 
							//&& (baseVariants[i][j].count[k] >= sqrt(refCount) 
							//&& baseVariants[i][j].count[k] >= refCount - 6 * sqrt(refCount))
							&& baseVariants[i][j].count[k] >= refCount * factor 
							&& k != nucToNum[s[j] - 'A'])
					{
						int id = candidateVariants.Size() ;
						struct _pair np ;
						np.a = i ;
						np.b = j ;	
						candidateVariants.PushBack(np) ;
						baseVariants[i][j].candidateId = id ;
						candidateVariantGroupId.PushBack(-1) ;
						//printf("%lf %lf\n", baseVariants[i][j].AllCountSum(), baseVariants[i][j].count[k]) ;
						break ;
					}
				}
			}
		}
	}

	void ExpandCandidateVariantsFromFragmentOverlap(char *read1, char *read2, std::vector<struct _fragmentOverlap> &fragmentAssignment, SimpleVector<struct _adjBaseVariantToBaseVariant> &adjVarToVar, std::vector<SimpleVector<int> > &seqCandidateAccuCount)
	{
		if (fragmentAssignment.size() <= 0)
			return ;

		int i, j, k ;

		SimpleVector<int> refPos, readPos ;
		SimpleVector<int> alignIdx ;
		SimpleVector<bool> validAssignment ; // check whether the overlap can be used for candidate variant expansion
		SimpleVector<char *> r ;
		int assignCnt = fragmentAssignment.size() ;

		refPos.ExpandTo(assignCnt) ;
		readPos.ExpandTo(assignCnt) ;
		r.ExpandTo(assignCnt) ;
		alignIdx.ExpandTo(assignCnt) ;
		validAssignment.ExpandTo(assignCnt) ;
		for (k = 0 ; k <= 1 ; ++k) // 0-read1, 1-read2
		{
			// Check whether there is variants in the alignment region
			if (k == 1 && !fragmentAssignment[0].hasMatePair)
				break ;

			for (i = 0 ; i < assignCnt ; ++i)
			{
				int seqIdx = fragmentAssignment[i].seqIdx ;
				struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
				if (containCandidateVar(o.seqStart, o.seqEnd, seqCandidateAccuCount[seqIdx])) ;
					break ;
			}

			if (i >= assignCnt) // no candidate variants
				continue ;
			
			char *read = read1 ;
			if (k == 1
					|| (k == 0 && fragmentAssignment[0].o1FromR2))
				read = read2 ;
			int len = strlen(read) ;
			//char *rc = (char *)malloc(sizeof(char) * (len + 1)) ;
			//refSet.ReverseComplement(rc, read, len) ;
			for (i = 0 ; i < assignCnt ; ++i) 
			{
				int seqIdx = fragmentAssignment[i].seqIdx ;
				struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
				/*if (o.strand == 1)
					r[i] = read ;
				else if (o.strand == -1)
					r[i] = rc ;*/
				refPos[i] = o.seqStart ;
				readPos[i] = o.readStart ;
			}	
			///free(rc) ;
			
			// They all should have the same start position in read position
			for (i = 1 ; i < assignCnt ; ++i)
			{
				if (readPos[i] != readPos[0])
					break ;
			}
			if (i < assignCnt)
				continue ;
			
			alignIdx.SetZero(0, assignCnt) ;
			for (j = 0 ; j < len ; ++j) // use the read pos as the anchor
			{
				// Expand the set of candidate variants
				int firstCandidateId = -1 ;
				int firstCandidateIdx = -1 ;
				for (i = 0 ; i < assignCnt ; ++i)
				{
					struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
					if (refPos[i] < refSet.GetSeqConsensusLen(o.seqIdx)) 
						validAssignment[i] = baseVariants[o.seqIdx][refPos[i]].IsGoodAssignment(o.matchCnt, o.similarity) ;
					else
						validAssignment[i] = false ;
				}

				for (i = 0 ; i < assignCnt ; ++i)
				{
					if (!validAssignment[i])
						continue ;
					struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
					if (refPos[i] < refSet.GetSeqConsensusLen(o.seqIdx) && baseVariants[o.seqIdx][refPos[i]].candidateId != -1)
					{
						firstCandidateId = baseVariants[o.seqIdx][refPos[i]].candidateId ;
						firstCandidateIdx = i ;
						break ;
					}
				}
				
				/*if (!strcmp("AGTGTCGTTAAATGTCCCCTCTCTGTGCAGAAGGAAGTGCTCAAACCTGACATCTGACCAACATTGCAGGATGACTGTCTCTTCTGATTTCACCAGGGGACCTGGGTGGGCCAGGAGGGAAGGTTTTCTGTGGACTCCTAGGAAGAGAGG", read1))
				{
					if (refPos[1] == 1052)
					{
						for (i = 0 ; i < assignCnt ; ++i)
						{
							struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
							printf("### %d %s %d %d. %d\n", fragmentAssignment[i].seqIdx, refSet.GetSeqName(fragmentAssignment[i].seqIdx), refPos[i], baseVariants[2][1052].candidateId, o.matchCnt) ;
						}
					}
				}*/

				if (firstCandidateId != -1)
				{
					// contains candidate varivants
					for (i = 0 ; i < assignCnt ; ++i)
					{
						if (!validAssignment[i])
							continue ;
						struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
						if (baseVariants[o.seqIdx][refPos[i]].candidateId == -1
								&& (o.align[alignIdx[i]] != -1 
									&& (o.align[ alignIdx[i] ] == EDIT_MATCH || o.align[ alignIdx[i] ] == EDIT_MISMATCH)))
						{
							int cid = candidateVariants.Size() ;
							struct _pair np ;
							np.a = o.seqIdx ;
							np.b = refPos[i] ;
							candidateVariants.PushBack(np) ;
							baseVariants[o.seqIdx][refPos[i]].candidateId = cid ;
							adjVarToVar[cid].varIdx = cid ;
							adjVarToVar[cid].rootCandidate = false ;
							adjVarToVar[cid].next = -1 ;
							candidateVariantGroupId.PushBack(-1) ;
							/*if (np.a == 8 && np.b == 703)
							{
								printf("strange %d %d %d %d %d: %s %s\n", o.seqIdx, refPos[i], cid, k, o.strand, read1, read2) ;

								char *r = read ;
								char *rc = NULL ;
								if (o.strand == -1)
								{
									rc = strdup(r) ;
									refSet.ReverseComplement(rc, r, len) ;
									r = rc ;
								}
								AlignAlgo::VisualizeAlignment(refSet.GetSeqConsensus(o.seqIdx) + o.seqStart, o.seqEnd - o.seqStart + 1, r + o.readStart, o.readEnd - o.readStart + 1, o.align) ;
								//if (rc != NULL)
								//	free(rc) ;
								
								struct _overlap tmpo = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[firstCandidateIdx]) ;
								r = read ;
								if (tmpo.strand == -1)
									r = rc ;
								printf("anchor %d %d %d %d\n", tmpo.seqIdx, refPos[firstCandidateIdx], k, o.strand) ;
								AlignAlgo::VisualizeAlignment(refSet.GetSeqConsensus(tmpo.seqIdx) + tmpo.seqStart, tmpo.seqEnd - tmpo.seqStart + 1, r + tmpo.readStart, tmpo.readEnd - tmpo.readStart + 1, tmpo.align) ;
								if (rc != NULL)
									free(rc) ;
							}*/
						}
						int cid = baseVariants[o.seqIdx][refPos[i]].candidateId ;
						//if (cid != -1)
						//	candidateVariantGroup[cid] = GetCandidateVariantGroup(firstCandidateId) ;
						if (cid != -1)
							candidateVariantGroupId[cid] = -1 ;
					}

					// Update the var to var abundance
					for (i = 0 ; i < assignCnt ; ++i)
					{	
						if (!validAssignment[i])
							continue ;
						int l ;
						for (l = 0 ; l < assignCnt ; ++l)
						{
							if (i == l || !validAssignment[l])
								continue ;

							struct _overlap oI = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
							struct _overlap oL = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[l]) ;
							int cidI = baseVariants[oI.seqIdx][refPos[i]].candidateId ;
							int cidL = baseVariants[oL.seqIdx][refPos[l]].candidateId ;
							if (cidI == -1 || cidL == -1)
								continue ;
							// add weight of i to j
							int p ;
							p = adjVarToVar[cidI].next ;
							while (p != -1)
							{
								if (adjVarToVar[p].varIdx == cidL)
								{
									++adjVarToVar[p].weight ;
									break ;
								}
								p = adjVarToVar[p].next ;
							}

							if (p == -1)
							{
								struct _adjBaseVariantToBaseVariant na ;
								na.varIdx = cidL ;
								na.weight = 1 ;
								na.rootCandidate = false ;
								na.next = adjVarToVar[cidI].next ;
								
								adjVarToVar[cidI].next = adjVarToVar.Size() ;
								adjVarToVar.PushBack(na) ;
							}
						}
					} // for i-assignCnt. update_var to var abundance
				}

				// Move to the next read position
				for (i = 0 ; i < assignCnt ; ++i)
				{
					struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
					char *align = o.align ;

					while (align[alignIdx[i]] != -1
							&& readPos[i] <= j)
					{
						int aidx = alignIdx[i] ;
						if (align[aidx] != EDIT_INSERT)
							++refPos[i] ;
						if (align[aidx] != EDIT_DELETE)
							++readPos[i] ;

						++alignIdx[i] ;
					}
				}
			} // for j on read position
		} // for k on read end selection
	}

	void BuildCandidateVariantGroup(int from, int tag, SimpleVector< struct _adjBaseVariantToBaseVariant > &adjVarToVar)
	{
		if (candidateVariantGroupId[from] != -1)
			return ;
		candidateVariantGroupId[from] = tag ;
		int p = adjVarToVar[from].next ;
		while (p != -1)
		{
			int to = adjVarToVar[p].varIdx ;
			int fromSeqIdx = candidateVariants[from].a ;
			int fromRefPos = candidateVariants[from].b ;
			int toSeqIdx = candidateVariants[to].a ;
			int toRefPos = candidateVariants[to].b ;
			
			//printf("%s %d %s %d %d %s %d %lf %lf\n", __func__, fromSeqIdx, refSet.GetSeqName(fromSeqIdx), fromRefPos, toSeqIdx, refSet.GetSeqName(toSeqIdx), toRefPos, adjVarToVar[p].weight, baseVariants[fromSeqIdx][fromRefPos].UnweightedCountSum()) ;
			if (adjVarToVar[p].weight >= baseVariants[fromSeqIdx][fromRefPos].UnweightedCountSum() * 0.15 || adjVarToVar[p].weight >= baseVariants[toSeqIdx][toRefPos].UnweightedCountSum() * 0.15)
				BuildCandidateVariantGroup(to, tag, adjVarToVar) ;

			p = adjVarToVar[p].next ;
		}
	}

	void BuildFragmentCandidateVarGraph(char *read1, char *read2, int fragIdx, std::vector<struct _fragmentOverlap> &fragmentAssignment, std::vector< SimpleVector<int> > &seqCandidateAccuCount, SimpleVector<struct _adjFragmentToBaseVariant> &adjFrag, SimpleVector<struct _adjBaseVariantToFragment> &adjVar)
	{
		int i, j, k ;
		int assignCnt = fragmentAssignment.size() ;
		if (assignCnt <= 0)
			return ;
		for (k = 0 ; k <= 1 ; ++k) // 0-read1, 1-read2
		{
			// Check whether there is variants in the alignment region
			if (k == 1 && !fragmentAssignment[0].hasMatePair)
				break ;

			for (i = 0 ; i < assignCnt ; ++i)
			{
				int seqIdx = fragmentAssignment[i].seqIdx ;
				struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
				if (containCandidateVar(o.seqStart, o.seqEnd, seqCandidateAccuCount[seqIdx])) ;
					break ;
			}

			if (i >= assignCnt) // no candidate variants
				continue ;
			
			char *read = read1 ;
			if (k == 1
					|| (k == 0 && fragmentAssignment[0].o1FromR2))
				read = read2 ;
			int len = strlen(read) ;
			char *rc = (char *)malloc(sizeof(char) * (len + 1)) ;
			refSet.ReverseComplement(rc, read, len) ;
			for (i = 0 ; i < assignCnt ; ++i) 
			{
				int seqIdx = fragmentAssignment[i].seqIdx ;
				struct _overlap o = SelectOverlapFromFragmentOverlap(k, fragmentAssignment[i]) ;
				char *r ;
				if (o.strand == 1)
					r = read ;
				else if (o.strand == -1)
					r = rc ;
				int refPos = o.seqStart ;
				int readPos = o.readStart ;
				char *align = o.align ;
				for (j = 0 ; align[j] != -1 ; ++j)
				{
					int cid = baseVariants[seqIdx][refPos].candidateId ;
					if (cid != -1 /*&& baseVariants[seqIdx][refPos].IsGoodAssignment(o.matchCnt, o.similarity)*/)
					{
						char var[5] ;
						var[0] = r[readPos] ;
						var[1] = '\0' ;

						// Check whether the edge has already been put.
						int p ;
					
						p = adjVar[cid].next ;
						while (p != -1)
						{
							if (adjVar[p].fragIdx == fragIdx 
									&& !strcmp(var, adjVar[p].nuc))
								break ;
							p = adjVar[p].next ;
						}
						
						if (p == -1)
						{
							// Add the edge
							struct _adjFragmentToBaseVariant nFragToBaseVar ;
							struct _adjBaseVariantToFragment nBaseVarToFrag ;

							strcpy(nFragToBaseVar.nuc, var) ;
							nFragToBaseVar.seqIdx = seqIdx ;
							nFragToBaseVar.refPos = refPos ;
							nFragToBaseVar.weight = 1 ;
							nFragToBaseVar.next = adjFrag[fragIdx].next ;
							adjFrag[fragIdx].next = adjFrag.Size() ;
							adjFrag.PushBack(nFragToBaseVar) ;

							strcpy(nBaseVarToFrag.nuc, var) ;
							nBaseVarToFrag.fragIdx = fragIdx ;
							nBaseVarToFrag.next = adjVar[cid].next ;
							adjVar[cid].next = adjVar.Size() ;
							adjVar.PushBack(nBaseVarToFrag) ;
						}
					}
					if (align[j] != EDIT_INSERT)
						++refPos ;
					if (align[j] != EDIT_DELETE)
						++readPos ;
				}
			}	// for j on read position
			free(rc) ;
		} // for k on read end selection
	}	
	
	void EnumerateVariants(int depth, SimpleVector<char> &choices, struct _enumVarResult &result, SimpleVector<int> &fragIds, SimpleVector<int> &vars, SimpleVector<struct _adjFragmentToBaseVariant> &adjFrag, SimpleVector<struct _adjBaseVariantToFragment> &adjVar)
	{
		int i ;
		if (depth == vars.Size())
		{
			SimpleVector<int> fragCovered ; 
			int fragCnt = fragIds.Size() ;
			int varCnt = vars.Size() ;
			int usedVarCnt = 0 ;	
			int maxFragIdx = 0 ;
			for (i = 0 ; i < fragCnt ; ++i)
				if (fragIds[i] > maxFragIdx)
					maxFragIdx = fragIds[i] ;
			++maxFragIdx ;
			fragCovered.ExpandTo(maxFragIdx) ;
			fragCovered.SetZero(0, maxFragIdx) ;
			for (i = 0 ; i < varCnt ; ++i)
			{
				int seqIdx = candidateVariants[vars[i]].a ;
				int refPos = candidateVariants[vars[i]].b ;
				if (varCnt <= 1 && seqCopy[ candidateVariants[vars[i]].a ] <= 1
						&& choices[i] != refSet.GetSeqConsensus(seqIdx)[refPos])
					continue ;
				int p = adjVar[ vars[i] ].next ;
				while (p != -1)
				{
					int fragIdx = adjVar[p].fragIdx ;
					if (fragIdx < maxFragIdx)
					{
						if (adjVar[p].nuc[0] == choices[i])
						{
							//if (candidateVariants[vars[0]].b==961)
							//	printf("%d %d\n", i, fragIdx);
							fragCovered[fragIdx] = 1 ;
						}
					}
					p = adjVar[p].next ;
				}
			}
			
			// We just want to test the contribution of the alternative nucleotide,
			// if there is some noise support the reference nuc, it is fine, we don't report those.
			// Only do this when the scenario is simple.
			for (i = 0 ; i < varCnt && varCnt <= 1 ; ++i) 
			{
				if (seqCopy[ candidateVariants[vars[i]].a ] != 1)
					continue ;
				int refContribution = 0 ;
				int altContribution = 0 ;
				int seqIdx = candidateVariants[vars[i]].a ;
				int refPos = candidateVariants[vars[i]].b ;

				if (choices[i] == refSet.GetSeqConsensus(seqIdx)[refPos])
					continue ;

				int p = adjVar[ vars[i] ].next ;
				while (p != -1)
				{
					if (adjVar[p].nuc[0] == choices[i])
						++altContribution ;
					else if (refSet.GetSeqConsensus(seqIdx)[refPos] == adjVar[p].nuc[0])
						++refContribution ;
					p = adjVar[p].next ;
				}	
				//if (candidateVariants[vars[0]].b==270)
				//printf("%d %s %d %d %d\n", i, refSet.GetSeqName(seqIdx), refPos, refContribution, altContribution) ;

				bool includeAlt = false ;
				if ( ((altContribution >= 2 && baseVariants[seqIdx][refPos].uniqCount[ nucToNum[choices[i] - 'A'] ] > 0) 
					  || (altContribution >= 10 )) &&
						altContribution > 0.15 * refContribution)
				{
					includeAlt = true ;
				}
				
				p = adjVar[ vars[i] ].next ;
				while (p != -1)
				{
					if (refSet.GetSeqConsensus(seqIdx)[refPos] == adjVar[p].nuc[0]
							|| (choices[i] == adjVar[p].nuc[0] && includeAlt))
					{
						int fragIdx = adjVar[p].fragIdx ;
						if (fragCovered[fragIdx] == 0)
							fragCovered[fragIdx] = 2 ;
					}
					p = adjVar[p].next ;
				}
			}

			double covered = 0 ;
			for (i = 0 ; i < fragCnt ; ++i)
			{
				if (fragCovered[fragIds[i]])
					++covered ;
			}
			
			for (i = 0 ; i < varCnt ; ++i)
			{
				int seqIdx = candidateVariants[ vars[i] ].a ;
				int refPos = candidateVariants[ vars[i] ].b ;
				if (refSet.GetSeqConsensus(seqIdx)[refPos] != choices[i])
					++usedVarCnt ; 
			}
			/*if (vars[0] == 2)
			{
				printf("%lf %d\n", covered, usedVarCnt) ;
			}*/
			/*if (candidateVariants[vars[0]].b==1007)
			{
				printf("%lf %c %c\n", covered, choices[0], choices[1]) ;
			}*/
			if (covered > result.bestCover
					|| (covered == result.bestCover && usedVarCnt < result.usedVarCnt))
			{
				result.bestCover = covered ;
				result.usedVarCnt = usedVarCnt ;
				result.bestEnumVariants = choices ;
				result.equalBestEnumVariants.Clear() ;
			}
			else if (covered == result.bestCover && usedVarCnt == result.usedVarCnt)
			{
				result.equalBestEnumVariants = choices ;
			}
			return ;	
		}
	
		for (i = 0 ; i < 4 ; ++i)
		{
			choices[depth] = numToNuc[i] ;
			EnumerateVariants(depth + 1, choices, result, fragIds, vars, adjFrag, adjVar) ;
		}
	}

	void SolveVariantGroup(SimpleVector<int> vars, SimpleVector<struct _adjFragmentToBaseVariant> &adjFrag, SimpleVector<struct _adjBaseVariantToFragment> &adjVar)
	{
		int i ;
		SimpleVector<char> choices ;
		int varCnt = vars.Size() ;
		struct _enumVarResult result ;
		SimpleVector<int> fragIds ;
		std::map<int, int> fragUsed ;
		std::map<int, int> seqIdxUsed ;	
		bool inExon = false;
		bool skip = false ;

		for (i = 0 ; i < varCnt ; ++i)
		{
			int seqIdx = candidateVariants[vars[i]].a ;
			int refPos = candidateVariants[vars[i]].b ;
			if (baseVariants[seqIdx][refPos].exon)
				inExon = true ;
			++seqIdxUsed[seqIdx] ;
			if (seqIdxUsed[seqIdx] > 1)
			{
				skip = true ;
				break ;
			}
		}
		if (skip || !inExon) // only compute for exons
			return ;
		
		/*for (i = 0 ; i < varCnt ; ++i)
		{
			int seqIdx = candidateVariants[vars[i]].a ;
			int refPos = candidateVariants[vars[i]].b ;
			printf("%d %d %s %d %d %c %lf %lf %lf %lf %lf %lf %lf %lf\n", i, seqIdx, refSet.GetSeqName(seqIdx), refPos, refSet.GetExonicPosition(seqIdx, refPos), refSet.GetSeqConsensus(seqIdx)[refPos],  
					baseVariants[seqIdx][refPos].count[0],
					baseVariants[seqIdx][refPos].count[1],
					baseVariants[seqIdx][refPos].count[2],
					baseVariants[seqIdx][refPos].count[3],
					baseVariants[seqIdx][refPos].uniqCount[0],
					baseVariants[seqIdx][refPos].uniqCount[1],
					baseVariants[seqIdx][refPos].uniqCount[2],
					baseVariants[seqIdx][refPos].uniqCount[3]
					) ;
		}*/
		choices.ExpandTo(varCnt) ;
		// Obtain related fragments
		for (i = 0 ; i < varCnt ; ++i)
		{
			int p = adjVar[ vars[i] ].next ;
			while (p != -1)
			{
				int fragIdx = adjVar[p].fragIdx ;
				if (fragUsed.find(fragIdx) == fragUsed.end())
				{
					int fragP = -1 ;
					/*if (varCnt > 1)
					{
						fragP = adjFrag[fragIdx].next ;
						while (fragP != -1)
						{
							int varIdx = baseVariants[adjFrag[fragP].seqIdx][adjFrag[fragP].refPos].candidateId ;
							if (candidateVariantGroupId[varIdx] != candidateVariantGroupId[vars[0]])
							{
								break ;
							}
							fragP = adjFrag[fragP].next ;
						}
					}*/
					if (fragP == -1)
					{
						fragUsed[fragIdx] = 1 ;
						fragIds.PushBack(fragIdx) ;
					}
				}
				p = adjVar[p].next ;
			}
		}
		
		result.bestCover = -1 ;	
		result.usedVarCnt = varCnt + 1 ;
		EnumerateVariants(0, choices, result, fragIds, vars, adjFrag, adjVar) ;

		// Process the final results.
		bool uniq = true ;
		if (result.equalBestEnumVariants.Size() > 0)
			uniq = false ;

		for (i = 0 ; i < varCnt ; ++i)
		{
			int seqIdx = candidateVariants[vars[i]].a ;
			int refPos = candidateVariants[vars[i]].b ;
			if (!baseVariants[seqIdx][refPos].exon)
				continue ;
			char refNuc = refSet.GetSeqConsensus(seqIdx)[refPos] ;
			char varNuc = result.bestEnumVariants[i] ;
			if (refNuc == varNuc)
				continue ;

			struct _variant nv ;
			nv.seqIdx = seqIdx ;
			nv.refStart = refPos ;
			nv.refEnd = refPos ;
			nv.ref[0] = refNuc ;
			nv.ref[1] = '\0' ;
			nv.var[0] = varNuc ;
			nv.var[1] = '\0' ;
			nv.allSupport = baseVariants[seqIdx][refPos].AllCountSum() ;
			nv.varSupport = baseVariants[seqIdx][refPos].count[ nucToNum[varNuc - 'A'] ] ;
			nv.varUniqSupport = baseVariants[seqIdx][refPos].uniqCount[ nucToNum[varNuc - 'A'] ] ;
			nv.varGroupId = candidateVariantGroupId[vars[i]] ;
			nv.outputGroupId = 0 ;
			if (uniq == false)
				nv.qual = 0 ;
			else
				nv.qual = 60 ;
			finalVariants.push_back(nv) ;
		}

		if (uniq == false)
		{
			for (i = 0 ; i < varCnt ; ++i)
			{
				int seqIdx = candidateVariants[vars[i]].a ;
				int refPos = candidateVariants[vars[i]].b ;
				if (!baseVariants[seqIdx][refPos].exon)
					continue ;
				char refNuc = refSet.GetSeqConsensus(seqIdx)[refPos] ;
				char varNuc = result.equalBestEnumVariants[i] ;
				if (refNuc == varNuc)
					continue ;

				struct _variant nv ;
				nv.seqIdx = seqIdx ;
				nv.refStart = refPos ;
				nv.refEnd = refPos ;
				nv.ref[0] = refNuc ;
				nv.ref[1] = '\0' ;
				nv.var[0] = varNuc ;
				nv.var[1] = '\0' ;
				nv.allSupport = baseVariants[seqIdx][refPos].AllCountSum() ;
				nv.varSupport = baseVariants[seqIdx][refPos].count[ nucToNum[varNuc - 'A'] ] ;
				nv.varUniqSupport = baseVariants[seqIdx][refPos].uniqCount[ nucToNum[varNuc - 'A'] ] ;
				nv.varGroupId = candidateVariantGroupId[vars[i]] ;
				nv.outputGroupId = 1 ;
				if (uniq == false)
					nv.qual = 0 ;
				else
					nv.qual = 60 ;
				finalVariants.push_back(nv) ;
			}

		}
	}

	void ComputeVariant(std::vector<char *> &read1, std::vector<char *> &read2, std::vector< std::vector<struct _fragmentOverlap> > &fragmentAssignments)
	{
		int fragCnt = fragmentAssignments.size() ;
		int seqCnt = refSet.Size() ;
		int i ;

		// Identify the preliminary set of candidate variants
		for (i = 0 ; i < fragCnt ; ++i)	
		{
			if (read2.size() > 0)
				UpdateBaseVariantFromFragmentOverlap(read1[i], read2[i], 1, fragmentAssignments[i]) ;
			else
				UpdateBaseVariantFromFragmentOverlap(read1[i], NULL, 1, fragmentAssignments[i]) ;
		}

		for (i = 0 ; i < fragCnt ; ++i)	
		{
			if (read2.size() > 0)
				UpdateBaseVariantFromFragmentOverlap(read1[i], read2[i], 0, fragmentAssignments[i]) ;
			else
				UpdateBaseVariantFromFragmentOverlap(read1[i], NULL, 0, fragmentAssignments[i]) ;
		}
		
		/*for (i = 0 ; i < seqCnt ; ++i)
		{
			int seqIdx = i;
			for (int refPos = 0 ; refPos < refSet.GetSeqConsensusLen(seqIdx) ; ++refPos )
			{
				printf("%d %s %d %d %c %lf %lf %lf %lf %lf %lf %lf %lf\n", seqIdx, refSet.GetSeqName(seqIdx), refPos, refSet.GetExonicPosition(seqIdx, refPos), refSet.GetSeqConsensus(seqIdx)[refPos],  
						baseVariants[seqIdx][refPos].count[0],
						baseVariants[seqIdx][refPos].count[1],
						baseVariants[seqIdx][refPos].count[2],
						baseVariants[seqIdx][refPos].count[3],
						baseVariants[seqIdx][refPos].uniqCount[0],
						baseVariants[seqIdx][refPos].uniqCount[1],
						baseVariants[seqIdx][refPos].uniqCount[2],
						baseVariants[seqIdx][refPos].uniqCount[3]
						) ;
			}
		}*/

		FindCandidateVariants() ;
		int candidateVarCnt = candidateVariants.Size() ;
		
		/*for (i = 0 ; i < candidateVariants.Size() ; ++i)
		{
			printf("Init: %d %s %d %d\n", i, refSet.GetSeqName(candidateVariants[i].a), candidateVariants[i].b,
					candidateVariantGroupId[i]) ;
		}*/
		// Identity the candidate variants on other sequences aligned with preliminary
		// candidate variants
		std::vector< SimpleVector<int> > seqCandidateVarAccuCount ; // useful to quickly determine whether there is a overlap of the read and candidate variations.
		SimpleVector<struct _adjBaseVariantToBaseVariant> adjVarToVar ; 
		seqCandidateVarAccuCount.resize(seqCnt) ;		
		int totalSeqLen = 0 ;
		for (i = 0 ; i < seqCnt ; ++i)
			totalSeqLen += refSet.GetSeqConsensusLen(i) ;
		adjVarToVar.Reserve(totalSeqLen + 2 * candidateVarCnt) ;
		adjVarToVar.Resize(totalSeqLen) ; // the first seqcount are used for the head nodes 

		for (i = 0 ; i < totalSeqLen ; ++i)
		{
			adjVarToVar[i].varIdx = -1 ;
			adjVarToVar[i].next = -1 ;
		}
		for (i = 0 ; i < candidateVarCnt ; ++i)
		{
			adjVarToVar[i].rootCandidate = true ;
			adjVarToVar[i].varIdx = i ;
			adjVarToVar[i].next = -1 ;
		}

		while (1) 
		{	
			int prevCandidateVarCnt = candidateVariants.Size() ;
			for (i = 0 ; i < prevCandidateVarCnt ; ++i) // reset the graph to recalculate the weight.
				adjVarToVar[i].next = -1 ;
			adjVarToVar.Resize(totalSeqLen) ;
		
			for (i = 0 ; i < seqCnt ; ++i)
				ComputeCandidateVarAccuCount(i, seqCandidateVarAccuCount[i]) ;

			for (i = 0 ; i < fragCnt ; ++i)
			{
				if (read2.size() > 0)
					ExpandCandidateVariantsFromFragmentOverlap(read1[i], read2[i], 
							fragmentAssignments[i], adjVarToVar, seqCandidateVarAccuCount ) ;
				else
					ExpandCandidateVariantsFromFragmentOverlap(read1[i], NULL, 
							fragmentAssignments[i], adjVarToVar, seqCandidateVarAccuCount ) ;
			}
			if (prevCandidateVarCnt == candidateVariants.Size())
				break ;
		}
		
		
		//for (i = 0 ; i < seqCnt ; ++i)
		//	ComputeCandidateVarAccuCount(i, seqCandidateVarAccuCount[i]) ;

		// Build the relation of fragments and variants
		candidateVarCnt = candidateVariants.Size() ;
		//struct _adjFragmentToBaseVariant *adjFrag = (struct _adjFragmentToBaseVariant*)malloc(sizeof(struct _adjFragmentToBaseVariant) * fragCnt) ;
		//struct _adjBaseVariantToFragment *adjVar = (struct _adjBaseVariantToFragment*)malloc(sizeof(_adjBaseVariantToFragment) * candidateVarCnt) ;
		int groupCnt = 0 ;
		for (i = 0 ; i < candidateVarCnt ; ++i)
		{
			if (adjVarToVar[i].rootCandidate && candidateVariantGroupId[i] == -1)
			{
				BuildCandidateVariantGroup(i, groupCnt, adjVarToVar) ;
				++groupCnt ;
			}
		}

		SimpleVector<struct _adjFragmentToBaseVariant> adjFrag ; 
		SimpleVector<struct _adjBaseVariantToFragment> adjVar ; 
		adjFrag.Reserve(fragCnt + candidateVarCnt) ;
		adjFrag.Resize(fragCnt) ;
		adjVar.Reserve(fragCnt + candidateVarCnt) ;
		adjVar.Resize(candidateVarCnt) ;

		/*for (i = 0 ; i < candidateVarCnt ; ++i)
		{
			printf("Expand: %d %d %s %d %d %d\n", i, candidateVariants[i].a, refSet.GetSeqName(candidateVariants[i].a), candidateVariants[i].b, adjVarToVar[i].rootCandidate, candidateVariantGroupId[i]) ;
		}*/
		
		for (i = 0 ; i < fragCnt ; ++i)
			adjFrag[i].next = -1 ;
		for (i = 0 ; i < candidateVarCnt ; ++i)
			adjVar[i].next = -1 ;

		for (i = 0 ; i < fragCnt ; ++i)
		{
			if (read2.size() > 0)
				BuildFragmentCandidateVarGraph(read1[i], read2[i], i, 
						fragmentAssignments[i], seqCandidateVarAccuCount, adjFrag, adjVar) ;	
			else
				BuildFragmentCandidateVarGraph(read1[i], NULL, i, 
						fragmentAssignments[i], seqCandidateVarAccuCount, adjFrag, adjVar) ;	
		}
		
		// Group variants
		std::vector< SimpleVector<int> > candidateVarGroup ;
		candidateVarGroup.resize(groupCnt) ;
		for (i = 0 ; i < candidateVarCnt ; ++i)
		{
			int gid = candidateVariantGroupId[i] ;
			if (gid == -1)
				continue ;
			candidateVarGroup[gid].PushBack(i) ;
		}

		// Solve each group
		int reducedGroupCnt = candidateVarGroup.size() ;
		for ( i = 0 ; i < reducedGroupCnt ; ++i)
		{
			SolveVariantGroup(candidateVarGroup[i], adjFrag, adjVar) ;
			//break ;
		}

		int finalVarCnt = finalVariants.size() ;
		for (i = 0 ; i < finalVarCnt ; ++i)
		{
			baseVariants[finalVariants[i].seqIdx][finalVariants[i].refStart].finalVariantIds.push_back(i) ;
		}
	}

	/*int GetSeqExonVariants(int seqIdx, std::vector<struct _variant> &variants)
	{
		int i, j, k ;
		const struct _seqWrapper &seq = seqs[seqIdx] ;
		k = 0 ;
		for (i = 0 ; i < seq.consensusLen ; ++i)
		{
			// major allele 
			if (!seq.isValidDiff[i].exon)
				continue ;
			double max = seq.baseVariants[i].count[0];
			int maxTag = 0 ;
			for (j = 1 ; j < 4 ; ++j)
			{
				if (seq.baseVariants[i].count[j] > max)
				{
					max = seq.baseVariants[i].count[j] ;
					maxTag = j ;
				}
			}
			//if (max == 0)
			//	printf("%s %d\n", seqs[seqIdx].name, i) ;
			if (numToNuc[maxTag] != seq.consensus[i] && max > 0)
			{
				// Variation happens
				struct _variant nv ;
				nv.seqIdx = seqIdx ;
				nv.refStart = k ;
				nv.refEnd = k ;
				nv.ref[0] = seq.consensus[i] ;
				nv.ref[1] = '\0' ;
				nv.var[0] = numToNuc[maxTag] ;
				nv.var[1] = '\0' ;
				nv.allSupport = seq.baseVariants[i].AllCountSum() ;
				nv.varSupport = max ;
				nv.varUniqSupport = seq.baseVariants[i].uniqCount[maxTag] ;
				variants.push_back(nv) ;
			}
			++k ;
		}
	}*/

	void ConvertVariantsToExonCoord()
	{
		return ;
		int i ;
		int varCnt = finalVariants.size() ;
		for (i = 0 ; i < varCnt ; ++i)
		{
			struct _variant &variant = finalVariants[i] ;
			variant.refStart = refSet.GetExonicPosition(variant.seqIdx, variant.refStart) ;
			variant.refEnd = refSet.GetExonicPosition(variant.seqIdx, variant.refEnd) ;
		}		
	}
	
	void OutputAlleleVCF(char *filename)
	{
		FILE *fp = fopen(filename, "w") ;
		int i ;
		int varCnt = finalVariants.size() ;	
		char buffer[10] = "PASS" ;	
		for (i = 0 ; i < varCnt ; ++i)
		{
			struct _variant &variant = finalVariants[i] ;
			if (variant.qual > 0)
				strcpy(buffer, "PASS") ;
			else
				strcpy(buffer, "FAIL") ;
			int exonRefStart = refSet.GetExonicPosition(variant.seqIdx, variant.refStart) ;
      char tmp[23] = "";
      tmp[22] = '\0'; 
      memcpy(tmp, refSet.GetSeqConsensus(variant.seqIdx) + variant.refStart - 10, 21);
			fprintf(fp, "%s %d . %s %s . %s %lf %lf %lf %d %d\n", 
					refSet.GetSeqName(variant.seqIdx), exonRefStart + 1, // the VCF file is 1-based
					variant.ref, variant.var, buffer, 
					variant.varSupport, variant.allSupport,
					variant.varUniqSupport, variant.refStart, variant.outputGroupId);
					//, tmp) ;
		}
		fclose(fp) ;
	}

	std::vector<struct _fragmentOverlap> AdjustFragmentAssignment(char *read1, char *read2, std::vector<struct _fragmentOverlap> &rawAssignments)
	{
		int i, j, k ;
		int assignCnt = rawAssignments.size() ;
		SimpleVector<double> changeScore ;
		changeScore.ExpandTo(assignCnt) ;
		changeScore.SetZero(0, assignCnt) ; 
		
		for (i = 0 ; i < assignCnt ; ++i)
		{
			for (k = 0 ; k < 2 ; ++k)
			{
				if (k == 1 && !rawAssignments[i].hasMatePair)
					continue ;
				char *read = read1 ;
				if (k == 1 || (k == 0 && rawAssignments[i].o1FromR2))
					read = read2 ;

				struct _overlap o = rawAssignments[i].overlap1 ;
				if (k == 1)
					o = rawAssignments[i].overlap2 ;
				if (o.align == NULL)
					continue ;

				char *r = read ;
				char *rc = NULL ;
				int readLen = strlen(r) ;
				if (o.strand == -1)
				{
					rc = strdup(read) ;
					refSet.ReverseComplement(rc, read, readLen) ;
					r = rc ;
				}
				
				char *align = o.align ;
				int readPos = o.readStart ;
				int refPos = o.seqStart ;
				int seqIdx = o.seqIdx ;
				for (j = 0 ; align[j] != -1 ; ++j)
				{
					if (align[j] == EDIT_MISMATCH)
					{
						int size = baseVariants[seqIdx][refPos].finalVariantIds.size() ;
						if (size > 0)
						{
							int l ;
							for (l = 0 ; l < size ; ++l)
							{
								int vid = baseVariants[seqIdx][refPos].finalVariantIds[l] ;
								if (finalVariants[vid].var[0] == r[readPos])
								{
									++changeScore[i] ;
									break ;
								}
							}
						}
					}
					
					if (align[j] != EDIT_INSERT)
						++refPos ;
					if (align[j] != EDIT_DELETE)
						++readPos ;
				}

				if (rc != NULL)
					free(rc) ;
			}
		}
		
		double maxScore = -1 ;
		for (i = 0 ; i < assignCnt ; ++i)
			if (changeScore[i] > maxScore)
				maxScore = changeScore[i] ;
		
		std::vector<struct _fragmentOverlap> ret ;
		for (i = 0 ; i < assignCnt ; ++i)
		{
			if (changeScore[i] == maxScore)
				ret.push_back( rawAssignments[i] ) ;
		}

		return ret ;
	}
} ;

#endif

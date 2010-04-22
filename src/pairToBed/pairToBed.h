/*****************************************************************************
  pairToBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "BamReader.h"
#include "BamWriter.h"
#include "BamAux.h"
using namespace BamTools;

#include "bedFile.h"
#include "bedFilePE.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;



/**************************************************
Helper function protoypes
**************************************************/
void IsCorrectMappingForBEDPE (const BamAlignment &bam, const RefVector &refs, BEDPE &a);



//************************************************
// Class methods and elements
//************************************************
class BedIntersectPE {

public:

	// constructor 
	BedIntersectPE(string bedAFilePE, string bedBFile, float overlapFraction, 
		string searchType, bool forceStrand, bool bamInput, bool bamOutput);

	// destructor
	~BedIntersectPE(void);

	void FindOverlaps(const BEDPE &, vector<BED> &hits1, vector<BED> &hits2, const string &type); 
	
	bool FindOneOrMoreOverlaps(const BEDPE &, const string &type); 

	void FindSpanningOverlaps(const BEDPE &a, vector<BED> &hits, const string &type); 
	bool FindOneOrMoreSpanningOverlaps(const BEDPE &a, const string &type);

	void IntersectBedPE(istream &bedInput);
	void IntersectBamPE(string bamFile);
	
	void DetermineBedPEInput();
	
private:

	string _bedAFilePE;
	string _bedBFile;
	float _overlapFraction;
	string _searchType;
	bool _forceStrand;
	bool _bamInput;
	bool _bamOutput;

	// instance of a paired-end bed file class.
	BedFilePE *_bedA;

	// instance of a bed file class.
	BedFile *_bedB;	

	inline 
	void ConvertBamToBedPE(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, BEDPE &a) {

		a.start1 = a.start2 = a.end1 = a.end2 = -1;
		a.chrom1 = a.chrom2 = a.strand1 = a.strand2 = ".";
		
		// end 1
		if (bam1.IsMapped()) {
			a.chrom1 = refs.at(bam1.RefID).RefName;
			a.start1 = bam1.Position;
			a.end1 = bam1.GetEndPosition();
			//a.end1 = bam.Position + bam.AlignedBases.size();
			a.strand1 = "+";
			if (bam1.IsReverseStrand()) a.strand1 = "-";	
		}
		
		// end 2
		if (bam2.IsMapped()) {
			a.chrom2 = refs.at(bam2.RefID).RefName;
			a.start2 = bam2.Position;
			a.end2 = bam2.GetEndPosition();
			//a.end1 = bam.Position + bam.AlignedBases.size();
			a.strand2 = "+";
			if (bam2.IsReverseStrand()) a.strand2 = "-";	
		}
		
		a.name = bam1.Name;		
		uint8_t edit1, edit2;
		bam1.GetEditDistance(edit1);
		bam2.GetEditDistance(edit2);
		a.score = edit1 + edit2;
	};
	
	inline
	void ProcessBamBlock (const vector<BamAlignment> &alignments, 
	                                      const RefVector &refs,
	                                      BamWriter &writer);
};

#endif /* PEINTERSECTBED_H */

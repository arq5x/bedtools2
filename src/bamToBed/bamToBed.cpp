/*****************************************************************************
  bamToBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "bedFile.h"
#include "version.h"
using namespace BamTools;

#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools bamtobed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void bamtobed_help(void);

void ConvertBamToBed(const string &bamFile, bool useEditDistance, const string &bamTag,
                     bool writeBed12, bool obeySplits, const string &color, 
                     bool useCigar,   bool useNovoalign, bool useBWA);
                     
void ConvertBamToBedpe(const string &bamFile, const bool &useEditDistance);

void PrintBed(const BamAlignment &bam, const RefVector &refs, bool useEditDistance, 
              const string &bamTag, bool obeySplits, bool useCigar, bool useNovoalign, bool useBWA);
              
void PrintBed12(const BamAlignment &bam, const RefVector &refs, bool useEditDistance, const string &bamTag, string color = "255,0,0");
void PrintBedPE(const BamAlignment &bam1, const BamAlignment &bam2,
                const RefVector &refs, bool useEditDistance);

void ParseCigarBed12(const vector<CigarOp> &cigar, vector<int> &blockStarts,
                     vector<int> &blockEnds, int &alignmentEnd);
string BuildCigarString(const vector<CigarOp> &cigar);

bool bamtobed_IsCorrectMappingForBEDPE (const BamAlignment &bam);


int bamtobed_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bamFile = "stdin";
    string color   = "255,0,0";
    string tag     = "";

    bool haveBam           = true;
    bool haveColor         = false;
    bool haveOtherTag      = false;
    bool writeBedPE        = false;
    bool writeBed12        = false;
    bool useEditDistance   = false;
    bool useAlignmentScore = false;
    bool useCigar          = false;
    bool obeySplits        = false;
    bool useNovoalign      = false;  // custom for Quinlan/Hall research
    bool useBWA            = false;  // custom for Quinlan/Hall research

    // check to see if we should print out some help

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bamtobed_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bamFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bedpe", 6, parameterLength)) {
                writeBedPE = true;
        }
        else if(PARAMETER_CHECK("-bed12", 6, parameterLength)) {
                writeBed12 = true;
        }
        else if(PARAMETER_CHECK("-split", 6, parameterLength)) {
                obeySplits = true;
        }
        else if(PARAMETER_CHECK("-ed", 3, parameterLength)) {
                useEditDistance = true;
        }
        else if(PARAMETER_CHECK("-cigar", 6, parameterLength)) {
                useCigar = true;
        }
        else if(PARAMETER_CHECK("-as", 3, parameterLength)) {
                useAlignmentScore = true;
        }
        else if(PARAMETER_CHECK("-novo", 5, parameterLength)) {
                useNovoalign = true;
        }
        else if(PARAMETER_CHECK("-bwa", 4, parameterLength)) {
                useBWA = true;
        }
        else if(PARAMETER_CHECK("-color", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveColor = true;
                color = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-tag", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveOtherTag = true;
                tag = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (haveBam == false) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveColor == true && writeBed12 == false) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot use color without BED12. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && obeySplits == true) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot use -ed with -splits.  Edit distances cannot be computed for each \'chunk\'." << endl << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && useCigar == true) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot use -cigar with -splits.  Not yet supported." << endl << "*****" << endl;
        showHelp = true;
    }
    if (useEditDistance == true && haveOtherTag == true) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot use -ed with -tag.  Choose one or the other." << endl << "*****" << endl;
        showHelp = true;
    }
    if (writeBedPE == true && haveOtherTag == true) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot use -tag with -bedpe." << endl << "*****" << endl;
        showHelp = true;
    }
    // if there are no problems, let's convert BAM to BED or BEDPE
    if (!showHelp) {
        if (writeBedPE == false)
            ConvertBamToBed(bamFile, useEditDistance, tag, writeBed12, obeySplits, color, useCigar, useNovoalign, useBWA);    // BED or "blocked BED"
        else
            ConvertBamToBedpe(bamFile, useEditDistance);                    // BEDPE
    }
    else {
        bamtobed_help();
    }
    return 0;
}


void bamtobed_help(void) {
    
    cerr << "\nTool:    bedtools bamtobed (aka bamToBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Converts BAM alignments to BED6 or BEDPE format." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-bedpe\t"    << "Write BEDPE format." << endl;
    cerr                    << "\t\t- Requires BAM to be grouped or sorted by query." << endl << endl;

    cerr << "\t-bed12\t"    << "Write \"blocked\" BED format (aka \"BED12\")." << endl << endl;
    cerr                    << "\t\thttp://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1" << endl << endl;

    cerr << "\t-split\t"    << "Report \"split\" BAM alignments as separate BED entries." << endl << endl;

    cerr << "\t-ed\t"       << "Use BAM edit distance (NM tag) for BED score." << endl;
    cerr                    << "\t\t- Default for BED is to use mapping quality." << endl;
    cerr                    << "\t\t- Default for BEDPE is to use the minimum of" << endl;
    cerr                    << "\t\t  the two mapping qualities for the pair." << endl;
    cerr                    << "\t\t- When -ed is used with -bedpe, the total edit" << endl;
    cerr                    << "\t\t  distance from the two mates is reported." << endl << endl;

    cerr << "\t-tag\t"      << "Use other NUMERIC BAM alignment tag for BED score." << endl;
    cerr                    << "\t\t- Default for BED is to use mapping quality." << endl;
    cerr                    << "\t\t  Disallowed with BEDPE output." << endl << endl;

    cerr << "\t-color\t"    << "An R,G,B string for the color used with BED12 format." << endl;
    cerr                    << "\t\tDefault is (255,0,0)." << endl << endl;

    cerr << "\t-cigar\t"    << "Add the CIGAR string to the BED entry as a 7th column." << endl << endl;


    // end the program here
    exit(1);
}


void ConvertBamToBed(const string &bamFile, bool useEditDistance, const string &bamTag,
                     bool writeBed12, bool obeySplits, const string &color, 
                     bool useCigar,   bool useNovoalign, bool useBWA) {
    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // rip through the BAM file and convert each mapped entry to BED
    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {
        if (bam.IsMapped() == true) {
            if (writeBed12 == false)                 // BED
                PrintBed(bam, refs, useEditDistance, bamTag, obeySplits, useCigar, useNovoalign, useBWA);
            else                                     //"blocked" BED
                PrintBed12(bam, refs, useEditDistance, bamTag, color);
        }
    }
    reader.Close();
}


/*
  Assumptions:
     1.  The BAM file is grouped/sorted by query name,
         not alignment position
*/
void ConvertBamToBedpe(const string &bamFile, const bool &useEditDistance) {
    // open the BAM file
    BamReader reader;
    reader.Open(bamFile);

    // get header & reference information
    string header = reader.GetHeaderText();
    RefVector refs = reader.GetReferenceData();

    // rip through the BAM file and convert each mapped entry to BEDPE
    BamAlignment bam1, bam2;
    while (reader.GetNextAlignment(bam1)) {
        
        reader.GetNextAlignment(bam2);        
        if (bam1.Name != bam2.Name) {
            while (bam1.Name != bam2.Name)
            {
                if (bam1.IsPaired()) 
                {
                    cerr << "*****WARNING: Query " << bam1.Name
                         << " is marked as paired, but it's mate does not occur"
                         << " next to it in your BAM file.  Skipping. " << endl;
                }
                bam1 = bam2;
                reader.GetNextAlignment(bam2);
            }
        }
        else if (bam1.IsPaired() && bam2.IsPaired()) {
            PrintBedPE(bam1, bam2, refs, useEditDistance);
        }
    }
    reader.Close();
}


string BuildCigarString(const vector<CigarOp> &cigar) {

    stringstream cigarString;
    for (size_t i = 0; i < cigar.size(); ++i) {
        switch (cigar[i].Type) {
            case ('M') :
            case ('I') :
            case ('D') :
            case ('N') :
            case ('S') :
            case ('H') :
            case ('P') :
                cigarString << cigar[i].Length << cigar[i].Type;
        }
    }
    return cigarString.str();
}


void PrintBed(const BamAlignment &bam,  const RefVector &refs, bool useEditDistance, 
              const string &bamTag, bool obeySplits, bool useCigar, bool useNovoalign, bool useBWA) {

    // set the strand
    string strand = "+";
    if (bam.IsReverseStrand() == true) strand = "-";

    // set the name of the feature based on the sequence
    string name = bam.Name;
    if (bam.IsFirstMate() == true)  name += "/1";
    if (bam.IsSecondMate() == true) name += "/2";

    // get the unpadded (parm = false) end position based on the CIGAR
    unsigned int alignmentEnd = bam.GetEndPosition(false, false);

    // report the entire BAM footprint as a single BED entry
    if (obeySplits == false) {
        
        if (!useNovoalign && !useBWA) {
            // report the alignment in BED6 format.
            if (useEditDistance == false && bamTag == "") {
                printf("%s\t%d\t%d\t\%s\t%d\t%s", refs.at(bam.RefID).RefName.c_str(), bam.Position,
                                              alignmentEnd, name.c_str(), bam.MapQuality, strand.c_str());
            }
            else if (useEditDistance == true && bamTag == "") {
                uint32_t editDistance;
                if (bam.GetTag("NM", editDistance)) {
                    printf("%s\t%d\t%d\t\%s\t%u\t%s", refs.at(bam.RefID).RefName.c_str(), bam.Position,
                                                  alignmentEnd, name.c_str(), editDistance, strand.c_str());
                }
                else {
                    cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
                    exit(1);
                }
            }
            else if (useEditDistance == false && bamTag != "") {
                int32_t tagValue;
                if (bam.GetTag(bamTag, tagValue)) {
                    printf("%s\t%d\t%d\t\%s\t%d\t%s", refs.at(bam.RefID).RefName.c_str(), bam.Position,
                                                  alignmentEnd, name.c_str(), tagValue, strand.c_str());
                }
                else {
                    cerr << "The requested tag (" << bamTag << ") was not found in the BAM file.  Exiting\n";
                    exit(1);
                }
            }

            // does the user want CIGAR as well?
            if (useCigar == false) {
                printf("\n");
            }
            else {
                string cigar = BuildCigarString(bam.CigarData);
                printf("\t%s\n", cigar.c_str());
            }
        }
        else if (useNovoalign && !useBWA) {
            // special BED format for Hydra using Novoalign.
            uint32_t editDistance;
            uint32_t numMappings;
            
            if (!bam.GetTag("NM", editDistance)) {
                cerr << "Unable to extract NM for query: " 
                     << bam.Name
                     << ".  Verify that your BAM was generated by Novoalign." << endl;
                exit(1);
            }
            if (!bam.GetTag("ZN", numMappings)) {
                // if ZN is missing, this means just one alignment was found.
                numMappings = 1;
            }
            printf("%s\t%d\t%d\t\%s\t%u_%u_%u\t%s\n", 
                   refs.at(bam.RefID).RefName.c_str(), bam.Position,
                   alignmentEnd, name.c_str(), 
                   bam.MapQuality, editDistance, numMappings, strand.c_str());
        }
        else if (!useNovoalign && useBWA) {
            // special BED format for Hydra using Novoalign.
            uint32_t editDistance;
            uint32_t numMappings;
            uint32_t x0 = 1;
            uint32_t x1 = 0;
            
            if (!bam.GetTag("NM", editDistance)) {
                cerr << "Unable to extract NM.  Verify that your BAM was generated by Novoalign." << endl;
                exit(1);
            }
            if (!bam.GetTag("X0", x0) && !bam.GetTag("X1", x1)) {
                // if ZN is missing, this means just one alignment was found.
                numMappings = 1;
            }
            else {
                numMappings = x0 + x1;
            }
            printf("%s\t%d\t%d\t\%s\t%u_%u_%u\t%s\n", 
                   refs.at(bam.RefID).RefName.c_str(), bam.Position,
                   alignmentEnd, name.c_str(), 
                   bam.MapQuality, editDistance, numMappings, strand.c_str());
        }
    }
    // Report each chunk of the BAM alignment as a discrete BED entry
    // For example 10M100N10M would be reported as two seprate BED entries of length 10
    else {
        // parse the CIGAR string and figure out the alignment blocks
        vector<BED> bedBlocks;
        string chrom = refs.at(bam.RefID).RefName;
        // extract the block starts and lengths from the CIGAR string
        GetBamBlocks(bam, chrom, bedBlocks, false, true);

        unsigned int i;
        for (i = 0; i < bedBlocks.size(); ++i) {
            BED curr = bedBlocks[i];
            printf("%s\t%d\t%d\t\%s\t%d\t%s\n", 
                   chrom.c_str(), curr.start, curr.end,
                   name.c_str(), bam.MapQuality, strand.c_str());
        }
    }
}


void PrintBed12(const BamAlignment &bam, const RefVector &refs, bool useEditDistance, const string &bamTag, string color) {

    // set the strand
    string strand = "+";
    if (bam.IsReverseStrand()) strand = "-";

    // set the name of the feature based on the sequence
    string name = bam.Name;
    if (bam.IsFirstMate()) name += "/1";
    if (bam.IsSecondMate()) name += "/2";

    // parse the CIGAR string and figure out the alignment blocks
    vector<BED> bedBlocks;
    string chrom = refs.at(bam.RefID).RefName;
    CHRPOS alignmentEnd = bam.GetEndPosition();
    // extract the block starts and lengths from the CIGAR string
    GetBamBlocks(bam, chrom, bedBlocks, false, true);

    // write BED6 portion
    if (useEditDistance == false && bamTag == "") {
        printf("%s\t%d\t%d\t\%s\t%d\t%s\t", refs.at(bam.RefID).RefName.c_str(), bam.Position,
            alignmentEnd, name.c_str(), bam.MapQuality, strand.c_str());
    }
    else if (useEditDistance == true && bamTag != "") {
        uint32_t editDistance;
        if (bam.GetTag("NM", editDistance)) {
            printf("%s\t%d\t%d\t\%s\t%u\t%s\t", refs.at(bam.RefID).RefName.c_str(), bam.Position,
                alignmentEnd, name.c_str(), editDistance, strand.c_str());
        }
        else {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }
    else if (useEditDistance == false && bamTag != "") {
        int32_t tagValue;
        if (bam.GetTag(bamTag, tagValue)) {
            printf("%s\t%d\t%d\t\%s\t%d\t%s\t", refs.at(bam.RefID).RefName.c_str(), bam.Position,
                                          alignmentEnd, name.c_str(), tagValue, strand.c_str());
        }
        else {
            cerr << "The requested tag (" << bamTag << ") was not found in the BAM file.  Exiting\n";
            exit(1);
        }
    }

    // write the colors, etc.
    printf("%d\t%d\t%s\t%d\t", bam.Position, alignmentEnd, color.c_str(), (int) bedBlocks.size());

    // now write the lengths portion
    unsigned int b;
    for (b = 0; b < bedBlocks.size() - 1; ++b) {
        printf("%d,", bedBlocks[b].end - bedBlocks[b].start);
    }
    printf("%d\t", bedBlocks[b].end - bedBlocks[b].start);

    // now write the starts portion
    for (b = 0; b < bedBlocks.size() - 1; ++b) {
        printf("%d,", bedBlocks[b].start - bam.Position);
    }
    printf("%d\n", bedBlocks[b].start - bam.Position);
}


void PrintBedPE(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, bool useEditDistance) {

    // initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int start1, start2, end1, end2;
    uint32_t editDistance1, editDistance2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    editDistance1 = editDistance2 = 0;
    uint16_t minMapQuality = 0;

    // extract relevant info for end 1
    if (bam1.IsMapped()) {
        chrom1 = refs.at(bam1.RefID).RefName;
        start1 = bam1.Position;
        end1   = bam1.GetEndPosition(false);
        strand1 = "+";
        if (bam1.IsReverseStrand()) strand1 = "-";

        // extract the edit distance from the NM tag
        // if possible. otherwise, complain.
        if (useEditDistance == true && bam1.GetTag("NM", editDistance1) == false) {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }

    // extract relevant info for end 2
    if (bam2.IsMapped()) {
        chrom2 = refs.at(bam2.RefID).RefName;
        start2 = bam2.Position;
        end2   = bam2.GetEndPosition(false);
        strand2 = "+";
        if (bam2.IsReverseStrand()) strand2 = "-";

        // extract the edit distance from the NM tag
        // if possible. otherwise, complain.
        if (useEditDistance == true && bam2.GetTag("NM", editDistance2) == false) {
            cerr << "The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting\n";
            exit(1);
        }
    }

    // swap the ends if necessary
    if ( chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)) ) {
        swap(chrom1, chrom2);
        swap(start1, start2);
        swap(end1, end2);
        swap(strand1, strand2);
    }

    // report BEDPE using min mapQuality
    if (useEditDistance == false) {
        // compute the minimum mapping quality b/w the two ends of the pair.
        if (bam1.IsMapped() == true && bam2.IsMapped() == true)
            minMapQuality = min(bam1.MapQuality, bam2.MapQuality);

        printf("%s\t%d\t%d\t\%s\t%d\t%d\t%s\t%d\t%s\t%s\n",
                chrom1.c_str(), start1, end1, chrom2.c_str(), start2, end2,
                bam1.Name.c_str(), minMapQuality, strand1.c_str(), strand2.c_str());
    }
    // report BEDPE using total edit distance
    else {
        uint16_t totalEditDistance = 0;
        if (bam1.IsMapped() == true && bam2.IsMapped() == true)
            totalEditDistance = editDistance1 + editDistance2;
        else if (bam1.IsMapped() == true)
            totalEditDistance = editDistance1;
        else if (bam2.IsMapped() == true)
            totalEditDistance = editDistance2;

        printf("%s\t%d\t%d\t\%s\t%d\t%d\t%s\t%d\t%s\t%s\n",
                chrom1.c_str(), start1, end1, chrom2.c_str(), start2, end2,
                bam1.Name.c_str(), totalEditDistance, strand1.c_str(), strand2.c_str());
    }
}


// deprecated.
bool bamtobed_IsCorrectMappingForBEDPE (const BamAlignment &bam) {

    if ( (bam.RefID == bam.MateRefID) && (bam.InsertSize > 0) ) {
        return true;
    }
    else if ( (bam.RefID == bam.MateRefID) && (bam.InsertSize == 0) && bam.IsFirstMate() ) {
        return true;
    }
    else if ( (bam.RefID != bam.MateRefID) && bam.IsFirstMate() ) {
        return true;
    }
    else return false;
}

/*****************************************************************************
  bedpeToBam.cpp

  (c) 2009 - Royden Clark, Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFilePE.h"
#include "GenomeFile.h"
#include "version.h"


#include "api/BamReader.h"
#include "api/BamAux.h"
#include "api/BamWriter.h"
using namespace BamTools;

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedpetobam"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void bedpetobam_help(void);
void ProcessBedPE(BedFilePE *bedpe, GenomeFile *genome,  int mapQual, bool uncompressedBam);
void ConvertBedPEToBam(const BEDPE &bedpe, BamAlignment &bam1,BamAlignment &bam2, map<string, int> &chromToId, int mapQual, int lineNum);

void bedpetobam_MakeBamHeader(const string &genomeFile, RefVector &refs, string &header, map<string, int> &chromToInt);
int  bedpetobam_reg2bin(int beg, int end);



int bedpetobam_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedpeFile = "stdin";
    string genomeFile;

    int mapQual = 255;

    bool haveBedPE         = true;
    bool haveGenome      = false;
    bool haveMapQual     = false;
   // bool isBED12         = false;
    bool uncompressedBam = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bedpetobam_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedpeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveGenome = true;
                genomeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-mapq", 5, parameterLength)) {
            haveMapQual = true;
            if ((i+1) < argc) {
                mapQual = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ubam", 5, parameterLength)) {
            uncompressedBam = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (!haveBedPE ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i (BEDPE) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (!haveGenome ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -g (genome) file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (mapQual < 0 || mapQual > 255) {
        cerr << endl << "*****" << endl << "*****ERROR: MAPQ must be in range [0,255]. " << endl << "*****" << endl;
        showHelp = true;
    }


    if (!showHelp) {
        BedFilePE *bedpe= new BedFilePE(bedpeFile);
        GenomeFile *genome = new GenomeFile(genomeFile);

       ProcessBedPE(bedpe, genome,  mapQual, uncompressedBam);
    }
    else {
        bedpetobam_help();
    }
    return 0;
}


void bedpetobam_help(void) {
    
    cerr << "\nTool:    bedtools bedpetobam (aka bedpeToBam)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Converts feature records to BAM format." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-mapq\t" << "Set the mappinq quality for the BAM records." << endl;
    cerr                    << "\t\t(INT) Default: 255" << endl << endl;

    cerr << "\t-ubam\t"     << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  BED files must be at least BED4 to create BAM (needs name field)." << endl << endl;


    // end the program here
    exit(1);
}

void ProcessBedPE(BedFilePE *bedpe, GenomeFile *genome,  int mapQual, bool uncompressedBam) {

    BamWriter *writer = new BamWriter();

    // build a BAM header from the genomeFile
    RefVector refs;
    string    bamHeader;
    map<string, int, std::less<string> > chromToId;
    bedpetobam_MakeBamHeader(genome->getGenomeFileName(), refs, bamHeader, chromToId);

    // set compression mode
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( uncompressedBam ) compressionMode = BamWriter::Uncompressed;
    writer->SetCompressionMode(compressionMode);
    // open a BAM and add the reference headers to the BAM file
    writer->Open("stdout", bamHeader, refs);


    // process each BED entry and convert to BAM
    BEDPE bedpeEntry, nullBedpe;
    int lineNum = 0;
    BedLineStatus bedpeStatus;
    // open the BED file for reading.
    bedpe->Open();
    while ((bedpeStatus = bedpe->GetNextBedPE(bedpeEntry, lineNum)) != BED_INVALID) {
        if (bedpeStatus == BED_VALID) {
            BamAlignment bamEntry1;
           	BamAlignment bamEntry2;

            if (bedpe->bedType >= 10) {
                ConvertBedPEToBam(bedpeEntry, bamEntry1, bamEntry2, chromToId,  mapQual, lineNum);
                writer->SaveAlignment(bamEntry1, true);
                writer->SaveAlignment(bamEntry2, true);

            }
            else {
                cerr << "Error: BEDPE entry without name found at line: " << lineNum << ".  Exiting!" << endl;
                exit (1);
            }
            bedpeEntry = nullBedpe;
        }
    }
    //close up
    bedpe->Close();
    writer->Close();
}


//void ConvertBedPEToBam(const BEDPE &bedpe, BamAlignment &bam1,BamAlignment &bam2, map<string, int, std::less<string> > &chromToId,
 //                    bool isBED12, int mapQual, int lineNum) {
void ConvertBedPEToBam(const BEDPE &bedpe, BamAlignment &bam1,BamAlignment &bam2, map<string, int, std::less<string> > &chromToId,
	                    int mapQual, int lineNum) {

    bam1.Name       = bedpe.name;
    bam1.Position   = bedpe.start1;
    bam1.Bin        = bedpetobam_reg2bin(bedpe.start1, bedpe.end1);
    bam2.Name       = bedpe.name;
    bam2.Position   = bedpe.start2;
    bam2.Bin        = bedpetobam_reg2bin(bedpe.start2, bedpe.end2);

    // hard-code the sequence and qualities.
    int bedpeLength1  = bedpe.end1 - bedpe.start1;
    int bedpeLength2  = bedpe.end2 - bedpe.start2;

    // set dummy seq and qual strings.  the input is BED,
    // so the sequence is inherently the same as it's
    // reference genome.
    // Thanks to James M. Ward for pointing this out.
    bam1.QueryBases = "";
    bam1.Qualities  = "";
 	bam2.QueryBases = "";
    bam2.Qualities  = "";

    // chrom and map quality
    bam1.RefID      = chromToId[bedpe.chrom1];
    bam1.MapQuality = mapQual;
    bam2.RefID      = chromToId[bedpe.chrom2];
    bam2.MapQuality = mapQual;

    // set the BAM FLAG
    bam1.AlignmentFlag = 0;
    bam2.AlignmentFlag = 0;
 
   if (bedpe.strand1 == "-"){
        bam1.SetIsReverseStrand(true);
		bam2.SetIsMateReverseStrand(true);
		
	}
	
   if (bedpe.strand2 == "-"){
        bam2.SetIsReverseStrand(true);
		bam1.SetIsMateReverseStrand(true);
	}
    bam1.MatePosition = bedpe.start2;
    
	if(chromToId[bedpe.chrom1] == chromToId[bedpe.chrom2]){
		bam1.InsertSize   = bedpe.start2-bedpe.start1;
		bam2.InsertSize   = bedpe.start2-bedpe.start1;
		if((bedpe.strand1 == "+") && (bedpe.strand2 == "-")){
			bam1.SetIsProperPair(true);
			bam2.SetIsProperPair(true);
		}
	}
	else{
		bam1.InsertSize   = 0;
		bam2.InsertSize   = 0;
		
	}
    bam1.MateRefID    = chromToId[bedpe.chrom2];
    bam2.MatePosition = bedpe.start1;
    bam2.MateRefID    = chromToId[bedpe.chrom1];
	
	bam1.SetIsFirstMate(true);
	bam2.SetIsSecondMate(true);
	bam1.SetIsPaired(true);
	bam2.SetIsPaired(true);
	
    bam1.CigarData.clear();
	bam2.CigarData.clear();

    CigarOp cOp1;
    cOp1.Type = 'M';
    cOp1.Length = bedpeLength1;
    bam1.CigarData.push_back(cOp1);
    CigarOp cOp2;
    cOp2.Type = 'M';
    cOp2.Length = bedpeLength2;
    bam2.CigarData.push_back(cOp2);
}


void bedpetobam_MakeBamHeader(const string &genomeFile, RefVector &refs, string &header,
                   map<string, int, std::less<string> > &chromToId) {

    // make a genome map of the genome file.
    GenomeFile genome(genomeFile);

    header += "@HD\tVN:1.0\tSO:unsorted\n";
    header += "@PG\tID:BEDTools_bedpeToBam\tVN:V";
    header += VERSION;
    header += "\n";

    int chromId = 0;
    vector<string> chromList = genome.getChromList();
    sort(chromList.begin(), chromList.end());

    // create a BAM header (@SQ) entry for each chrom in the BEDTools genome file.
    vector<string>::const_iterator genomeItr  = chromList.begin();
    vector<string>::const_iterator genomeEnd  = chromList.end();
    for (; genomeItr != genomeEnd; ++genomeItr) {
        chromToId[*genomeItr] = chromId;
        chromId++;

        // add to the header text
        int size = genome.getChromSize(*genomeItr);
        string chromLine = "@SQ\tSN:" + *genomeItr + "\tAS:" + genomeFile + "\tLN:" + ToString(size) + "\n";
        header += chromLine;

        // create a chrom entry and add it to
        // the RefVector
        RefData chrom;
        chrom.RefName            = *genomeItr;
        chrom.RefLength          = size;
        refs.push_back(chrom);
    }
}


/* Taken directly from the SAMTools spec
calculate bin given an alignment in [beg,end) (zero-based, half-close, half-open) */
int bedpetobam_reg2bin(int beg, int end) {
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}



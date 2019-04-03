/*****************************************************************************
  bedFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef BEDFILE_H
#define BEDFILE_H

// "local" includes
#include "BedtoolsTypes.h"
#include "gzstream.h"
#include "lineFileUtilities.h"
#include "fileType.h"

// standard includes
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <stdint.h>
#include <cstdio>
//#include <tr1/unordered_map>  // Experimental.
using namespace std;


//*************************************************
// Data type tydedef
//*************************************************
typedef uint16_t BINLEVEL;
typedef uint32_t BIN;
typedef uint16_t USHORT;
typedef uint32_t UINT;

//*************************************************
// Genome binning constants
//*************************************************
const BINLEVEL _binLevels = 8;

const BIN _binOffsetsExtended[] = {262144+32678+4096+512+64+8+1, 32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

const USHORT _binFirstShift = 14;       /* How much to shift to get to finest bin. */
const USHORT _binNextShift  = 3;        /* How much to shift to get to next larger bin. */

const BIN _numBins = (1 << (_binNextShift * _binLevels)) / ((1 << _binNextShift) - 1);


//*************************************************
// Common data structures
//*************************************************

struct DEPTH {
    UINT starts;
    UINT ends;
};


/*
    Structure for regular BED records
*/
struct BED {

    // Regular BED fields
    string chrom;
    CHRPOS start;
    CHRPOS end;
    string name;
    string score;
    string strand;
    // all of the original fields in the record
    vector<string> fields;
    // indices of the "other" fields
    vector<uint16_t> other_idxs;
    // is this a zero length feature: i.e., start == end
    bool   zeroLength;
    double weight;

public:
    // constructors

    // Null
    BED()
    : chrom(""),
      start(0),
      end(0),
      name(""),
      score(""),
      strand(""),
      fields(),
      other_idxs(),
      zeroLength(false),
      weight(0.0)
    {}

    // BED3
    BED(string chrom, CHRPOS start, CHRPOS end)
    : chrom(chrom),
      start(start),
      end(end),
      name(""),
      score(""),
      strand(""),
      fields(),
      other_idxs(),
      zeroLength(false),
      weight(0.0)
    {}

    // BED4
    BED(string chrom, CHRPOS start, CHRPOS end, string strand)
    : chrom(chrom),
      start(start),
      end(end),
      name(""),
      score(""),
      strand(strand),
      fields(),
      other_idxs(),
      zeroLength(false),
      weight(0.0)
    {}

    // BED6
    BED(string chrom, CHRPOS start, CHRPOS end, string name,
        string score, string strand)
    : chrom(chrom),
      start(start),
      end(end),
      name(name),
      score(score),
      strand(strand),
      fields(),
      other_idxs(),
      zeroLength(false),
      weight(0.0)
    {}

    // BEDALL
    BED(string chrom, CHRPOS start, CHRPOS end, string name,
        string score, string strand, vector<string> fields, 
        vector<uint16_t> other_idxs)
    : chrom(chrom),
      start(start),
      end(end),
      name(name),
      score(score),
      strand(strand),
      fields(fields),
      other_idxs(other_idxs),
      zeroLength(false),
      weight(0.0)
    {}
    
    CHRPOS size() const {
        return end-start;
    }

}; // BED


/*
    Structure for each end of a paired BED record
    mate points to the other end.
*/
struct MATE {
    BED bed;
    int lineNum;
    MATE *mate;
};


/*
    Structure for regular BED COVERAGE records
*/
struct BEDCOV {

    string chrom;
    
    // Regular BED fields
    CHRPOS start;
    CHRPOS end;
    string name;
    string score;
    string strand;

    // all of the original fields in the record
    vector<string> fields;
    // indices of the "other" fields
    vector<uint16_t> other_idxs;
    // is this a zero length feature: i.e., start == end
    bool   zeroLength;
    
    // Additional fields specific to computing coverage
    map<unsigned int, DEPTH> depthMap;
    unsigned int count;
    CHRPOS minOverlapStart;
    
    
    public:
    // constructors
    // Null
    BEDCOV()
    : chrom(""),
      start(0),
      end(0),
      name(""),
      score(""),
      strand(""),
      fields(),
      other_idxs(),
      zeroLength(false),
      depthMap(),
      count(0),
      minOverlapStart(0)
    {}
};


/*
    Structure for BED COVERAGE records having lists of
    multiple coverages
*/
struct BEDCOVLIST {

    // Regular BED fields
    string chrom;
    CHRPOS start;
    CHRPOS end;
    string name;
    string score;
    string strand;

    // all of the original fields in the record
    vector<string> fields;
    // indices of the "other" fields
    vector<uint16_t> other_idxs;
    // is this a zero length feature: i.e., start == end
    bool   zeroLength;

    // Additional fields specific to computing coverage
    vector< map<CHRPOS, DEPTH> > depthMapList;
    vector<unsigned int> counts;
    vector<CHRPOS> minOverlapStarts;
    
    
    public:
    // constructors
    // Null
    BEDCOVLIST()
    : chrom(""),
      start(0),
      end(0),
      name(""),
      score(""),
      strand(""),
      fields(),
      other_idxs(),
      zeroLength(false),
      depthMapList(),
      counts(0),
      minOverlapStarts(0)
    {}
};


// enum to flag the state of a given line in a BED file.
enum BedLineStatus
{
    BED_INVALID = -1,
    BED_HEADER  = 0,
    BED_BLANK   = 1,
    BED_VALID   = 2
};

// enum to indicate the type of file we are dealing with
enum FileType
{
    BED_FILETYPE,
    GFF_FILETYPE,
    VCF_FILETYPE
};

//*************************************************
// Data structure typedefs
//*************************************************
typedef vector<BED>    bedVector;
typedef vector<BEDCOV> bedCovVector;
typedef vector<MATE> mateVector;
typedef vector<BEDCOVLIST> bedCovListVector;

typedef map<BIN, bedVector> binsToBeds;
typedef map<BIN, bedCovVector> binsToBedCovs;
typedef map<BIN, mateVector> binsToMates;
typedef map<BIN, bedCovListVector> binsToBedCovLists;

typedef map<string, binsToBeds>    masterBedMap;
typedef map<string, binsToBedCovs> masterBedCovMap;
typedef map<string, binsToMates> masterMateMap;
typedef map<string, binsToBedCovLists> masterBedCovListMap;
typedef map<string, bedVector>     masterBedMapNoBin;


// return the genome "bin" for a feature with this start and end
inline
BIN getBin(CHRPOS start, CHRPOS end) {
    --end;
    start >>= _binFirstShift;
    end   >>= _binFirstShift;

    for (short i = 0; i < _binLevels; ++i) {
        if (start == end) return _binOffsetsExtended[i] + start;
        start >>= _binNextShift;
        end   >>= _binNextShift;
    }
    cerr << "start " << start << ", end " << end 
         << " out of range in findBin (max is 512M)" 
         << endl;
    return 0;
}

/****************************************************
// isInteger(s): Tests if string s is a valid integer
*****************************************************/
inline bool isInteger(const std::string& s) {
    int len = s.length();
    for (int i = 0; i < len; i++) {
        if (!std::isdigit(s[i])) return false;
    }
    return true;
}


// return the amount of overlap between two features.  Negative if none and the
// number of negative bases is the distance between the two.
inline
CHRPOS overlaps(CHRPOS aS, CHRPOS aE, CHRPOS bS, CHRPOS bE) {
    return min(aE, bE) - max(aS, bS);
}

// is A after (to the right of) B?
inline
bool after(const BED &a, const BED &b) {
    return (a.start >= b.end);
}


// Ancillary functions
void splitBedIntoBlocks(const BED &bed, bedVector &bedBlocks);


// BED Sorting Methods
bool sortByChrom(const BED &a, const BED &b);
bool sortByStart(const BED &a, const BED &b);
bool sortBySizeAsc(const BED &a, const BED &b);
bool sortBySizeDesc(const BED &a, const BED &b);
bool sortByScoreAsc(const BED &a, const BED &b);
bool sortByScoreDesc(const BED &a, const BED &b);
bool byChromThenStart(BED const &a, BED const &b);



//************************************************
// BedFile Class methods and elements
//************************************************
class BedFile {

public:

    // Constructor
    BedFile(string &);
    BedFile(void);

    // Destructor
    ~BedFile(void);

    /********* File management ********/
    // Open a BED file for reading (creates an istream pointer)
    void Open(void);

    // Close an opened BED file.
    void Close(void);

    // are the any intervals left in the file?
    bool Empty(void);

    // Rewind the pointer back to the beginning of the file
    void Rewind(void);

    // Jump to a specific byte in the file
    void Seek(unsigned long offset);

    // dump the header, which is collected as part of Open()
    void PrintHeader(void);

    // Get the next BED entry in an opened BED file.
    bool GetNextBed (BED &bed, bool forceSorted = false);
    
    // Returns the next MERGED (i.e., non-overlapping) interval in 
    // an opened BED file
    // NOTE: assumes input file is sorted by chrom then start
    bool GetNextMergedBed(BED &merged_bed);

    // load a BED file into a map keyed by chrom, then bin. value is 
    // vector of BEDs
    void loadBedFileIntoMap();
    void loadBedFileIntoMergedMap();

    // load a BED entry into and existing map
    void addBEDIntoMap(BED bedEntry);

    // load a BED file into a map keyed by chrom, then bin. value is 
    // vector of BEDCOVs
    void loadBedCovFileIntoMap();

    // load a BED file into a map keyed by chrom, then bin. value is 
    // vector of BEDCOVLISTs
    void loadBedCovListFileIntoMap();

    // load a BED file into a map keyed by chrom. value is vector of BEDs
    void loadBedFileIntoMapNoBin();

    // load a BED file into a vector of BEDs
    void loadBedFileIntoVector();

    // load a BED file into a vector ordered in decreasing order by size
    void assignWeightsBasedOnSize();

    BED * sizeWeightedSearch(double val);

    // Given a chrom, start, end and strand for a single feature,
    // search for all overlapping features in another BED file.
    // Searches through each relevant genome bin on the same chromosome
    // as the single feature. Note: Adapted from kent source "binKeeperFind"
    void allHits(string chrom, CHRPOS start, CHRPOS end, string strand, 
                 vector<BED> &hits, bool sameStrand, bool diffStrand, 
                 float overlapFraction, bool reciprocal);

    // return true if at least one overlap was found.  otherwise, return false.
    bool anyHits(string chrom, CHRPOS start, CHRPOS end, string strand,
                bool sameStrand, bool diffStrand, float overlapFraction, bool reciprocal);


    // Given a chrom, start, end and strand for a single feature,
    // increment a the number of hits for each feature in B file
    // that the feature overlaps
    void countHits(const BED &a, bool sameStrand = false, 
                   bool diffStrand = false, bool countsOnly = false);

    // same as above, but has special logic that processes a set of
    // BED "blocks" from a single entry so as to avoid over-counting
    // each "block" of a single BAM/BED12 as distinct coverage.  That is,
    // if one read has four block, we only want to count the coverage as
    // coming from one read, not four.
    void countSplitHits(const vector<BED> &bedBlock, bool sameStrand = false, 
                        bool diffStrand = false, bool countsOnly = false);

    // Given a chrom, start, end and strand for a single feature,
    // increment a the number of hits for each feature in B file
    // that the feature overlaps
    void countListHits(const BED &a, int index, 
                       bool sameStrand, bool diffStrand);
    
    
    // return the total length of all the intervals in the file.
    // use with GetNextBed()
    unsigned long getTotalLength(void);

    // return the total _flattened_ length of all the intervals in the file.
    // use with GetNextMergedBed()
    unsigned long getTotalFlattenedLength(void);
    

    // the bedfile with which this instance is associated
    string bedFile;
    unsigned int bedType;  // 3-6, 12 for BED
                           // 9 for GFF
    bool isBed12;          // is it file of true blocked BED12 records?
    bool isZeroBased;

    // Main data structures used by BEDTools
    masterBedCovMap      bedCovMap;
    masterBedCovListMap  bedCovListMap;
    masterBedMap         bedMap;
    bedVector            bedList;
    masterBedMapNoBin    bedMapNoBin;
    
    BedLineStatus _status;
    int _lineNum;

private:

    // data
    bool _isGff;
    bool _isVcf;
    bool _typeIsKnown;        // do we know the type?   (i.e., BED, GFF, VCF)
    FileType   _fileType;     // what is the file type? (BED? GFF? VCF?)
    istream   *_bedStream;
    string _bedLine;

    BED _nullBed;
    string _header;
    bool _firstLine;
    vector<string> _bedFields;
    unsigned int _numFields;
    CHRPOS _merged_start;
    CHRPOS _merged_end;
    string _merged_chrom;
    CHRPOS _prev_start;
    string _prev_chrom;
    unsigned long _total_length;
    unsigned long _total_flattened_length;

    void setZeroBased(bool zeroBased);
    void setGff (bool isGff);
    void setVcf (bool isVcf);
    void setFileType (FileType type);
    void setBedType (int colNums);
    void setBed12 (bool isBed12);

    /************ Private utilities ***********************/
    void GetHeader(void);

    /******************************************************
    Private definitions to circumvent linker issues with
    templated member functions.
    *******************************************************/

    /*
        parseLine: converts a lineVector into either BED or BEDCOV (templated, hence in header to avoid linker issues.)
    */
    template <typename T>
    inline BedLineStatus parseLine (T &bed, const vector<string> &fields) {
        
        // clear out the data from the last line.
        bed = _nullBed;
        // bail out if we have a blank line
        if (_numFields == 0) {
            return BED_BLANK;
        }
        // bail out if we have a comment line
        if ( (fields[0].find("#")       == 0) ||
             (fields[0].find("browser") == 0) ||
             (fields[0].find("track")   == 0) 
           )
        {
            return BED_HEADER;
        }

        if (_numFields >= 3) {
            // line parsing for all lines after the first non-header line
            if (_typeIsKnown == true) {
                switch(_fileType) {
                    case BED_FILETYPE:
                        if (parseBedLine(bed, fields) == true) 
                            return BED_VALID;
                    case VCF_FILETYPE:
                        if (parseVcfLine(bed, fields) == true)
                        {
                            return BED_VALID;
                        }
                    case GFF_FILETYPE:
                        if (parseGffLine(bed, fields) == true) 
                                return BED_VALID;
                    default:
                        printf("ERROR: file type encountered. Exiting\n");
                        exit(1);
                }
            }
            // line parsing for first non-header line: figure out file contents
            else {
                // it's BED format if columns 2 and 3 are integers
                if (isInteger(fields[1]) && isInteger(fields[2])) {
                    setGff(false);
                    setZeroBased(true);
                    setFileType(BED_FILETYPE);
                    // we now expect numFields columns in each line
                    setBedType(_numFields);
                    
                    // test to see if the file has true blocked BED12 records
                    if (_numFields == 12) {
                        CHRPOS cdsStart = stoll(fields[6].c_str());
                        CHRPOS cdsEnd   = stoll(fields[7].c_str());
                        int numExons = atoi(fields[9].c_str());

                        if (cdsStart > 0 && cdsEnd > 0&& numExons > 0 &&
                            fields[10].find(",") == 0 &&
                            fields[11].find(",") == 0)
                        {
                            setBed12(true);
                        }
                        else setBed12(false);
                    }
                    if (parseBedLine(bed, fields) == true) 
                        return BED_VALID;
                }
                // it's VCF, assuming the second column is numeric and 
                // there are at least 8 fields.
                else if (isInteger(fields[1]) && _numFields >= 8) {
                    setGff(false);
                    setVcf(true);
                    setZeroBased(false);
                    setFileType(VCF_FILETYPE);
                    // we now expect numFields columns in each line
                    setBedType(_numFields);
                    if (parseVcfLine(bed, fields) == true) 
                        return BED_VALID;
                }
                // it's GFF, assuming columns columns 4 and 5 are numeric 
                // and we have 9 fields total.
                else if ((_numFields >= 8) && 
                          isInteger(fields[3]) && 
                          isInteger(fields[4])) 
                {
                    setGff(true);
                    setZeroBased(false);
                    setFileType(GFF_FILETYPE);
                    // we now expect numFields columns in each line
                    setBedType(_numFields);
                    if (parseGffLine(bed, fields) == true) 
                    {
                        return BED_VALID;
                    }
                }
                else {
                    cerr << "Unexpected file format.  "
                         << "Please use tab-delimited BED, GFF, or VCF. " 
                         << "Perhaps you have non-integer starts or ends "
                         << "at line " 
                         << _lineNum 
                         << "?" 
                         << endl;
                    exit(1);
                }
            }
        }
        else {
            cerr << "It looks as though you have less than 3 columns at line " 
                 << _lineNum 
                 << " in file " << bedFile << ".  Are you sure your files are tab-delimited?" 
                 << endl;
            exit(1);
        }
        // default
        return BED_INVALID;
    }


    /*
        parseBedLine: converts a lineVector into either BED or BEDCOV (templated, hence in header to avoid linker issues.)
    */
    template <typename T>
    inline bool parseBedLine (T &bed, const vector<string> &fields) 
    {
        // process as long as the number of fields in this
        // line matches what we expect for this file.
        if (_numFields == this->bedType) {
            bed.fields = fields;
            bed.chrom = fields[0];

            CHRPOS i;
            i = stoll(fields[1].c_str());
            if (i<0) {
                 cerr << "Error: malformed BED entry at line " 
                      << _lineNum 
                      << ". Start Coordinate detected that is < 0. Exiting." 
                      << endl;
                 exit(1);
            }
            bed.start = i;
            i = stoll(fields[2].c_str());
            if (i<0) {
                cerr << "Error: malformed BED entry at line " 
                     << _lineNum 
                     << ". End Coordinate detected that is < 0. Exiting." 
                     << endl;
                exit(1);
            }
            bed.end = i;
            
            // handle starts == end (e.g., insertions in reference genome)
            if (bed.start == bed.end) {
                bed.start--;
                bed.end++;
                bed.zeroLength = true;
            }
            
            if (this->bedType == 4) {
                bed.name = fields[3];
            }
            else if (this->bedType == 5) {
                bed.name = fields[3];
                bed.score = fields[4];
            }
            else if (this->bedType == 6) {
                bed.name = fields[3];
                bed.score = fields[4];
                bed.strand = fields[5];
            }
            else if (this->bedType > 6) {
                bed.name = fields[3];
                bed.score = fields[4];
                bed.strand = fields[5];
                for (unsigned int i = 6; i < fields.size(); ++i) {
                    bed.other_idxs.push_back(i);
                }
            }
            else if (this->bedType != 3) {
                cerr << "Error: unexpected number of fields at line: " 
                     << _lineNum
                     << ".  Verify that your files are TAB-delimited. "     
                     << "Exiting..." 
                     << endl;
                exit(1);
            }

            // sanity checks.
            if (bed.start <= bed.end) {
                return true;
            }
            else {
                cerr << "Error: malformed BED entry at line " 
                     << _lineNum 
                     << ". Start was greater than end. Exiting." 
                     << endl;
                exit(1);
            }
        }
        else if (_numFields == 1) {
            cerr << "Only one BED field detected: " 
                 << _lineNum 
                 << ".  Verify that your files are TAB-delimited.  Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields != this->bedType) && (_numFields != 0)) {
            cerr << "Differing number of BED fields encountered at line: " 
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields < 3) && (_numFields != 0)) {
            cerr << "TAB delimited BED file with at least 3 fields"
                 << " (chrom, start, end) is required at line: "
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        return false;
    }


    /*
        parseVcfLine: converts a lineVector into either BED or BEDCOV (templated, hence in header to avoid linker issues.)
    */
    template <typename T>
    inline bool parseVcfLine (T &bed, const vector<string> &fields) 
    {
        if (_numFields >= this->bedType) {
            bed.fields = fields;
            bed.chrom  = fields[0];
            // VCF is one-based
            bed.start  = stoll(fields[1].c_str()) - 1;  
            // VCF 4.0 stores the size of the affected REF allele.
            bed.end    = bed.start + fields[3].size(); 
            bed.strand = "+";
            // construct the name from the ref and alt alleles.
            // if it's an annotated variant, add the rsId as well.
            bed.name   = fields[3] + "/" + fields[4];
            if (fields[2] != ".") {
                bed.name += "_" + fields[2];
            }

            if (this->bedType > 2) {
                for (unsigned int i = 2; i < _numFields; ++i)
                    bed.other_idxs.push_back(i);
            }

            if ((bed.start <= bed.end) && (bed.start >= 0) && ((bed.end) >= 0)) {
                return true;
            }
            else if (bed.start > bed.end) {
                cerr << "Error: malformed VCF entry at line " 
                    << _lineNum 
                    << ". Start was greater than end. Exiting." 
                    << endl;
                exit(1);
            }
            else if ( (bed.start < 0) || ((bed.end) < 0) ) {
                cerr << "Error: malformed VCF entry at line " 
                     << _lineNum << ". Coordinate detected that is < 0. "
                     << "Exiting." 
                     << endl;
                exit(1);
            }
        }
        else if (_numFields == 1) {
            cerr << "Only one VCF field detected: " 
                 << _lineNum 
                 << ".  Verify that your files are TAB-delimited. "
                 << "Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields != this->bedType) && (_numFields != 0)) {
            cerr << "Differing number of VCF fields encountered at line: " 
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields < 2) && (_numFields != 0)) {
            cerr << "TAB delimited VCF file with at least 2 fields "
                 << "(chrom, pos) is required at line: "
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        return false;
    }



    /*
        parseGffLine: converts a lineVector into either BED or BEDCOV (templated, hence in header to avoid linker issues.)
    */
    template <typename T>
    inline bool parseGffLine (T &bed, const vector<string> &fields) 
    {
        if (_numFields == this->bedType) {
            bed.fields = fields;
            if (this->bedType >= 8 && _isGff) {
                bed.chrom = fields[0];
                if (isInteger(fields[3]))
                    bed.start  = stoll(fields[3].c_str());
                if (isInteger(fields[4]))
                    bed.end  = stoll(fields[4].c_str());
                bed.name   = fields[2];
                bed.score  = fields[5];
                bed.strand = fields[6].c_str();
                // add GFF "source". unused in BED
                bed.other_idxs.push_back(1); 
                // add GFF "fname". unused in BED
                bed.other_idxs.push_back(7);  
                // handle the optional 9th field.
                if (this->bedType == 9)
                    // add GFF "group". unused in BED
                    bed.other_idxs.push_back(8); 
                bed.start--;
            }
            else {
                cerr << "Error: unexpected number of fields at line: " 
                    << _lineNum 
                    << ".  Verify that your files are TAB-delimited and that "
                    << "your GFF file has 8 or 9 fields.  Exiting..." 
                    << endl;
                exit(1);
            }
            if (bed.start > bed.end) {
                cerr << "Error: malformed GFF entry at line " 
                     << _lineNum 
                     << ". Start was greater than end. Exiting." 
                     << endl;
                exit(1);
            }
            if ( (bed.start < 0) || ((bed.end) < 0) ) {
                cerr << "Error: malformed GFF entry at line " 
                     << _lineNum 
                     << ". Coordinate detected that is < 1. Exiting." 
                     << endl;
                exit(1);
            }
            return true;
        }
        else if (_numFields == 1) {
            cerr << "Only one GFF field detected: " 
                 << _lineNum 
                 << ".  Verify that your files are TAB-delimited.  Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields != this->bedType) && (_numFields != 0)) {
            cerr << "Differing number of GFF fields encountered at line: " 
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        else if ((_numFields < 8) && (_numFields != 0)) {
            cerr << "TAB delimited GFF file with 8 or 9 fields is required"
                 << " at line: "
                 << _lineNum 
                 << ".  Exiting..." 
                 << endl;
            exit(1);
        }
        return false;
    }


public:

    /*
        reportBedTab

        Writes the _original_ BED entry with a TAB
        at the end of the line.
        Works for BED3 - BED6.
    */
    template <typename T>
    inline void reportBedTab(const T &bed) {
        // if it is azeroLength feature, we need to
        // correct the start and end coords to what they were
        // in the original file
        CHRPOS start = bed.start;
        CHRPOS end   = bed.end;
        if (bed.zeroLength) {
            if (_isGff == false)
                start++;
            end--;
        }
        
        // BED
        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t", bed.chrom.c_str(), start, end);
            }
            else if (this->bedType == 4) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t",
                    bed.chrom.c_str(), start, end, bed.name.c_str());
            }
            else if (this->bedType == 5) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, 
                    bed.name.c_str(), bed.score.c_str());
            }
            else if (this->bedType == 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, 
                    bed.name.c_str(), bed.score.c_str(), bed.strand.c_str());
            }
            else if (this->bedType > 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());
                
                vector<uint16_t>::const_iterator 
                    othIt  = bed.other_idxs.begin();
                vector<uint16_t>::const_iterator 
                    othEnd = bed.other_idxs.end();
                for ( ; othIt != othEnd; ++othIt) {
                    printf("%s\t", bed.fields[*othIt].c_str());
                }
            }
        }
        // VCF
        else if (_isGff == false && _isVcf == true) {
            printf ("%s\t%" PRId_CHRPOS "\t", bed.chrom.c_str(), start+1);

            vector<uint16_t>::const_iterator othIt  = bed.other_idxs.begin();
            vector<uint16_t>::const_iterator othEnd = bed.other_idxs.end();
            for ( ; othIt != othEnd; ++othIt) {
                printf("%s\t", bed.fields[*othIt].c_str());
            }
        }
        // GFF
        else if (_isGff == true) {
            // "GFF-8"
            if (this->bedType == 8) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end, bed.score.c_str(), 
                    bed.strand.c_str(), bed.fields[bed.other_idxs[1]].c_str());
            }
            // "GFF-9"
            else if (this->bedType == 9) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str(), 
                    bed.fields[bed.other_idxs[2]].c_str());
            }
        }
    }



    /*
        reportToFileBedNewLine

        Writes the _original_ BED entry with a NEWLINE
        at the end of the line.
        Works for BED3 - BED6.
    */
    template <typename T>
    inline void reportToFileBedNewLine(FILE* out,const T &bed) {
        
        // if it is azeroLength feature, we need to
        // correct the start and end coords to what they were
        // in the original file
        CHRPOS start = bed.start;
        CHRPOS end   = bed.end;
        if (bed.zeroLength) {
            if (_isGff == false)
                start++;
            end--;
        }
        //BED
        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                fprintf(out,"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\n", bed.chrom.c_str(), start, end);
            }
            else if (this->bedType == 4) {
                fprintf(out,"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\n",
                    bed.chrom.c_str(), start, end, bed.name.c_str());
            }
            else if (this->bedType == 5) {
                fprintf(out,"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\n",
                    bed.chrom.c_str(), start, end, 
                    bed.name.c_str(), bed.score.c_str());
            }
            else if (this->bedType == 6) {
                fprintf(out,"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());
            }
            else if (this->bedType > 6) {
                fprintf(out,"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());

                vector<uint16_t>::const_iterator 
                    othIt  = bed.other_idxs.begin();
                vector<uint16_t>::const_iterator 
                    othEnd = bed.other_idxs.end();
                for ( ; othIt != othEnd; ++othIt) {
                    fprintf(out,"\t%s", bed.fields[*othIt].c_str());
                }
                fprintf(out,"\n");
            }
        }
        // VCF
        else if (_isGff == false && _isVcf == true) {
            fprintf(out,"%s\t%" PRId_CHRPOS, bed.chrom.c_str(), start+1);

            vector<uint16_t>::const_iterator othIt  = bed.other_idxs.begin();
            vector<uint16_t>::const_iterator othEnd = bed.other_idxs.end();
            for ( ; othIt != othEnd; ++othIt) {
                fprintf(out,"\t%s", bed.fields[*othIt].c_str());
            }
            fprintf(out,"\n");
        }
        // GFF
        else if (_isGff == true) {
            // "GFF-8"
            if (this->bedType == 8) {
                fprintf (out,"%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str());
            }
            // "GFF-9"
            else if (this->bedType == 9) {
                fprintf (out,"%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str(), 
                    bed.fields[bed.other_idxs[2]].c_str());
            }
        }
    }

    /* specialized version of reportBedNewLine printing to stdout */
     template <typename T>
     inline void reportBedNewLine(const T &bed) {
     reportToFileBedNewLine(stdout,bed);
     }

    /*
        reportBedRangeNewLine

        Writes a custom start->end for a BED entry
        with a NEWLINE at the end of the line.

        Works for BED3 - BED6.
    */
    template <typename T>
    inline void reportBedRangeTab(const T &bed, CHRPOS start, CHRPOS end) {
        // if it is azeroLength feature, we need to
        // correct the start and end coords to what they were
        // in the original file
        if (bed.zeroLength) {
            start = bed.start + 1;
            end   = bed.end - 1;
        }
        // BED
        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t", bed.chrom.c_str(), start, end);
            }
            else if (this->bedType == 4) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t",
                    bed.chrom.c_str(), start, end, bed.name.c_str());
            }
            else if (this->bedType == 5) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, 
                    bed.name.c_str(), bed.score.c_str());
            }
            else if (this->bedType == 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());
            }
            else if (this->bedType > 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());
                
                vector<uint16_t>::const_iterator 
                    othIt  = bed.other_idxs.begin();
                vector<uint16_t>::const_iterator 
                    othEnd = bed.other_idxs.end();
                for ( ; othIt != othEnd; ++othIt) {
                    printf("%s\t", bed.fields[*othIt].c_str());
                }
            }
        }
        // VCF
        else if (_isGff == false && _isVcf == true) {
            printf ("%s\t%" PRId_CHRPOS "\t", bed.chrom.c_str(), bed.start+1);
            vector<uint16_t>::const_iterator othIt  = bed.other_idxs.begin();
            vector<uint16_t>::const_iterator othEnd = bed.other_idxs.end();
            for ( ; othIt != othEnd; ++othIt) {
                printf("%s\t", bed.fields[*othIt].c_str());
            }
        }
        // GFF
        else if (_isGff == true) {
            // "GFF-8"
            if (this->bedType == 8) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str());
            }
            // "GFF-9"
            else if (this->bedType == 9) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\t",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str(), 
                    bed.fields[bed.other_idxs[2]].c_str());
            }
        }
    }



    /*
        reportBedRangeTab

        Writes a custom start->end for a BED entry
        with a TAB at the end of the line.

        Works for BED3 - BED6.
    */
    template <typename T>
    inline void reportBedRangeNewLine(const T &bed, CHRPOS start, CHRPOS end) {
        
        // if it is azeroLength feature, we need to
        // correct the start and end coords to what they were
        // in the original file
        if (bed.zeroLength) {
            start = bed.start + 1;
            end   = bed.end - 1;
        }
        // BED
        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\n", bed.chrom.c_str(), start, end);
            }
            else if (this->bedType == 4) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\n",
                    bed.chrom.c_str(), start, end, bed.name.c_str());
            }
            else if (this->bedType == 5) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\n",
                    bed.chrom.c_str(), start, end, 
                    bed.name.c_str(), bed.score.c_str());
            }
            else if (this->bedType == 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());
            }
            else if (this->bedType > 6) {
                printf ("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s",
                    bed.chrom.c_str(), start, end, bed.name.c_str(),
                    bed.score.c_str(), bed.strand.c_str());

                vector<uint16_t>::const_iterator 
                    othIt  = bed.other_idxs.begin();
                vector<uint16_t>::const_iterator 
                    othEnd = bed.other_idxs.end();
                for ( ; othIt != othEnd; ++othIt) {
                    printf("\t%s", bed.fields[*othIt].c_str());
                }
                printf("\n");
            }
        }
        // VCF
        else if (_isGff == false && _isVcf == true) {
            printf ("%s\t%" PRId_CHRPOS, bed.chrom.c_str(), bed.start+1);
            vector<uint16_t>::const_iterator othIt  = bed.other_idxs.begin();
            vector<uint16_t>::const_iterator othEnd = bed.other_idxs.end();
            for ( ; othIt != othEnd; ++othIt) {
                printf("\t%s", bed.fields[*othIt].c_str());
            }
            printf("\n");
        }
        // GFF
        else if (_isGff == true) {
            // "GFF-8"
            if (this->bedType == 8) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),                                          bed.fields[bed.other_idxs[1]].c_str());
            }
            // "GFF-9"
            else if (this->bedType == 9) {
                printf ("%s\t%s\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\n",
                    bed.chrom.c_str(), bed.fields[bed.other_idxs[0]].c_str(),
                    bed.name.c_str(), start+1, end,
                    bed.score.c_str(), bed.strand.c_str(),
                    bed.fields[bed.other_idxs[1]].c_str(), 
                    bed.fields[bed.other_idxs[2]].c_str());
            }
        }
    }


    /*
        reportNullBedTab
    */
    void reportNullBedTab() {

        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                printf (".\t-1\t-1\t");
            }
            else if (this->bedType == 4) {
                printf (".\t-1\t-1\t.\t");
            }
            else if (this->bedType == 5) {
                printf (".\t-1\t-1\t.\t-1\t");
            }
            else if (this->bedType == 6) {
                printf (".\t-1\t-1\t.\t-1\t.\t");
            }
            else if (this->bedType > 6) {
                printf (".\t-1\t-1\t.\t-1\t.\t");
                for (unsigned int i = 6; i < this->bedType; ++i) {
                    printf(".\t");
                }
            }
        }
        else if (_isGff == true && _isVcf == false) {
            if (this->bedType == 8) {
                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t");
            }
            else if (this->bedType == 9) {
                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t.\t");
            }
        }
    }


    /*
        reportNullBedTab
    */
    void reportNullBedNewLine() {

        if (_isGff == false && _isVcf == false) {
            if (this->bedType == 3) {
                printf (".\t-1\t-1\n");
            }
            else if (this->bedType == 4) {
                printf (".\t-1\t-1\t.\n");
            }
            else if (this->bedType == 5) {
                printf (".\t-1\t-1\t.\t-1\n");
            }
            else if (this->bedType == 6) {
                printf (".\t-1\t-1\t.\t-1\t.\n");
            }
            else if (this->bedType > 6) {
                printf (".\t-1\t-1\t.\t-1\t.");
                for (unsigned int i = 6; i < this->bedType; ++i) {
                    printf("\t.");
                }
                printf("\n");
            }
        }
        else if (_isGff == true && _isVcf == false) {
            if (this->bedType == 8) {
                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\n");
            }
            else if (this->bedType == 9) {
                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t.\n");
            }
        }
    }


};

#endif /* BEDFILE_H */

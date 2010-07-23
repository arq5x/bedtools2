// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 15 July 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

// C++ includes
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>

// BamTools includes
#include "BGZF.h"
#include "BamReader.h"
#include "BamIndex.h"
using namespace BamTools;
using namespace std;

struct BamReader::BamReaderPrivate {

    // -------------------------------
    // structs, enums, typedefs
    // -------------------------------
    enum RegionState { BEFORE_REGION = 0
                     , WITHIN_REGION
                     , AFTER_REGION
                     };

    // -------------------------------
    // data members
    // -------------------------------

    // general file data
    BgzfData  mBGZF;
    string    HeaderText;
    //BamIndex  Index;
    BamIndex* NewIndex;
    RefVector References;
    bool      IsIndexLoaded;
    int64_t   AlignmentsBeginOffset;
    string    Filename;
    string    IndexFilename;
    
    // system data
    bool IsBigEndian;

    // user-specified region values
    BamRegion Region;
    bool IsLeftBoundSpecified;
    bool IsRightBoundSpecified;
    
    bool IsRegionSpecified;
    int  CurrentRefID;
    int  CurrentLeft;

    // parent BamReader
    BamReader* Parent;
    
    // BAM character constants
    const char* DNA_LOOKUP;
    const char* CIGAR_LOOKUP;

    // -------------------------------
    // constructor & destructor
    // -------------------------------
    BamReaderPrivate(BamReader* parent);
    ~BamReaderPrivate(void);

    // -------------------------------
    // "public" interface
    // -------------------------------

    // file operations
    void Close(void);
    bool Jump(int refID, int position = 0);
    bool Open(const string& filename, const string& indexFilename = "");
    bool Rewind(void);
    bool SetRegion(const BamRegion& region);

    // access alignment data
    bool GetNextAlignment(BamAlignment& bAlignment);
    bool GetNextAlignmentCore(BamAlignment& bAlignment);

    // access auxiliary data
    int GetReferenceID(const string& refName) const;

    // index operations
    bool CreateIndex(bool useDefaultIndex);

    // -------------------------------
    // internal methods
    // -------------------------------

    // *** reading alignments and auxiliary data *** //

    // fills out character data for BamAlignment data
    bool BuildCharData(BamAlignment& bAlignment);
    // checks to see if alignment overlaps current region
    RegionState IsOverlap(BamAlignment& bAlignment);
    // retrieves header text from BAM file
    void LoadHeaderData(void);
    // retrieves BAM alignment under file pointer
    bool LoadNextAlignment(BamAlignment& bAlignment);
    // builds reference data structure from BAM file
    void LoadReferenceData(void);

    // *** index file handling *** //

    // clear out inernal index data structure
    void ClearIndex(void);
    // loads index from BAM index file
    bool LoadIndex(void);
};

// -----------------------------------------------------
// BamReader implementation (wrapper around BRPrivate)
// -----------------------------------------------------
// constructor
BamReader::BamReader(void) {
    d = new BamReaderPrivate(this);
}

// destructor
BamReader::~BamReader(void) {
    delete d;
    d = 0;
}

// file operations
void BamReader::Close(void) { d->Close(); }
bool BamReader::IsOpen(void) const { return d->mBGZF.IsOpen; }
bool BamReader::Jump(int refID, int position) { 
    d->Region.LeftRefID = refID;
    d->Region.LeftPosition = position;
    d->IsLeftBoundSpecified = true;
    d->IsRightBoundSpecified = false;
    return d->Jump(refID, position); 
}
bool BamReader::Open(const string& filename, const string& indexFilename) { return d->Open(filename, indexFilename); }
bool BamReader::Rewind(void) { return d->Rewind(); }
bool BamReader::SetRegion(const BamRegion& region) { return d->SetRegion(region); }
bool BamReader::SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound) {
    return d->SetRegion( BamRegion(leftRefID, leftBound, rightRefID, rightBound) );
}

// access alignment data
bool BamReader::GetNextAlignment(BamAlignment& bAlignment) { return d->GetNextAlignment(bAlignment); }
bool BamReader::GetNextAlignmentCore(BamAlignment& bAlignment) { return d->GetNextAlignmentCore(bAlignment); }

// access auxiliary data
const string BamReader::GetHeaderText(void) const { return d->HeaderText; }
int BamReader::GetReferenceCount(void) const { return d->References.size(); }
const RefVector& BamReader::GetReferenceData(void) const { return d->References; }
int BamReader::GetReferenceID(const string& refName) const { return d->GetReferenceID(refName); }
const std::string BamReader::GetFilename(void) const { return d->Filename; }

// index operations
bool BamReader::CreateIndex(bool useDefaultIndex) { return d->CreateIndex(useDefaultIndex); }

// -----------------------------------------------------
// BamReaderPrivate implementation
// -----------------------------------------------------

// constructor
BamReader::BamReaderPrivate::BamReaderPrivate(BamReader* parent)
    : NewIndex(0)
    , IsIndexLoaded(false)
    , AlignmentsBeginOffset(0)
    , IsLeftBoundSpecified(false)
    , IsRightBoundSpecified(false)
    , IsRegionSpecified(false)
    , CurrentRefID(0)
    , CurrentLeft(0)
    , Parent(parent)
    , DNA_LOOKUP("=ACMGRSVTWYHKDBN")
    , CIGAR_LOOKUP("MIDNSHP")
{ 
    IsBigEndian = SystemIsBigEndian();
}

// destructor
BamReader::BamReaderPrivate::~BamReaderPrivate(void) {
    Close();
}

bool BamReader::BamReaderPrivate::BuildCharData(BamAlignment& bAlignment) {
  
    // calculate character lengths/offsets
    const unsigned int dataLength      = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    const unsigned int cigarDataOffset = bAlignment.SupportData.QueryNameLength;
    const unsigned int seqDataOffset   = bAlignment.SupportData.QueryNameLength + (bAlignment.SupportData.NumCigarOperations * 4);
    const unsigned int qualDataOffset  = seqDataOffset + (bAlignment.SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset   = qualDataOffset + bAlignment.SupportData.QuerySequenceLength;
    const unsigned int tagDataLength   = dataLength - tagDataOffset;
      
    // set up char buffers
    const char*     allCharData = bAlignment.SupportData.AllCharData.data();
          uint32_t* cigarData   = (uint32_t*)(allCharData + cigarDataOffset);
    const char*     seqData     = ((const char*)allCharData) + seqDataOffset;
    const char*     qualData    = ((const char*)allCharData) + qualDataOffset;
          char*     tagData     = ((char*)allCharData) + tagDataOffset;
  
    // store alignment name (depends on null char as terminator)
    bAlignment.Name.assign((const char*)(allCharData));    
          
    // save CigarOps 
    CigarOp op;
    bAlignment.CigarData.clear();
    bAlignment.CigarData.reserve(bAlignment.SupportData.NumCigarOperations);
    for (unsigned int i = 0; i < bAlignment.SupportData.NumCigarOperations; ++i) {

        // swap if necessary
        if ( IsBigEndian ) { SwapEndian_32(cigarData[i]); }
      
        // build CigarOp structure
        op.Length = (cigarData[i] >> BAM_CIGAR_SHIFT);
        op.Type   = CIGAR_LOOKUP[ (cigarData[i] & BAM_CIGAR_MASK) ];

        // save CigarOp
        bAlignment.CigarData.push_back(op);
    }
          
          
    // save query sequence
    bAlignment.QueryBases.clear();
    bAlignment.QueryBases.reserve(bAlignment.SupportData.QuerySequenceLength);
    for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
        char singleBase = DNA_LOOKUP[ ( ( seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
        bAlignment.QueryBases.append(1, singleBase);
    }
  
    // save qualities, converting from numeric QV to 'FASTQ-style' ASCII character
    bAlignment.Qualities.clear();
    bAlignment.Qualities.reserve(bAlignment.SupportData.QuerySequenceLength);
    for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
        char singleQuality = (char)(qualData[i]+33);
        bAlignment.Qualities.append(1, singleQuality);
    }
    
    // if QueryBases is empty (and this is a allowed case)
    if ( bAlignment.QueryBases.empty() ) 
        bAlignment.AlignedBases = bAlignment.QueryBases;
    
    // if QueryBases contains data, then build AlignedBases using CIGAR data
    else {
    
        // resize AlignedBases
        bAlignment.AlignedBases.clear();
        bAlignment.AlignedBases.reserve(bAlignment.SupportData.QuerySequenceLength);
      
        // iterate over CigarOps
        int k = 0;
        vector<CigarOp>::const_iterator cigarIter = bAlignment.CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = bAlignment.CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            
            const CigarOp& op = (*cigarIter);
            switch(op.Type) {
              
                case ('M') :
                case ('I') :
                    bAlignment.AlignedBases.append(bAlignment.QueryBases.substr(k, op.Length)); // for 'M', 'I' - write bases
                    // fall through
                
                case ('S') :
                    k += op.Length;                                     // for 'S' - soft clip, skip over query bases
                    break;
                    
                case ('D') :
                    bAlignment.AlignedBases.append(op.Length, '-');     // for 'D' - write gap character
                    break;
                    
                case ('P') :
                    bAlignment.AlignedBases.append( op.Length, '*' );   // for 'P' - write padding character
                    break;
                    
                case ('N') :
                    bAlignment.AlignedBases.append( op.Length, 'N' );  // for 'N' - write N's, skip bases in original query sequence
                    break;
                    
                case ('H') :
                    break;  // for 'H' - hard clip, do nothing to AlignedBases, move to next op
                    
                default:
                    printf("ERROR: Invalid Cigar op type\n"); // shouldn't get here
                    exit(1);
            }
        }
    }
 
    // -----------------------
    // Added: 3-25-2010 DB
    // Fixed: endian-correctness for tag data
    // -----------------------
    if ( IsBigEndian ) {
        int i = 0;
        while ( (unsigned int)i < tagDataLength ) {
          
            i += 2; // skip tag type (e.g. "RG", "NM", etc)
            uint8_t type = toupper(tagData[i]);     // lower & upper case letters have same meaning 
            ++i;                                    // skip value type
    
            switch (type) {
                
                case('A') :
                case('C') : 
                    ++i;
                    break;

                case('S') : 
                    SwapEndian_16p(&tagData[i]); 
                    i += sizeof(uint16_t);
                    break;
                    
                case('F') :
                case('I') : 
                    SwapEndian_32p(&tagData[i]);
                    i += sizeof(uint32_t);
                    break;
                
                case('D') : 
                    SwapEndian_64p(&tagData[i]);
                    i += sizeof(uint64_t);
                    break;
                
                case('H') :
                case('Z') : 
                    while (tagData[i]) { ++i; }
                    ++i; // increment one more for null terminator
                    break;
                
                default : 
                    printf("ERROR: Invalid tag value type\n"); // shouldn't get here
                    exit(1);
            }
        }
    }
    
    // store TagData
    bAlignment.TagData.clear();
    bAlignment.TagData.resize(tagDataLength);
    memcpy((char*)bAlignment.TagData.data(), tagData, tagDataLength);
    
    // clear the core-only flag
    bAlignment.SupportData.HasCoreOnly = false;
    
    // return success
    return true;
}

// clear index data structure
void BamReader::BamReaderPrivate::ClearIndex(void) {
    delete NewIndex;
    NewIndex = 0;
}

// closes the BAM file
void BamReader::BamReaderPrivate::Close(void) {
    
    // close BGZF file stream
    mBGZF.Close();
    
    // clear out index data
    ClearIndex();
    
    // clear out header data
    HeaderText.clear();
    
    // clear out region flags
    IsLeftBoundSpecified  = false;
    IsRightBoundSpecified = false;
    IsRegionSpecified     = false;
}

// create BAM index from BAM file (keep structure in memory) and write to default index output file
bool BamReader::BamReaderPrivate::CreateIndex(bool useDefaultIndex) {

    // clear out prior index data
    ClearIndex();
    
    // create default index
    if ( useDefaultIndex )
        NewIndex = new BamDefaultIndex(&mBGZF, Parent, IsBigEndian);
    // create BamTools 'custom' index
    else
        NewIndex = new BamToolsIndex(&mBGZF, Parent, IsBigEndian);
    
    bool ok = true;
    ok &= NewIndex->Build();
    ok &= NewIndex->Write(Filename); 
    
    // return success/fail
    return ok;
}

// get next alignment (from specified region, if given)
bool BamReader::BamReaderPrivate::GetNextAlignment(BamAlignment& bAlignment) {

    // if valid alignment found, attempt to parse char data, and return success/failure
    if ( GetNextAlignmentCore(bAlignment) )
        return BuildCharData(bAlignment);
    
    // no valid alignment found
    else
        return false;
}

// retrieves next available alignment core data (returns success/fail)
// ** DOES NOT parse any character data (read name, bases, qualities, tag data)
//    these can be accessed, if necessary, from the supportData 
// useful for operations requiring ONLY positional or other alignment-related information
bool BamReader::BamReaderPrivate::GetNextAlignmentCore(BamAlignment& bAlignment) {

    // if valid alignment available
    if ( LoadNextAlignment(bAlignment) ) {

        // set core-only flag
        bAlignment.SupportData.HasCoreOnly = true;
      
        // if region not specified, return success
        if ( !IsLeftBoundSpecified ) return true;

        // determine region state (before, within, after)
        BamReader::BamReaderPrivate::RegionState state = IsOverlap(bAlignment);
      
        // if alignment lies after region, return false
        if ( state == AFTER_REGION ) 
            return false;

        while ( state != WITHIN_REGION ) {
            // if no valid alignment available (likely EOF) return failure
            if ( !LoadNextAlignment(bAlignment) ) return false;
            // if alignment lies after region, return false (no available read within region)
            state = IsOverlap(bAlignment);
            if ( state == AFTER_REGION) return false;
            
        }

        // return success (alignment found that overlaps region)
        return true;
    }

    // no valid alignment
    else
        return false;
}

// returns RefID for given RefName (returns References.size() if not found)
int BamReader::BamReaderPrivate::GetReferenceID(const string& refName) const {

    // retrieve names from reference data
    vector<string> refNames;
    RefVector::const_iterator refIter = References.begin();
    RefVector::const_iterator refEnd  = References.end();
    for ( ; refIter != refEnd; ++refIter) {
        refNames.push_back( (*refIter).RefName );
    }

    // return 'index-of' refName ( if not found, returns refNames.size() )
    return distance(refNames.begin(), find(refNames.begin(), refNames.end(), refName));
}

// returns region state - whether alignment ends before, overlaps, or starts after currently specified region
// this *internal* method should ONLY called when (at least) IsLeftBoundSpecified == true
BamReader::BamReaderPrivate::RegionState BamReader::BamReaderPrivate::IsOverlap(BamAlignment& bAlignment) {
    
    // --------------------------------------------------
    // check alignment start against right bound cutoff
    
    // if full region of interest was given
    if ( IsRightBoundSpecified ) {
      
        // read starts on right bound reference, but AFTER right bound position
        if ( bAlignment.RefID == Region.RightRefID && bAlignment.Position > Region.RightPosition )
            return AFTER_REGION;
      
        // if read starts on reference AFTER right bound, return false
        if ( bAlignment.RefID > Region.RightRefID ) 
            return AFTER_REGION;
    }
  
    // --------------------------------------------------------
    // no right bound given OR read starts before right bound
    // so, check if it overlaps left bound 
  
    // if read starts on left bound reference AND after left boundary, return success
    if ( bAlignment.RefID == Region.LeftRefID && bAlignment.Position >= Region.LeftPosition)
        return WITHIN_REGION;
  
    // if read is on any reference sequence before left bound, return false
    if ( bAlignment.RefID < Region.LeftRefID )
        return BEFORE_REGION;

    // --------------------------------------------------------
    // read is on left bound reference, but starts before left bound position

    // if it overlaps, return WITHIN_REGION
    if ( bAlignment.GetEndPosition() >= Region.LeftPosition )
        return WITHIN_REGION;
    // else begins before left bound position
    else
        return BEFORE_REGION;
}

// jumps to specified region(refID, leftBound) in BAM file, returns success/fail
bool BamReader::BamReaderPrivate::Jump(int refID, int position) {

    // -----------------------------------------------------------------------
    // check for existing index 
    if ( NewIndex == 0 ) return false; 
    // see if reference has alignments
    if ( !NewIndex->HasAlignments(refID) ) return false; 
    // make sure position is valid
    if ( position > References.at(refID).RefLength ) return false;
    
    // determine possible offsets
    vector<int64_t> offsets;
    if ( !NewIndex->GetOffsets(Region, IsRightBoundSpecified, offsets) ) {
        printf("ERROR: Could not jump: unable to calculate offset for specified region.\n");
        return false;
    }
      
    // iterate through offsets
    BamAlignment bAlignment;
    bool result = true;
    for ( vector<int64_t>::const_iterator o = offsets.begin(); o != offsets.end(); ++o) {
        
        // attempt seek & load first available alignment
        result &= mBGZF.Seek(*o);
        LoadNextAlignment(bAlignment);
        
        // if this alignment corresponds to desired position
        // return success of seeking back to 'current offset'
        if ( (bAlignment.RefID == refID && bAlignment.Position + bAlignment.Length > position) || (bAlignment.RefID > refID) )
            return mBGZF.Seek(*o);
    }
    
    return result;
}

// load BAM header data
void BamReader::BamReaderPrivate::LoadHeaderData(void) {

    // check to see if proper BAM header
    char buffer[4];
    if (mBGZF.Read(buffer, 4) != 4) {
        printf("Could not read header type\n");
        exit(1);
    }

    if (strncmp(buffer, "BAM\001", 4)) {
        printf("wrong header type!\n");
        exit(1);
    }

    // get BAM header text length
    mBGZF.Read(buffer, 4);
    unsigned int headerTextLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(headerTextLength); }
    
    // get BAM header text
    char* headerText = (char*)calloc(headerTextLength + 1, 1);
    mBGZF.Read(headerText, headerTextLength);
    HeaderText = (string)((const char*)headerText);

    // clean up calloc-ed temp variable
    free(headerText);
}

// load existing index data from BAM index file (".bai"), return success/fail
bool BamReader::BamReaderPrivate::LoadIndex(void) {

    // clear out any existing index data
    ClearIndex();

    // skip if index file empty
    if ( IndexFilename.empty() )
        return false;

    // check supplied filename for index type
    size_t defaultExtensionFound = IndexFilename.find(".bai");
    size_t customExtensionFound  = IndexFilename.find(".bti");
    
    // if SAM/BAM default (".bai")
    if ( defaultExtensionFound != string::npos )
        NewIndex = new BamDefaultIndex(&mBGZF, Parent, IsBigEndian);
    
    // if BamTools custom index (".bti")
    else if ( customExtensionFound != string::npos )
        NewIndex = new BamToolsIndex(&mBGZF, Parent, IsBigEndian);
    
    // else unknown
    else {
        printf("ERROR: Unknown index file extension.\n");
        return false;
    }
    
    // return success of loading index data
    return NewIndex->Load(IndexFilename);
}

// populates BamAlignment with alignment data under file pointer, returns success/fail
bool BamReader::BamReaderPrivate::LoadNextAlignment(BamAlignment& bAlignment) {

    // read in the 'block length' value, make sure it's not zero
    char buffer[4];
    mBGZF.Read(buffer, 4);
    bAlignment.SupportData.BlockLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(bAlignment.SupportData.BlockLength); }
    if ( bAlignment.SupportData.BlockLength == 0 ) { return false; }

    // read in core alignment data, make sure the right size of data was read
    char x[BAM_CORE_SIZE];
    if ( mBGZF.Read(x, BAM_CORE_SIZE) != BAM_CORE_SIZE ) { return false; }

    if ( IsBigEndian ) {
        for ( int i = 0; i < BAM_CORE_SIZE; i+=sizeof(uint32_t) ) { 
          SwapEndian_32p(&x[i]); 
        }
    }
    
    // set BamAlignment 'core' and 'support' data
    bAlignment.RefID    = BgzfData::UnpackSignedInt(&x[0]);  
    bAlignment.Position = BgzfData::UnpackSignedInt(&x[4]);
    
    unsigned int tempValue = BgzfData::UnpackUnsignedInt(&x[8]);
    bAlignment.Bin        = tempValue >> 16;
    bAlignment.MapQuality = tempValue >> 8 & 0xff;
    bAlignment.SupportData.QueryNameLength = tempValue & 0xff;

    tempValue = BgzfData::UnpackUnsignedInt(&x[12]);
    bAlignment.AlignmentFlag = tempValue >> 16;
    bAlignment.SupportData.NumCigarOperations = tempValue & 0xffff;

    bAlignment.SupportData.QuerySequenceLength = BgzfData::UnpackUnsignedInt(&x[16]);
    bAlignment.MateRefID    = BgzfData::UnpackSignedInt(&x[20]);
    bAlignment.MatePosition = BgzfData::UnpackSignedInt(&x[24]);
    bAlignment.InsertSize   = BgzfData::UnpackSignedInt(&x[28]);
    
    // set BamAlignment length
    bAlignment.Length = bAlignment.SupportData.QuerySequenceLength;
    
    // read in character data - make sure proper data size was read
    bool readCharDataOK = false;
    const unsigned int dataLength = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    char* allCharData = (char*)calloc(sizeof(char), dataLength);
    
    if ( mBGZF.Read(allCharData, dataLength) == (signed int)dataLength) { 
      
        // store 'allCharData' in supportData structure
        bAlignment.SupportData.AllCharData.assign((const char*)allCharData, dataLength);
        
        // set success flag
        readCharDataOK = true;
    }

    free(allCharData);
    return readCharDataOK;
}

// loads reference data from BAM file
void BamReader::BamReaderPrivate::LoadReferenceData(void) {

    // get number of reference sequences
    char buffer[4];
    mBGZF.Read(buffer, 4);
    unsigned int numberRefSeqs = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(numberRefSeqs); }
    if (numberRefSeqs == 0) { return; }
    References.reserve((int)numberRefSeqs);

    // iterate over all references in header
    for (unsigned int i = 0; i != numberRefSeqs; ++i) {

        // get length of reference name
        mBGZF.Read(buffer, 4);
        unsigned int refNameLength = BgzfData::UnpackUnsignedInt(buffer);
        if ( IsBigEndian ) { SwapEndian_32(refNameLength); }
        char* refName = (char*)calloc(refNameLength, 1);

        // get reference name and reference sequence length
        mBGZF.Read(refName, refNameLength);
        mBGZF.Read(buffer, 4);
        int refLength = BgzfData::UnpackSignedInt(buffer);
        if ( IsBigEndian ) { SwapEndian_32(refLength); }

        // store data for reference
        RefData aReference;
        aReference.RefName   = (string)((const char*)refName);
        aReference.RefLength = refLength;
        References.push_back(aReference);

        // clean up calloc-ed temp variable
        free(refName);
    }
}

// opens BAM file (and index)
bool BamReader::BamReaderPrivate::Open(const string& filename, const string& indexFilename) {

    Filename = filename;
    IndexFilename = indexFilename;

    // open the BGZF file for reading, return false on failure
    if ( !mBGZF.Open(filename, "rb") ) 
        return false;
    
    // retrieve header text & reference data
    LoadHeaderData();
    LoadReferenceData();

    // store file offset of first alignment
    AlignmentsBeginOffset = mBGZF.Tell();

    // open index file & load index data (if exists)
    if ( !IndexFilename.empty() )
        LoadIndex();
    
    // return success
    return true;
}

// returns BAM file pointer to beginning of alignment data
bool BamReader::BamReaderPrivate::Rewind(void) {
   
    // rewind to first alignment
    if ( !mBGZF.Seek(AlignmentsBeginOffset) ) return false;
  
    // retrieve first alignment data
    BamAlignment al;
    if ( !LoadNextAlignment(al) ) return false;
      
    // reset default region info using first alignment in file
    Region.LeftRefID      = al.RefID;
    Region.LeftPosition   = al.Position;
    Region.RightRefID     = -1;
    Region.RightPosition  = -1;
    IsLeftBoundSpecified  = false;
    IsRightBoundSpecified = false; 

    // rewind back to before first alignment
    // return success/fail of seek
    return mBGZF.Seek(AlignmentsBeginOffset);
}

// sets a region of interest (with left & right bound reference/position)
// attempts a Jump() to left bound as well
// returns success/failure of Jump()
bool BamReader::BamReaderPrivate::SetRegion(const BamRegion& region) {
    
    // save region of interest
    Region = region;
    
    // set flags
    if ( region.LeftRefID >= 0 && region.LeftPosition >= 0 ) 
        IsLeftBoundSpecified = true;
    if ( region.RightRefID >= 0 && region.RightPosition >= 0 ) 
        IsRightBoundSpecified = true;
    
    // attempt jump to beginning of region, return success/fail of Jump()
    return Jump( Region.LeftRefID, Region.LeftPosition );
}

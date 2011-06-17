// ***************************************************************************
// BamReader_p.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 10 May 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#include <api/BamConstants.h>
#include <api/BamReader.h>
#include <api/internal/BamHeader_p.h>
#include <api/internal/BamRandomAccessController_p.h>
#include <api/internal/BamReader_p.h>
#include <api/internal/BamStandardIndex_p.h>
#include <api/internal/BamToolsIndex_p.h>
#include <api/internal/BgzfStream_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
using namespace std;

// constructor
BamReaderPrivate::BamReaderPrivate(BamReader* parent)
    : m_alignmentsBeginOffset(0)
    , m_parent(parent)
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

// destructor
BamReaderPrivate::~BamReaderPrivate(void) {
    Close();
}

// closes the BAM file
void BamReaderPrivate::Close(void) {

    // clear header & reference data
    m_references.clear();
    m_header.Clear();

    // close internal
    m_randomAccessController.Close();
    m_stream.Close();

    // clear filename
    m_filename.clear();
}

// creates an index file of requested type on current BAM file
bool BamReaderPrivate::CreateIndex(const BamIndex::IndexType& type) {
    if ( !IsOpen() ) return false;
    return m_randomAccessController.CreateIndex(this, type);
}

// return path & filename of current BAM file
const string BamReaderPrivate::Filename(void) const {
    return m_filename;
}

// return header data as std::string
string BamReaderPrivate::GetHeaderText(void) const {
    return m_header.ToString();
}

// return header data as SamHeader object
SamHeader BamReaderPrivate::GetSamHeader(void) const {
    return m_header.ToSamHeader();
}

// get next alignment (with character data fully parsed)
bool BamReaderPrivate::GetNextAlignment(BamAlignment& alignment) {

    // if valid alignment found
    if ( GetNextAlignmentCore(alignment) ) {

        // store alignment's "source" filename
        alignment.Filename = m_filename;

        // return success/failure of parsing char data
        return alignment.BuildCharData();
    }

    // no valid alignment found
    return false;
}

// retrieves next available alignment core data (returns success/fail)
// ** DOES NOT populate any character data fields (read name, bases, qualities, tag data, filename)
//    these can be accessed, if necessary, from the supportData
// useful for operations requiring ONLY positional or other alignment-related information
bool BamReaderPrivate::GetNextAlignmentCore(BamAlignment& alignment) {

    // skip if region is set but has no alignments
    if ( m_randomAccessController.HasRegion() &&
         !m_randomAccessController.RegionHasAlignments() )
    {
        return false;
    }

    // if can't read next alignment
    if ( !LoadNextAlignment(alignment) )
        return false;

    // check alignment's region-overlap state
    BamRandomAccessController::RegionState state = m_randomAccessController.AlignmentState(alignment);

    // if alignment starts after region, no need to keep reading
    if ( state == BamRandomAccessController::AfterRegion )
        return false;

    // read until overlap is found
    while ( state != BamRandomAccessController::OverlapsRegion ) {

        // if can't read next alignment
        if ( !LoadNextAlignment(alignment) )
            return false;

        // check alignment's region-overlap state
        state = m_randomAccessController.AlignmentState(alignment);

        // if alignment starts after region, no need to keep reading
        if ( state == BamRandomAccessController::AfterRegion )
            return false;
    }

    // if we get here, we found the next 'valid' alignment
    // (e.g. overlaps current region if one was set, simply the next alignment if not)
    alignment.SupportData.HasCoreOnly = true;
    return true;
}

int BamReaderPrivate::GetReferenceCount(void) const {
    return m_references.size();
}

const RefVector& BamReaderPrivate::GetReferenceData(void) const {
    return m_references;
}

// returns RefID for given RefName (returns References.size() if not found)
int BamReaderPrivate::GetReferenceID(const string& refName) const {

    // retrieve names from reference data
    vector<string> refNames;
    RefVector::const_iterator refIter = m_references.begin();
    RefVector::const_iterator refEnd  = m_references.end();
    for ( ; refIter != refEnd; ++refIter)
        refNames.push_back( (*refIter).RefName );

    // return 'index-of' refName (or -1 if not found)
    int index = distance(refNames.begin(), find(refNames.begin(), refNames.end(), refName));
    if ( index == (int)m_references.size() ) return -1;
    else return index;
}

bool BamReaderPrivate::HasIndex(void) const {
    return m_randomAccessController.HasIndex();
}

bool BamReaderPrivate::IsOpen(void) const {
    return m_stream.IsOpen;
}

// load BAM header data
bool BamReaderPrivate::LoadHeaderData(void) {
     return m_header.Load(&m_stream);
}

// populates BamAlignment with alignment data under file pointer, returns success/fail
bool BamReaderPrivate::LoadNextAlignment(BamAlignment& alignment) {

    // read in the 'block length' value, make sure it's not zero
    char buffer[sizeof(uint32_t)];
    m_stream.Read(buffer, sizeof(uint32_t));
    alignment.SupportData.BlockLength = BamTools::UnpackUnsignedInt(buffer);
    if ( m_isBigEndian ) BamTools::SwapEndian_32(alignment.SupportData.BlockLength);
    if ( alignment.SupportData.BlockLength == 0 ) return false;

    // read in core alignment data, make sure the right size of data was read
    char x[Constants::BAM_CORE_SIZE];
    if ( m_stream.Read(x, Constants::BAM_CORE_SIZE) != Constants::BAM_CORE_SIZE )
        return false;

    // swap core endian-ness if necessary
    if ( m_isBigEndian ) {
        for ( int i = 0; i < Constants::BAM_CORE_SIZE; i+=sizeof(uint32_t) )
            BamTools::SwapEndian_32p(&x[i]);
    }

    // set BamAlignment 'core' and 'support' data
    alignment.RefID    = BamTools::UnpackSignedInt(&x[0]);
    alignment.Position = BamTools::UnpackSignedInt(&x[4]);

    unsigned int tempValue = BamTools::UnpackUnsignedInt(&x[8]);
    alignment.Bin        = tempValue >> 16;
    alignment.MapQuality = tempValue >> 8 & 0xff;
    alignment.SupportData.QueryNameLength = tempValue & 0xff;

    tempValue = BamTools::UnpackUnsignedInt(&x[12]);
    alignment.AlignmentFlag = tempValue >> 16;
    alignment.SupportData.NumCigarOperations = tempValue & 0xffff;

    alignment.SupportData.QuerySequenceLength = BamTools::UnpackUnsignedInt(&x[16]);
    alignment.MateRefID    = BamTools::UnpackSignedInt(&x[20]);
    alignment.MatePosition = BamTools::UnpackSignedInt(&x[24]);
    alignment.InsertSize   = BamTools::UnpackSignedInt(&x[28]);

    // set BamAlignment length
    alignment.Length = alignment.SupportData.QuerySequenceLength;

    // read in character data - make sure proper data size was read
    bool readCharDataOK = false;
    const unsigned int dataLength = alignment.SupportData.BlockLength - Constants::BAM_CORE_SIZE;
    char* allCharData = (char*)calloc(sizeof(char), dataLength);

    if ( m_stream.Read(allCharData, dataLength) == (signed int)dataLength ) {

        // store 'allCharData' in supportData structure
        alignment.SupportData.AllCharData.assign((const char*)allCharData, dataLength);

        // set success flag
        readCharDataOK = true;

        // save CIGAR ops
        // need to calculate this here so that  BamAlignment::GetEndPosition() performs correctly,
        // even when GetNextAlignmentCore() is called
        const unsigned int cigarDataOffset = alignment.SupportData.QueryNameLength;
        uint32_t* cigarData = (uint32_t*)(allCharData + cigarDataOffset);
        CigarOp op;
        alignment.CigarData.clear();
        alignment.CigarData.reserve(alignment.SupportData.NumCigarOperations);
        for ( unsigned int i = 0; i < alignment.SupportData.NumCigarOperations; ++i ) {

            // swap endian-ness if necessary
            if ( m_isBigEndian ) BamTools::SwapEndian_32(cigarData[i]);

            // build CigarOp structure
            op.Length = (cigarData[i] >> Constants::BAM_CIGAR_SHIFT);
            op.Type   = Constants::BAM_CIGAR_LOOKUP[ (cigarData[i] & Constants::BAM_CIGAR_MASK) ];

            // save CigarOp
            alignment.CigarData.push_back(op);
        }
    }

    // clean up & return parsing success/failure
    free(allCharData);
    return readCharDataOK;
}

// loads reference data from BAM file
bool BamReaderPrivate::LoadReferenceData(void) {

    // get number of reference sequences
    char buffer[sizeof(uint32_t)];
    m_stream.Read(buffer, sizeof(uint32_t));
    uint32_t numberRefSeqs = BamTools::UnpackUnsignedInt(buffer);
    if ( m_isBigEndian ) BamTools::SwapEndian_32(numberRefSeqs);
    m_references.reserve((int)numberRefSeqs);

    // iterate over all references in header
    for ( unsigned int i = 0; i != numberRefSeqs; ++i ) {

        // get length of reference name
        m_stream.Read(buffer, sizeof(uint32_t));
        uint32_t refNameLength = BamTools::UnpackUnsignedInt(buffer);
        if ( m_isBigEndian ) BamTools::SwapEndian_32(refNameLength);
        char* refName = (char*)calloc(refNameLength, 1);

        // get reference name and reference sequence length
        m_stream.Read(refName, refNameLength);
        m_stream.Read(buffer, sizeof(int32_t));
        int32_t refLength = BamTools::UnpackSignedInt(buffer);
        if ( m_isBigEndian ) BamTools::SwapEndian_32(refLength);

        // store data for reference
        RefData aReference;
        aReference.RefName   = (string)((const char*)refName);
        aReference.RefLength = refLength;
        m_references.push_back(aReference);

        // clean up calloc-ed temp variable
        free(refName);
    }

    // return success
    return true;
}

bool BamReaderPrivate::LocateIndex(const BamIndex::IndexType& preferredType) {
    return m_randomAccessController.LocateIndex(this, preferredType);
}

// opens BAM file (and index)
bool BamReaderPrivate::Open(const string& filename) {

    // close current BAM file if open
    if ( m_stream.IsOpen )
        Close();

    // attempt to open BgzfStream for reading
    if ( !m_stream.Open(filename, "rb") ) {
        cerr << "BamReader ERROR: Could not open BGZF stream for " << filename << endl;
        return false;
    }

    // attempt to load header data
    if ( !LoadHeaderData() ) {
        cerr << "BamReader ERROR: Could not load header data for " << filename << endl;
        Close();
        return false;
    }

    // attempt to load reference data
    if ( !LoadReferenceData() ) {
        cerr << "BamReader ERROR: Could not load reference data for " << filename << endl;
        Close();
        return false;
    }

    // if all OK, store filename & offset of first alignment
    m_filename = filename;
    m_alignmentsBeginOffset = m_stream.Tell();

    // return success
    return true;
}

bool BamReaderPrivate::OpenIndex(const std::string& indexFilename) {
    return m_randomAccessController.OpenIndex(indexFilename, this);
}

// returns BAM file pointer to beginning of alignment data
bool BamReaderPrivate::Rewind(void) {

    // attempt rewind to first alignment
    if ( !m_stream.Seek(m_alignmentsBeginOffset) )
        return false;

    // verify that we can read first alignment
    BamAlignment al;
    if ( !LoadNextAlignment(al) )
        return false;

    // reset region
    m_randomAccessController.ClearRegion();

    // rewind back to beginning of first alignment
    // return success/fail of seek
    return m_stream.Seek(m_alignmentsBeginOffset);
}

bool BamReaderPrivate::Seek(const int64_t& position) {
    return m_stream.Seek(position);
}

void BamReaderPrivate::SetIndex(BamIndex* index) {
    m_randomAccessController.SetIndex(index);
}

// change the index caching behavior
void BamReaderPrivate::SetIndexCacheMode(const BamIndex::IndexCacheMode& mode) {
    m_randomAccessController.SetIndexCacheMode(mode);
}

// sets current region & attempts to jump to it
// returns success/failure
bool BamReaderPrivate::SetRegion(const BamRegion& region) {
    return m_randomAccessController.SetRegion(this, region, m_references.size());
}

int64_t BamReaderPrivate::Tell(void) const {
    return m_stream.Tell();
}

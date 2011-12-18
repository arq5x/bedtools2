// ***************************************************************************
// BamReader_p.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 28 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#include "api/BamConstants.h"
#include "api/BamReader.h"
#include "api/IBamIODevice.h"
#include "api/internal/bam/BamHeader_p.h"
#include "api/internal/bam/BamRandomAccessController_p.h"
#include "api/internal/bam/BamReader_p.h"
#include "api/internal/index/BamStandardIndex_p.h"
#include "api/internal/index/BamToolsIndex_p.h"
#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <cassert>
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
bool BamReaderPrivate::Close(void) {

    // clear BAM metadata
    m_references.clear();
    m_header.Clear();

    // clear filename
    m_filename.clear();

    // close random access controller
    m_randomAccessController.Close();

    // if stream is open, attempt close
    if ( IsOpen() ) {
        try {
            m_stream.Close();
        } catch ( BamException& e ) {
            const string streamError = e.what();
            const string message = string("encountered error closing BAM file: \n\t") + streamError;
            SetErrorString("BamReader::Close", message);
            return false;
        }
    }

    // return success
    return true;
}

// creates an index file of requested type on current BAM file
bool BamReaderPrivate::CreateIndex(const BamIndex::IndexType& type) {

    // skip if BAM file not open
    if ( !IsOpen() ) {
        SetErrorString("BamReader::CreateIndex", "cannot create index on unopened BAM file");
        return false;
    }

    // attempt to create index
    if ( m_randomAccessController.CreateIndex(this, type) )
        return true;
    else {
        const string bracError = m_randomAccessController.GetErrorString();
        const string message = string("could not create index: \n\t") + bracError;
        SetErrorString("BamReader::CreateIndex", message);
        return false;
    }
}

// return path & filename of current BAM file
const string BamReaderPrivate::Filename(void) const {
    return m_filename;
}

string BamReaderPrivate::GetErrorString(void) const {
    return m_errorString;
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
        if ( alignment.BuildCharData() )
            return true;
        else {
            const string alError = alignment.GetErrorString();
            const string message = string("could not populate alignment data: \n\t") + alError;
            SetErrorString("BamReader::GetNextAlignment", message);
            return false;
        }
    }

    // no valid alignment found
    return false;
}

// retrieves next available alignment core data (returns success/fail)
// ** DOES NOT populate any character data fields (read name, bases, qualities, tag data, filename)
//    these can be accessed, if necessary, from the supportData
// useful for operations requiring ONLY positional or other alignment-related information
bool BamReaderPrivate::GetNextAlignmentCore(BamAlignment& alignment) {

    // skip if stream not opened
    if ( !m_stream.IsOpen() )
        return false;

    try {

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

    } catch ( BamException& e ) {
        const string streamError = e.what();
        const string message = string("encountered error reading BAM alignment: \n\t") + streamError;
        SetErrorString("BamReader::GetNextAlignmentCore", message);
        return false;
    }
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
    return m_stream.IsOpen();
}

// load BAM header data
void BamReaderPrivate::LoadHeaderData(void) {
    m_header.Load(&m_stream);
}

// populates BamAlignment with alignment data under file pointer, returns success/fail
bool BamReaderPrivate::LoadNextAlignment(BamAlignment& alignment) {

    // read in the 'block length' value, make sure it's not zero
    char buffer[sizeof(uint32_t)];
    m_stream.Read(buffer, sizeof(uint32_t));
    alignment.SupportData.BlockLength = BamTools::UnpackUnsignedInt(buffer);
    if ( m_isBigEndian ) BamTools::SwapEndian_32(alignment.SupportData.BlockLength);
    if ( alignment.SupportData.BlockLength == 0 )
        return false;

    // read in core alignment data, make sure the right size of data was read
    char x[Constants::BAM_CORE_SIZE];
    if ( m_stream.Read(x, Constants::BAM_CORE_SIZE) != Constants::BAM_CORE_SIZE )
        return false;

    // swap core endian-ness if necessary
    if ( m_isBigEndian ) {
        for ( unsigned int i = 0; i < Constants::BAM_CORE_SIZE; i+=sizeof(uint32_t) )
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
    RaiiBuffer allCharData(dataLength);

    if ( m_stream.Read(allCharData.Buffer, dataLength) == dataLength ) {

        // store 'allCharData' in supportData structure
        alignment.SupportData.AllCharData.assign((const char*)allCharData.Buffer, dataLength);

        // set success flag
        readCharDataOK = true;

        // save CIGAR ops
        // need to calculate this here so that  BamAlignment::GetEndPosition() performs correctly,
        // even when GetNextAlignmentCore() is called
        const unsigned int cigarDataOffset = alignment.SupportData.QueryNameLength;
        uint32_t* cigarData = (uint32_t*)(allCharData.Buffer + cigarDataOffset);
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

    // return success/failure
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
        RaiiBuffer refName(refNameLength);

        // get reference name and reference sequence length
        m_stream.Read(refName.Buffer, refNameLength);
        m_stream.Read(buffer, sizeof(int32_t));
        int32_t refLength = BamTools::UnpackSignedInt(buffer);
        if ( m_isBigEndian ) BamTools::SwapEndian_32(refLength);

        // store data for reference
        RefData aReference;
        aReference.RefName   = (string)((const char*)refName.Buffer);
        aReference.RefLength = refLength;
        m_references.push_back(aReference);
    }

    // return success
    return true;
}

bool BamReaderPrivate::LocateIndex(const BamIndex::IndexType& preferredType) {

    if ( m_randomAccessController.LocateIndex(this, preferredType) )
        return true;
    else {
        const string bracError = m_randomAccessController.GetErrorString();
        const string message = string("could not locate index: \n\t") + bracError;
        SetErrorString("BamReader::LocateIndex", message);
        return false;
    }
}

// opens BAM file (and index)
bool BamReaderPrivate::Open(const string& filename) {

    try {

        // make sure we're starting with fresh state
        Close();

        // open BgzfStream
        m_stream.Open(filename, IBamIODevice::ReadOnly);

        // load BAM metadata
        LoadHeaderData();
        LoadReferenceData();

        // store filename & offset of first alignment
        m_filename = filename;
        m_alignmentsBeginOffset = m_stream.Tell();

        // return success
        return true;

    } catch ( BamException& e ) {
        const string error = e.what();
        const string message = string("could not open file: ") + filename +
                               "\n\t" + error;
        SetErrorString("BamReader::Open", message);
        return false;
    }
}

bool BamReaderPrivate::OpenIndex(const std::string& indexFilename) {

    if ( m_randomAccessController.OpenIndex(indexFilename, this) )
        return true;
    else {
        const string bracError = m_randomAccessController.GetErrorString();
        const string message = string("could not open index: \n\t") + bracError;
        SetErrorString("BamReader::OpenIndex", message);
        return false;
    }
}

// returns BAM file pointer to beginning of alignment data
bool BamReaderPrivate::Rewind(void) {

    // reset region
    m_randomAccessController.ClearRegion();

    // return status of seeking back to first alignment
    if ( Seek(m_alignmentsBeginOffset) )
        return true;
    else {
        const string currentError = m_errorString;
        const string message = string("could not rewind: \n\t") + currentError;
        SetErrorString("BamReader::Rewind", message);
        return false;
    }
}

bool BamReaderPrivate::Seek(const int64_t& position) {

    // skip if BAM file not open
    if ( !IsOpen() ) {
        SetErrorString("BamReader::Seek", "cannot seek on unopened BAM file");
        return false;
    }

    try {
        m_stream.Seek(position);
        return true;
    }
    catch ( BamException& e ) {
        const string streamError = e.what();
        const string message = string("could not seek in BAM file: \n\t") + streamError;
        SetErrorString("BamReader::Seek", message);
        return false;
    }
}

void BamReaderPrivate::SetErrorString(const string& where, const string& what) {
    static const string SEPARATOR = ": ";
    m_errorString = where + SEPARATOR + what;
}

void BamReaderPrivate::SetIndex(BamIndex* index) {
    m_randomAccessController.SetIndex(index);
}

// sets current region & attempts to jump to it
// returns success/failure
bool BamReaderPrivate::SetRegion(const BamRegion& region) {

    if ( m_randomAccessController.SetRegion(region, m_references.size()) )
        return true;
    else {
        const string bracError = m_randomAccessController.GetErrorString();
        const string message = string("could not set region: \n\t") + bracError;
        SetErrorString("BamReader::SetRegion", message);
        return false;
    }
}

int64_t BamReaderPrivate::Tell(void) const {
    return m_stream.Tell();
}

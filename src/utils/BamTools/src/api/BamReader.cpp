// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides read access to BAM files.
// ***************************************************************************

#include "api/BamReader.h"
#include "api/internal/bam/BamReader_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
using namespace std;

/*! \class BamTools::BamReader
    \brief Provides read access to BAM files.
*/

/*! \fn BamReader::BamReader(void)
    \brief constructor
*/
BamReader::BamReader(void)
    : d(new BamReaderPrivate(this))
{ }

/*! \fn BamReader::~BamReader(void)
    \brief destructor
*/
BamReader::~BamReader(void) {
    delete d;
    d = 0;
}

/*! \fn bool BamReader::Close(void)
    \brief Closes the current BAM file.

    Also clears out all header and reference data.

    \return \c true if file closed OK
    \sa IsOpen(), Open()
*/
bool BamReader::Close(void) {
    return d->Close();
}

/*! \fn bool BamReader::CreateIndex(const BamIndex::IndexType& type)
    \brief Creates an index file for current BAM file.

    \param[in] type file format to create, see BamIndex::IndexType for available formats
    \return \c true if index created OK
    \sa LocateIndex(), OpenIndex()
*/
bool BamReader::CreateIndex(const BamIndex::IndexType& type) {
    return d->CreateIndex(type);
}

/*! \fn std::string BamReader::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
string BamReader::GetErrorString(void) const {
    return d->GetErrorString();
}

/*! \fn const std::string BamReader::GetFilename(void) const
    \brief Returns name of current BAM file.

    Retrieved filename will contain whatever was passed via Open().
    If you need full directory paths here, be sure to include them
    when you open the BAM file.

    \returns name of open BAM file. If no file is open, returns an empty string.
    \sa IsOpen()
*/
const std::string BamReader::GetFilename(void) const {
    return d->Filename();
}

/*! \fn SamHeader BamReader::GetHeader(void) const
    \brief Returns SAM header data.

    Header data is wrapped in a SamHeader object that can be conveniently queried & modified.

    \note Modifying the retrieved SamHeader object does NOT affect the
    current BAM file. This file has been opened in a read-only mode.
    However, your modified SamHeader object can be used in conjunction with
    BamWriter to generate a new BAM file with the appropriate header information.

    \returns header data object
    \sa GetHeaderText()
*/
SamHeader BamReader::GetHeader(void) const {
    return d->GetSamHeader();
}

/*! \fn std::string BamReader::GetHeaderText(void) const
    \brief Returns SAM header data, as SAM-formatted text.

    \note Modifying the retrieved text does NOT affect the current
    BAM file. This file has been opened in a read-only mode. However,
    your modified header text can be used in conjunction with BamWriter
    to generate a new BAM file with the appropriate header information.

    \returns SAM-formatted header text
    \sa GetHeader()
*/
std::string BamReader::GetHeaderText(void) const {
    return d->GetHeaderText();
}

/*! \fn bool BamReader::GetNextAlignment(BamAlignment& alignment)
    \brief Retrieves next available alignment.

    Attempts to read the next alignment record from BAM file, and checks to see
    if it overlaps the current region. If no region is currently set, then the
    next alignment available is always considered valid.

    If a region has been set, via Jump() or SetRegion(), an alignment is only
    considered valid if it overlaps the region. If the actual 'next' alignment record
    in the BAM file does not overlap this region, then this function will read sequentially
    through the file until the next alignment that overlaps this region is found.
    Once the region has been exhausted (i.e. the next alignment loaded is beyond the region),
    the function aborts and returns \c false. In this case, there is no point to continue
    reading, assuming properly sorted alignments.

    This function fully populates all of the alignment's available data fields,
    including the string data fields (read name, bases, qualities, tags, filename).
    If only positional data (refID, position, CIGAR ops, alignment flags, etc.)
    are required, consider using GetNextAlignmentCore() for a significant
    performance boost.

    \param[out] alignment destination for alignment record data
    \returns \c true if a valid alignment was found
*/
bool BamReader::GetNextAlignment(BamAlignment& alignment) {
    return d->GetNextAlignment(alignment);
}

/*! \fn bool BamReader::GetNextAlignmentCore(BamAlignment& alignment)
    \brief Retrieves next available alignment, without populating the alignment's string data fields.

    Equivalent to GetNextAlignment() with respect to what is a valid overlapping alignment.

    However, this method does NOT populate the alignment's string data fields
    (read name, bases, qualities, tags, filename). This provides a boost in speed
    when these fields are not required for every alignment. These fields can be
    populated 'lazily' (as needed) by calling BamAlignment::BuildCharData() later.

    \param[out] alignment destination for alignment record data
    \returns \c true if a valid alignment was found
    \sa SetRegion()
*/
bool BamReader::GetNextAlignmentCore(BamAlignment& alignment) {
    return d->GetNextAlignmentCore(alignment);
}

/*! \fn int BamReader::GetReferenceCount(void) const
    \brief Returns number of reference sequences.
*/
int BamReader::GetReferenceCount(void) const {
    return d->GetReferenceCount();
}

/*! \fn const RefVector& BamReader::GetReferenceData(void) const
    \brief Returns all reference sequence entries.
    \sa RefData
*/
const RefVector& BamReader::GetReferenceData(void) const {
    return d->GetReferenceData();
}

/*! \fn int BamReader::GetReferenceID(const std::string& refName) const
    \brief Returns the ID of the reference with this name.

    If \a refName is not found, returns -1.

    \param[in] refName name of reference to look up
*/
int BamReader::GetReferenceID(const std::string& refName) const {
    return d->GetReferenceID(refName);
}

/*! \fn bool BamReader::HasIndex(void) const
    \brief Returns \c true if index data is available.
*/
bool BamReader::HasIndex(void) const {
    return d->HasIndex();
}

/*! \fn bool BamReader::IsOpen(void) const
    \brief Returns \c true if a BAM file is open for reading.
*/
bool BamReader::IsOpen(void) const {
    return d->IsOpen();
}

/*! \fn bool BamReader::Jump(int refID, int position)
    \brief Performs a random-access jump within BAM file.

    This is a convenience method, equivalent to calling SetRegion()
    with only a left boundary specified.

    \param[in] refID    left-bound reference ID
    \param[in] position left-bound position

    \returns \c true if jump was successful
    \sa HasIndex()
*/
bool BamReader::Jump(int refID, int position) {
    return d->SetRegion( BamRegion(refID, position) );
}

/*! \fn bool BamReader::LocateIndex(const BamIndex::IndexType& preferredType)
    \brief Looks in BAM file's directory for a matching index file.

    Use this function when you need an index file, and perhaps have a
    preferred index format, but do not depend heavily on which format
    actually gets loaded at runtime.

    This function will defer to your \a preferredType whenever possible.
    However, if an index file of \a preferredType can not be found, then
    it will look for any other index file that corresponds to this BAM file.

    If you want precise control over which index file is loaded, use OpenIndex()
    with the desired index filename. If that function returns false, you can use
    CreateIndex() to then build an index of the exact requested format.

    \param[in] preferredType desired index file format, see BamIndex::IndexType for available formats

    \returns \c true if (any) index file could be found
*/
bool BamReader::LocateIndex(const BamIndex::IndexType& preferredType) {
    return d->LocateIndex(preferredType);
}

/*! \fn bool BamReader::Open(const std::string& filename)
    \brief Opens a BAM file.

    If BamReader is already opened on another file, this function closes
    that file, then attempts to open requested \a filename.

    \param[in] filename name of BAM file to open

    \returns \c true if BAM file was opened successfully
    \sa Close(), IsOpen(), OpenIndex()
*/
bool BamReader::Open(const std::string& filename) {
    return d->Open(filename);
}

/*! \fn bool BamReader::OpenIndex(const std::string& indexFilename)
    \brief Opens a BAM index file.

    \param[in] indexFilename name of BAM index file to open

    \returns \c true if BAM index file was opened & data loaded successfully
    \sa LocateIndex(), Open(), SetIndex()
*/
bool BamReader::OpenIndex(const std::string& indexFilename) {
    return d->OpenIndex(indexFilename);
}

/*! \fn bool BamReader::Rewind(void)
    \brief Returns the internal file pointer to the first alignment record.

    Useful for performing multiple sequential passes through a BAM file.
    Calling this function clears any prior region that may have been set.

    \note This function sets the file pointer to first alignment record
    in the BAM file, NOT the beginning of the file.

    \returns \c true if rewind operation was successful
    \sa Jump(), SetRegion()
*/
bool BamReader::Rewind(void) {
    return d->Rewind();
}

/*! \fn void BamReader::SetIndex(BamIndex* index)
    \brief Sets a custom BamIndex on this reader.

    Only necessary for custom BamIndex subclasses. Most clients should
    never have to use this function.

    Example:
    \code
        BamReader reader;
        reader.SetIndex(new MyCustomBamIndex);
    \endcode

    \note BamReader takes ownership of \a index - i.e. the BamReader will
    take care of deleting it when the reader is destructed, when the current
    BAM file is closed, or when a new index is requested.

    \param[in] index custom BamIndex subclass created by client
    \sa CreateIndex(), LocateIndex(), OpenIndex()
*/
void BamReader::SetIndex(BamIndex* index) {
    d->SetIndex(index);
}

/*! \fn bool BamReader::SetRegion(const BamRegion& region)
    \brief Sets a target region of interest

    Requires that index data be available. Attempts a random-access
    jump in the BAM file, near \a region left boundary position.

    Subsequent calls to GetNextAlignment() or GetNextAlignmentCore()
    will only return \c true when alignments can be found that overlap
    this \a region.

    A \a region with no right boundary is considered open-ended, meaning
    that all alignments that lie downstream of the left boundary are
    considered valid, continuing to the end of the BAM file.

    \warning BamRegion now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] region desired region-of-interest to activate

    \returns \c true if reader was able to jump successfully to the region's left boundary
    \sa HasIndex(), Jump()
*/
bool BamReader::SetRegion(const BamRegion& region) {
    return d->SetRegion(region);
}

/*! \fn bool BamReader::SetRegion(const int& leftRefID,
                                  const int& leftPosition,
                                  const int& rightRefID,
                                  const int& rightPosition)
    \brief Sets a target region of interest.

    This is an overloaded function.

    \warning This function expects a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] leftRefID     referenceID of region's left boundary
    \param[in] leftPosition  position of region's left boundary
    \param[in] rightRefID    reference ID of region's right boundary
    \param[in] rightPosition position of region's right boundary

    \returns \c true if reader was able to jump successfully to the region's left boundary
    \sa HasIndex(), Jump()
*/
bool BamReader::SetRegion(const int& leftRefID,
                          const int& leftBound,
                          const int& rightRefID,
                          const int& rightBound)
{
    return d->SetRegion( BamRegion(leftRefID, leftBound, rightRefID, rightBound) );
}

// ***************************************************************************
// BamMultiReader.cpp (c) 2010 Erik Garrison, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Convenience class for reading multiple BAM files.
//
// This functionality allows applications to work on very large sets of files
// without requiring intermediate merge, sort, and index steps for each file
// subset. It also improves the performance of our merge system as it
// precludes the need to sort merged files.
// ***************************************************************************

#include "api/BamMultiReader.h"
#include "api/internal/bam/BamMultiReader_p.h"
using namespace BamTools;

#include <string>
#include <vector>
using namespace std;

/*! \class BamTools::BamMultiReader
    \brief Convenience class for reading multiple BAM files.
*/

/*! \fn BamMultiReader::BamMultiReader(void)
    \brief constructor
*/
BamMultiReader::BamMultiReader(void)
    : d(new Internal::BamMultiReaderPrivate)
{ }

/*! \fn BamMultiReader::~BamMultiReader(void)
    \brief destructor
*/
BamMultiReader::~BamMultiReader(void) {
    delete d;
    d = 0;
}

/*! \fn void BamMultiReader::Close(void)
    \brief Closes all open BAM files.

    Also clears out all header and reference data.

    \sa CloseFile(), IsOpen(), Open(), BamReader::Close()
*/
bool BamMultiReader::Close(void) {
    return d->Close();
}

/*! \fn void BamMultiReader::CloseFile(const std::string& filename)
    \brief Closes requested BAM file.

    Leaves any other file(s) open, along with header and reference data.

    \param[in] filename name of specific BAM file to close

    \sa Close(), IsOpen(), Open(), BamReader::Close()
*/
bool BamMultiReader::CloseFile(const std::string& filename) {
    return d->CloseFile(filename);
}

/*! \fn bool BamMultiReader::CreateIndexes(const BamIndex::IndexType& type)
    \brief Creates index files for the current BAM files.

    \param[in] type file format to create, see BamIndex::IndexType for available formats
    \return \c true if index files created OK
    \sa LocateIndexes(), OpenIndexes(), BamReader::CreateIndex()
*/
bool BamMultiReader::CreateIndexes(const BamIndex::IndexType& type) {
    return d->CreateIndexes(type);
}

/*! \fn const std::vector<std::string> BamMultiReader::Filenames(void) const
    \brief Returns list of filenames for all open BAM files.

    Retrieved filenames will contain whatever was passed via Open().
    If you need full directory paths here, be sure to include them
    when you open the BAM files.

    \returns names of open BAM files. If no files are open, returns an empty vector.
    \sa IsOpen(), BamReader::GetFilename()
*/
const std::vector<std::string> BamMultiReader::Filenames(void) const {
    return d->Filenames();
}

/*! \fn std::string BamMultiReader::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamMultiReader::GetErrorString(void) const {
    return d->GetErrorString();
}

/*! \fn SamHeader BamMultiReader::GetHeader(void) const
    \brief Returns unified SAM-format header for all files

    \note Modifying the retrieved text does NOT affect the current
    BAM files. These files have been opened in a read-only mode. However,
    your modified header text can be used in conjunction with BamWriter
    to generate a new BAM file with the appropriate header information.

    \returns header data wrapped in SamHeader object
    \sa GetHeaderText(), BamReader::GetHeader()
*/
SamHeader BamMultiReader::GetHeader(void) const {
    return d->GetHeader();
}

/*! \fn std::string BamMultiReader::GetHeaderText(void) const
    \brief Returns unified SAM-format header text for all files

    \note Modifying the retrieved text does NOT affect the current
    BAM files. These files have been opened in a read-only mode. However,
    your modified header text can be used in conjunction with BamWriter
    to generate a new BAM file with the appropriate header information.

    \returns SAM-formatted header text
    \sa GetHeader(), BamReader::GetHeaderText()
*/
std::string BamMultiReader::GetHeaderText(void) const {
    return d->GetHeaderText();
}

/*! \fn bool BamMultiReader::GetNextAlignment(BamAlignment& alignment)
    \brief Retrieves next available alignment.

    Equivalent to BamReader::GetNextAlignment() with respect to what is a valid
    overlapping alignment and what data gets populated.

    This method takes care of determining which alignment actually is 'next'
    across multiple files, depending on their sort order.

    \param[out] alignment destination for alignment record data
    \returns \c true if a valid alignment was found
    \sa GetNextAlignmentCore(), SetRegion(), BamReader::GetNextAlignment()
*/
bool BamMultiReader::GetNextAlignment(BamAlignment& nextAlignment) {
    return d->GetNextAlignment(nextAlignment);
}

/*! \fn bool BamMultiReader::GetNextAlignmentCore(BamAlignment& alignment)
    \brief Retrieves next available alignment.

    Equivalent to BamReader::GetNextAlignmentCore() with respect to what is a valid
    overlapping alignment and what data gets populated.

    This method takes care of determining which alignment actually is 'next'
    across multiple files, depending on their sort order.

    \param[out] alignment destination for alignment record data
    \returns \c true if a valid alignment was found
    \sa GetNextAlignment(), SetRegion(), BamReader::GetNextAlignmentCore()
*/
bool BamMultiReader::GetNextAlignmentCore(BamAlignment& nextAlignment) {
    return d->GetNextAlignmentCore(nextAlignment);
}

/*! \fn int BamMultiReader::GetReferenceCount(void) const
    \brief Returns number of reference sequences.
    \sa BamReader::GetReferenceCount()
*/
int BamMultiReader::GetReferenceCount(void) const {
    return d->GetReferenceCount();
}

/*! \fn const RefVector& BamMultiReader::GetReferenceData(void) const
    \brief Returns all reference sequence entries.
    \sa RefData, BamReader::GetReferenceData()
*/
const BamTools::RefVector BamMultiReader::GetReferenceData(void) const {
    return d->GetReferenceData();
}

/*! \fn int BamMultiReader::GetReferenceID(const std::string& refName) const
    \brief Returns the ID of the reference with this name.

    If \a refName is not found, returns -1.

    \param[in] refName name of reference to look up
    \sa BamReader::GetReferenceID()
*/
int BamMultiReader::GetReferenceID(const std::string& refName) const {
    return d->GetReferenceID(refName);
}

/*! \fn bool BamMultiReader::HasIndexes(void) const
    \brief Returns \c true if all BAM files have index data available.
    \sa BamReader::HasIndex()
*/
bool BamMultiReader::HasIndexes(void) const {
    return d->HasIndexes();
}

/*! \fn bool BamMultiReader::HasOpenReaders(void) const
    \brief Returns \c true if there are any open BAM files.
*/
bool BamMultiReader::HasOpenReaders(void) const {
    return d->HasOpenReaders();
}

/*! \fn bool BamMultiReader::Jump(int refID, int position)
    \brief Performs a random-access jump within current BAM files.

    This is a convenience method, equivalent to calling SetRegion()
    with only a left boundary specified.

    \param[in] refID    ID of reference to jump to
    \param[in] position (0-based) left boundary

    \returns \c true if jump was successful
    \sa HasIndex(), BamReader::Jump()
*/

bool BamMultiReader::Jump(int refID, int position) {
    return d->Jump(refID, position);
}

/*! \fn bool BamMultiReader::LocateIndexes(const BamIndex::IndexType& preferredType)
    \brief Looks for index files that match current BAM files.

    Use this function when you need index files, and perhaps have a
    preferred index format, but do not depend heavily on which indexes
    actually get loaded at runtime.

    For each BAM file, this function will defer to your \a preferredType
    whenever possible. However, if an index file of \a preferredType can
    not be found, then it will look for any other index file that matches
    that BAM file.

    An example case would look this:
    \code
        BamMultiReader reader;

        // do setup...

        // ensure that all files have an index
        if ( !reader.LocateIndexes() )      // opens any existing index files that match our BAM files
            reader.CreateIndexes();         // creates index files for any BAM files that still lack one

        // do interesting stuff using random-access...

    \endcode

    If you want precise control over which index files are loaded, use OpenIndexes()
    with the desired index filenames. If that function returns false, you can use
    CreateIndexes() to then build index files of the exact requested format.

    \param[in] preferredType desired index file format, see BamIndex::IndexType for available formats
    \returns \c true if index files could be found for \b ALL open BAM files
    \sa BamReader::LocateIndex()
*/
bool BamMultiReader::LocateIndexes(const BamIndex::IndexType& preferredType) {
    return d->LocateIndexes(preferredType);
}

/*! \fn bool BamMultiReader::Open(const std::vector<std::string>& filenames)
    \brief Opens BAM files.

    \note Opening BAM files will invalidate any current region set on the multireader.
    All file pointers will be returned to the beginning of the alignment data. Follow
    this with Jump() or SetRegion() to establish a region of interest.

    \param[in] filenames list of BAM filenames to open
    \returns \c true if BAM files were opened successfully
    \sa Close(), HasOpenReaders(), OpenFile(), OpenIndexes(), BamReader::Open()
*/
bool BamMultiReader::Open(const std::vector<std::string>& filenames) {
    return d->Open(filenames);
}

/*! \fn bool BamMultiReader::OpenFile(const std::string& filename)
    \brief Opens a single BAM file.

    Adds another BAM file to multireader "on-the-fly".

    \note Opening a BAM file will invalidate any current region set on the multireader.
    All file pointers will be returned to the beginning of the alignment data. Follow
    this with Jump() or SetRegion() to establish a region of interest.

    \param[in] filename BAM filename to open
    \returns \c true if BAM file was opened successfully
    \sa Close(), HasOpenReaders(), Open(), OpenIndexes(), BamReader::Open()
*/
bool BamMultiReader::OpenFile(const std::string& filename) {
    return d->OpenFile(filename);
}

/*! \fn bool BamMultiReader::OpenIndexes(const std::vector<std::string>& indexFilenames)
    \brief Opens index files for current BAM files.

    \note Currently assumes that index filenames match the order (and number) of
    BAM files passed to Open().

    \param[in] indexFilenames list of BAM index file names
    \returns \c true if BAM index file was opened & data loaded successfully
    \sa LocateIndex(), Open(), SetIndex(), BamReader::OpenIndex()
*/
bool BamMultiReader::OpenIndexes(const std::vector<std::string>& indexFilenames) {
    return d->OpenIndexes(indexFilenames);
}

/*! \fn bool BamMultiReader::Rewind(void)
    \brief Returns the internal file pointers to the beginning of alignment records.

    Useful for performing multiple sequential passes through BAM files.
    Calling this function clears any prior region that may have been set.

    \returns \c true if rewind operation was successful
    \sa Jump(), SetRegion(), BamReader::Rewind()
*/
bool BamMultiReader::Rewind(void) {
    return d->Rewind();
}

/*! \fn bool BamMultiReader::SetRegion(const BamRegion& region)
    \brief Sets a target region of interest

    Equivalent to calling BamReader::SetRegion() on all open BAM files.

    \warning BamRegion now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] region desired region-of-interest to activate
    \returns \c true if ALL readers set the region successfully
    \sa HasIndexes(), Jump(), BamReader::SetRegion()
*/
bool BamMultiReader::SetRegion(const BamRegion& region) {
    return d->SetRegion(region);
}

/*! \fn bool BamMultiReader::SetRegion(const int& leftRefID,
                                       const int& leftPosition,
                                       const int& rightRefID,
                                       const int& rightPosition)
    \brief Sets a target region of interest

    This is an overloaded function. Equivalent to calling BamReader::SetRegion() on all open BAM files.

    \warning This function now expects a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] leftRefID     referenceID of region's left boundary
    \param[in] leftPosition  position of region's left boundary
    \param[in] rightRefID    reference ID of region's right boundary
    \param[in] rightPosition position of region's right boundary

    \returns \c true if ALL readers set the region successfully
    \sa HasIndexes(), Jump(), BamReader::SetRegion()
*/
bool BamMultiReader::SetRegion(const int& leftRefID,
                               const int& leftPosition,
                               const int& rightRefID,
                               const int& rightPosition)
{
    return d->SetRegion( BamRegion(leftRefID, leftPosition, rightRefID, rightPosition) );
}

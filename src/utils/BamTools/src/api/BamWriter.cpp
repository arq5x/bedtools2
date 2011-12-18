// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "api/SamHeader.h"
#include "api/internal/bam/BamWriter_p.h"
using namespace BamTools;
using namespace BamTools::Internal;
using namespace std;

/*! \class BamTools::BamWriter
    \brief Provides write access for generating BAM files.
*/
/*! \enum BamTools::BamWriter::CompressionMode
    \brief This enum describes the compression behaviors for output BAM files.
*/
/*! \var BamWriter::CompressionMode BamWriter::Compressed
    \brief Use normal BAM compression
*/
/*! \var BamWriter::CompressionMode BamWriter::Uncompressed
    \brief Disable BAM compression

    Useful in situations where the BAM data is streamed (e.g. piping).
    It would be wasteful to compress, and then immediately decompress
    the data.
*/

/*! \fn BamWriter::BamWriter(void)
    \brief constructor
*/
BamWriter::BamWriter(void)
    : d(new BamWriterPrivate)
{ }

/*! \fn BamWriter::~BamWriter(void)
    \brief destructor
*/
BamWriter::~BamWriter(void) {
    delete d;
    d = 0;
}

/*! \fn BamWriter::Close(void)
    \brief Closes the current BAM file.
    \sa Open()
*/
void BamWriter::Close(void) {
    d->Close();
}

/*! \fn std::string BamWriter::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamWriter::GetErrorString(void) const {
    return d->GetErrorString();
}

/*! \fn bool BamWriter::IsOpen(void) const
    \brief Returns \c true if BAM file is open for writing.
    \sa Open()
*/
bool BamWriter::IsOpen(void) const {
    return d->IsOpen();
}

/*! \fn bool BamWriter::Open(const std::string& filename,
                             const std::string& samHeaderText,
                             const RefVector& referenceSequences)
    \brief Opens a BAM file for writing.

    Will overwrite the BAM file if it already exists.

    \param[in] filename           name of output BAM file
    \param[in] samHeaderText      header data, as SAM-formatted string
    \param[in] referenceSequences list of reference entries

    \return \c true if opened successfully
    \sa Close(), IsOpen(), BamReader::GetHeaderText(), BamReader::GetReferenceData()
*/
bool BamWriter::Open(const std::string& filename,
                     const std::string& samHeaderText,
                     const RefVector& referenceSequences)
{
    return d->Open(filename, samHeaderText, referenceSequences);
}

/*! \fn bool BamWriter::Open(const std::string& filename,
                             const SamHeader& samHeader,
                             const RefVector& referenceSequences)
    \brief Opens a BAM file for writing.

    This is an overloaded function.

    Will overwrite the BAM file if it already exists.

    \param[in] filename           name of output BAM file
    \param[in] samHeader          header data, wrapped in SamHeader object
    \param[in] referenceSequences list of reference entries

    \return \c true if opened successfully
    \sa Close(), IsOpen(), BamReader::GetHeader(), BamReader::GetReferenceData()
*/
bool BamWriter::Open(const std::string& filename,
                     const SamHeader& samHeader,
                     const RefVector& referenceSequences)
{
    return d->Open(filename, samHeader.ToString(), referenceSequences);
}

/*! \fn void BamWriter::SaveAlignment(const BamAlignment& alignment)
    \brief Saves an alignment to the BAM file.

    \param[in] alignment BamAlignment record to save
    \sa BamReader::GetNextAlignment(), BamReader::GetNextAlignmentCore()
*/
bool BamWriter::SaveAlignment(const BamAlignment& alignment) {
    return d->SaveAlignment(alignment);
}

/*! \fn void BamWriter::SetCompressionMode(const BamWriter::CompressionMode& compressionMode)
    \brief Sets the output compression mode.

    Default mode is BamWriter::Compressed.

    \note Changing the compression mode is disabled on open files (i.e. the request will
    be ignored). Be sure to call this function before opening the BAM file.

    \code
        BamWriter writer;
        writer.SetCompressionMode(BamWriter::Uncompressed);
        writer.Open( ... );
        // ...
    \endcode

    \param[in] compressionMode desired output compression behavior
    \sa IsOpen(), Open()
*/
void BamWriter::SetCompressionMode(const BamWriter::CompressionMode& compressionMode) {
    d->SetWriteCompressed( compressionMode == BamWriter::Compressed );
}

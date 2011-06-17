// ***************************************************************************
// BamWriter.h (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 4 March 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_H
#define BAMWRITER_H

#include <api/api_global.h>
#include <api/BamAux.h>
#include <string>

namespace BamTools {

class BamAlignment;
class SamHeader;

namespace Internal {
    class BamWriterPrivate;
} // namespace Internal

class API_EXPORT BamWriter {

    public: enum CompressionMode { Compressed = 0
                                 , Uncompressed
                                 };

    // ctor & dtor
    public:
        BamWriter(void);
        ~BamWriter(void);

    // public interface
    public:
        //  closes the current BAM file
        void Close(void);
        // returns true if BAM file is open for writing
        bool IsOpen(void) const;
        // opens a BAM file for writing
        bool Open(const std::string& filename, 
                  const std::string& samHeaderText,
                  const RefVector& referenceSequences);
        // opens a BAM file for writing
        bool Open(const std::string& filename,
                  const SamHeader& samHeader,
                  const RefVector& referenceSequences);
        // saves the alignment to the alignment archive
        void SaveAlignment(const BamAlignment& alignment);
        // sets the output compression mode
        void SetCompressionMode(const CompressionMode& compressionMode);

    // private implementation
    private:
        Internal::BamWriterPrivate* d;
};

} // namespace BamTools

#endif // BAMWRITER_H

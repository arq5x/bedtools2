// ***************************************************************************
// BamWriter_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 24 February 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_P_H
#define BAMWRITER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include <api/BamAux.h>
#include <api/internal/BgzfStream_p.h>
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

class BamWriterPrivate {

    // ctor & dtor
    public:
        BamWriterPrivate(void);
        ~BamWriterPrivate(void);

    // interface methods
    public:
        void Close(void);
        bool IsOpen(void) const;
        bool Open(const std::string& filename,
                  const std::string& samHeaderText,
                  const BamTools::RefVector& referenceSequences);
        void SaveAlignment(const BamAlignment& al);
        void SetWriteCompressed(bool ok);

    // 'internal' methods
    public:
        unsigned int CalculateMinimumBin(const int begin, int end) const;
        void CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string& packedCigar);
        void EncodeQuerySequence(const std::string& query, std::string& encodedQuery);
        void WriteMagicNumber(void);
        void WriteReferences(const BamTools::RefVector& referenceSequences);
        void WriteSamHeaderText(const std::string& samHeaderText);

    // data members
    private:
        BgzfStream m_stream;
        bool m_isBigEndian;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMWRITER_P_H

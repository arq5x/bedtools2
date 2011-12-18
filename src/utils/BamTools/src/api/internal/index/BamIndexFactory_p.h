// ***************************************************************************
// BamIndexFactory_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides interface for generating BamIndex implementations
// ***************************************************************************

#ifndef BAMINDEX_FACTORY_P_H
#define BAMINDEX_FACTORY_P_H

#include "api/BamIndex.h"
#include <string>

namespace BamTools {
namespace Internal {

class BamIndexFactory {

    // static interface methods
    public:
        // creates a new BamIndex object, depending on extension of @indexFilename
        static BamIndex* CreateIndexFromFilename(const std::string& indexFilename,
                                                 BamReaderPrivate* reader);
        // creates a new BamIndex object, of requested @type
        static BamIndex* CreateIndexOfType(const BamIndex::IndexType& type,
                                           BamReaderPrivate* reader);
        // returns name of existing index file that corresponds to @bamFilename
        // will defer to @preferredType if possible
        // if @preferredType not found, will attempt to load any supported index type
        // returns empty string if no index file (of any type) is found
        static const std::string FindIndexFilename(const std::string& bamFilename,
                                                   const BamIndex::IndexType& preferredType);

    // internal methods
    public:
        // generates index filename from BAM filename (depending on requested type)
        // if type is unknown, returns empty string
        static const std::string CreateIndexFilename(const std::string& bamFilename,
                                                     const BamIndex::IndexType& type);
        // retrieves file extension (including '.')
        static const std::string FileExtension(const std::string& filename);
};

} // namespace Internal
} // namespace BamTools

#endif // BAMINDEX_FACTORY_P_H

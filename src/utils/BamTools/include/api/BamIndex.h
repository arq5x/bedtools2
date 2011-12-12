// ***************************************************************************
// BamIndex.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides basic BAM index interface
// ***************************************************************************

#ifndef BAM_INDEX_H
#define BAM_INDEX_H

#include <api/api_global.h>
#include <api/BamAux.h>
#include <string>

namespace BamTools {

namespace Internal {
    class BamReaderPrivate;
} // namespace Internal

/*! \class BamTools::BamIndex
    \brief Provides methods for generating & loading BAM index files.

    This class straddles the line between public API and internal
    implementation detail. Most client code should never have to use this
    class directly.

    It is exposed to the public API to allow advanced users to implement
    their own custom indexing schemes.

    More documentation on methods & enums coming soon.
*/

class API_EXPORT BamIndex {

    // enums
    public:
        // specify index-caching behavior
        enum IndexCacheMode { FullIndexCaching = 0 // store entire index file contents in memory
                            , LimitedIndexCaching  // store only index data for current reference
                            , NoIndexCaching       // do not store any index data between jumps
                            };

        // list of supported BamIndex types
        enum IndexType { BAMTOOLS = 0
                       , STANDARD
                       };
  
    // ctor & dtor
    public:
        BamIndex(Internal::BamReaderPrivate* reader) : m_reader(reader) { }
        virtual ~BamIndex(void) { }
        
    // index interface
    public:
        // builds index from associated BAM file & writes out to index file
        virtual bool Create(void) =0; // creates index file from BAM file
        // returns whether reference has alignments or no
        virtual bool HasAlignments(const int& referenceID) const =0;
        // attempts to use index data to jump to @region, returns success/fail
        // a "successful" jump indicates no error, but not whether this region has data
        //   * thus, the method sets a flag to indicate whether there are alignments
        //     available after the jump position
        virtual bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion) =0;
        // loads existing data from file into memory
        virtual bool Load(const std::string& filename) =0;
        // change the index caching behavior
        virtual void SetCacheMode(const BamIndex::IndexCacheMode& mode) =0;

    // data members
    protected:
        Internal::BamReaderPrivate* m_reader; // copy, not ownedprivate:
};

} // namespace BamTools

#endif // BAM_INDEX_H

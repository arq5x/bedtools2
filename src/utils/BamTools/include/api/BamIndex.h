// ***************************************************************************
// BamIndex.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides basic BAM index interface
// ***************************************************************************

#ifndef BAM_INDEX_H
#define BAM_INDEX_H

#include "api/api_global.h"
#include "api/BamAux.h"
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
*/

class API_EXPORT BamIndex {

    // enums
    public:

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
        virtual bool Create(void) =0;

        // returns a human-readable description of the last error encountered
        std::string GetErrorString(void) { return m_errorString; }

        // returns whether reference has alignments or no
        virtual bool HasAlignments(const int& referenceID) const =0;

        // attempts to use index data to jump to @region, returns success/fail
        // a "successful" jump indicates no error, but not whether this region has data
        //   * thus, the method sets a flag to indicate whether there are alignments
        //     available after the jump position
        virtual bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion) =0;

        // loads existing data from file into memory
        virtual bool Load(const std::string& filename) =0;

        // returns the 'type' enum for derived index format
        virtual BamIndex::IndexType Type(void) const =0;

    //! \cond

    // internal methods
    protected:
        void SetErrorString(const std::string& where, const std::string& what) const {
            m_errorString = where + ": " + what;
        }

    // data members
    protected:
        Internal::BamReaderPrivate* m_reader; // copy, not owned
        mutable std::string m_errorString;

    //! \endcond
};

} // namespace BamTools

#endif // BAM_INDEX_H

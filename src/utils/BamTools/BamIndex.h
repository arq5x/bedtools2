// ***************************************************************************
// BamIndex.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 17 August 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index functionality - both for the default (standardized) BAM 
// index format (.bai) as well as a BamTools-specific (nonstandard) index 
// format (.bti).
// ***************************************************************************

#ifndef BAM_INDEX_H
#define BAM_INDEX_H

#include <string>
#include <vector>
#include "BamAux.h"

namespace BamTools {

class BamReader;
class BgzfData;
  
// --------------------------------------------------  
// BamIndex base class
class BamIndex {

    public:
        BamIndex(BamTools::BgzfData*  bgzf, 
                 BamTools::BamReader* reader,
                 bool isBigEndian);
        virtual ~BamIndex(void) { }

    public:
        // creates index data (in-memory) from current reader data
        virtual bool Build(void) =0;
        // calculates offset(s) for a given region
        virtual bool GetOffsets(const BamTools::BamRegion& region, const bool isRightBoundSpecified, std::vector<int64_t>& offsets) =0;
        // loads existing data from file into memory
        virtual bool Load(const std::string& filename)  =0;
        // returns whether reference has alignments or no
        virtual bool HasAlignments(const int& referenceID); 
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        virtual bool Write(const std::string& bamFilename) =0;
        
    protected:
        BamTools::BgzfData*  m_BGZF;
        BamTools::BamReader* m_reader;
        BamTools::RefVector  m_references;
        bool m_isBigEndian;
};

// --------------------------------------------------
// BamDefaultIndex class
// 
// implements default (per SAM/BAM spec) index file ops
class BamDefaultIndex : public BamIndex {

  
    // ctor & dtor
    public:
        BamDefaultIndex(BamTools::BgzfData*  bgzf, 
                        BamTools::BamReader* reader,
                        bool isBigEndian);
        ~BamDefaultIndex(void);
        
    // interface (implements BamIndex virtual methods)
    public:
        // creates index data (in-memory) from current reader data
        bool Build(void);
        // calculates offset(s) for a given region
        bool GetOffsets(const BamTools::BamRegion& region, const bool isRightBoundSpecified, std::vector<int64_t>& offsets);
         // loads existing data from file into memory
        bool Load(const std::string& filename);
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        bool Write(const std::string& bamFilename);
      
    // internal implementation
    private:
        struct BamDefaultIndexPrivate;
        BamDefaultIndexPrivate* d;
};

// --------------------------------------------------
// BamToolsIndex class
//
// implements BamTools-specific index file ops
class BamToolsIndex : public BamIndex {

    // ctor & dtor
    public:
        BamToolsIndex(BamTools::BgzfData*  bgzf, 
                      BamTools::BamReader* reader,
                      bool isBigEndian);
        ~BamToolsIndex(void);
        
    // interface (implements BamIndex virtual methods)
    public:
        // creates index data (in-memory) from current reader data
        bool Build(void);
        // calculates offset(s) for a given region
        bool GetOffsets(const BamTools::BamRegion& region, const bool isRightBoundSpecified, std::vector<int64_t>& offsets);
         // loads existing data from file into memory
        bool Load(const std::string& filename);
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        bool Write(const std::string& bamFilename);
    
    // internal implementation
    private:
        struct BamToolsIndexPrivate;
        BamToolsIndexPrivate* d;
};

} // namespace BamTools

#endif // BAM_INDEX_H
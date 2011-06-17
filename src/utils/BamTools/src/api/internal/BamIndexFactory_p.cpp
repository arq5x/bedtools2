// ***************************************************************************
// BamIndexFactory_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides interface for generating BamIndex implementations
// ***************************************************************************

#include <api/BamAux.h>
#include <api/internal/BamIndexFactory_p.h>
#include <api/internal/BamStandardIndex_p.h>
#include <api/internal/BamToolsIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
using namespace std;

// generates index filename from BAM filename (depending on requested type)
// if type is unknown, returns empty string
const string BamIndexFactory::CreateIndexFilename(const string& bamFilename,
                                                  const BamIndex::IndexType& type)
{
    switch ( type ) {
        case ( BamIndex::STANDARD ) : return ( bamFilename + BamStandardIndex::Extension() );
        case ( BamIndex::BAMTOOLS ) : return ( bamFilename + BamToolsIndex::Extension() );
        default :
            cerr << "BamIndexFactory ERROR: unknown index type" << type << endl;
            return string();
    }
}

// creates a new BamIndex object, depending on extension of @indexFilename
BamIndex* BamIndexFactory::CreateIndexFromFilename(const string& indexFilename, BamReaderPrivate* reader) {

    // if file doesn't exist, return null index
    if ( !BamTools::FileExists(indexFilename) )
        return 0;

    // get file extension from index filename, including dot (".EXT")
    // if can't get file extension, return null index
    const string extension = FileExtension(indexFilename);
    if ( extension.empty() )
        return 0;

    // create index based on extension
    if      ( extension == BamStandardIndex::Extension() ) return new BamStandardIndex(reader);
    else if ( extension == BamToolsIndex::Extension()    ) return new BamToolsIndex(reader);
    else
        return 0;
}

// creates a new BamIndex, object of requested @type
BamIndex* BamIndexFactory::CreateIndexOfType(const BamIndex::IndexType& type,
                                             BamReaderPrivate* reader)
{
    switch ( type ) {
        case ( BamIndex::STANDARD ) : return new BamStandardIndex(reader);
        case ( BamIndex::BAMTOOLS ) : return new BamToolsIndex(reader);
        default :
            cerr << "BamIndexFactory ERROR: unknown index type " << type << endl;
            return 0;
    }
}

// retrieves file extension (including '.')
const string BamIndexFactory::FileExtension(const string& filename) {

    // if filename cannot contain valid path + extension, return empty string
    if ( filename.empty() || filename.length() <= 4 )
        return string();

    // look for last dot in filename
    size_t lastDotPosition = filename.find_last_of('.');

    // if none found, return empty string
    if ( lastDotPosition == string::npos )
        return string();

    // return substring from last dot position
    return filename.substr(lastDotPosition);
}

// returns name of existing index file that corresponds to @bamFilename
// will defer to @preferredType if possible, if not will attempt to load any supported type
// returns empty string if not found
const string BamIndexFactory::FindIndexFilename(const string& bamFilename,
                                                const BamIndex::IndexType& preferredType)
{
    // try to find index of preferred type first
    // return index filename if found
    string indexFilename = CreateIndexFilename(bamFilename, preferredType);
    if ( !indexFilename.empty() && BamTools::FileExists(indexFilename) )
        return indexFilename;

    // couldn't find preferred type, try the other supported types
    // return index filename if found
    if ( preferredType != BamIndex::STANDARD ) {
        indexFilename = CreateIndexFilename(bamFilename, BamIndex::STANDARD);
        if ( !indexFilename.empty() && BamTools::FileExists(indexFilename) )
            return indexFilename;
    }
    if ( preferredType != BamIndex::BAMTOOLS ) {
        indexFilename = CreateIndexFilename(bamFilename, BamIndex::BAMTOOLS);
        if ( !indexFilename.empty() && BamTools::FileExists(indexFilename) )
            return indexFilename;
    }

    // otherwise couldn't find any index matching this filename
    return string();
}

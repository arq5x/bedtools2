// ***************************************************************************
// SamProgram.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM header program records.
// ***************************************************************************

#ifndef SAM_PROGRAM_H
#define SAM_PROGRAM_H

#include "api/api_global.h"
#include <string>

namespace BamTools {

class SamProgramChain;

struct API_EXPORT SamProgram {

    // ctor & dtor
    SamProgram(void);
    SamProgram(const std::string& id);
    SamProgram(const SamProgram& other);
    ~SamProgram(void);

    // query/modify entire program record
    void Clear(void);                      // clears all data fields

    // convenience query methods
    bool HasCommandLine(void) const;       // returns true if program record has a command line entry
    bool HasID(void) const;                // returns true if program record has an ID
    bool HasName(void) const;              // returns true if program record has a name
    bool HasPreviousProgramID(void) const; // returns true if program record has a 'previous program ID'
    bool HasVersion(void) const;           // returns true if program record has a version

    // data members
    std::string CommandLine;               // CL:<CommandLine>
    std::string ID;                        // ID:<ID>          *Required for valid SAM header*
    std::string Name;                      // PN:<Name>
    std::string PreviousProgramID;         // PP:<PreviousProgramID>
    std::string Version;                   // VN:<Version>

    // internal (non-standard) methods & fields
    private:
        bool HasNextProgramID(void) const;
        std::string NextProgramID;
        friend class BamTools::SamProgramChain;
};

/*! \fn bool operator==(const SamProgram& lhs, const SamProgram& rhs)
    \brief tests equality by comparing program IDs
*/
API_EXPORT inline bool operator==(const SamProgram& lhs, const SamProgram& rhs) {
    return lhs.ID == rhs.ID;
}

} // namespace BamTools

#endif // SAM_PROGRAM_H

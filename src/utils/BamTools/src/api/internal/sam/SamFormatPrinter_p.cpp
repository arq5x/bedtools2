// ***************************************************************************
// SamFormatPrinter.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for printing formatted SAM header to string
// ***************************************************************************

#include "api/SamConstants.h"
#include "api/SamHeader.h"
#include "api/internal/sam/SamFormatPrinter_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

static inline
const string FormatTag(const string& tag, const string& value) {
    return string(Constants::SAM_TAB + tag + Constants::SAM_COLON + value);
}

// ---------------------------------
// SamFormatPrinter implementation
// ---------------------------------

SamFormatPrinter::SamFormatPrinter(const SamHeader& header)
    : m_header(header)
{ }

SamFormatPrinter::~SamFormatPrinter(void) { }

const string SamFormatPrinter::ToString(void) const {

    // clear out stream
    stringstream out("");

    // generate formatted header text
    PrintHD(out);
    PrintSQ(out);
    PrintRG(out);
    PrintPG(out);
    PrintCO(out);

    // return result
    return out.str();
}

void SamFormatPrinter::PrintHD(std::stringstream& out) const {

    // if header has @HD data
    if ( m_header.HasVersion() ) {

        // @HD VN:<Version>
        out << Constants::SAM_HD_BEGIN_TOKEN
            << FormatTag(Constants::SAM_HD_VERSION_TAG, m_header.Version);

        // SO:<SortOrder>
        if ( m_header.HasSortOrder() )
            out << FormatTag(Constants::SAM_HD_SORTORDER_TAG, m_header.SortOrder);

        // GO:<GroupOrder>
        if ( m_header.HasGroupOrder() )
            out << FormatTag(Constants::SAM_HD_GROUPORDER_TAG, m_header.GroupOrder);

        // newline
        out << endl;
    }
}

void SamFormatPrinter::PrintSQ(std::stringstream& out) const {

    // iterate over sequence entries
    SamSequenceConstIterator seqIter = m_header.Sequences.ConstBegin();
    SamSequenceConstIterator seqEnd  = m_header.Sequences.ConstEnd();
    for ( ; seqIter != seqEnd; ++seqIter ) {
        const SamSequence& seq = (*seqIter);

        // @SQ SN:<Name> LN:<Length>
        out << Constants::SAM_SQ_BEGIN_TOKEN
            << FormatTag(Constants::SAM_SQ_NAME_TAG, seq.Name)
            << FormatTag(Constants::SAM_SQ_LENGTH_TAG, seq.Length);

        // AS:<AssemblyID>
        if ( seq.HasAssemblyID() )
            out << FormatTag(Constants::SAM_SQ_ASSEMBLYID_TAG, seq.AssemblyID);

        // M5:<Checksum>
        if ( seq.HasChecksum() )
            out << FormatTag(Constants::SAM_SQ_CHECKSUM_TAG, seq.Checksum);

        // SP:<Species>
        if ( seq.HasSpecies() )
            out << FormatTag(Constants::SAM_SQ_SPECIES_TAG, seq.Species);

        // UR:<URI>
        if ( seq.HasURI() )
            out << FormatTag(Constants::SAM_SQ_URI_TAG, seq.URI);

        // newline
        out << endl;
    }
}

void SamFormatPrinter::PrintRG(std::stringstream& out) const {

    // iterate over read group entries
    SamReadGroupConstIterator rgIter = m_header.ReadGroups.ConstBegin();
    SamReadGroupConstIterator rgEnd  = m_header.ReadGroups.ConstEnd();
    for ( ; rgIter != rgEnd; ++rgIter ) {
        const SamReadGroup& rg = (*rgIter);

        // @RG ID:<ID>
        out << Constants::SAM_RG_BEGIN_TOKEN
            << FormatTag(Constants::SAM_RG_ID_TAG, rg.ID);

        // CN:<SequencingCenter>
        if ( rg.HasSequencingCenter() )
            out << FormatTag(Constants::SAM_RG_SEQCENTER_TAG, rg.SequencingCenter);

        // DS:<Description>
        if ( rg.HasDescription() )
            out << FormatTag(Constants::SAM_RG_DESCRIPTION_TAG, rg.Description);

        // DT:<ProductionDate>
        if ( rg.HasProductionDate() )
            out << FormatTag(Constants::SAM_RG_PRODUCTIONDATE_TAG, rg.ProductionDate);

        // FO:<FlowOrder>
        if ( rg.HasFlowOrder() )
            out << FormatTag(Constants::SAM_RG_FLOWORDER_TAG, rg.FlowOrder);

        // KS:<KeySequence>
        if ( rg.HasKeySequence() )
            out << FormatTag(Constants::SAM_RG_KEYSEQUENCE_TAG, rg.KeySequence);

        // LB:<Library>
        if ( rg.HasLibrary() )
            out << FormatTag(Constants::SAM_RG_LIBRARY_TAG, rg.Library);

        // PG:<Program>
        if ( rg.HasProgram() )
            out << FormatTag(Constants::SAM_RG_PROGRAM_TAG, rg.Program);

        // PI:<PredictedInsertSize>
        if ( rg.HasPredictedInsertSize() )
            out << FormatTag(Constants::SAM_RG_PREDICTEDINSERTSIZE_TAG, rg.PredictedInsertSize);

        // PL:<SequencingTechnology>
        if ( rg.HasSequencingTechnology() )
            out << FormatTag(Constants::SAM_RG_SEQTECHNOLOGY_TAG, rg.SequencingTechnology);

        // PU:<PlatformUnit>
        if ( rg.HasPlatformUnit() )
            out << FormatTag(Constants::SAM_RG_PLATFORMUNIT_TAG, rg.PlatformUnit);

        // SM:<Sample>
        if ( rg.HasSample() )
            out << FormatTag(Constants::SAM_RG_SAMPLE_TAG, rg.Sample);

        // newline
        out << endl;
    }
}

void SamFormatPrinter::PrintPG(std::stringstream& out) const {

    // iterate over program record entries
    SamProgramConstIterator pgIter = m_header.Programs.ConstBegin();
    SamProgramConstIterator pgEnd  = m_header.Programs.ConstEnd();
    for ( ; pgIter != pgEnd; ++pgIter ) {
        const SamProgram& pg = (*pgIter);

        // @PG ID:<ID>
        out << Constants::SAM_PG_BEGIN_TOKEN
            << FormatTag(Constants::SAM_PG_ID_TAG, pg.ID);

        // PN:<Name>
        if ( pg.HasName() )
            out << FormatTag(Constants::SAM_PG_NAME_TAG, pg.Name);

        // CL:<CommandLine>
        if ( pg.HasCommandLine() )
            out << FormatTag(Constants::SAM_PG_COMMANDLINE_TAG, pg.CommandLine);

        // PP:<PreviousProgramID>
        if ( pg.HasPreviousProgramID() )
            out << FormatTag(Constants::SAM_PG_PREVIOUSPROGRAM_TAG, pg.PreviousProgramID);

        // VN:<Version>
        if ( pg.HasVersion() )
            out << FormatTag(Constants::SAM_PG_VERSION_TAG, pg.Version);

        // newline
        out << endl;
    }
}

void SamFormatPrinter::PrintCO(std::stringstream& out) const {

    // iterate over comments
    vector<string>::const_iterator commentIter = m_header.Comments.begin();
    vector<string>::const_iterator commentEnd  = m_header.Comments.end();
    for ( ; commentIter != commentEnd; ++commentIter ) {

        // @CO <Comment>
        out << Constants::SAM_CO_BEGIN_TOKEN
            << Constants::SAM_TAB
            << (*commentIter)
            << endl;
    }
}

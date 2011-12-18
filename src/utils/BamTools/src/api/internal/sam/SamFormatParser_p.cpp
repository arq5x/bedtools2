// ***************************************************************************
// SamFormatParser.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for parsing SAM header text into SamHeader object
// ***************************************************************************

#include "api/SamConstants.h"
#include "api/SamHeader.h"
#include "api/internal/sam/SamFormatParser_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

SamFormatParser::SamFormatParser(SamHeader& header)
    : m_header(header)
{ }

SamFormatParser::~SamFormatParser(void) { }

void SamFormatParser::Parse(const string& headerText) {

    // clear header's prior contents
    m_header.Clear();

    // empty header is OK, but skip processing
    if ( headerText.empty() )
        return;

    // other wise parse SAM lines
    istringstream headerStream(headerText);
    string headerLine("");
    while ( getline(headerStream, headerLine) )
         ParseSamLine(headerLine);
}

void SamFormatParser::ParseSamLine(const string& line) {

    // skip if line is not long enough to contain true values
    if ( line.length() < 5 ) return;

    // determine token at beginning of line
    const string firstToken = line.substr(0,3);
    const string restOfLine = line.substr(4);
    if      ( firstToken == Constants::SAM_HD_BEGIN_TOKEN) ParseHDLine(restOfLine);
    else if ( firstToken == Constants::SAM_SQ_BEGIN_TOKEN) ParseSQLine(restOfLine);
    else if ( firstToken == Constants::SAM_RG_BEGIN_TOKEN) ParseRGLine(restOfLine);
    else if ( firstToken == Constants::SAM_PG_BEGIN_TOKEN) ParsePGLine(restOfLine);
    else if ( firstToken == Constants::SAM_CO_BEGIN_TOKEN) ParseCOLine(restOfLine);
}

void SamFormatParser::ParseHDLine(const string& line) {

    // split HD lines into tokens
    vector<string> tokens = Split(line, Constants::SAM_TAB);

    // iterate over tokens
    vector<string>::const_iterator tokenIter = tokens.begin();
    vector<string>::const_iterator tokenEnd  = tokens.end();
    for ( ; tokenIter != tokenEnd; ++tokenIter ) {

        // get tag/value
        const string tokenTag = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set header contents
        if      ( tokenTag == Constants::SAM_HD_VERSION_TAG    ) m_header.Version    = tokenValue;
        else if ( tokenTag == Constants::SAM_HD_SORTORDER_TAG  ) m_header.SortOrder  = tokenValue;
        else if ( tokenTag == Constants::SAM_HD_GROUPORDER_TAG ) m_header.GroupOrder = tokenValue;
    }

    // check for required tags
    if ( !m_header.HasVersion() )
        throw BamException("SamFormatParser::ParseHDLine", "@HD line is missing VN tag");
}

void SamFormatParser::ParseSQLine(const string& line) {

    SamSequence seq;

    // split SQ line into tokens
    vector<string> tokens = Split(line, Constants::SAM_TAB);

    // iterate over tokens
    vector<string>::const_iterator tokenIter = tokens.begin();
    vector<string>::const_iterator tokenEnd  = tokens.end();
    for ( ; tokenIter != tokenEnd; ++tokenIter ) {

        // get tag/value
        const string tokenTag = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set sequence contents
        if      ( tokenTag == Constants::SAM_SQ_NAME_TAG       ) seq.Name = tokenValue;
        else if ( tokenTag == Constants::SAM_SQ_LENGTH_TAG     ) seq.Length = tokenValue;
        else if ( tokenTag == Constants::SAM_SQ_ASSEMBLYID_TAG ) seq.AssemblyID = tokenValue;
        else if ( tokenTag == Constants::SAM_SQ_CHECKSUM_TAG   ) seq.Checksum = tokenValue;
        else if ( tokenTag == Constants::SAM_SQ_SPECIES_TAG    ) seq.Species = tokenValue;
        else if ( tokenTag == Constants::SAM_SQ_URI_TAG        ) seq.URI = tokenValue;
    }

    // check for required tags
    if ( !seq.HasName() )
        throw BamException("SamFormatParser::ParseSQLine", "@SQ line is missing SN tag");
    if ( !seq.HasLength() )
        throw BamException("SamFormatParser::ParseSQLine", "@SQ line is missing LN tag");

    // store SAM sequence entry
    m_header.Sequences.Add(seq);
}

void SamFormatParser::ParseRGLine(const string& line) {

    SamReadGroup rg;

    // split string into tokens
    vector<string> tokens = Split(line, Constants::SAM_TAB);

    // iterate over tokens
    vector<string>::const_iterator tokenIter = tokens.begin();
    vector<string>::const_iterator tokenEnd  = tokens.end();
    for ( ; tokenIter != tokenEnd; ++tokenIter ) {

        // get token tag/value
        const string tokenTag = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set read group contents
        if      ( tokenTag == Constants::SAM_RG_ID_TAG                  ) rg.ID = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_DESCRIPTION_TAG         ) rg.Description = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_FLOWORDER_TAG           ) rg.FlowOrder = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_KEYSEQUENCE_TAG         ) rg.KeySequence = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_LIBRARY_TAG             ) rg.Library = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_PLATFORMUNIT_TAG        ) rg.PlatformUnit = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_PREDICTEDINSERTSIZE_TAG ) rg.PredictedInsertSize = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_PRODUCTIONDATE_TAG      ) rg.ProductionDate = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_PROGRAM_TAG             ) rg.Program = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_SAMPLE_TAG              ) rg.Sample = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_SEQCENTER_TAG           ) rg.SequencingCenter = tokenValue;
        else if ( tokenTag == Constants::SAM_RG_SEQTECHNOLOGY_TAG       ) rg.SequencingTechnology = tokenValue;
    }

    // check for required tags
    if ( !rg.HasID() )
        throw BamException("SamFormatParser::ParseRGLine", "@RG line is missing ID tag");

    // store SAM read group entry
    m_header.ReadGroups.Add(rg);
}

void SamFormatParser::ParsePGLine(const string& line) {

    SamProgram pg;

    // split string into tokens
    vector<string> tokens = Split(line, Constants::SAM_TAB);

    // iterate over tokens
    vector<string>::const_iterator tokenIter = tokens.begin();
    vector<string>::const_iterator tokenEnd  = tokens.end();
    for ( ; tokenIter != tokenEnd; ++tokenIter ) {

        // get token tag/value
        const string tokenTag = (*tokenIter).substr(0,2);
        const string tokenValue = (*tokenIter).substr(3);

        // set program record contents
        if      ( tokenTag == Constants::SAM_PG_ID_TAG              ) pg.ID = tokenValue;
        else if ( tokenTag == Constants::SAM_PG_NAME_TAG            ) pg.Name = tokenValue;
        else if ( tokenTag == Constants::SAM_PG_COMMANDLINE_TAG     ) pg.CommandLine = tokenValue;
        else if ( tokenTag == Constants::SAM_PG_PREVIOUSPROGRAM_TAG ) pg.PreviousProgramID = tokenValue;
        else if ( tokenTag == Constants::SAM_PG_VERSION_TAG         ) pg.Version = tokenValue;
    }

    // check for required tags
    if ( !pg.HasID() )
        throw BamException("SamFormatParser::ParsePGLine", "@PG line is missing ID tag");

    // store SAM program entry
    m_header.Programs.Add(pg);
}

void SamFormatParser::ParseCOLine(const string& line) {
    // simply add line to comments list
    m_header.Comments.push_back(line);
}

const vector<string> SamFormatParser::Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
        tokens.push_back(token);
    return tokens;
}

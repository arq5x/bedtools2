// ***************************************************************************
// SamSequenceDictionary.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 16 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides methods for operating on a collection of SamSequence entries.
// *************************************************************************

#include "api/SamSequenceDictionary.h"
using namespace BamTools;

#include <iostream>
using namespace std;

/*! \class BamTools::SamSequenceDictionary
    \brief Container of SamSequence entries.

    Provides methods for operating on a collection of SamSequence entries.
*/

/*! \fn SamSequenceDictionary::SamSequenceDictionary(void)
    \brief constructor
*/
SamSequenceDictionary::SamSequenceDictionary(void) { }

/*! \fn SamSequenceDictionary::SamSequenceDictionary(const SamSequenceDictionary& other)
    \brief copy constructor
*/
SamSequenceDictionary::SamSequenceDictionary(const SamSequenceDictionary& other)
    : m_data(other.m_data)
    , m_lookupData(other.m_lookupData)
{ }

/*! \fn SamSequenceDictionary::~SamSequenceDictionary(void)
    \brief destructor
*/
SamSequenceDictionary::~SamSequenceDictionary(void) { }

/*! \fn void SamSequenceDictionary::Add(const SamSequence& sequence)
    \brief Appends a sequence to the dictionary.

    Duplicate entries are silently discarded.

    \param[in] sequence entry to be added
*/
void SamSequenceDictionary::Add(const SamSequence& sequence) {
    if ( IsEmpty() || !Contains(sequence) ) {
        m_data.push_back(sequence);
        m_lookupData[sequence.Name] = m_data.size() - 1;
    }
}

/*! \fn void SamSequenceDictionary::Add(const std::string& name, const int& length)
    \brief Appends a sequence to the dictionary.

    This is an overloaded function.

    \param[in] name name of sequence entry to be added
    \param[in] length length of sequence entry to be added
    \sa Add()
*/
void SamSequenceDictionary::Add(const std::string& name, const int& length) {
    Add( SamSequence(name, length) );
}

/*! \fn void SamSequenceDictionary::Add(const SamSequenceDictionary& sequences)
    \brief Appends another sequence dictionary to this one

    This is an overloaded function.

    \param[in] sequences sequence dictionary to be appended
    \sa Add()
*/
void SamSequenceDictionary::Add(const SamSequenceDictionary& sequences) {
    SamSequenceConstIterator seqIter = sequences.ConstBegin();
    SamSequenceConstIterator seqEnd  = sequences.ConstEnd();
    for ( ; seqIter != seqEnd; ++seqIter )
        Add(*seqIter);
}

/*! \fn void SamSequenceDictionary::Add(const std::vector<SamSequence>& sequences)
    \brief Appends multiple sequences to the dictionary.

    This is an overloaded function.

    \param[in] sequences entries to be added
    \sa Add()
*/
void SamSequenceDictionary::Add(const std::vector<SamSequence>& sequences) {
    vector<SamSequence>::const_iterator seqIter = sequences.begin();
    vector<SamSequence>::const_iterator seqEnd  = sequences.end();
    for ( ; seqIter!= seqEnd; ++seqIter )
        Add(*seqIter);
}

/*! \fn void SamSequenceDictionary::Add(const std::map<std::string, int>& sequenceMap)
    \brief Appends multiple sequences to the dictionary.

    This is an overloaded function.

    \param[in] sequenceMap map of sequence entries (name => length) to be added
    \sa Add()
*/
void SamSequenceDictionary::Add(const std::map<std::string, int>& sequenceMap) {
    map<string, int>::const_iterator seqIter = sequenceMap.begin();
    map<string, int>::const_iterator seqEnd  = sequenceMap.end();
    for ( ; seqIter != seqEnd; ++seqIter ) {
        const string& name = (*seqIter).first;
        const int& length = (*seqIter).second;
        Add( SamSequence(name, length) );
    }
}

/*! \fn SamSequenceIterator SamSequenceDictionary::Begin(void)
    \return an STL iterator pointing to the first sequence
    \sa ConstBegin(), End()
*/
SamSequenceIterator SamSequenceDictionary::Begin(void) {
    return m_data.begin();
}

/*! \fn SamSequenceConstIterator SamSequenceDictionary::Begin(void) const
    \return an STL const_iterator pointing to the first sequence

    This is an overloaded function.

    \sa ConstBegin(), End()
*/
SamSequenceConstIterator SamSequenceDictionary::Begin(void) const {
    return m_data.begin();
}

/*! \fn void SamSequenceDictionary::Clear(void)
    \brief Clears all sequence entries.
*/
void SamSequenceDictionary::Clear(void) {
    m_data.clear();
    m_lookupData.clear();
}

/*! \fn SamSequenceConstIterator SamSequenceDictionary::ConstBegin(void) const
    \return an STL const_iterator pointing to the first sequence
    \sa Begin(), ConstEnd()
*/
SamSequenceConstIterator SamSequenceDictionary::ConstBegin(void) const {
    return m_data.begin();
}

/*! \fn SamSequenceConstIterator SamSequenceDictionary::ConstEnd(void) const
    \return an STL const_iterator pointing to the imaginary entry after the last sequence
    \sa End(), ConstBegin()
*/
SamSequenceConstIterator SamSequenceDictionary::ConstEnd(void) const {
    return m_data.end();
}

/*! \fn bool SamSequenceDictionary::Contains(const std::string& sequenceName) const
    \brief Returns true if dictionary contains sequence.

    \param[in] sequenceName search for sequence matching this name
    \return \c true if dictionary contains a sequence with this name
*/
bool SamSequenceDictionary::Contains(const std::string& sequenceName) const {
    return ( m_lookupData.find(sequenceName) != m_lookupData.end() );
}

/*! \fn bool SamSequenceDictionary::Contains(const SamSequence& sequence) const
    \brief Returns true if dictionary contains sequence (matches on name).

    This is an overloaded function.

    \param[in] sequence search for this sequence
    \return \c true if dictionary contains sequence (matching on name)
*/
bool SamSequenceDictionary::Contains(const SamSequence& sequence) const {
    return Contains(sequence.Name);
}

/*! \fn SamSequenceIterator SamSequenceDictionary::End(void)
    \return an STL iterator pointing to the imaginary entry after the last sequence
    \sa Begin(), ConstEnd()
*/
SamSequenceIterator SamSequenceDictionary::End(void) {
    return m_data.end();
}

/*! \fn SamSequenceConstIterator SamSequenceDictionary::End(void) const
    \return an STL const_iterator pointing to the imaginary entry after the last sequence

    This is an overloaded function.

    \sa Begin(), ConstEnd()
*/
SamSequenceConstIterator SamSequenceDictionary::End(void) const {
    return m_data.end();
}

/*! \fn bool SamSequenceDictionary::IsEmpty(void) const
    \brief Returns \c true if dictionary contains no sequences
    \sa Size()
*/
bool SamSequenceDictionary::IsEmpty(void) const {
    return m_data.empty();
}

/*! \fn void SamSequenceDictionary::Remove(const SamSequence& sequence)
    \brief Removes sequence from dictionary, if found (matches on name).

    This is an overloaded function.

    \param[in] sequence SamSequence to remove (matching on name)
*/
void SamSequenceDictionary::Remove(const SamSequence& sequence) {
    Remove(sequence.Name);
}

/*! \fn void SamSequenceDictionary::Remove(const std::string& sequenceName)
    \brief Removes sequence from dictionary, if found.

    \param[in] sequenceName name of sequence to remove
    \sa Remove()
*/
void SamSequenceDictionary::Remove(const std::string& sequenceName) {

    // skip if empty dictionary or if name unknown
    if ( IsEmpty() || !Contains(sequenceName) )
        return;

    // update 'lookup index' for every entry after @sequenceName
    const size_t indexToRemove = m_lookupData[sequenceName];
    const size_t numEntries = m_data.size();
    for ( size_t i = indexToRemove+1; i < numEntries; ++i ) {
        const SamSequence& sq = m_data.at(i);
        --m_lookupData[sq.Name];
    }

    // erase entry from containers
    m_data.erase( Begin() + indexToRemove );
    m_lookupData.erase(sequenceName);
}

/*! \fn void SamSequenceDictionary::Remove(const std::vector<SamSequence>& sequences)
    \brief Removes multiple sequences from dictionary.

    This is an overloaded function.

    \param[in] sequences sequences to remove
    \sa Remove()
*/
void SamSequenceDictionary::Remove(const std::vector<SamSequence>& sequences) {
    vector<SamSequence>::const_iterator rgIter = sequences.begin();
    vector<SamSequence>::const_iterator rgEnd  = sequences.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Remove(*rgIter);
}

/*! \fn void SamSequenceDictionary::Remove(const std::vector<std::string>& sequenceNames)
    \brief Removes multiple sequences from dictionary.

    This is an overloaded function.

    \param[in] sequenceNames names of the sequences to remove
    \sa Remove()
*/
void SamSequenceDictionary::Remove(const std::vector<std::string>& sequenceNames) {
    vector<string>::const_iterator rgIter = sequenceNames.begin();
    vector<string>::const_iterator rgEnd  = sequenceNames.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Remove(*rgIter);
}

/*! \fn int SamSequenceDictionary::Size(void) const
    \brief Returns number of sequences in dictionary.
    \sa IsEmpty()
*/
int SamSequenceDictionary::Size(void) const {
    return m_data.size();
}

/*! \fn SamSequence& SamSequenceDictionary::operator[](const std::string& sequenceName)
    \brief Retrieves the modifiable SamSequence that matches \a sequenceName.

    \note If the dictionary contains no sequence matching this name, this function inserts
    a new one with this name (length:0), and returns a reference to it. If you want to avoid
    this insertion behavior, check the result of Contains() before using this operator.

    \param[in] sequenceName name of sequence to retrieve
    \return a modifiable reference to the SamSequence associated with the name
*/
SamSequence& SamSequenceDictionary::operator[](const std::string& sequenceName) {

    if ( !Contains(sequenceName) ) {
        SamSequence seq(sequenceName, 0);
        m_data.push_back(seq);
        m_lookupData[sequenceName] = m_data.size() - 1;
    }

    const size_t index = m_lookupData[sequenceName];
    return m_data.at(index);
}

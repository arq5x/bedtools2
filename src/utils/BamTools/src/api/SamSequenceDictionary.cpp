// ***************************************************************************
// SamSequenceDictionary.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides methods for operating on a collection of SamSequence entries.
// *************************************************************************

#include <api/SamSequenceDictionary.h>
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
{ }

/*! \fn SamSequenceDictionary::~SamSequenceDictionary(void)
    \brief destructor
*/
SamSequenceDictionary::~SamSequenceDictionary(void) { }

/*! \fn void SamSequenceDictionary::Add(const SamSequence& sequence)
    \brief Adds a sequence to the dictionary.

    Duplicate entries are silently discarded.

    \param sequence entry to be added
*/
void SamSequenceDictionary::Add(const SamSequence& sequence) {

    // TODO: report error on attempted duplicate?

    if ( IsEmpty() || !Contains(sequence) )
        m_data.push_back(sequence);
}

/*! \fn void SamSequenceDictionary::Add(const std::string& name, const int& length)
    \brief Adds a sequence to the dictionary.

    This is an overloaded function.

    \param name name of sequence entry to be added
    \param length length of sequence entry to be added
    \sa Add()
*/
void SamSequenceDictionary::Add(const std::string& name, const int& length) {
    Add( SamSequence(name, length) );
}

/*! \fn void SamSequenceDictionary::Add(const std::vector<SamSequence>& sequences)
    \brief Adds multiple sequences to the dictionary.

    This is an overloaded function.

    \param sequences entries to be added
    \sa Add()
*/
void SamSequenceDictionary::Add(const std::vector<SamSequence>& sequences) {
    vector<SamSequence>::const_iterator seqIter = sequences.begin();
    vector<SamSequence>::const_iterator seqEnd  = sequences.end();
    for ( ; seqIter!= seqEnd; ++seqIter )
        Add(*seqIter);
}

/*! \fn void SamSequenceDictionary::Add(const std::map<std::string, int>& sequenceMap)
    \brief Adds multiple sequences to the dictionary.

    This is an overloaded function.

    \param sequenceMap map of sequence entries (name => length) to be added
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
    \param sequenceName search for sequence matching this name
    \return \c true if dictionary contains a sequence with this name
*/
bool SamSequenceDictionary::Contains(const std::string& sequenceName) const {
    return ( IndexOf(sequenceName) != (int)m_data.size() );
}

/*! \fn bool SamSequenceDictionary::Contains(const SamSequence& sequence) const
    \brief Returns true if dictionary contains sequence (matches on name).

    This is an overloaded function.

    \param sequence search for this sequence
    \return \c true if dictionary contains sequence (matching on name)
*/
bool SamSequenceDictionary::Contains(const SamSequence& sequence) const {
    return ( IndexOf(sequence.Name) != (int)m_data.size() );
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

/*! \fn int SamSequenceDictionary::IndexOf(const std::string& name) const
    \internal
    \return index of sequence if found (matching on name).  Otherwise, returns vector::size() (invalid index).
*/
int SamSequenceDictionary::IndexOf(const std::string& name) const {
    SamSequenceConstIterator begin = ConstBegin();
    SamSequenceConstIterator iter  = begin;
    SamSequenceConstIterator end   = ConstEnd();
    for ( ; iter != end; ++iter ) {
        const SamSequence& currentSeq = (*iter);
        if ( currentSeq.Name == name )
            break;
    }
    return distance( begin, iter );
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

    \param sequence SamSequence to remove (matching on name)
*/
void SamSequenceDictionary::Remove(const SamSequence& sequence) {
    Remove( sequence.Name );
}

/*! \fn void SamSequenceDictionary::Remove(const std::string& sequenceName)
    \brief Removes sequence from dictionary, if found.

    \param sequenceName name of sequence to remove
    \sa Remove()
*/
void SamSequenceDictionary::Remove(const std::string& sequenceName) {
    if ( Contains(sequenceName) )
        m_data.erase( m_data.begin() + IndexOf(sequenceName) );
}

/*! \fn void SamSequenceDictionary::Remove(const std::vector<SamSequence>& sequences)
    \brief Removes multiple sequences from dictionary.

    This is an overloaded function.

    \param sequences sequences to remove
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

    \param sequenceNames names of the sequences to remove
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

    NOTE - If the dictionary contains no sequence matching this name, this function inserts
    a new one with this name (length:0), and returns a reference to it.

    If you want to avoid this insertion behavior, check the result of Contains() before
    using this operator.

    \param sequenceName name of sequence to retrieve
    \return a modifiable reference to the SamSequence associated with the name
*/
SamSequence& SamSequenceDictionary::operator[](const std::string& sequenceName) {

    // look up sequence ID
    int index = IndexOf(sequenceName);

    // if found, return sequence at index
    if ( index != (int)m_data.size() )
        return m_data[index];

    // otherwise, append new sequence and return reference
    else {
        m_data.push_back( SamSequence(sequenceName, 0) );
        return m_data.back();
    }
}

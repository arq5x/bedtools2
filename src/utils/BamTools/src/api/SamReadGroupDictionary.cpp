// ***************************************************************************
// SamReadGroupDictionary.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 16 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides methods for operating on a collection of SamReadGroup entries.
// ***************************************************************************

#include "api/SamReadGroupDictionary.h"
using namespace BamTools;

#include <iostream>
using namespace std;

/*! \class BamTools::SamReadGroupDictionary
    \brief Container of SamReadGroup entries.

    Provides methods for operating on a collection of SamReadGroup entries.
*/

/*! \fn SamReadGroupDictionary::SamReadGroupDictionary(void)
    \brief constructor
*/
SamReadGroupDictionary::SamReadGroupDictionary(void) { }

/*! \fn SamReadGroupDictionary::SamReadGroupDictionary(const SamReadGroupDictionary& other)
    \brief copy constructor
*/
SamReadGroupDictionary::SamReadGroupDictionary(const SamReadGroupDictionary& other)
    : m_data(other.m_data)
    , m_lookupData(other.m_lookupData)
{ }

/*! \fn SamReadGroupDictionary::~SamReadGroupDictionary(void)
    \brief destructor
*/
SamReadGroupDictionary::~SamReadGroupDictionary(void) { }

/*! \fn void SamReadGroupDictionary::Add(const SamReadGroup& readGroup)
    \brief Appends a read group to the dictionary.

    Duplicate entries are silently discarded.

    \param[in] readGroup entry to be added
*/
void SamReadGroupDictionary::Add(const SamReadGroup& readGroup) {
    if ( IsEmpty() || !Contains(readGroup) ) {
        m_data.push_back(readGroup);
        m_lookupData[readGroup.ID] = m_data.size() - 1;
    }
}

/*! \fn void SamReadGroupDictionary::Add(const std::string& readGroupId)
    \brief Appends a read group to the dictionary.

    This is an overloaded function.

    \param[in] readGroupId ID of read group to be added
    \sa Add()
*/
void SamReadGroupDictionary::Add(const std::string& readGroupId) {
    Add( SamReadGroup(readGroupId) );
}

/*! \fn void SamReadGroupDictionary::Add(const SamReadGroupDictionary& readGroups)
    \brief Appends another read group dictionary to this one.

    This is an overloaded function.

    \param[in] readGroups entries to be added
    \sa Add()
*/
void SamReadGroupDictionary::Add(const SamReadGroupDictionary& readGroups) {
    SamReadGroupConstIterator rgIter = readGroups.ConstBegin();
    SamReadGroupConstIterator rgEnd  = readGroups.ConstEnd();
    for ( ; rgIter != rgEnd; ++rgIter )
        Add(*rgIter);
}

/*! \fn void SamReadGroupDictionary::Add(const std::vector<SamReadGroup>& readGroups)
    \brief Appends multiple read groups to the dictionary.

    This is an overloaded function.

    \param[in] readGroups entries to be added
    \sa Add()
*/
void SamReadGroupDictionary::Add(const std::vector<SamReadGroup>& readGroups) {
    vector<SamReadGroup>::const_iterator rgIter = readGroups.begin();
    vector<SamReadGroup>::const_iterator rgEnd  = readGroups.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Add(*rgIter);
}

/*! \fn void SamReadGroupDictionary::Add(const std::vector<std::string>& readGroupIds)
    \brief Appends multiple read groups to the dictionary.

    This is an overloaded function.

    \param[in] readGroupIds IDs of read groups to be added
    \sa Add()
*/
void SamReadGroupDictionary::Add(const std::vector<std::string>& readGroupIds) {
    vector<string>::const_iterator rgIter = readGroupIds.begin();
    vector<string>::const_iterator rgEnd  = readGroupIds.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Add(*rgIter);
}

/*! \fn SamReadGroupIterator SamReadGroupDictionary::Begin(void)
    \return an STL iterator pointing to the first read group
    \sa ConstBegin(), End()
*/
SamReadGroupIterator SamReadGroupDictionary::Begin(void) {
    return m_data.begin();
}

/*! \fn SamReadGroupConstIterator SamReadGroupDictionary::Begin(void) const
    \return an STL const_iterator pointing to the first read group

    This is an overloaded function.

    \sa ConstBegin(), End()
*/
SamReadGroupConstIterator SamReadGroupDictionary::Begin(void) const {
    return m_data.begin();
}

/*! \fn void SamReadGroupDictionary::Clear(void)
    \brief Clears all read group entries.
*/
void SamReadGroupDictionary::Clear(void) {
    m_data.clear();
    m_lookupData.clear();
}

/*! \fn SamReadGroupConstIterator SamReadGroupDictionary::ConstBegin(void) const
    \return an STL const_iterator pointing to the first read group
    \sa Begin(), ConstEnd()
*/
SamReadGroupConstIterator SamReadGroupDictionary::ConstBegin(void) const {
    return m_data.begin();
}

/*! \fn SamReadGroupConstIterator SamReadGroupDictionary::ConstEnd(void) const
    \return an STL const_iterator pointing to the imaginary entry after the last read group
    \sa ConstBegin(), End()
*/
SamReadGroupConstIterator SamReadGroupDictionary::ConstEnd(void) const {
    return m_data.end();
}

/*! \fn bool SamReadGroupDictionary::Contains(const std::string& readGroupId) const
    \brief Returns true if dictionary contains read group.

    \param[in] readGroupId search for read group matching this ID
    \return \c true if dictionary contains a read group with this ID
*/
bool SamReadGroupDictionary::Contains(const std::string& readGroupId) const {
    return ( m_lookupData.find(readGroupId) != m_lookupData.end() );
}

/*! \fn bool SamReadGroupDictionary::Contains(const SamReadGroup& readGroup) const
    \brief Returns true if dictionary contains read group (matching on ID).

    This is an overloaded function.

    \param[in] readGroup search for this read group
    \return \c true if dictionary contains read group (matching on ID).
*/
bool SamReadGroupDictionary::Contains(const SamReadGroup& readGroup) const {
    return Contains(readGroup.ID);
}

/*! \fn SamReadGroupIterator SamReadGroupDictionary::End(void)
    \return an STL iterator pointing to the imaginary entry after the last read group
    \sa Begin(), ConstEnd()
*/
SamReadGroupIterator SamReadGroupDictionary::End(void) {
    return m_data.end();
}

/*! \fn SamReadGroupConstIterator SamReadGroupDictionary::End(void) const
    \return an STL const_iterator pointing to the imaginary entry after the last read group

    This is an overloaded function.

    \sa Begin(), ConstEnd()
*/
SamReadGroupConstIterator SamReadGroupDictionary::End(void) const {
    return m_data.end();
}

/*! \fn bool SamReadGroupDictionary::IsEmpty(void) const
    \brief Returns \c true if dictionary contains no read groups
    \sa Size()
*/
bool SamReadGroupDictionary::IsEmpty(void) const {
    return m_data.empty();
}

/*! \fn void SamReadGroupDictionary::Remove(const SamReadGroup& readGroup)
    \brief Removes read group from dictionary, if found (matching on ID).

    This is an overloaded function.

    \param[in] readGroup read group to remove (matches on ID)
*/
void SamReadGroupDictionary::Remove(const SamReadGroup& readGroup) {
    Remove(readGroup.ID);
}

/*! \fn void SamReadGroupDictionary::Remove(const std::string& readGroupId)
    \brief Removes read group from dictionary, if found.

    \param[in] readGroupId ID of read group to remove
    \sa Remove()
*/
void SamReadGroupDictionary::Remove(const std::string& readGroupId) {

    // skip if empty dictionary or if ID unknown
    if ( IsEmpty() || !Contains(readGroupId) )
        return;

    // update 'lookup index' for every entry after @readGroupId
    const size_t indexToRemove = m_lookupData[readGroupId];
    const size_t numEntries = m_data.size();
    for ( size_t i = indexToRemove+1; i < numEntries; ++i ) {
        const SamReadGroup& rg = m_data.at(i);
        --m_lookupData[rg.ID];
    }

    // erase entry from containers
    m_data.erase( Begin() + indexToRemove );
    m_lookupData.erase(readGroupId);
}

/*! \fn void SamReadGroupDictionary::Remove(const std::vector<SamReadGroup>& readGroups)
    \brief Removes multiple read groups from dictionary (matching on ID).

    This is an overloaded function.

    \param[in] readGroups read groups to remove
    \sa Remove()
*/
void SamReadGroupDictionary::Remove(const std::vector<SamReadGroup>& readGroups) {
    vector<SamReadGroup>::const_iterator rgIter = readGroups.begin();
    vector<SamReadGroup>::const_iterator rgEnd  = readGroups.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Remove(*rgIter);
}

/*! \fn void SamReadGroupDictionary::Remove(const std::vector<std::string>& readGroupIds)
    \brief Removes multiple read groups from dictionary.

    This is an overloaded function.

    \param[in] readGroupIds IDs of the read groups to remove
    \sa Remove()
*/
void SamReadGroupDictionary::Remove(const std::vector<std::string>& readGroupIds) {
    vector<string>::const_iterator rgIter = readGroupIds.begin();
    vector<string>::const_iterator rgEnd  = readGroupIds.end();
    for ( ; rgIter!= rgEnd; ++rgIter )
        Remove(*rgIter);
}

/*! \fn int SamReadGroupDictionary::Size(void) const
    \brief Returns number of read groups in dictionary.
    \sa IsEmpty()
*/
int SamReadGroupDictionary::Size(void) const {
    return m_data.size();
}

/*! \fn SamReadGroup& SamReadGroupDictionary::operator[](const std::string& readGroupId)
    \brief Retrieves the modifiable SamReadGroup that matches \a readGroupId.

    \note If the dictionary contains no read group matching this ID, this function inserts
    a new one with this ID, and returns a reference to it. If you want to avoid this insertion
    behavior, check the result of Contains() before using this operator.

    \param[in] readGroupId ID of read group to retrieve
    \return a modifiable reference to the SamReadGroup associated with the ID
*/
SamReadGroup& SamReadGroupDictionary::operator[](const std::string& readGroupId) {

    if ( !Contains(readGroupId) ) {
        SamReadGroup rg(readGroupId);
        m_data.push_back(rg);
        m_lookupData[readGroupId] = m_data.size() - 1;
    }

    const size_t index = m_lookupData[readGroupId];
    return m_data.at(index);
}

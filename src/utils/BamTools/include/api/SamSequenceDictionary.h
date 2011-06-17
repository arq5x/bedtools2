// ***************************************************************************
// SamSequenceDictionary.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 April 2011
// ---------------------------------------------------------------------------
// Provides methods for operating on a collection of SamSequence entries.
// ***************************************************************************

#ifndef SAM_SEQUENCE_DICTIONARY_H
#define SAM_SEQUENCE_DICTIONARY_H

#include <api/api_global.h>
#include <api/SamSequence.h>
#include <string>
#include <map>
#include <vector>

namespace BamTools {

typedef std::vector<SamSequence>             SamSequenceContainer;
typedef SamSequenceContainer::iterator       SamSequenceIterator;
typedef SamSequenceContainer::const_iterator SamSequenceConstIterator;

class API_EXPORT SamSequenceDictionary {

    // ctor & dtor
    public:
        SamSequenceDictionary(void);
        SamSequenceDictionary(const SamSequenceDictionary& other);
        ~SamSequenceDictionary(void);

    // query/modify sequence data
    public:
        // adds a sequence
        void Add(const SamSequence& sequence);
        void Add(const std::string& name, const int& length);

        // adds multiple sequences
        void Add(const std::vector<SamSequence>& sequences);
        void Add(const std::map<std::string, int>& sequenceMap);

        // clears all sequence entries
        void Clear(void);

        // returns true if dictionary contains this sequence
        bool Contains(const SamSequence& sequence) const;
        bool Contains(const std::string& sequenceName) const;

        // returns true if dictionary is empty
        bool IsEmpty(void) const;

        // removes sequence, if found
        void Remove(const SamSequence& sequence);
        void Remove(const std::string& sequenceName);

        // removes multiple sequences
        void Remove(const std::vector<SamSequence>& sequences);
        void Remove(const std::vector<std::string>& sequenceNames);

        // returns number of sequences in dictionary
        int Size(void) const;

        // retrieves a modifiable reference to the SamSequence object associated with this name
        SamSequence& operator[](const std::string& sequenceName);

    // retrieve STL-compatible iterators
    public:
        SamSequenceIterator      Begin(void);               // returns iterator to begin()
        SamSequenceConstIterator Begin(void) const;         // returns const_iterator to begin()
        SamSequenceConstIterator ConstBegin(void) const;    // returns const_iterator to begin()
        SamSequenceIterator      End(void);                 // returns iterator to end()
        SamSequenceConstIterator End(void) const;           // returns const_iterator to end()
        SamSequenceConstIterator ConstEnd(void) const;      // returns const_iterator to end()

    // internal methods
    private:
        int IndexOf(const std::string& name) const;

    // data members
    private:
        SamSequenceContainer m_data;
};

} // namespace BamTools

#endif // SAM_SEQUENCE_DICTIONARY_H


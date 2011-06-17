// ***************************************************************************
// BamMultiMerger_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 March 2011 (DB)
// ---------------------------------------------------------------------------
// Provides merging functionality for BamMultiReader.  At this point, supports
// sorting results by (refId, position) or by read name.
// ***************************************************************************

#ifndef BAMMULTIMERGER_P_H
#define BAMMULTIMERGER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <map>
#include <queue>
#include <string>
#include <utility>

namespace BamTools {
namespace Internal {

typedef std::pair<BamReader*, BamAlignment*> ReaderAlignment;

// generic MultiMerger interface
class IBamMultiMerger {

    public:
        IBamMultiMerger(void) { }
        virtual ~IBamMultiMerger(void) { }

    public:
        virtual void Add(const ReaderAlignment& value) =0;
        virtual void Clear(void) =0;
        virtual const ReaderAlignment& First(void) const =0;
        virtual bool IsEmpty(void) const =0;
        virtual void Remove(BamReader* reader) =0;
        virtual int Size(void) const =0;
        virtual ReaderAlignment TakeFirst(void) =0;
};

// IBamMultiMerger implementation - sorted on BamAlignment: (RefId, Position)
class PositionMultiMerger : public IBamMultiMerger {

    public:
        PositionMultiMerger(void) : IBamMultiMerger() { }
        ~PositionMultiMerger(void) { }

    public:
        void Add(const ReaderAlignment& value);
        void Clear(void);
        const ReaderAlignment& First(void) const;
        bool IsEmpty(void) const;
        void Remove(BamReader* reader);
        int Size(void) const;
        ReaderAlignment TakeFirst(void);

    private:
        typedef std::pair<int, int>           KeyType;
        typedef ReaderAlignment               ValueType;
        typedef std::pair<KeyType, ValueType> ElementType;

        typedef std::multimap<KeyType, ValueType> ContainerType;
        typedef ContainerType::iterator           DataIterator;
        typedef ContainerType::const_iterator     DataConstIterator;

        ContainerType m_data;
};

// IBamMultiMerger implementation - sorted on BamAlignment: Name
class ReadNameMultiMerger : public IBamMultiMerger {

    public:
        ReadNameMultiMerger(void) : IBamMultiMerger() { }
        ~ReadNameMultiMerger(void) { }

    public:
        void Add(const ReaderAlignment& value);
        void Clear(void);
        const ReaderAlignment& First(void) const;
        bool IsEmpty(void) const;
        void Remove(BamReader* reader);
        int Size(void) const;
        ReaderAlignment TakeFirst(void);

    private:
        typedef std::string                   KeyType;
        typedef ReaderAlignment               ValueType;
        typedef std::pair<KeyType, ValueType> ElementType;

        typedef std::multimap<KeyType, ValueType> ContainerType;
        typedef ContainerType::iterator           DataIterator;
        typedef ContainerType::const_iterator     DataConstIterator;

        ContainerType m_data;
};

// IBamMultiMerger implementation - unsorted BAM file(s)
class UnsortedMultiMerger : public IBamMultiMerger {

    public:
        UnsortedMultiMerger(void) : IBamMultiMerger() { }
        ~UnsortedMultiMerger(void) { }

    public:
        void Add(const ReaderAlignment& value);
        void Clear(void);
        const ReaderAlignment& First(void) const;
        bool IsEmpty(void) const;
        void Remove(BamReader* reader);
        int Size(void) const;
        ReaderAlignment TakeFirst(void);

    private:
        typedef ReaderAlignment ElementType;
        typedef std::vector<ReaderAlignment>  ContainerType;
        typedef ContainerType::iterator       DataIterator;
        typedef ContainerType::const_iterator DataConstIterator;

        ContainerType m_data;
};

// ------------------------------------------
// PositionMultiMerger implementation

inline void PositionMultiMerger::Add(const ReaderAlignment& value) {
    const KeyType key( value.second->RefID, value.second->Position );
    m_data.insert( ElementType(key, value) );
}

inline void PositionMultiMerger::Clear(void) {
    m_data.clear();
}

inline const ReaderAlignment& PositionMultiMerger::First(void) const {
    const ElementType& entry = (*m_data.begin());
    return entry.second;
}

inline bool PositionMultiMerger::IsEmpty(void) const {
    return m_data.empty();
}

inline void PositionMultiMerger::Remove(BamReader* reader) {

    if ( reader == 0 ) return;
    const std::string filenameToRemove = reader->GetFilename();

    // iterate over readers in cache
    DataIterator dataIter = m_data.begin();
    DataIterator dataEnd  = m_data.end();
    for ( ; dataIter != dataEnd; ++dataIter ) {
        const ValueType& entry = (*dataIter).second;
        const BamReader* entryReader = entry.first;
        if ( entryReader == 0 ) continue;

        // remove iterator on match
        if ( entryReader->GetFilename() == filenameToRemove ) {
            m_data.erase(dataIter);
            return;
        }
    }
}

inline int PositionMultiMerger::Size(void) const {
    return m_data.size();
}

inline ReaderAlignment PositionMultiMerger::TakeFirst(void) {
    DataIterator first = m_data.begin();
    ReaderAlignment next = (*first).second;
    m_data.erase(first);
    return next;
}

// ------------------------------------------
// ReadNameMultiMerger implementation

inline void ReadNameMultiMerger::Add(const ReaderAlignment& value) {
    const KeyType key(value.second->Name);
    m_data.insert( ElementType(key, value) );
}

inline void ReadNameMultiMerger::Clear(void) {
    m_data.clear();
}

inline const ReaderAlignment& ReadNameMultiMerger::First(void) const {
    const ElementType& entry = (*m_data.begin());
    return entry.second;
}

inline bool ReadNameMultiMerger::IsEmpty(void) const {
    return m_data.empty();
}

inline void ReadNameMultiMerger::Remove(BamReader* reader) {

    if ( reader == 0 ) return;
    const std::string filenameToRemove = reader->GetFilename();

    // iterate over readers in cache
    DataIterator dataIter = m_data.begin();
    DataIterator dataEnd  = m_data.end();
    for ( ; dataIter != dataEnd; ++dataIter ) {
        const ValueType& entry = (*dataIter).second;
        const BamReader* entryReader = entry.first;
        if ( entryReader == 0 ) continue;

        // remove iterator on match
        if ( entryReader->GetFilename() == filenameToRemove ) {
            m_data.erase(dataIter);
            return;
        }
    }

}

inline int ReadNameMultiMerger::Size(void) const {
    return m_data.size();
}

inline ReaderAlignment ReadNameMultiMerger::TakeFirst(void) {
    DataIterator first = m_data.begin();
    ReaderAlignment next = (*first).second;
    m_data.erase(first);
    return next;
}

// ------------------------------------------
// UnsortedMultiMerger implementation

inline void UnsortedMultiMerger::Add(const ReaderAlignment& value) {
    m_data.push_back(value);
}

inline void UnsortedMultiMerger::Clear(void) {
    for (size_t i = 0; i < m_data.size(); ++i )
        m_data.pop_back();
}

inline const ReaderAlignment& UnsortedMultiMerger::First(void) const {
    return m_data.front();
}

inline bool UnsortedMultiMerger::IsEmpty(void) const {
    return m_data.empty();
}

inline void UnsortedMultiMerger::Remove(BamReader* reader) {

    if ( reader == 0 ) return;
    const std::string filenameToRemove = reader->GetFilename();

    // iterate over readers in cache
    DataIterator dataIter = m_data.begin();
    DataIterator dataEnd  = m_data.end();
    for ( ; dataIter != dataEnd; ++dataIter ) {
        const BamReader* entryReader = (*dataIter).first;
        if ( entryReader == 0 ) continue;

        // remove iterator on match
        if ( entryReader->GetFilename() == filenameToRemove ) {
            m_data.erase(dataIter);
            return;
        }
    }
}

inline int UnsortedMultiMerger::Size(void) const {
    return m_data.size();
}

inline ReaderAlignment UnsortedMultiMerger::TakeFirst(void) {
    ReaderAlignment first = m_data.front();
    m_data.erase( m_data.begin() );
    return first;
}

} // namespace Internal
} // namespace BamTools

#endif // BAMMULTIMERGER_P_H

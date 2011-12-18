// ***************************************************************************
// BamMultiMerger_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
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

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/algorithms/Sort.h"
#include <deque>
#include <functional>
#include <set>
#include <string>

namespace BamTools {
namespace Internal {

struct MergeItem {

    // data members
    BamReader*    Reader;
    BamAlignment* Alignment;

    // ctors & dtor
    MergeItem(BamReader* reader = 0,
              BamAlignment* alignment = 0)
        : Reader(reader)
        , Alignment(alignment)
    { }

    MergeItem(const MergeItem& other)
        : Reader(other.Reader)
        , Alignment(other.Alignment)
    { }

    ~MergeItem(void) { }
};

template<typename Compare>
struct MergeItemSorter : public std::binary_function<MergeItem, MergeItem, bool> {

    public:
        MergeItemSorter(const Compare& comp = Compare())
            : m_comp(comp)
        { }

        bool operator()(const MergeItem& lhs, const MergeItem& rhs) {
            const BamAlignment& l = *lhs.Alignment;
            const BamAlignment& r = *rhs.Alignment;
            return m_comp(l,r);
        }

    private:
        Compare m_comp;
};

// pure ABC so we can just work polymorphically with any specific merger implementation
class IMultiMerger {

    public:
        IMultiMerger(void) { }
        virtual ~IMultiMerger(void) { }
    public:
        virtual void Add(MergeItem item) =0;
        virtual void Clear(void) =0;
        virtual const MergeItem& First(void) const =0;
        virtual bool IsEmpty(void) const =0;
        virtual void Remove(BamReader* reader) =0;
        virtual int Size(void) const =0;
        virtual MergeItem TakeFirst(void) =0;
};

// general merger
template<typename Compare>
class MultiMerger : public IMultiMerger {

    public:
        typedef Compare                      CompareType;
        typedef MergeItemSorter<CompareType> MergeType;

    public:
        explicit MultiMerger(const Compare& comp = Compare())
            : IMultiMerger()
            , m_data( MergeType(comp) )
        { }
        ~MultiMerger(void) { }

    public:
        void Add(MergeItem item);
        void Clear(void);
        const MergeItem& First(void) const;
        bool IsEmpty(void) const;
        void Remove(BamReader* reader);
        int Size(void) const;
        MergeItem TakeFirst(void);

    private:
        typedef MergeItem                              ValueType;
        typedef std::multiset<ValueType, MergeType>    ContainerType;
        typedef typename ContainerType::iterator       DataIterator;
        typedef typename ContainerType::const_iterator DataConstIterator;
        ContainerType m_data;
};

template <typename Compare>
inline void MultiMerger<Compare>::Add(MergeItem item) {

    // N.B. - any future custom Compare types must define this method
    //        see algorithms/Sort.h

    if ( CompareType::UsesCharData() )
        item.Alignment->BuildCharData();
    m_data.insert(item);
}

template <typename Compare>
inline void MultiMerger<Compare>::Clear(void) {
    m_data.clear();
}

template <typename Compare>
inline const MergeItem& MultiMerger<Compare>::First(void) const {
    const ValueType& entry = (*m_data.begin());
    return entry;
}

template <typename Compare>
inline bool MultiMerger<Compare>::IsEmpty(void) const {
    return m_data.empty();
}
template <typename Compare>
inline void MultiMerger<Compare>::Remove(BamReader* reader) {

    if ( reader == 0 ) return;
    const std::string& filenameToRemove = reader->GetFilename();

    // iterate over readers in cache
    DataIterator dataIter = m_data.begin();
    DataIterator dataEnd  = m_data.end();
    for ( ; dataIter != dataEnd; ++dataIter ) {
        const MergeItem& item = (*dataIter);
        const BamReader* itemReader = item.Reader;
        if ( itemReader == 0 ) continue;

        // remove iterator on match
        if ( itemReader->GetFilename() == filenameToRemove ) {
            m_data.erase(dataIter);
            return;
        }
    }
}
template <typename Compare>
inline int MultiMerger<Compare>::Size(void) const {
    return m_data.size();
}

template <typename Compare>
inline MergeItem MultiMerger<Compare>::TakeFirst(void) {
    DataIterator firstIter = m_data.begin();
    MergeItem    firstItem = (*firstIter);
    m_data.erase(firstIter);
    return firstItem;
}

// unsorted "merger"
template<>
class MultiMerger<Algorithms::Sort::Unsorted> : public IMultiMerger {

    public:
        explicit MultiMerger(const Algorithms::Sort::Unsorted& comp = Algorithms::Sort::Unsorted())
            : IMultiMerger()
        { }
        ~MultiMerger(void) { }

    public:
        void Add(MergeItem item);
        void Clear(void);
        const MergeItem& First(void) const;
        bool IsEmpty(void) const;
        void Remove(BamReader* reader);
        int Size(void) const;
        MergeItem TakeFirst(void);

    private:
        typedef MergeItem                     ValueType;
        typedef std::deque<ValueType>         ContainerType;
        typedef ContainerType::iterator       DataIterator;
        typedef ContainerType::const_iterator DataConstIterator;
        ContainerType m_data;
};

inline
void MultiMerger<Algorithms::Sort::Unsorted>::Add(MergeItem item) {
    m_data.push_back(item);
}

inline
void MultiMerger<Algorithms::Sort::Unsorted>::Clear(void) {
    m_data.clear();
}

inline
const MergeItem& MultiMerger<Algorithms::Sort::Unsorted>::First(void) const {
    return m_data.front();
}

inline
bool MultiMerger<Algorithms::Sort::Unsorted>::IsEmpty(void) const {
    return m_data.empty();
}

inline
void MultiMerger<Algorithms::Sort::Unsorted>::Remove(BamReader* reader) {

    if ( reader == 0 ) return;
    const std::string filenameToRemove = reader->GetFilename();

    // iterate over readers in cache
    DataIterator dataIter = m_data.begin();
    DataIterator dataEnd  = m_data.end();
    for ( ; dataIter != dataEnd; ++dataIter ) {
        const MergeItem& item = (*dataIter);
        const BamReader* itemReader = item.Reader;
        if ( itemReader == 0 ) continue;

        // remove iterator on match
        if ( itemReader->GetFilename() == filenameToRemove ) {
            m_data.erase(dataIter);
            return;
        }
    }
}

inline
int MultiMerger<Algorithms::Sort::Unsorted>::Size(void) const {
    return m_data.size();
}

inline
MergeItem MultiMerger<Algorithms::Sort::Unsorted>::TakeFirst(void) {
    MergeItem firstItem = m_data.front();
    m_data.pop_front();
    return firstItem;
}

} // namespace Internal
} // namespace BamTools

#endif // BAMMULTIMERGER_P_H

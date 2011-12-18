// ***************************************************************************
// Sort.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides sorting functionality.
// ***************************************************************************

#ifndef ALGORITHMS_SORT_H
#define ALGORITHMS_SORT_H

#include "api/api_global.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include <cassert>
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

namespace BamTools {
namespace Algorithms {

/*! \struct BamTools::Algorithms::Sort
    \brief Provides classes & methods related to sorting BamAlignments
*/
struct API_EXPORT Sort {

    //! Provides explicit values for specifying desired sort ordering
    enum Order { AscendingOrder = 0
               , DescendingOrder
               };

    /*! \fn template<typename ElemType> static inline bool sort_helper(const Sort::Order& order, const ElemType& lhs, const ElemType& rhs)
        \internal

        Determines necessary STL function object depending on requested Sort::Order
    */
    template<typename ElemType>
    static inline bool sort_helper(const Sort::Order& order, const ElemType& lhs, const ElemType& rhs) {
        switch ( order ) {
            case ( Sort::AscendingOrder  ) : { std::less<ElemType> comp;    return comp(lhs, rhs); }
            case ( Sort::DescendingOrder ) : { std::greater<ElemType> comp; return comp(lhs, rhs); }
            default : BT_ASSERT_UNREACHABLE;
        }
        return false; // <-- unreachable
    }

    //! Base class for our sorting function objects
    typedef std::binary_function<BamAlignment, BamAlignment, bool> AlignmentSortBase;

    /*! \struct BamTools::Algorithms::Sort::ByName
        \brief Function object for comparing alignments by name

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by name, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByName() );
            std::sort( a.begin(), a.end(), Sort::ByName(Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByName(Sort::DescendingOrder) );
        \endcode
    */
    struct ByName : public AlignmentSortBase {

        // ctor
        ByName(const Sort::Order& order = Sort::AscendingOrder)
            : m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) {
            return sort_helper(m_order, lhs.Name, rhs.Name);
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return true; }

        // data members
        private:
            const Sort::Order& m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::ByPosition
        \brief Function object for comparing alignments by position

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by position, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByPosition() );
            std::sort( a.begin(), a.end(), Sort::ByPosition(Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByPosition(Sort::DescendingOrder) );
        \endcode
    */
    struct ByPosition : public AlignmentSortBase {

        // ctor
        ByPosition(const Sort::Order& order = Sort::AscendingOrder)
            : m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) {

            // force unmapped aligmnents to end
            if ( lhs.RefID == -1 ) return false;
            if ( rhs.RefID == -1 ) return true;

            // if on same reference, sort on position
            if ( lhs.RefID == rhs.RefID )
                return sort_helper(m_order, lhs.Position, rhs.Position);

            // otherwise sort on reference ID
            return sort_helper(m_order, lhs.RefID, rhs.RefID);
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return false; }

        // data members
        private:
            Sort::Order m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::ByTag
        \brief Function object for comparing alignments by tag value

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by edit distance, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM") );
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM", Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM", Sort::DescendingOrder) );
        \endcode
    */
    template<typename T>
    struct ByTag : public AlignmentSortBase {

        // ctor
        ByTag(const std::string& tag,
              const Sort::Order& order = Sort::AscendingOrder)
            : m_tag(tag)
            , m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) {

            // force alignments without tag to end
            T lhsTagValue;
            T rhsTagValue;
            if ( !lhs.GetTag(m_tag, lhsTagValue) ) return false;
            if ( !rhs.GetTag(m_tag, rhsTagValue) ) return true;

            // otherwise compare on tag values
            return sort_helper(m_order, lhsTagValue, rhsTagValue);
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return true; }

        // data members
        private:
            std::string m_tag;
            Sort::Order m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::Unsorted
        \brief Placeholder function object

        This function object exists purely to allow for dropping a "do not care" ordering
        into methods, containers, etc that are designed to work with the other sorting objects.

        \code
            std::set<BamAlignment, Sort::ByName>;   // STL set, ordered on alignment name
            std::set<BamAlignment, Sort::Unsorted>; // STL set, unsorted (but probably insertion order)
        \endcode
    */
    struct Unsorted : public AlignmentSortBase {

        // comparison function
        inline bool operator()(const BamTools::BamAlignment&, const BamTools::BamAlignment&) {
            return false;   // returning false tends to retain insertion order
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return false; }
    };

    /*! Sorts a std::vector of alignments (in-place), using the provided compare function.

        \code
            std::vector<BamAlignemnt> a;
            // populate data

            // sort our alignment list by edit distance
            Sort::SortAlignments(a, Sort::ByTag<int>("NM"));
        \endcode

        \param[in,out] data vector of alignments to be sorted
        \param[in]     comp comparison function object
    */
    template<typename Compare>
    static inline void SortAlignments(std::vector<BamAlignment>& data,
                                      const Compare& comp = Compare())
    {
        std::sort(data.begin(), data.end(), comp);
    }

    /*! Returns a sorted copy of the input alignments, using the provided compare function.

        \code
            std::vector<BamAlignemnt> a;
            // populate data

            // get a copy of our original data, sorted by edit distance (descending order)
            std::vector<BamAligment> sortedData;
            sortedData = Sort::SortAlignments(a, Sort::ByTag<int>("NM", Sort::DescendingOrder));
        \endcode

        \param[in] input vector of alignments to be sorted
        \param[in] comp  comparison function object
        \return sorted copy of the input data
    */
    template<typename Compare>
    static inline std::vector<BamAlignment> SortAlignments(const std::vector<BamAlignment>& input,
                                                           const Compare& comp = Compare())
    {
        std::vector<BamAlignment> output(input);
        SortAlignments(output, comp);
        return output;
    }

    /*! Reads a region of alignments from a position-sorted BAM file,
        then sorts by the provided compare function

        \code
            BamReader reader;
            // open BAM file & index file

            BamRegion region;
            // define a region of interest (i.e. a exon or some other feature)

            // get all alignments covering that region, sorted by read group name
            std::vector<BamAlignments> a;
            a = Sort::GetSortedRegion(reader, region, Sort::ByTag<std::string>("RG"));
        \endcode

        \param[in] reader BamReader opened on desired BAM file
        \param[in] region desired region-of-interest
        \param[in] comp   comparison function object
        \return sorted vector of the region's alignments
    */
    template<typename Compare>
    static std::vector<BamAlignment> GetSortedRegion(BamReader& reader,
                                                     const BamRegion& region,
                                                     const Compare& comp = Compare())
    {
        // return empty container if unable to find region
        if ( !reader.IsOpen() )          return std::vector<BamAlignment>();
        if ( !reader.SetRegion(region) ) return std::vector<BamAlignment>();

        // iterate through region, grabbing alignments
        BamAlignment al;
        std::vector<BamAlignment> results;
        while ( reader.GetNextAlignmentCore(al) )
            results.push_back(al);

        // sort & return alignments
        SortAlignments(results, comp);
        return results;
    }

    /*! Reads a region of alignments from position-sorted BAM files,
        then sorts by the provided compare function

        \code
            BamMultiReader reader;
            // open BAM files & index files

            BamRegion region;
            // define a region of interest (i.e. a exon or some other feature)

            // get all alignments covering that region, sorted by read group name
            std::vector<BamAlignments> a;
            a = Sort::GetSortedRegion(reader, region, Sort::ByTag<std::string>("RG"));
        \endcode

        \param[in] reader BamMultiReader opened on desired BAM files
        \param[in] region desired region-of-interest
        \param[in] comp   comparison function object
        \return sorted vector of the region's alignments
    */
    template<typename Compare>
    static std::vector<BamAlignment> GetSortedRegion(BamMultiReader& reader,
                                                     const BamRegion& region,
                                                     const Compare& comp = Compare())
    {
        // return empty container if unable to find region
        if ( !reader.HasOpenReaders() )  return std::vector<BamAlignment>();
        if ( !reader.SetRegion(region) ) return std::vector<BamAlignment>();

        // iterate through region, grabbing alignments
        BamAlignment al;
        std::vector<BamAlignment> results;
        while ( reader.GetNextAlignmentCore(al) )
            results.push_back(al);

        // sort & return alignments
        SortAlignments(results, comp);
        return results;
    }
};

} // namespace Algorithms
} // namespace BamTools

#endif // ALGORITHMS_SORT_H

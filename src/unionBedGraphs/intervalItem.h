/*****************************************************************************
  intervalItem.h

  (c) 2010 - Assaf Gordon
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef INTERVALITEM_H
#define INTERVALITEM_H

#include <string>
#include <queue>

enum COORDINATE_TYPE {
    START,
    END
};

/*
   An interval item in the priority queue.

   An IntervalItem can mark either a START position or an END position.
 */
class IntervalItem
{
private:
    IntervalItem();

public:
    int source_index;           // which source BedGraph file this came from
    COORDINATE_TYPE coord_type; // is this the start or the end position?
    CHRPOS coord;
    std::string depth;

    IntervalItem(int _index, COORDINATE_TYPE _type, CHRPOS _coord, std::string _depth) :
        source_index(_index),
        coord_type(_type),
        coord(_coord),
        depth(_depth)
    {}

    IntervalItem(const IntervalItem &other) :
        source_index(other.source_index),
        coord_type(other.coord_type),
        coord(other.coord),
        depth(other.depth)
    {}

    bool operator< ( const IntervalItem& other ) const
    {
        return this->coord > other.coord;
    }
};

// our priority queue
typedef std::priority_queue<IntervalItem> INTERVALS_PRIORITY_QUEUE;

#endif

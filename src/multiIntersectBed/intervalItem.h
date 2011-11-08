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


public:
    int source_index;           // which source BedGraph file this came from
    COORDINATE_TYPE coord_type; // is this the start or the end position?
    CHRPOS coord;

    IntervalItem () :
       source_index(-1),
       coord_type(START),
       coord(0)
    {}

    IntervalItem(int _index, COORDINATE_TYPE _type, CHRPOS _coord) :
        source_index(_index),
        coord_type(_type),
        coord(_coord)
    {}

    IntervalItem(const IntervalItem &other) :
        source_index(other.source_index),
        coord_type(other.coord_type),
        coord(other.coord)
    {}

    bool operator< ( const IntervalItem& other ) const
    {
        if (this->coord > other.coord) {
            return true;
        }
        else if (this->coord < other.coord) {
            return false;
        }
        // prefer ENDs to come before STARTS if the same coordinate.
        // needed for stability of the -cluster algorithm.
        else {
            if (this->coord_type == END) {
                return false;
            }
            else {
                return true;
            }
        }
    }
};

// our priority queue
typedef std::priority_queue<IntervalItem> INTERVALS_PRIORITY_QUEUE;

#endif

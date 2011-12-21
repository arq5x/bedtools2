/*****************************************************************************
  Point.h

  (c) 2010 - Assaf Gordon, Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef POINT_H
#define POINT_H

#include <string>
#include <queue>

enum COORDINATE_TYPE {
    START,
    END
};


/*
   An Point can mark either a START position or an END position.
 */
class Point {
public:
    int source_index;           // which source BedGraph file this came from
    COORDINATE_TYPE coord_type; // is this the start or the end position?
    CHRPOS coord;

    Point () :
       source_index(-1),
       coord_type(START),
       coord(0)
    {}

    Point(int _index, COORDINATE_TYPE _type, CHRPOS _coord) :
        source_index(_index),
        coord_type(_type),
        coord(_coord)
    {}

    Point(const Point &other) :
        source_index(other.source_index),
        coord_type(other.coord_type),
        coord(other.coord)
    {}

    bool operator< ( const Point& other ) const
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


/*
   A specialized point for BEDGRAPH files.
 */
class PointWithDepth {
  public:
      int source_index;           // which source BedGraph file this came from
      COORDINATE_TYPE coord_type; // is this the start or the end position?
      CHRPOS coord;
      std::string depth;

      PointWithDepth(int _index, COORDINATE_TYPE _type, CHRPOS _coord, std::string _depth) :
          source_index(_index),
          coord_type(_type),
          coord(_coord),
          depth(_depth)
      {}

      PointWithDepth(const PointWithDepth &other) :
          source_index(other.source_index),
          coord_type(other.coord_type),
          coord(other.coord),
          depth(other.depth)
      {}
      
      bool operator< ( const PointWithDepth& other ) const
      {
          return this->coord > other.coord;
      }
};


// our priority queue
typedef std::priority_queue<Point> POINT_PQUEUE;
typedef std::priority_queue<PointWithDepth> POINTWITHDEPTH_PQUEUE;

#endif

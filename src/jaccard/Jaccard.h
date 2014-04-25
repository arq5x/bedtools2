/*****************************************************************************
  jaccard.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef JACCARD_H
#define JACCARD_H

#include "ContextJaccard.h"

class BlockMgr;

class Jaccard {

public:

    Jaccard(ContextJaccard *context);
    ~Jaccard();

    bool calculate();

private:
    ContextJaccard *_context;
    BlockMgr *_blockMgr;
    unsigned long _intersectionVal;
    unsigned long _unionVal;
    int _numIntersections;

    bool getIntersectionAndUnion();
    unsigned long getTotalIntersection(RecordKeyList *hits);
};

#endif /* JACCARD_H */

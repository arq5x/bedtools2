/*
 * newClosestFile.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#ifndef NEWCLOSESTFILE_H_
#define NEWCLOSESTFILE_H_

#include "ContextClosest.h"

using namespace std;

class RecordOutputMgr;

class ClosestFile {

public:
    ClosestFile(ContextClosest *context);
    ~ClosestFile(void);

    bool getClosest();

private:
    ContextClosest *_context;
    RecordOutputMgr *_recordOutputMgr;
};


#endif /* NEWCLOSESTFILE_H_ */

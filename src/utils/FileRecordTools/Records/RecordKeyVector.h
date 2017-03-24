/*
 * RecordKeyVector.h
 *
 *  Created on: Dec 14, 2012
 *      Author: nek3d
 */

#ifndef KEYVECTOR_H_
#define KEYVECTOR_H_


#include "Record.h"
#include <vector>

using namespace std;


class RecordKeyVector {
public:
	typedef Record * elemType;
	typedef vector<Record *> vecType;
	typedef vecType::iterator iterator_type;
	RecordKeyVector();
    RecordKeyVector(elemType item);
    RecordKeyVector(elemType item, vecType *vector);
    ~RecordKeyVector();

    const RecordKeyVector &operator=(const RecordKeyVector &other);
    iterator_type begin();
    iterator_type next();
    iterator_type end();
    void erase();
    size_t size() const ;
    bool empty() const ; //only checks whether list is empty. Doesn't check key.
    void push_back(elemType item);
    elemType getKey() const ;
    void setKey(elemType key);

    //setVectorNoCopy will make our list share the nodes of another
    //list, not copy them.
    void setVector(vecType *vec);

    // WARNING! clearVector will NOT delete records pointed to by list nodes. Caller must do that separately, since the RecordKeyVector
    // does not have it's own RecordMgr member.
    void clearVector();

    void clearAll() {
    	setKey(NULL);
    	clearVector();
    }
    bool allClear() { return (_key == NULL && empty()); }

    void sortVector();
    void swap(RecordKeyVector &other);


private:
    elemType _key;
    vecType *_recVec;
    int _currPos;
    bool _mustDeleteVec;
};


#endif /* RECORDKEYLIST_H_ */

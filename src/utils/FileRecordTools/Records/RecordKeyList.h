/*
 * RecordKeyList.h
 *
 *  Created on: Dec 14, 2012
 *      Author: nek3d
 */

#ifndef KEYLIST_H_
#define KEYLIST_H_

#include "Record.h"
#include "RecordList.h"

class RecordKeyList {
public:
	typedef Record * elemType;
	typedef RecordList listType;
	typedef const RecordListNode *const_iterator_type;
	RecordKeyList();
    RecordKeyList(elemType item);
    RecordKeyList(elemType item, listType &list);
    ~RecordKeyList();

    RecordKeyList &operator=(RecordKeyList &other);
    const_iterator_type begin();
    const_iterator_type next();
    const_iterator_type end();
    size_t size() const ;
    bool empty() const ; //only checks whether list is empty. Doesn't check key.
    void push_back(elemType item);
    elemType getKey() const ;
    void setKey(elemType key);

    //setListNoCopy will make our list share the nodes of another
    //list, not copy them.
    void setListNoCopy(listType &list);

    // WARNING! clearList will NOT delete records pointed to by list nodes. Caller must do that separately, since the RecordKeyList
    // does not have it's own RecordMgr member.
    void clearList();

    void clearAll() {
    	setKey(NULL);
    	clearList();
    }
    bool allClear() { return (_key == NULL && empty()); }

    void sort() { _list.sort(); }


private:
    elemType _key;
    listType _list;
};


#endif /* RECORDKEYLIST_H_ */

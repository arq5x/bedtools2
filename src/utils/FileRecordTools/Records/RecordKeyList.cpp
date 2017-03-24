#include "RecordKeyList.h"

RecordKeyList::RecordKeyList()
: _key(NULL)
{
}

RecordKeyList::RecordKeyList(Record * item)
: _key(item)
{
}

RecordKeyList::RecordKeyList(Record * item, listType &list)
: _key(item)
{
	_list = list;
}

RecordKeyList::~RecordKeyList() {
}

RecordKeyList &RecordKeyList::operator=(RecordKeyList &other)
{
	setKey(other._key);
	_list = other._list;
	return *this;
}

RecordKeyList::const_iterator_type RecordKeyList::begin()  {
	return _list.begin();
}

RecordKeyList::const_iterator_type RecordKeyList::next()  {
	return _list.next();
}


RecordKeyList::const_iterator_type RecordKeyList::end() {
	return _list.end();
}

size_t RecordKeyList::size() const {
	return _list.size();
}

bool RecordKeyList::empty() const {
	return _list.empty();
}

void RecordKeyList::push_back(elemType item) {
	_list.push_back(item);
}

Record *RecordKeyList::getKey() const {
	return _key;
}

void RecordKeyList::setKey(elemType key) {
	_key = key;
}

void RecordKeyList::setListNoCopy(listType &list) {
	_list.assignNoCopy(list);
}

void RecordKeyList::clearList() {
	_list.clear();
}

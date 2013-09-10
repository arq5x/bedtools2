#include "RecordKeyList.h"

RecordKeyList::RecordKeyList()
: _key(NULL)
{
}

RecordKeyList::RecordKeyList(const Record * item)
: _key(item)
{
}

RecordKeyList::RecordKeyList(const Record * item, const listType &list)
: _key(item)
{
	_list = list;
}

RecordKeyList::~RecordKeyList() {
}

const RecordKeyList &RecordKeyList::operator=(const RecordKeyList &other)
{
	setKey(other._key);
	_list = other._list;
	return *this;
}

const RecordKeyList::const_iterator_type RecordKeyList::begin()  {
	return _list.begin();
}

const RecordKeyList::const_iterator_type RecordKeyList::next()  {
	return _list.next();
}


const RecordKeyList::const_iterator_type RecordKeyList::end() {
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

const Record *RecordKeyList::getKey() const {
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

#include "SingleLineDelimTransferBuffer.h"
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

SingleLineDelimTransferBuffer::SingleLineDelimTransferBuffer(int numFields, char delimChar)
: _numFields(numFields),
  _delimChar(delimChar)
{
	_fields = new char *[_numFields];
	for (int i=0; i < numFields; i++) {
		_fields[i] = new char[MAX_FIELD_SIZE];
		memset(_fields[i], 0, MAX_FIELD_SIZE);
	}
}

SingleLineDelimTransferBuffer::~SingleLineDelimTransferBuffer()
{
	for (int i=0; i < _numFields; i++) {
		delete [] _fields[i];
		_fields[i] = NULL;
	}
	delete [] _fields;
	_fields = NULL;
}

bool SingleLineDelimTransferBuffer::initFromInput(const char *inBuffer)
{
	clear();
	if (strlen(inBuffer) == 0) {
		return false;
	}
	const char *oldCursor = (const char *)inBuffer;
	const char *newCursor = (const char *)inBuffer;
	int fieldNum = 0;
	while (oldCursor[0] != '\0' && oldCursor[0] != '\n') {
		while (newCursor[0] != _delimChar && newCursor[0] != '\0') {
			newCursor++;
		}
		if (newCursor != oldCursor) {
			memcpy(_fields[fieldNum], oldCursor, newCursor - oldCursor);
			fieldNum++;
		}
		if (newCursor[0] != '\0') {
			oldCursor = newCursor +1;
			newCursor = oldCursor;
		} else {
			break;
		}
	}
	if (fieldNum == _numFields) {
		return true;
	}
	return false;
}

const char *SingleLineDelimTransferBuffer::getField(int fieldNum) const
{
	if (fieldNum < _numFields) {
		return _fields[fieldNum];
	}
	cerr << "Error: Requested field " << fieldNum << "from TransferBuffer with only " << _numFields << " fields." <<endl;
	exit(1);
}

void SingleLineDelimTransferBuffer::clear(void)
{
	for (int i=0; i < _numFields; i++) {
		memset(_fields[i], 0, MAX_FIELD_SIZE);
	}
}

/*
 * Tokenizer.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: nek3d
 */
#include "Tokenizer.h"
#include <cstring>
#include <cstdio>
Tokenizer::Tokenizer()
: _numExpectedElems(0),
  _keepFinalIncElem(USE_NOW),
  _numValidElems(0) {

	_elems.resize(INITIAL_NUM_ELEMS, NULL);
	for (int i = 0; i < INITIAL_NUM_ELEMS; i++) {
		_elems[i] = new QuickString();
	}
}

Tokenizer::~Tokenizer() {
	resize(0); //easy way to delete elems without repeating code.
}

void Tokenizer::setNumExpectedItems(int newSize) {
	_numExpectedElems = newSize;
	resize(newSize);
}

int Tokenizer::tokenize(const QuickString &str, char delimiter, bool eofHit, bool isCompressed) {

	int strLen = (int)str.size();

	int startPos = 0;
	int currPos = 0;

	int currIdx = 0;

	while (startPos < strLen) {
		while (str[currPos] != delimiter && currPos < strLen) {
			currPos++;
		}
		if (currPos > startPos) {
			if ((currPos == strLen && _keepFinalIncElem != USE_NOW) &&
					( (!(delimiter == '\n' && eofHit)) || isCompressed)) {
				//we found an incomplete final element.
				// if we're ignoring incomplete elems, do nothing with it.
				currIdx--; //make sure it's not included in the final count of valid elems.

			} else {
				QuickString *newStr = fetchElem(currIdx);
				newStr->assign(str.c_str() + startPos, min(currPos, strLen) - startPos);

				// If splitting lines, strip any white space from the end of the line
				// including DOS newline characters and excess tabs.
				if (delimiter == '\n') {
					int lastPos = newStr->size();
					while (isspace(newStr->at(lastPos-1))) lastPos--;
					newStr->resize(lastPos);
				}
			}
		}
		startPos = currPos +1;
		currPos = startPos;
		currIdx++;
	}
	_numValidElems = currIdx;
	return currIdx;
}

void Tokenizer::setKeepFinalIncompleteElem(lastElemCode code) {
	_keepFinalIncElem = code;
}

QuickString *Tokenizer::fetchElem(int idx)
{
	if (idx >= (int)_elems.size()) {
		resize(idx +1);
	}
	return _elems[idx];
}


void Tokenizer::resize(int newSize) {
	int oldSize = (int)_elems.size();

	if (newSize > oldSize) { //need to add items.
		_elems.resize(newSize);
		for (int i=oldSize; i < newSize; i++) {
			_elems[i] = new QuickString();
		}
	} else if (oldSize > newSize) {
		//need to remove items.
		for (int i = oldSize - 1; i >= newSize; i--) {
			delete _elems[i];
			_elems[i] = NULL;
		}
		_elems.resize(newSize);
	}
	//if oldSize is the same as newSize, do nothing.
}


#include "PlusFields.h"
#include "SingleLineDelimTextFileReader.h"

PlusFields::PlusFields()
: _numFields(0), _split(NULL)
{
}

PlusFields::PlusFields(const PlusFields &other)
: _text(other._text), _numFields(other._numFields), _lens(other._lens), _split(NULL),
  _numOffsetFields(other._numOffsetFields)
{
}

PlusFields &PlusFields::operator=(const PlusFields &other)
{
	if (this != &other) {
		_text = other._text;
		_numFields = other._numFields;
		_lens = other._lens;
		_numOffsetFields = other._numOffsetFields;
		delete _split;
		_split = NULL;
	}
	return *this;
}

PlusFields::~PlusFields()
{
	delete _split;
}

bool PlusFields::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	size_t numFields = fileReader->getNumFields() - _numOffsetFields;

	_numFields = numFields;
	if (_split != NULL) {
		delete _split;
		_split = NULL;
	}

	if (numFields == 0) {
		_text.clear();
		_lens.clear();
		return true;
	}
	const char *line = fileReader->getLineData();
	const int *delims = fileReader->getDelimPositions() + _numOffsetFields;
	int startPos = delims[0] + 1;
	_text.assign(line + startPos, delims[numFields] - startPos);
	_lens.clear();
	if (numFields > 1 && line[delims[1]] != '\t') {
		//output joins columns with tabs; keep the column lengths to split them again
		for (size_t i=0; i < numFields; i++) {
			_lens.push_back(delims[i+1] - delims[i] - 1);
			if (i > 0) _text[delims[i] - startPos] = '\t';
		}
	}
	return true;
}

void PlusFields::clear() {
	_text.clear();
	_lens.clear();
	if (_split != NULL) {
		delete _split;
		_split = NULL;
	}
}

const string &PlusFields::getField(int fieldNum) const
{
	if (_split == NULL) {
		_split = new vector<string>();
		_split->reserve(_numFields);
		size_t pos = 0;
		for (size_t i=0; i < _numFields; i++) {
			size_t end;
			if (_lens.empty()) {
				end = _text.find('\t', pos);
				if (end == string::npos) end = _text.size();
			} else {
				end = pos + _lens[i];
			}
			_split->push_back(_text.substr(pos, end - pos));
			pos = end + 1;
		}
	}
	return (*_split)[fieldNum - _numOffsetFields - 1];
}

void PlusFields::printFields(string &outBuf) const {
	outBuf.append(_text);
}


#include "Record.h"
#include <cstdio>

Record::Record()
: _chrId(-1),
  _startPos(-1),
  _endPos(-1),
  _strandVal(UNKNOWN),
  _zeroLength(false),
  _isUnmapped(false),
  _isMateUnmapped(false)
{
}

Record::~Record() {
}

const Record &Record::operator=(const Record &other)
{
	_chrName = other._chrName;
	_chrId = other._chrId;
	_startPos = other._startPos;
	_endPos = other._endPos;
	_strand = other._strand;
	_strandVal = other._strandVal;
	_name = other._name;
	return *this;
}

void Record::clear() {
	_chrName.clear();
	_chrId = -1;
	_startPos = -1;
	_endPos = -1;
	_name.clear();
	_score.clear();
	_strand.clear();
	_strandVal = UNKNOWN;
	_startPosStr.clear();
	_endPosStr.clear();
	_zeroLength = false;
	_isUnmapped = false;
	_isMateUnmapped = false;
}

bool Record::operator < (const Record &other) const
{

	if (!sameChrom(&other)) {
		return chromBefore(&other);
	}
	if (_startPos != other._startPos) {
		return _startPos < other._startPos;
	}
	if (_endPos != other._endPos) {
		return _endPos < other._endPos;
	}
	return false;
}

bool Record::operator > (const Record &other) const {
	if (!sameChrom(&other)) {
		return chromAfter(&other);
	}
	if (_startPos != other._startPos) {
		return _startPos > other._startPos;
	}
	if (_endPos != other._endPos) {
		return _endPos > other._endPos;
	}
	return false;
}

bool Record::sameChrom(const Record *other) const {
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName == other->_chrName) : (_chrId == other->_chrId);
}

bool Record::chromBefore(const Record *other) const
{
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName < other->_chrName) : (_chrId < other->_chrId);
}

bool Record::chromAfter(const Record *other) const
{
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName > other->_chrName) : (_chrId > other->_chrId);
}



bool Record::after(const Record *other) const
{
	return (_chrId == other->_chrId && _startPos >= other->_endPos);
}

bool Record::intersects(const Record *record,
			bool wantSameStrand, bool wantDiffStrand, float overlapFraction, bool reciprocal) const
{
	//must be on same chromosome
	if (!sameChrom(record)) {
		return false;
	}
	return sameChromIntersects(record, wantSameStrand, wantDiffStrand, overlapFraction, reciprocal);
}

bool Record::sameChromIntersects(const Record *record,
		bool wantSameStrand, bool wantDiffStrand, float overlapFraction, bool reciprocal) const
{
	// Special: For records that are unmapped, intersect should automatically return false
	if (_isUnmapped || record->isUnmapped()) {
		return false;
	}

	//If user requested hits only on same strand, or only on different strands,
	//rule out different strandedness first.
	//If the strand is unknown in either case, then queries regarding strandedness
	//can not be answered, so we return false;
	bool isSameStrand = (_strandVal == record->_strandVal && _strandVal != UNKNOWN);
	bool isDiffStrand = ( _strandVal != UNKNOWN && record->_strandVal != UNKNOWN && _strandVal != record->_strandVal);

	if (wantSameStrand && !isSameStrand) {
		return false; //want same, but they're not same.
	}
	if (wantDiffStrand && !isDiffStrand) {
		return false; //want different, but they're not different.
	}

	int otherStart = record->getStartPos();
	int otherEnd = record->getEndPos();

	int maxStart = max(_startPos, otherStart);
	int minEnd = min(_endPos, otherEnd);

	//rule out all cases of no intersection at all
	if (minEnd < maxStart) {
		return false;
	}


	int overlapBases = minEnd - maxStart;
	int len = _endPos - _startPos;
	int otherLen = otherEnd - otherStart;

	if ((float)overlapBases / (float)len < overlapFraction) {
		return false;
	}

	//at this point, we've ruled out strandedness, non-intersection,
	//and query coverage (overlapFraction). The only thing left is
	//database coverage.
	if (!reciprocal) {
		return true;
	}

	if ((float)overlapBases / (float)otherLen >= overlapFraction) {
		return true;
	}

	return false;
}

bool Record::coordsValid() {
	if (_startPos < 0 || _endPos < 0 || _endPos < _startPos) {
		return false;
	}
	adjustZeroLength();
	return true;
}

void Record::adjustZeroLength()
{
	if (_startPos == _endPos) {
		_zeroLength = true;
		_startPos--;
		_endPos++;
	}
}

void Record::undoZeroLength()
{
	if (_zeroLength) {
		_startPos++;
		_endPos--;
		_zeroLength = false;
	}
}

ostream &operator << (ostream &out, const Record &record)
{
	QuickString outBuf;
	record.print(outBuf);
	out << outBuf;
	return out;
}

const QuickString &Record::getField(int fieldNum) const
{
    cerr << endl << "*****" << endl
         << "*****ERROR: requested column " << fieldNum <<
         " , but record only has fields 1 - " << getNumFields() << ". Exiting." << endl
          << endl << "*****" << endl;
    exit(1);
}

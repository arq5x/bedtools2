#include "Record.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstdio>

Record::Record()
: _fileIdx(-1),
  _chrId(-1),
  _startPos(-1),
  _endPos(-1),
  _strandVal(UNKNOWN),
  _posMode(POS_NONE),
  _origStart(0),
  _origEnd(0),
  _nsMode(NS_EMPTY),
  _zeroLength(false),
  _isUnmapped(false),
  _isMateUnmapped(false),
  _isValidHit(true),
  _frm(NULL)
{}

static void chrposToDecimal(CHRPOS v, string &out) {
	char tmp[24];
	int n = 0;
	bool neg = v < 0;
	unsigned long long u = neg ? 0ULL - (unsigned long long)v : (unsigned long long)v;
	do {
		tmp[23 - n++] = '0' + (char)(u % 10);
		u /= 10;
	} while (u);
	if (neg) tmp[23 - n++] = '-';
	out.assign(tmp + 24 - n, n);
}

void Record::materializePosStrs() const {
	PosStrBox *b = _posBox.make();
	if (_posMode == POS_PLAIN) {
		chrposToDecimal(_origStart, b->start);
		chrposToDecimal(_origEnd, b->end);
	} else {
		b->start.clear();
		b->end.clear();
	}
	_posMode = POS_STRINGS;
}

void Record::nsToFields() const {
	NSBox *b = _nsBox.make();
	if (_nsMode == NS_TEXT) {
		size_t t1 = _nsText.find('\t');
		size_t t2 = _nsText.find('\t', t1 + 1);
		b->name.assign(_nsText, 0, t1);
		b->score.assign(_nsText, t1 + 1, t2 - t1 - 1);
		b->strand.assign(_nsText, t2 + 1, string::npos);
		_nsText.clear();
	}
	_nsMode = NS_FIELDS;
}

bool Record::initNameScoreStrand(SingleLineDelimTextFileReader *fileReader) {
	const char *line = fileReader->getLineData();
	const int *d = fileReader->getDelimPositions();
	if (line[d[4]] != '\t' || line[d[5]] != '\t') return false;
	setNameScoreStrandText(line + d[3] + 1, d[6] - d[3] - 1);
	char c = line[d[5] + 1];
	bool one = (d[6] - d[5] == 2);
	_strandVal = (one && c == '+') ? FORWARD : ((one && c == '-') ? REVERSE : UNKNOWN);
	return true;
}

void Record::setNameScoreStrandText(const char *p, size_t len) {
	if (_nsMode == NS_FIELDS) {
		NSBox *b = _nsBox.get();
		b->name.clear();
		b->score.clear();
		b->strand.clear();
	}
	_nsText.assign(p, len);
	_nsMode = NS_TEXT;
}

void Record::appendNameScoreStrand(string &buf) const {
	if (_nsMode == NS_TEXT) {
		buf.append(_nsText);
	} else if (_nsMode == NS_FIELDS) {
		const NSBox *b = _nsBox.get();
		append_cstr(buf, b->name);
		buf.push_back('\t');
		append_cstr(buf, b->score);
		buf.push_back('\t');
		append_cstr(buf, b->strand);
	} else {
		buf.append("\t\t", 2);
	}
}

void Record::setFastCoords(const ChrRef &chr, CHRPOS s, CHRPOS e, int fileIdx) {
	_fileIdx = fileIdx;
	_chrName = chr;
	_startPos = s;
	_endPos = e;
	_posMode = POS_PLAIN;
	_origStart = s;
	_origEnd = e;
	_zeroLength = false;
	adjustZeroLength();
}

// p..p+len is "name\tscore\tstrand" of a tab-delimited line.
void Record::setFastNameScoreStrand(const char *p, size_t len) {
	setNameScoreStrandText(p, len);
	const char *q = p + len;
	while (q > p && q[-1] != '\t') q--;
	_strandVal = (p + len - q == 1) ? (*q == '+' ? FORWARD : (*q == '-' ? REVERSE : UNKNOWN)) : UNKNOWN;
}

Record::~Record()
{
	//_frm->deleteRecord(this);
}

const Record &Record::operator=(const Record &other)
{
	_fileIdx = other._fileIdx;
	_chrName = other._chrName;
	_chrId = other._chrId;
	_startPos = other._startPos;
	_endPos = other._endPos;
	strandRef() = other.strandRef();
	_strandVal = other._strandVal;
	nameRef() = other.nameRef();
	return *this;
}

void Record::clear() {
	_fileIdx = -1;
	_chrName.clear();
	_chrId = -1;
	_startPos = -1;
	_endPos = -1;
	if (_nsMode == NS_FIELDS) {
		NSBox *b = _nsBox.get();
		b->name.clear();
		b->score.clear();
		b->strand.clear();
	} else if (_nsMode == NS_TEXT) {
		_nsText.clear();
	}
	_nsMode = NS_EMPTY;
	_strandVal = UNKNOWN;
	_posMode = POS_NONE;
	_zeroLength = false;
	_isUnmapped = false;
	_isMateUnmapped = false;
	_isValidHit = false;
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

bool Record::lessThan(const Record *other) const
{
	if (!sameChrom(other)) {
		return chromBefore(other);
	}
	if (_startPos != other->_startPos) {
		return _startPos < other->_startPos;
	}
	if (_endPos != other->_endPos) {
		return _endPos < other->_endPos;
	}
	return false;
}

bool Record::greaterThan(const Record *other) const
{
	if (!sameChrom(other)) {
		return chromAfter(other);
	}
	if (_startPos != other->_startPos) {
		return _startPos > other->_startPos;
	}
	if (_endPos != other->_endPos) {
		return _endPos > other->_endPos;
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
	return (sameChrom(other) && _startPos >= other->_endPos);
}

bool Record::intersects(const Record *record,
                        bool sameStrand,
                        bool diffStrand,
                        float overlapFractionA,
                        float overlapFractionB,
                        bool reciprocalFraction,
                        bool eitherFraction,
                        bool obeySplits) const
{
	//must be on same chromosome
	if (!sameChrom(record)) {
		return false;
	}
	return sameChromIntersects(record,
                               sameStrand, diffStrand,
                               overlapFractionA, overlapFractionB,
                               reciprocalFraction, eitherFraction, obeySplits);
}

bool Record::sameChromIntersects(const Record *record,
                                 bool sameStrand,
                                 bool diffStrand,
                                 float overlapFractionA,
                                 float overlapFractionB,
                                 bool reciprocalFraction,
                                 bool eitherFraction,
                                 bool obeySplits) const
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

	if (sameStrand && !isSameStrand) {
		return false; //want same, but they're not same.
	}
	if (diffStrand && !isDiffStrand) {
		return false; //want different, but they're not different.
	}
	CHRPOS otherStart = record->getStartPos();
	CHRPOS otherEnd = record->getEndPos();

	bool otherZeroLen = (otherStart - otherEnd == 0);
	CHRPOS maxStart = max(_startPos, otherStart);
	CHRPOS minEnd = min(_endPos, otherEnd);

	bool localZeroLen = (_endPos - _startPos == 0);
	//rule out all cases of no intersection at all
	if (minEnd < maxStart) {
		return false;
	}

	if ((overlapFractionA == 0.0) && (overlapFractionB == 0.0))
    {
		//don't care about amount of overlap.
		//however, if minEnd and maxStart are equal, and
		//neither record is zeroLen, return false.
		if (minEnd == maxStart && !otherZeroLen && !localZeroLen) {
			return false;
		}
		return true;
	}

	CHRPOS overlapBases = minEnd - maxStart;
	CHRPOS aLen = _endPos - _startPos;
	CHRPOS bLen = otherEnd - otherStart;

    float overlapA = (float)overlapBases / (float)aLen;
    float overlapB = (float)overlapBases / (float)bLen;

    bool sufficentFractionA = (overlapA >= overlapFractionA);
    bool sufficentFractionB = (overlapB >= overlapFractionB);

    // if there is any overlap at all, we need to consider if when
    // -split is involved, as overlaps can span multiple
    // database records and query chunks, so the final calculatiom
    // of the overlap fractions must be handled in findBlockedOverlaps()
    if ((overlapBases > 0) && (obeySplits == true))
    {
    	return true;
    }
    
    if (!eitherFraction)
    {
        if (sufficentFractionA && sufficentFractionB) { return true; }
        return false;
    }
    else {
        if (sufficentFractionA || sufficentFractionB) { return true; }
        return false;
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
	string outBuf;
	record.print(outBuf);
	out << outBuf;
	return out;
}

const string &Record::getField(int fieldNum) const
{
    cerr << endl << "*****" << endl
         << "*****ERROR: requested column " << fieldNum <<
         " , but record only has fields 1 - " << getNumFields() << ". Exiting." << endl
          << endl << "*****" << endl;
    exit(1);
}


bool Record::hasChrInChromName() const {
	const char *str = _chrName.c_str();
	//if the chrom name has at least 3 characters,
	//and the first 3 are c, h, r, case-insensitive,
	//return true. Otherwise, return false.
	return ((_chrName.size() >= 3) &&
			(str[0] == 'c' || str[0] == 'C') &&
			(str[1] == 'h' || str[1] == 'H') &&
			(str[2] == 'r' || str[2] == 'R'));
}

bool Record::hasLeadingZeroInChromName(bool chrKnown) const {
	const char *str = _chrName.c_str();
	// only check to see if the digit zero follows the occurance of
	// "chr", case insensitive.
	return (_chrName.size() >= 4 && str[3] == '0' && (chrKnown || hasChrInChromName()));
}

void Record::print(FILE *fp, bool newline) const {
	string buf;
	print(buf);
	fprintf(fp, "%s", buf.c_str());
	if(newline) fprintf(fp, "\n");
}

CHRPOS Record::getLength(bool obeySplits) const {
	//only bed12 and BAM need to check splits
	return _endPos - _startPos;
}

void Record::setFileRecordManager(FileRecordMgr *frm)
{
	_frm = frm;
}

FileRecordMgr * Record::getFileRecordManager()
{
	return _frm;
}

/*
 * Record.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef RECORD_H_
#define RECORD_H_

#include <string>
#include <unordered_set>
#include <ostream>
#include <cstdarg>
#include "BedtoolsTypes.h"
#include "FreeList.h"
#include "string.h"
#include "FileRecordTypeChecker.h"


using namespace std;

// A chromosome name shared between records: one pointer into a table of interned strings.
class ChrRef {
public:
	ChrRef() : _p(&table().emptyStr) {}
	operator const string &() const { return *_p; }
	const string &str() const { return *_p; }
	bool empty() const { return _p->empty(); }
	size_t size() const { return _p->size(); }
	const char *c_str() const { return _p->c_str(); }
	ChrRef &operator=(const string &s) { assign(s.data(), s.size()); return *this; }
	ChrRef &operator=(const char *s) { assign(s, strlen(s)); return *this; }
	void clear() { _p = &table().emptyStr; }
	void assign(const char *p, size_t n) {
		Table &t = table();
		if (t.last && t.last->size() == n && memcmp(t.last->data(), p, n) == 0) { _p = t.last; return; }
		_p = t.last = &*t.set.insert(string(p, n)).first;
	}
	bool operator==(const ChrRef &o) const { return _p == o._p; }
	bool operator<(const ChrRef &o) const { return *_p < *o._p; }
	bool operator>(const ChrRef &o) const { return *_p > *o._p; }
private:
	struct Table {
		unordered_set<string> set;
		string emptyStr;
		const string *last;
		Table() : last(NULL) {}
	};
	static Table &table() { static Table t; return t; }
	const string *_p;
};
inline ostream &operator<<(ostream &o, const ChrRef &c) { return o << c.str(); }

class FileRecordMgr;
class FileReader;
class ChromIdLookup;

static inline const char* buffer_printf(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
static inline const char* buffer_printf(const char* fmt, ...) {
	static char static_buffer[1024];
	static char* dynamic_buffer = NULL;
	static size_t dynamic_buffer_size;

	char* current_buffer = dynamic_buffer ? dynamic_buffer : static_buffer;
	size_t current_buffer_size =  dynamic_buffer ? dynamic_buffer_size : sizeof(static_buffer);

	va_list ap;

	va_start(ap, fmt);

	for(;;) {
		va_list ap_copy;
		va_copy(ap_copy, ap);
		int required_size = vsnprintf(current_buffer, current_buffer_size, fmt, ap_copy);
		if(required_size < (int)current_buffer_size)
			break;

		free(dynamic_buffer);
		dynamic_buffer = current_buffer = (char*)malloc(dynamic_buffer_size = required_size + 1);
		current_buffer_size = required_size + 1;
	}

	va_end(ap);

	return current_buffer;
}


// Appends v as "%" PRId_CHRPOS would print it.
static inline void append_chrpos(string& buf, CHRPOS v) {
	char tmp[24];
	char* e = tmp + sizeof(tmp);
	char* p = e;
	uint64_t u = v < 0 ? 0 - (uint64_t)v : (uint64_t)v;
	do {
		*--p = (char)('0' + u % 10);
		u /= 10;
	} while (u);
	if (v < 0) *--p = '-';
	buf.append(p, e - p);
}

static inline void append_cstr(string& buf, const string& s) {
	buf.append(s.c_str());
}

class Record {
public:
	friend class RecordMgr;
	friend class RecordOutputMgr;

	friend class FreeList<Record>;

	virtual ~Record(); //by making the destructor protected, only the friend class(es) can actually delete Record objects, or objects derived from Record.

	typedef enum { FORWARD, REVERSE, UNKNOWN } strandType;
	Record();
	virtual bool initFromFile(FileReader *) =0;
	virtual void clear();
	virtual void print(string &) const {}
	virtual void print(string &, CHRPOS, CHRPOS ) const {}
	virtual void print(string &, const string &, const string &) const {}
	virtual void print(FILE *fp, bool newline = false) const;
	virtual void printNull(string &) const {}
	friend ostream &operator << (ostream &out, const Record &record);

	const Record & operator=(const Record &);

	virtual bool isZeroBased() const {return true;};

	virtual void setValid(const bool valid)  { _isValidHit = valid; }
	virtual bool isValid() const { return _isValidHit; }

	virtual const string &getChrName() const { return _chrName.str(); }
	virtual void setChrName(const string &chr) { _chrName = chr; }
	virtual void setChrName(const char *chr) { _chrName = chr; }

	virtual int getFileIdx() const { return _fileIdx; }
	virtual void setFileIdx(int fileIdx) { _fileIdx = fileIdx; }

	virtual int getChromId() const { return _chrId; }
	virtual void setChromId(int id) { _chrId = id; }

	virtual CHRPOS getStartPos() const { return _startPos; }
	virtual void setStartPos(CHRPOS startPos) { _startPos = startPos; }
	virtual const string &getStartPosStr() const { return startStr(); }
	virtual void setStartPosStr(const string &str) { startStr() = str; _posMode = POS_STRINGS; }

	virtual CHRPOS getEndPos() const { return _endPos; }
	virtual void setEndPos(CHRPOS endPos) { _endPos = endPos; }
	virtual const string &getEndPosStr() const { return endStr(); }
	virtual void setEndPosStr(const string &str) { endStr() = str; _posMode = POS_STRINGS; }

	virtual bool getZeroLength() const { return _zeroLength; }
	virtual void setZeroLength(bool val) { _zeroLength = val; }

	virtual const string &getStrand() const { return _strand; }
	virtual void setStrand(const string &val) {
		_strand = val;
		_strandVal = (val == "+" ? FORWARD : (val == "-" ? REVERSE : UNKNOWN));
	}
	virtual void setStrand(char val) { _strand = val;
		_strandVal = (val == '+' ? FORWARD : (val == '-' ? REVERSE : UNKNOWN));
	}
	virtual void adjustStrandVal() {
		_strandVal = (_strand == "+" ? FORWARD : (_strand == "-" ? REVERSE : UNKNOWN));
	}

	virtual strandType getStrandVal() const {return _strandVal; }

	virtual const string &getName() const { return _name; }
	virtual void setName(const string &name) { _name = name; }
	virtual void setName(const char *chr) { _name = chr; }

	virtual const string &getScore() const { return _score; }
	virtual void setScore(const string &score) { _score = score; }
	virtual void setScore(const char *chr) { _score = chr; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  = 0;

	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::UNKNOWN_RECORD_TYPE; }

	virtual bool coordsValid(); //test that no coords negative, end not less than start, check zeroLength (see below).

	//Some files can have insertions of the form 2,2. If found this should translate to cover the base before and after,
	//thus meaning the startPos is decremented and the endPos is incremented. This method will find and handle that case.
	//Don't adjust the startPosStr and endPosStr strings because they aren't used in
	//calculation. They're only used in output, and it would be slower to change them
	//and then change them back.
	virtual void adjustZeroLength();
	virtual void undoZeroLength(); //change it back just before output;
	virtual bool isZeroLength() const { return _zeroLength; }

	// "Unmapped" only applies to BamRecord, but for design reasons, it has to be here,
	// because we want to short circuit the intersects method if either record is an unmapped
	// Bam record.
	bool isUnmapped() const { return _isUnmapped; }
	bool isMateUnmapped() const { return _isMateUnmapped; }
	virtual void printUnmapped(string &) const {}



	virtual bool operator < (const Record &other) const;
	virtual bool operator > (const Record &other) const;
	virtual bool lessThan(const Record *other) const;
	virtual bool greaterThan(const Record *other) const;

	//is this on the same chromosome as another record?
	bool sameChrom(const Record *other) const;
	bool chromBefore(const Record *other) const;
	bool chromAfter(const Record *other) const;

	//is this record after the other one?
	virtual bool after(const Record *other) const;

	//does this record intersect with another record?
	virtual bool intersects(const Record *otherRecord,
							bool sameStrand,
							bool diffStrand,
							float overlapFractionA,
							float overlapFractionB,
							bool reciprocalFraction,
						    bool eitherFraction,
						    bool obeySplits) const;

	// *** WARNING !!! ** sameChromIntersects is a faster version of the intersects method,
	// BUT the caller MUST ensure that the records are on the same
	//chromosome. If you're not absolutely sure, use the regular intersects method.
	virtual bool sameChromIntersects(const Record *otherRecord,
									 bool sameStrand,
									 bool diffStrand,
									 float overlapFractionA,
									 float overlapFractionB,
									 bool reciprocalFraction,
									 bool eitherFraction,
									 bool obeySplits) const;

//	virtual static bool isNumericField(int fieldNum) const = 0;

	bool hasChrInChromName() const;
	bool hasLeadingZeroInChromName(bool chrKnown = false) const;
	virtual CHRPOS getLength(bool obeySplits) const;

	void setFileRecordManager(FileRecordMgr *frm);
	FileRecordMgr * getFileRecordManager();


protected:

	int _fileIdx; //associated file the record came from
	ChrRef _chrName;
	int _chrId;
	CHRPOS _startPos;
	CHRPOS _endPos;
	//It is actually faster to also store the start and end positions as their original strings than to
	//have to convert their integer representations back to strings when printing them.
	// The strings live in a side object, built on first use: a record parsed with plain
	// decimal coordinates keeps only the two numbers (_origStart/_origEnd) until asked.
	enum { POS_NONE = 0, POS_PLAIN = 1, POS_STRINGS = 2 };
	struct PosStrBox {
		string start;
		string end;
	};
	class PosStrPtr {
	public:
		PosStrPtr() : _p(NULL) {}
		PosStrPtr(const PosStrPtr &o) : _p(o._p ? new PosStrBox(*o._p) : NULL) {}
		PosStrPtr &operator=(const PosStrPtr &o) {
			if (this != &o) { PosStrBox *n = o._p ? new PosStrBox(*o._p) : NULL; delete _p; _p = n; }
			return *this;
		}
		~PosStrPtr() { delete _p; }
		PosStrBox *get() const { return _p; }
		PosStrBox *make() { if (!_p) _p = new PosStrBox(); return _p; }
	private:
		PosStrBox *_p;
	};
	mutable PosStrPtr _posBox;
	mutable char _posMode;
	CHRPOS _origStart;
	CHRPOS _origEnd;
	void materializePosStrs() const;
	string &startStr() const { if (_posMode != POS_STRINGS) materializePosStrs(); return _posBox.get()->start; }
	string &endStr() const { if (_posMode != POS_STRINGS) materializePosStrs(); return _posBox.get()->end; }
	string _name;
	string _score;
	string _strand;
	strandType _strandVal;
	bool _zeroLength;
	bool _isUnmapped;
	bool _isMateUnmapped;
	bool _isValidHit;
	FileRecordMgr *_frm;
};

class RecordPtrSortAscFunctor {
public:
	bool operator()(const Record *rec1, const Record *rec2) const { return *rec1 < *rec2; }
};

class RecordPtrSortDescFunctor {
public:
	bool operator()(const Record *rec1, const Record *rec2) const { return *rec1 > *rec2; }
};
#endif /* RECORD_H_ */

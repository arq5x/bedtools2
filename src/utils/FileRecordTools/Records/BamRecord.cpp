
/*
 * BamRecord.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: nek3d
 */

#include "BamRecord.h"
#include "BamFileReader.h"
#include "RecordKeyVector.h"

BamRecord::BamRecord()
: _bamChromId(-1)
{

}

BamRecord::~BamRecord()
{

}

const BamRecord &BamRecord::operator=(const BamRecord &other)
{
	Bed6Interval::operator=(other);
	_bamAlignment = other._bamAlignment;
	_bamChromId = other._bamChromId;

	_mateChrName = other._mateChrName;
	_matePos = other._matePos;
	_cigarStr = other._cigarStr;
	_insertSize = other._insertSize;
	_queryBases = other._queryBases;
	_qualities = other._qualities;

	return *this;
}


bool BamRecord::initFromFile(FileReader *fileReader)
{
        BamFileReader *bamFileReader = static_cast<BamFileReader*>(fileReader);
        return initFromFile(bamFileReader);
}

bool BamRecord::initFromFile(BamFileReader *bamFileReader)
{
	setFileIdx(bamFileReader->getFileIdx());
	_bamAlignment = bamFileReader->getAlignment();
	bamFileReader->getChrName(_chrName);

	_bamChromId = bamFileReader->getCurrChromdId();
	_startPos = bamFileReader->getStartPos();
        // NOTE: casting to int here, but bam currently only supports int32_t for position
	int2str((int)_startPos, _startPosStr);
	_endPos = bamFileReader->getEndPos();
	int2str((int)_endPos, _endPosStr);
	bamFileReader->getName(_name);
	bamFileReader->getScore(_score);
	char strandChar = bamFileReader->getStrand();
	setStrand(strandChar);

	 _isUnmapped = !bamFileReader->getAlignment().IsMapped();
	_isMateUnmapped = !bamFileReader->getAlignment().IsMateMapped();


	// Get the cigar data into a string.
	buildCigarStr();

	bamFileReader->getMateChrName(_mateChrName);
	int2str((CHRPOS)_bamAlignment.MatePosition, _matePos);
	int2str((CHRPOS)_bamAlignment.InsertSize, _insertSize);
	_queryBases = _bamAlignment.QueryBases;
	_qualities = _bamAlignment.Qualities;

	return true;
}

void BamRecord::clear()
{
	Bed6Interval::clear();
	_bamChromId = -1;

	_cigarStr.clear();
	_mateChrName.clear();
	_matePos.clear();
	_insertSize.clear();
	_queryBases.clear();
	_qualities.clear();

    //Clear the BamAlignment object. Sadly, it does not have a clear() method,
    //so we have to do each member manually.
    _bamAlignment.Name.clear();
    _bamAlignment.Length = 0;
    _bamAlignment.QueryBases.clear();
    _bamAlignment.AlignedBases.clear();
    _bamAlignment.Qualities.clear();
    _bamAlignment.TagData.clear();
    _bamAlignment.RefID = -1;
    _bamAlignment.Position = -1;
    _bamAlignment.Bin = 0;
    _bamAlignment.MapQuality = 0;
    _bamAlignment.AlignmentFlag = 0;
    _bamAlignment.CigarData.clear();
    _bamAlignment.MateRefID = -1;
    _bamAlignment.MatePosition = -1;
    _bamAlignment.InsertSize = -1;
    _bamAlignment.Filename.clear();

    _bamAlignment.SupportData.AllCharData.clear();
    _bamAlignment.SupportData.BlockLength = 0;
    _bamAlignment.SupportData.NumCigarOperations = 0;
    _bamAlignment.SupportData.QueryNameLength = 0;
    _bamAlignment.SupportData.QuerySequenceLength = 0;
    _bamAlignment.SupportData.HasCoreOnly = false;

    _bamAlignment.ErrorString.clear();

}

void BamRecord::print(string &outBuf, RecordKeyVector *keyList) const
{
	Bed6Interval::print(outBuf);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::print(string &outBuf, int start, int end, RecordKeyVector *keyList) const
{
	Bed6Interval::print(outBuf, start, end);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::print(string &outBuf, const string & start, const string & end, RecordKeyVector *keyList) const
{
	Bed6Interval::print(outBuf, start, end);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::printNull(string &outBuf) const
{
        Bed6Interval::printNull(outBuf);
        outBuf.append("\t.\t.\t.\t.\t.\t.", 12);
}

void BamRecord::printRemainingBamFields(string &outBuf, RecordKeyVector *keyList) const
{
	ostringstream s;
    s << "\t";
    s << _startPosStr;
    s << "\t";
    s << _endPos;
    s << "\t0,0,0";
    s << "\t";

    int numBlocks = (int)keyList->size();

    if (numBlocks > 0) {
            s << numBlocks;

            vector<int> blockLengths;
            vector<int> blockStarts;
            for (RecordKeyVector::iterator_type iter = keyList->begin(); iter != keyList->end(); iter = keyList->next()) {
                    const Record *block = *iter;
                    blockLengths.push_back((int)(block->getEndPos() - block->getStartPos()));
                    blockStarts.push_back((int)(block->getStartPos() - _startPos));
            }

            s << "\t";
            for (int i=0; i < (int)blockLengths.size(); i++) {
                    s << blockLengths[i];
                    s << ',';
            }
            s << "\t";
            for (int i=0; i < (int)blockStarts.size(); i++) {
                    s <<  blockStarts[i];
                    s << ',';
            }
    }
    else {
            s << "1\t0,\t0,";
    }
    outBuf.append(s.str());
}

void BamRecord::printUnmapped(string &outBuf) const {
	ostringstream s;
    s << (_chrName.empty() ? "." : _chrName);
    s << "\t-1\t-1\t";
    s << (_name.empty() ? "." : _name);
    s << "\t";
    s << (_score.empty() ? "." : _score);
    s << "\t.\t-1\t-1\t-1\t0,0,0\t0\t.\t."; // dot for strand, -1 for blockStarts and blockEnd
    outBuf.append(s.str());
}

const string &BamRecord::getField(int fieldNum) const
{
	// Unlike other records, Bam records should return the fields specified in
	// in the BAM spec. In order, they are:
//	1 - QNAME 	Query template NAME
//	2 - FLAG	bitwise Flag
//	3 - RNAME	Reference sequence NAME
//	4 - POS		1-based left most mapping POSition
//	5 - MAPQ	MAPping Quality
//	6 - CIGAR	CIGAR String
//	7 - RNEXT	Ref name of the mate/NEXT read
//	8 - PNEXT	Position of the mate/NEXT read
//	9 - TLEN	observed Template LENgth
//	10 - SEQ	segement SEQuence
//	11 - QUAL	ASCII of Phred-scaled base QUALity+33

	switch (fieldNum) {
	case 1:
		return _name;
		break;
	case 2:
		// TBD - right now, there isn't a direct way to get the flag field.
		// Context will have to error out if they try.
		break;
	case 3:
		return _chrName;
		break;
	case 4:
		return _startPosStr;
		break;
	case 5:
		return _score;
		break;
	case 6:
		return getCigarStr();
		break;
	case 7:
		return _mateChrName;
		break;
	case 8:
		return _matePos;
		break;
	case 9:
		return _insertSize;
		break;
	case 10:
		return _queryBases;
		break;
	case 11:
		return _qualities;
		break;
	default:
		cerr << "***** Error: requested invalid field number " << fieldNum << " from BAM file." << endl;
		exit(1);
		break;
	}

	// control can never reach here, but compiler is complaining about
	// reaching end of non-void function, so this dummy statement will
	// keep it happy.
	return _name;
}


bool BamRecord::isNumericField(int fieldNum) {

	//TBD: As with getField, this isn't defined for BAM.
	return (fieldNum > 6 ? false : Bed6Interval::isNumericField(fieldNum));
}

void BamRecord::buildCigarStr() {

	const vector<BamTools::CigarOp> &cigarData = _bamAlignment.CigarData;
	size_t cigarVecLen = cigarData.size();
	_cigarStr.clear();

	static char* buffer = NULL;
	static size_t buffer_size = 0;

	size_t new_buffer_size = buffer_size;

	while(new_buffer_size < cigarVecLen * 11)
		new_buffer_size = new_buffer_size == 0 ? 128 : new_buffer_size * 2;

	if(new_buffer_size != buffer_size)
	{
		free(buffer);
		buffer = (char*)malloc(new_buffer_size + 1);
		buffer_size = new_buffer_size;
	}


	if(cigarVecLen)
	{
		char* ptr = buffer + buffer_size - 1;
		ptr[1] = 0;

		for (size_t i = cigarVecLen ; i > 0; i--) {
			for(uint32_t val = cigarData[i - 1].Length; val > 0; val /= 10)
				*ptr-- = (val % 10) + '0';
			if(cigarData[i - 1].Length == 0)
				*ptr-- = '0';
			*ptr-- = cigarData[i - 1].Type;
		}

		_cigarStr.assign(ptr + 1, buffer + buffer_size - 1 - ptr);
	}
}


CHRPOS BamRecord::getLength(bool obeySplits) const {
	//only bed12 and BAM need to check splits
	if (!obeySplits) {
		return _endPos - _startPos;
	} else {
		int length = 0;
    	//parse the Cigar Data.
		const vector<BamTools::CigarOp> &cigarData = _bamAlignment.CigarData;
	    for (int i=0; i < (int)cigarData.size(); i++) {
	    	char op = cigarData[i].Type;
	    	if (op == 'M' || op == 'N' || op == 'I') {
	    		length += (int)cigarData[i].Length;
	    	}
	    }
	    return (CHRPOS)length;
	}
}



#include "BamFileReader.h"
#include "ParseTools.h"
#include <cstdio>
BamFileReader::BamFileReader()
:  _bamReader(NULL),
   _eof(false),
   _useTags(true)
{

}

BamFileReader::~BamFileReader()
{
}

bool BamFileReader::open()
{
    _bamHeader = _bamReader->GetHeaderText();
    _references = _bamReader->GetReferenceData();

	return true;
}

bool BamFileReader::isOpen() const
{
	return _bamReader->IsOpen();
}

void BamFileReader::close()
{
//	_bamReader->Close();
}

bool BamFileReader::readEntry()
{
	if (_useTags) {
		if (_bamReader->GetNextAlignment(_bamAlignment)) {
			return true;
		}
	} else {
		if (_bamReader->GetNextAlignmentCore(_bamAlignment)) {
			return true;
		}
	}
	//hit end of file
	_eof = true;
	return false;
}

void BamFileReader::getChrName(string &str) const
{
	int refId = _bamAlignment.RefID;
	if (refId < 0) {
		return;
	}
	str =  _references[refId].RefName;
}

int BamFileReader::getBamChrId() const
{
	return _bamAlignment.RefID;
}

int BamFileReader::getStartPos() const
{
	return _bamAlignment.Position;
}

int BamFileReader::getEndPos() const
{
	return _bamAlignment.GetEndPosition(false, false);
}

void BamFileReader::getName(string &str) const
{
	if (!_useTags) {
		str = _bamAlignment.SupportData.AllCharData.c_str();
	} else {
		str = _bamAlignment.Name;
	}
    if (_bamAlignment.IsFirstMate()) {
    	str += "/1";
    }
    else if (_bamAlignment.IsSecondMate()) {
    	str += "/2";
    }
}

void BamFileReader::getScore(string &str) const
{
	int2str((int)_bamAlignment.MapQuality, str);
}

char BamFileReader::getStrand() const
{
	if (_bamAlignment.IsReverseStrand()) {
		return '-';
	}
	return '+';
}

void BamFileReader::getMateChrName(string &str) const
{
	int refId = _bamAlignment.MateRefID;
	if (refId < 0) {
		return;
	}
	str =  _references[refId].RefName;
}


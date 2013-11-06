#include "BamRecord.h"
#include "BamFileReader.h"
#include "RecordKeyList.h"

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
	return *this;
}


bool BamRecord::initFromFile(FileReader *fileReader)
{
	BamFileReader *bamFileReader = static_cast<BamFileReader*>(fileReader);
	return initFromFile(bamFileReader);

}

bool BamRecord::initFromFile(BamFileReader *bamFileReader)
{
	bamFileReader->getChrName(_chrName);

	_bamChromId = bamFileReader->getCurrChromdId();
	_startPos = bamFileReader->getStartPos();
	int2str(_startPos, _startPosStr);
	_endPos = bamFileReader->getEndPos();
	int2str(_endPos, _endPosStr);
	bamFileReader->getName(_name);
	bamFileReader->getScore(_score);
	char strandChar = bamFileReader->getStrand();
	setStrand(strandChar);

	_bamAlignment = bamFileReader->getAlignment();
	_isUnmapped = !_bamAlignment.IsMapped();
	_isMateUnmapped = !_bamAlignment.IsMateMapped();
	return true;
}

void BamRecord::clear()
{
	Bed6Interval::clear();
	_bamChromId = -1;


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

void BamRecord::print(QuickString &outBuf, RecordKeyList *keyList) const
{
	Bed6Interval::print(outBuf);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::print(QuickString &outBuf, int start, int end, RecordKeyList *keyList) const
{
	Bed6Interval::print(outBuf, start, end);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::print(QuickString &outBuf, const QuickString & start, const QuickString & end, RecordKeyList *keyList) const
{
	Bed6Interval::print(outBuf, start, end);
    printRemainingBamFields(outBuf, keyList);
}

void BamRecord::printNull(QuickString &outBuf) const
{
	Bed6Interval::printNull(outBuf);
	outBuf.append("\t.\t.\t.\t.\t.\t.", 12);
}

void BamRecord::printRemainingBamFields(QuickString &outBuf, RecordKeyList *keyList) const
{
	outBuf.append('\t');
	outBuf.append(_bamAlignment.Position);
	outBuf.append('\t');
	outBuf.append(_endPos);
	outBuf.append("\t0,0,0", 6);
	outBuf.append('\t');

	int numBlocks = (int)keyList->size();

	if (numBlocks > 0) {
		outBuf.append(numBlocks);

		vector<int> blockLengths;
		vector<int> blockStarts;
		for (RecordKeyList::const_iterator_type iter = keyList->begin(); iter != keyList->end(); iter = keyList->next()) {
			const Record *block = iter->value();
			blockLengths.push_back(block->getEndPos() - block->getStartPos());
			blockStarts.push_back(block->getStartPos() - _bamAlignment.Position);
		}

		outBuf.append('\t');
		for (int i=0; i < (int)blockLengths.size(); i++) {
			outBuf.append(blockLengths[i]);
			outBuf.append(',');
		}
		outBuf.append('\t');
		for (int i=0; i < (int)blockStarts.size(); i++) {
			outBuf.append( blockStarts[i]);
			outBuf.append(',');
		}
	}
	else {
		outBuf.append("1\t0,\t0,");
	}
}

void BamRecord::printUnmapped(QuickString &outBuf) const {
	outBuf.append(_chrName.empty() ? "." : _chrName);
	outBuf.append("\t-1\t-1\t");
	outBuf.append(_name.empty() ? "." : _name);
	outBuf.append('\t');
	outBuf.append(_score.empty() ? "." : _score);
	outBuf.append("\t.\t-1\t-1\t-1\t0,0,0\t0\t.\t."); // dot for strand, -1 for blockStarts and blockEnd
}

#include "complementFile.h"
#include "NewGenomeFile.h"

ComplementFile::ComplementFile(ContextComplement *context)
    : ToolBase(context),
      _genomeFile(context->getGenomeFile()),
      _currStartPos(0),
      _onlyChromsWithBedRecords(context->getOnlyChromsWithBedRecords()),
      _outputMgr(NULL),
      _chromList(_genomeFile->getChromList()),
      _currPosInGenomeList(-1) {}

ComplementFile::~ComplementFile() {}

bool ComplementFile::init() {
  _frm = static_cast<FileRecordMergeMgr *>(upCast(_context)->getFile(0));
  return true;
}

bool ComplementFile::findNext(RecordKeyVector &hits) {
  while (!_frm->eof()) {
    _frm->getNextRecord(&hits);
    if (hits.getKey() == NULL) continue;
    return true;
  }
  return false;
}

void ComplementFile::processHits(RecordOutputMgr *outputMgr,
                                 RecordKeyVector &hits) {
  _outputMgr = outputMgr;
  const Record *rec = hits.getKey();

  // test for chrom change.
  const string &newChrom = rec->getChrName();
  if (_currChrom != newChrom) {
    outPutLastRecordInPrevChrom();

    // if record's chrom doesn't exist in the genome file, do
    // nothing
    if (!fastForward(newChrom)) return;

    // we've switched to a new chromosome that is in both the DB
    // and genome file.
    _currStartPos = 0;
    _currChrom = newChrom;
    _outRecord.setChrName(newChrom);
  }

  // warn if the record's interval is beyond the
  // length of the chromosome
  checkCoordinatesAgainstChromLength(rec);

  // safe guard against the first record for the chrom
  // starting with 0.
  if (rec->getStartPos() != 0) {
    CHRPOS endPos = rec->getStartPos();
    printRecord(endPos);
  }
  _currStartPos = rec->getEndPos();
}

void ComplementFile::checkCoordinatesAgainstChromLength(const Record *rec) {
  CHRPOS maxChromSize = _genomeFile->getChromSize(_currChrom);
  if (rec->getStartPos() > maxChromSize || rec->getEndPos() > maxChromSize) {
    cerr << "***** WARNING: " << rec->getChrName() << ":" << rec->getStartPos()
         << "-" << rec->getEndPos() << " exceeds the length of chromosome ("
         << _currChrom << ")" << endl;
  }
}

void ComplementFile::cleanupHits(RecordKeyVector &hits) {
  _frm->deleteMergedRecord(hits);
}

void ComplementFile::giveFinalReport(RecordOutputMgr *outputMgr) {
  _outputMgr = outputMgr;
  outPutLastRecordInPrevChrom();
  fastForward("");
}

void ComplementFile::outPutLastRecordInPrevChrom() {
  const string &chrom = _outRecord.getChrName();

  // do nothing if triggered by first record in DB. At this point,
  // there was no prev chrom, so nothing is stored in the output Record yet.
  if (chrom.empty()) return;
  CHRPOS maxChromSize = _genomeFile->getChromSize(chrom);
  if (_currStartPos >= maxChromSize)
    return;  // chrom already covered and reported.

  printRecord(maxChromSize);
}

bool ComplementFile::fastForward(const string &newChrom) {
  if (!newChrom.empty() && !_genomeFile->hasChrom(newChrom)) return false;

  int i = _currPosInGenomeList + 1;
  while (i < (int)_chromList.size() && _chromList[i] != newChrom) {
    _outRecord.setChrName(_chromList[i]);
    _currStartPos = 0;
    // by default, print all chroms in genomefile regardless of records
    // in the input BED file. If the _onlyChromsWithBedRecords is True,
    // however, we skip chroms without data.
    if (!_onlyChromsWithBedRecords) {
      CHRPOS endPos = _genomeFile->getChromSize(_chromList[i]);
      printRecord(endPos);
    }
    i++;
  }
  if (newChrom.empty()) return true;

  if (i == (int)_chromList.size()) {
    // reached end but didn't find new chrom. Genome and DB are not sorted in
    // same order.
    cerr << "***** ERROR: genome file and input file are not sorted in same "
            "order. Exiting..."
         << endl;
    exit(1);
    // this is where we'd return false if we weren't exiting.
  }
  _currChrom = newChrom;
  _currPosInGenomeList = i;
  return true;
}

void ComplementFile::printRecord(CHRPOS endPos) {
  _outRecord.setStartPos(_currStartPos);
  stringstream startStr;
  startStr << _currStartPos;
  _outRecord.setStartPosStr(startStr.str());

  _outRecord.setEndPos(endPos);
  stringstream endStr;
  endStr << endPos;
  _outRecord.setEndPosStr(endStr.str());

  _outputMgr->printRecord(&_outRecord);
  _outputMgr->newline();
}

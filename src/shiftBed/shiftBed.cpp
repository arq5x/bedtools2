/*****************************************************************************
  shiftBed.cpp

  (c) 2016 - David Richardson
  European Molecular Biology Laboratory, European Bioinformatics Institute
  davidr@ebi.ac.uk

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "shiftBed.h"

BedShift::BedShift(string &bedFile, string &genomeFile, float shiftMinus,
                   float shiftPlus, bool fractional, bool printHeader) {

  _bedFile = bedFile;
  _genomeFile = genomeFile;
  _shiftMinus = shiftMinus;
  _shiftPlus = shiftPlus;
  _fractional = fractional;
  _printHeader = printHeader;

  _bed = new BedFile(bedFile);
  _genome = new GenomeFile(genomeFile);

  // get going, shift it around.
  ShiftBed();
}

BedShift::~BedShift(void) {}

void BedShift::ShiftBed() {

  BED bedEntry; // used to store the current BED line from the BED file.

  _bed->Open();
  // report header first if asked.
  if (_printHeader == true) {
    _bed->PrintHeader();
  }
  while (_bed->GetNextBed(bedEntry)) {
    if (_bed->_status == BED_VALID) {
      AddShift(bedEntry);
      _bed->reportBedNewLine(bedEntry);
    }
  }
  _bed->Close();
}

void BedShift::AddShift(BED &bed) {

  CHRPOS chromSize = (CHRPOS)_genome->getChromSize(bed.chrom);

  double shift;

  if (bed.strand == "-") {
    shift = _shiftMinus;
  } else {
    shift = _shiftPlus;
  }
  if (_fractional == true)
    shift = shift * (double)bed.size();

  if ((bed.start + shift) < 0)
    bed.start = 0;
  else if ((bed.start + shift) > (chromSize - 1))
    bed.start = (chromSize - 1);
  else
    bed.start = bed.start + shift;

  if ((bed.end + shift) <= 0)
    bed.end = 1;
  else if ((bed.end + shift) > chromSize)
    bed.end = chromSize;
  else
    bed.end = bed.end + shift;
}

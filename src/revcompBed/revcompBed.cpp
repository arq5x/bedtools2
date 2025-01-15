#include "revcompBed.h"

#include "lineFileUtilities.h"

BedRevcomp::BedRevcomp(const std::string &bedFile, const std::string &genomeFile, bool printHeader):
 _bedFile(bedFile),
 _genomeFile(genomeFile),
 _printHeader(printHeader),
 _bed(_bedFile),
 _genome(genomeFile)
{
  Revcomp();
}

void BedRevcomp::Revcomp() {

  BED bedEntry; // used to store the current BED line from the BED file.

  _bed.Open();
  // report header first if asked.
  if (_printHeader == true) {
    _bed.PrintHeader();
  }
  while (_bed.GetNextBed(bedEntry)) {
    if (_bed._status == BED_VALID) {
      Revcomp(bedEntry);
      _bed.reportBedNewLine(bedEntry);
    }
  }
  _bed.Close();
}

namespace {
// like std::clamp from C++17
CHRPOS clamp(const CHRPOS& a, const CHRPOS& b, const CHRPOS& c) {
  return std::min(std::max(a, b), c);
}
}

void BedRevcomp::Revcomp(BED &bed) {
  CHRPOS chromSize = _genome.getChromSize(bed.chrom);

  auto rcStart = clamp(0, chromSize - bed.end, chromSize);
  auto rcEnd   = clamp(0, chromSize - bed.start, chromSize);

  bed.end = rcEnd;
  bed.start = rcStart;
}

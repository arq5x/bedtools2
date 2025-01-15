#pragma once

#include "bedFile.h"
#include "GenomeFile.h"

#include <string>

class BedRevcomp {

public:

    BedRevcomp(const std::string &bedFile, const std::string &genomeFile, bool printHeader);

private:

    std::string _bedFile;
    std::string _genomeFile;

    bool   _printHeader;

    BedFile _bed;
    GenomeFile _genome;

    void Revcomp();
    void Revcomp(BED &bed);
};

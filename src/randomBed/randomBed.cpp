/*****************************************************************************
  randomBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "randomBed.h"
#include "Random.h"


BedRandom::BedRandom(string &genomeFile, uint32_t numToGenerate, int seed,
                       bool haveSeed, CHRPOS length) {
    _genomeFile        = genomeFile;
    _numToGenerate     = numToGenerate;
    _length            = length;
    _haveSeed          = haveSeed;
    _seed              = seed;

    // use the supplied seed for the random
    // number generation if given.  else,
    // roll our own.
    if (!_haveSeed) {
        // thanks to Rob Long for the tip.
        _seed = (unsigned)time(0)+(unsigned)getpid();
    }
    rand_set_seed(_seed);
    Generate();
}


BedRandom::~BedRandom(void) 
{}


void BedRandom::Generate() 
{
    _genome  = new GenomeFile(_genomeFile);
    CHRPOS genomeSize  = _genome->getGenomeSize();
    
    string chrom;
    CHRPOS start;
    CHRPOS end;
    char strand;
    CHRPOS chromSize;
    uint32_t numGenerated = 0;
    while (numGenerated < _numToGenerate)
    {
        do 
        {
            CHRPOS randStart = rand_range(genomeSize);
            // use the above randomStart (e.g., for human 0..3.1billion) 
            // to identify the chrom and start on that chrom.
            pair<string, CHRPOS> location = _genome->projectOnGenome(randStart);
            chrom     = location.first;
            start     = location.second;
            end       = start + _length;
            chromSize = _genome->getChromSize(location.first);
        // keep looking if we have exceeded the end of the chrom.
        } while (end > chromSize);
        numGenerated++;
        // flip a coin for strand
        strand = (rand_range(2) >= 1)? '+' : '-';
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%d\t%" PRId_CHRPOS "\t%c\n",
            chrom.c_str(), start, end, numGenerated, end-start, strand);
    }
}



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


BedRandom::BedRandom(string &genomeFile, uint32_t numToGenerate, int seed,
                       bool haveSeed, uint32_t length) {
    _genomeFile        = genomeFile;
    _numToGenerate     = numToGenerate;
    _length            = length;
    _haveSeed          = haveSeed;
    _seed              = seed;

    // use the supplied seed for the random
    // number generation if given.  else,
    // roll our own.
    if (_haveSeed) {
        _seed = seed;
        srand(seed);
    }
    else {
        // thanks to Rob Long for the tip.
        _seed = (unsigned)time(0)+(unsigned)getpid();
        srand(_seed);
    }
    Generate();
}


BedRandom::~BedRandom(void) 
{}


void BedRandom::Generate() 
{
    _genome  = new GenomeFile(_genomeFile);
    uint32_t genomeSize  = _genome->getGenomeSize();
    
    string chrom;
    uint32_t start;
    uint32_t end;
    char strand;
    uint32_t chromSize;
    uint32_t numGenerated = 0;
    while (numGenerated < _numToGenerate)
    {
        do 
        {
            // we need to combine two consective calls to rand()
            // because RAND_MAX is 2^31 (2147483648), whereas
            // mammalian genomes are obviously much larger.
            uint32_t randStart = ((rand() << 31) | rand()) % genomeSize;
            // use the above randomStart (e.g., for human 0..3.1billion) 
            // to identify the chrom and start on that chrom.
            pair<string, int> location = _genome->projectOnGenome(randStart);
            chrom     = location.first;
            start     = location.second;
            end       = start + _length;
            chromSize = _genome->getChromSize(location.first);
        // keep looking if we have exceeded the end of the chrom.
        } while (end > chromSize);
        numGenerated++;
        // flip a coin for strand
        (rand() / double(RAND_MAX)) > 0.5 ? strand = '+' : strand = '-';
        printf("%s\t%d\t%d\t%d\t%d\t%c\n", chrom.c_str(), start, end, numGenerated, end-start, strand);
    }
}



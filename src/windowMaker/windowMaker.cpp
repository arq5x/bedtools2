/*****************************************************************************
  windowMaker.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "windowMaker.h"

WindowMaker::WindowMaker(string &genomeFile, uint32_t size, uint32_t step) 
    : _genomeFile(genomeFile)
    , _size(size)
    , _step(step)
{
    _genome = new GenomeFile(genomeFile);
    MakeWindows();
}

WindowMaker::~WindowMaker(void) {}


void WindowMaker::MakeWindows() {

    // get a list of the chroms in the user's genome
    vector<string> chromList =  _genome->getChromList();

    // process each chrom in the genome
    for (size_t c = 0; c < chromList.size(); ++c) {
        string chrom = chromList[c];
        uint32_t chrom_size = _genome->getChromSize(chrom);
        
    	for (uint32_t start = 0; start <= chrom_size; start += _step) {
    		if ((start + _size) <= chrom_size) {
    			cout << chrom << "\t" << start << "\t" << start + _size << endl;
    		}
    		else if (start < chrom_size) {
    			cout << chrom << "\t" << start << "\t" << chrom_size << endl;
    		}
    	}
    }
}


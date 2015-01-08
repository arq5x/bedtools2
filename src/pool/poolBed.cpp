/*****************************************************************************
  poolBed.h

  (c) 2015 - Pierre Lindenbaum PhD
  @yokofakun http://plindenbaum.blogspot.com
  Univ. Nantes, France

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <math.h>
#include "lineFileUtilities.h"
#include "poolBed.h"


class BedPoolItems
	{
	public:
		vector<BED> items;
		double _nbases;
		BedPoolItems():_nbases(-1)
			{
			}
		~BedPoolItems()
			{
			}
		double length(const BED* bedEntry)
		    {
	        if(bedEntry==NULL && this->_nbases>=0) return _nbases;
		    double nBases=0.0;      
	        for(size_t i=0;i< items.size();++i)
	            {
	            nBases += items[i].size();
	            }
		    if( bedEntry!=NULL) 
		        {
		        nBases += bedEntry->size();
		        }
		    else
		        {
		        this->_nbases = nBases;
		        }
		    return nBases;
		    }
	};	



BedPool::BedPool(string &bedFile, string& outfileprefix,size_t num_split):
    _bedFile(bedFile),
    _outfileprefix(outfileprefix),
    _num_split(num_split),
    _bed(0)
    {
    _bed = new BedFile(bedFile);
    doWork();
    }


BedPool::~BedPool(void) {
    if(_bed!=0) delete _bed;
    }


void BedPool::doWork() {
    size_t i,nLine(0);
    vector<BedPoolItems*> pools;
    BED bedEntry;
    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
           //pools can be created
           if( pools.size() < this->_num_split )
           		{
           		BedPoolItems* newpool=new BedPoolItems;
           		pools.push_back(newpool);
           		//pools.back().bf->bedType =  _bed->bedType ;
           		//pools.back().bf->addBEDIntoMap(bedEntry);
           		newpool->items.push_back(bedEntry);
           		}
           	//only one pool ?
           	else if( this->_num_split == 1 )
           	    {
           	    pools.front()->items.push_back(bedEntry);
           	    }
           	else
           	    {
           	    size_t index_insert = 0UL;
           	    vector<double> stdevs;
           	    double lowest_stddev=-1;
           	    /** try to insert the bedEntry in any of the pool */
           	    for(size_t  try_index_insert = 0;
           	            try_index_insert < pools.size();
           	            ++try_index_insert
           	            )
           	            {
           	            vector<double> pool_sizes;
           	            double mean = 0;
           	            /* get the mean size of the pool */
                   	    for(i=0;i< pools.size() ; ++i)
                   	        {
                   	        double num_bases = pools[i]->length(i==try_index_insert?&bedEntry:NULL);
                   	        pool_sizes.push_back( num_bases  );
                   	        mean += num_bases;
                   	        }
                   	    mean = mean / pools.size() ;
                   	    /* get the stddev of the pool */
                   	    double stddev=0.0;
                   	    for(i=0;i< pools.size() ; ++i)
                   	        {
                   	        stddev += fabs( pool_sizes[i] - mean);
                   	        }
                   	    /* this position in BED index is better = stddev of size lowest */  
                   	    if( try_index_insert == 0 || lowest_stddev> stddev )
                   	        {
                   	        lowest_stddev = stddev ;
                   	        index_insert = try_index_insert;
                   	        }
           	            }
           	    pools[index_insert]->items.push_back(bedEntry);    
           	    pools[index_insert]->_nbases=-1;/* reset size cache */
           	    }
        }
    }
    for(i=0;i< pools.size() ; ++i)
        {
        string filename(this->_outfileprefix);
        char tmp[10];
        sprintf(tmp,"%05d",(i+1));
        filename.append(".").append(tmp).append(".bed");
        FILE* out = fopen(filename.c_str(),"w");
        if(out==NULL)
            {
            fprintf(stderr,"Cannot open \"%s\". %s\n",
                filename.c_str(),
                strerror(errno)
                );
            exit(EXIT_FAILURE);
            }
        BedFile bf;
        bf.bedType = this->_bed->bedType;
        for(size_t j=0;j< pools[i]->items.size();++j)
            {
            bf.reportToFileBedNewLine(out,pools[i]->items[j]);
            }
        fclose(out);
        }

    _bed->Close();
    }


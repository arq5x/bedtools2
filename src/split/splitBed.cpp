/*****************************************************************************
  splitBed.h

  (c) 2015 - Pierre Lindenbaum PhD
  @yokofakun http://plindenbaum.blogspot.com
  Univ. Nantes, France

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include "lineFileUtilities.h"
#include "version.h"
#include "splitBed.h"

// define our program name
#define PROGRAM_NAME "bedtools split"

#define SAVE_BEANS(splits) do {\
    for(size_t _i = 0; _i < splits.size();++_i) {\
            BedSplitItems* bean=splits[_i];\
            saveBean(bean,_i);\
            delete bean;\
            } } while(0)

class Random
    {
    private:
         unsigned int seedp;
    public:
        Random():seedp(0)
            {
            }
        Random(unsigned int seed):seedp(seed)
            {
            }
        
         int nextInt()
            {
            return rand_r(&seedp);
            }
         int nextInt(int max)
            {
            return nextInt()%max;
            }
         int nextInt(int beg,int end)
            {
            return beg+nextInt(end-beg);
            }
         uint64_t nextULong()
            {
            uint64_t num;
            num = nextInt();
            num = (num << 32) | nextInt();
            return num;
            }
         uint64_t nextULong(uint64_t max)
            {
             return nextULong()%max;
            }
         uint64_t nextULong(uint64_t beg,uint64_t end)
            {
            return beg+nextULong(end-beg);
            }
         double nextDouble()
            {
            return nextInt()/(double)RAND_MAX;
            }
    };


class BedSplitItems
	{
	public:
		vector<BED*> items;
		double nbases;
		BedSplitItems():nbases(0.0)
			{
			}
		~BedSplitItems()
			{
			}
		void add(BED* entry)
		    {
		    items.push_back(entry);
		    nbases += (double)entry->size();
		    }
		vector<BED*>::size_type size() const
		    {
		    return items.size();
		    }
	};	



BedSplit::BedSplit():outfileprefix("_split"),num_beans(0)
    {
    }


BedSplit::~BedSplit(void) {
    }

int BedSplit::main(int argc,char** argv)
    {
    char* algorithm = NULL;
    for(;;) {
        static struct option long_options[] =
            {
            {"help",      no_argument, 0, 'h'},
            {"version",   no_argument, 0, 'v'},
            {"prefix",    required_argument, 0, 'p'},
            {"input",    required_argument, 0, 'i'},
            {"bed",    required_argument, 0, 'i'},
            {"number",    required_argument, 0, 'n'},
            {"algorithm",    required_argument, 0, 'a'},
            {0, 0, 0, 0}
            };
        int option_index = 0;
        int c = getopt_long (argc, argv, "vhp:i:n:a:",
            long_options,
            &option_index
            );
        if (c == -1) break;
        switch (c)
            {
            case 'v': cout << VERSION << endl; return EXIT_SUCCESS;
            case 'h': usage(cout); return EXIT_SUCCESS;
            case 'p': outfileprefix.assign(optarg); break;
            case 'i': bedFileName.assign(optarg); break;
            case 'a': algorithm = optarg ; break;
            case 'n': num_beans = (unsigned int)atoi(optarg); break;
            case '?': cerr << "Unknown option -"<< (char)optopt<< ".\n"; return EXIT_FAILURE;
            default: cerr << "Bad input" << endl; return EXIT_FAILURE;
            }
        }
    if(optind+1==argc && bedFileName.empty())
        {
        bedFileName.assign(argv[optind++]);
        }
    else if(optind==argc && bedFileName.empty())
        {
        bedFileName.assign("-");
        }
    if(optind!=argc)
        {
        cerr << "Illegal number of arguments.\n" << endl;
        return EXIT_FAILURE;
        }
    if(outfileprefix.empty())
        {
        cerr << "Error: output file prefix is empty.\n" << endl;
        return EXIT_FAILURE;
        }
    if(num_beans<=0)
        {
        cerr << "Error: num_beans==0.\n" << endl;
        usage(cerr);
        return EXIT_FAILURE;
        }
    if(algorithm==NULL || strcmp(algorithm,"size")==0 )
        {
        return doEuristicSplitOnTotalSize();
        }
    else if(strcmp(algorithm,"simple")==0)
        {
        return doSimpleSplit();
        }
    else 
        {
        cerr << "Unknow split algorithm " << algorithm << endl;
        return EXIT_FAILURE;
        }
    }


 void BedSplit::saveBean(void* data,size_t file_index)
    {
    BedSplitItems* bean=(BedSplitItems*)data;
    string filename(this->outfileprefix);
    char tmp[10];
    sprintf(tmp,"%05d",(file_index+1));
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
    bf.bedType = this->bedType;
    
    for(size_t j=0;j< bean->items.size();++j)
        {
        bf.reportToFileBedNewLine(out,*(bean->items[j]));
        }
    fflush(out);
    fclose(out);
    cout << filename << "\t"
    	<< bean->nbases << "\t"
    	<< bean->items.size()
    	<< endl;
    }

    
/* load whole bed in memory */
void BedSplit::loadBed()
  {
  BED bedEntry;
  BedFile bed(this->bedFileName);
  bed.Open();
  while(bed.GetNextBed(bedEntry))
    {
    if (bed._status != BED_VALID) continue;
    this->items.push_back(bedEntry);
    }
  bed.Close();
  this->bedType = bed.bedType;
  if(this->items.empty())
    {
    cerr << "[WARNING] no BED record found in "<< this->bedFileName << endl;
    }
  }


int BedSplit::doSimpleSplit()
    {
    vector<BedSplitItems*> splits;
    loadBed();
    for(size_t i = 0; i< this->items.size();++i)
            {
            BED* bedEntry = &(this->items[i]);
            if( splits.size() < num_beans )
           		{
           		BedSplitItems* bean=new BedSplitItems;
           		splits.push_back(bean);
           		bean->add(bedEntry);
           		}
           	else 
           	    {
           	    splits[ i % this->num_beans]->add(bedEntry);
           	    }
            }
    SAVE_BEANS(splits);
    return 0;
    }

int BedSplit::doEuristicSplitOnTotalSize()
    {
    double total_bases = 0.0;
    vector<BedSplitItems*> splits;
    loadBed();
   
    
    //biggest first
    std::sort(items.begin(),items.end(),sortBySizeDesc);
    
    
    for(size_t item_index = 0; item_index< this->items.size();++item_index)
        {
        BED* bedEntry = &(this->items[item_index]);
        if( splits.size() < num_beans )
       		{
       		BedSplitItems* bean=new BedSplitItems;
       		splits.push_back(bean);
       		bean->add(bedEntry);
       		total_bases += bedEntry->size();
       		}
       	else if( num_beans == 1 )
       	    {
       	    splits.front()->add(bedEntry);
       	    total_bases += bedEntry->size();
       	    }
       	else /* what is the best place to insert this record ?*/
       	    {
       	    size_t index_insert = 0UL;
       	    vector<double> stdevs;
       	    double lowest_stddev=-1;
       	    /* expected mean number of bases in each bean */
       	    double mean = (total_bases + bedEntry->size() )/ (double)num_beans;
       	    
       	    /** try to insert the bedEntry in any of the split */
       	    for(size_t  try_index_insert = 0;
       	            try_index_insert < splits.size();
       	            ++try_index_insert
       	            )
       	            {
               	    /* get the stddev of the split */
               	    double stddev=0.0;
               	    for(size_t i=0;i< splits.size() ; ++i)
               	        {
               	        double split_size = splits[i]->nbases + ( i==try_index_insert ? bedEntry->size() : 0 );
               	        stddev += fabs( split_size - mean);
               	        }
               	    /* this position in BED index is better = stddev of size lowest */  
               	    if( try_index_insert == 0 || lowest_stddev> stddev )
               	        {
               	        lowest_stddev = stddev ;
               	        index_insert = try_index_insert;
               	        }
       	            }
       	    splits[index_insert]->add(bedEntry);       	    
       	    total_bases += bedEntry->size();
       	    }
        }
    SAVE_BEANS(splits);
    return 0;
    }

void BedSplit::usage(std::ostream& out)
    {
    out << "\nTool:    bedtools split" << endl;
    out << "Version: " << VERSION << "\n";    
    out << "Summary: Split a Bed file." << endl << endl;
    out << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed> -n number-of-files" << endl << endl;
    out << "Options: " << endl;
    out << "\t-i|--input (file) BED input." << endl;
    out << "\t-n|--number (int)\tNumber of beans to create." << endl;
    out << "\t-p|--prefix (string)\toutput BED file prefix." << endl;
    out << "\t-a|--algorithm (string) algorithm used to split data." << endl;
    out << "\t\t* size (default): uses a heuristic algorithm to group the items so\n";
    out << "\t\t  all beans contain the ~ same number of bases\n";
    out << "\t\t* simple : slit items as they come\n";
    out << "\t\t  all beans contain the ~ same number of bases\n";
    out << "\t-h|--help\tprint help (this screen)." << endl;
    out << "\t-v|--version\tprint version." << endl;
    out << "\n\nThis programs stores the BED items in memory." << endl;
    out << endl;
    }


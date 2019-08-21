/*****************************************************************************
  splitBed.cpp

  (c) 2015 - Pierre Lindenbaum PhD
  @yokofakun http://plindenbaum.blogspot.com
  Univ. Nantes, France

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <cmath>
#include <climits>
#include <sstream>
#include <iomanip>
#include <getopt.h>
#include "lineFileUtilities.h"
#include "version.h"
#include "splitBed.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools split"

#define SAVE_BEDITEMS(splits) do {\
    for(size_t _i = 0; _i < splits.size();++_i) {\
            BedSplitItems* bedItems=splits[_i];\
            saveBedItems(bedItems,_i);\
            delete bedItems;\
            } } while(0)

#define  REPORT_CHUNK(FILENAME,NBASES,COUNT) \
    cout << FILENAME << "\t";\
    if( NBASES >= LONG_MAX) \
    	{\
    	cout << NBASES << "\t";\
    	}\
    else\
    	{\
    	cout << (long)(NBASES) << "\t";\
    	}\
    cout << COUNT << endl

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



BedSplit::BedSplit():outfileprefix("_split"),num_chuncks(0)
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
            case 'n': num_chuncks = (unsigned int)atoi(optarg); break;
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
    if(num_chuncks<=0)
        {
        cerr << "Error: num_chunks==0.\n" << endl;
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
        cerr << "Unknown split algorithm " << algorithm << endl;
        return EXIT_FAILURE;
        }
    }

std::FILE* BedSplit::saveFileChunk(std::string& filename,size_t file_index)
   {
    ostringstream name;
    name << this->outfileprefix << '.'
         << setfill('0') << setw(5) << file_index+1 << ".bed";
    filename = name.str();
    FILE* out = fopen(filename.c_str(),"w");
    
    if(out==NULL)
        {
        fprintf(stderr,"Cannot open \"%s\". %s\n",
            filename.c_str(),
            strerror(errno)
            );
        exit(EXIT_FAILURE);
        }
    return out;
    }

 void BedSplit::saveBedItems(void* data,size_t file_index)
    {
    BedSplitItems* bedItems=(BedSplitItems*)data;
    string filename;
    FILE* out = saveFileChunk(filename,file_index);
    BedFile bf;
    bf.bedType = this->bedType;
    
    for(size_t j=0;j< bedItems->items.size();++j)
        {
        bf.reportToFileBedNewLine(out,*(bedItems->items[j]));
        }
    fflush(out);
    fclose(out);
    REPORT_CHUNK(filename, bedItems->nbases,bedItems->items.size());
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

class SimpleSplitInfo
    {
    public:
        FILE* out;
        std::string filename;
        size_t count;
        double nbases;
        SimpleSplitInfo():out(0),count(0),nbases(0)
         {
         }
    };

int BedSplit::doSimpleSplit()
    {
    if( this->num_chuncks + 3 < FOPEN_MAX)//3 for security ?, not too much...
    	{
    	size_t i=0,count=0;
    	vector<SimpleSplitInfo*> outputs;
    	
    	
    	BED bedEntry;
        BedFile bed(this->bedFileName);
        bed.Open();
        while(bed.GetNextBed(bedEntry))
            {
            if (bed._status != BED_VALID) continue;
            if(count==0) this->bedType = bed.bedType;
            size_t file_index= count % this->num_chuncks;
            if(file_index == outputs.size())
                {
                SimpleSplitInfo* newinfo = new SimpleSplitInfo;
                outputs.push_back(newinfo);
                newinfo->out = saveFileChunk(newinfo->filename,file_index);
                } 
            SimpleSplitInfo* info = outputs[file_index];
            info->nbases += bedEntry.size();
            info->count++;
            bed.reportToFileBedNewLine(info->out,bedEntry);
            ++count;
            }
        
        
        for(i=0;i< outputs.size();++i)
            {
            SimpleSplitInfo* info = outputs[i];
            fflush(info->out);
            fclose(info->out);
            REPORT_CHUNK(info->filename, info->nbases,info->count);
            delete info;
            }
        bed.Close();
    	}
    else  {
        vector<BedSplitItems*> splits;
        loadBed();
        for(size_t i = 0; i< this->items.size();++i)
            {
            BED* bedEntry = &(this->items[i]);
            if( splits.size() < this->num_chuncks )
           		{
           		BedSplitItems* bedItems=new BedSplitItems;
           		splits.push_back(bedItems);
           		bedItems->add(bedEntry);
           		}
           	else 
           	    {
           	    splits[ i % this->num_chuncks]->add(bedEntry);
           	    }
            }
        SAVE_BEDITEMS(splits);
	    }
    return 0;
    }

int BedSplit::doEuristicSplitOnTotalSize()
    {
    long double total_bases = 0.0;
    vector<BedSplitItems*> splits;
    loadBed();
   
    
    //biggest first
    std::sort(items.begin(),items.end(),sortBySizeDesc);
    
    
    for(size_t item_index = 0; item_index< this->items.size();++item_index)
        {
        BED* bedEntry = &(this->items[item_index]);
        if( splits.size() < num_chuncks )
       		{
       		BedSplitItems* bedItems=new BedSplitItems;
       		splits.push_back(bedItems);
       		bedItems->add(bedEntry);
       		total_bases += bedEntry->size();
       		}
       	else if( num_chuncks == 1 )
       	    {
       	    splits.front()->add(bedEntry);
       	    total_bases += bedEntry->size();
       	    }
       	else /* what is the best place to insert this record ?*/
       	    {
       	    size_t index_insert = 0UL;
       	    vector<double> stdevs;
       	    double lowest_stddev=-1;
       	    /* expected mean number of bases in each bedItems */
       	    double mean = (total_bases + bedEntry->size() )/ (double)num_chuncks;
       	    
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
    SAVE_BEDITEMS(splits);
    return 0;
    }

void BedSplit::usage(std::ostream& out)
    {
    out << "\nTool:    bedtools split" << endl;
    out << "Version: " << VERSION << "\n";    
    out << "Summary: Split a Bed file." << endl << endl;
    out << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed> -n number-of-files" << endl << endl;
    out << "Options: " << endl;
    out << "\t-i|--input (file)\tBED input file (req'd)." << endl;
    out << "\t-n|--number (int)\tNumber of files to create (req'd)." << endl;
    out << "\t-p|--prefix (string)\tOutput BED file prefix." << endl;
    out << "\t-a|--algorithm (string) Algorithm used to split data." << endl;
    out << "\t\t* size (default): uses a heuristic algorithm to group the items" << endl;
    out << "\t\t  so all files contain the ~ same number of bases" << endl;
    out << "\t\t* simple : route records such that each split file has" << endl;
    out << "\t\t  approximately equal records (like Unix split)." << endl << endl;
    out << "\t-h|--help\t\tPrint help (this screen)." << endl;
    out << "\t-v|--version\t\tPrint version." << endl;
    out << "\n\nNote: This programs stores the input BED records in memory." << endl;
    out << endl;
    }


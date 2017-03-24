#include "ContextIntersect.h"
#include "FileRecordMgr.h"
#include <iostream>
#include <cstdio>
#include "RecordKeyList.h"
#include "NewChromsweep.h"
//#include "DualQueue.h"
#include "ParseTools.h"
#include <sstream>
#include <iomanip>
//#include "FormatGuess.h"
#include <map>

#include "PushBackStreamBuf.h"
#include "InflateStreamBuf.h"
#include "InputStreamMgr.h"
#include "BufferedStreamMgr.h"
#include "api/internal/io/BgzfStream_p.h"

//void doSweep(const Context *context);
//void testDualQueue(Context *context);
//
//void test_streams();


using namespace std;

int nek_sandbox1_main2(int argc,char** argv);

int nek_sandbox1_main(int argc,char** argv)
{

	multimap<int, char> myMap;
	myMap.insert(pair<int, char>(1, 'a'));
	myMap.insert(pair<int, char>(2, 'b'));
	myMap.insert(pair<int, char>(3, 'c'));
	myMap.insert(pair<int, char>(3, 'x'));
	myMap.insert(pair<int, char>(3, 'y'));
	myMap.insert(pair<int, char>(3, 'z'));
	myMap.insert(pair<int, char>(4, 'd'));

	cout << "Multimap contains:" << endl;
	for (multimap<int, char>::iterator iter = myMap.begin(); iter != myMap.end(); iter++) {
		cout << iter-> first << ":" << iter->second << endl;
	}
	cout << "Size is " << myMap.size() << "." << endl;

//	for (int i=0; i < 4000; i++) {
//		cout << "# This is line " << i << " of a file with a large header." << endl;
//	}
//	return 0;

//	if (argc < 2) {
//		cerr << "Error: Need one input file. Use \"-\" for stdin." << endl;
//	}
//	ifstream inFileStream(argv[1]);
//	static const int BUF_SIZE = 8192;
//	BamTools::Internal::BgzfStream bgStream;
//	bgStream.OpenStream(&inFileStream, BamTools::IBamIODevice::ReadOnly);
//	char sLine[BUF_SIZE];
//	while (bgStream.IsOpen()) {
//		memset(sLine, 0, BUF_SIZE);
//		bgStream.Read(sLine, BUF_SIZE-1);
//		if((int)strlen(sLine) < BUF_SIZE-1) {
//			bgStream.Close();
//		}
//		printf("%s", sLine);
//	}
//	return 0;
//	string filename(argv[1]);
//	istream *inputStream = NULL;
//	if (filename  == "-") {
//		inputStream = &cin;
//	} else {
//		inputStream = new ifstream(filename.c_str());
//	}
//
//	BamTools::BamReader _bamReader;
//	try {
//		_bamReader.OpenStream(inputStream);
//	}
//	catch (...) {
//		fprintf(stderr, "ERROR: Unable to open BAM file from standard input.\n");
//		exit(1);
//	}
////	try {
////		_bamReader.Open(argv[1]);
////	}
////	catch (...) {
////		fprintf(stderr, "ERROR: Unable to open BAM file %s\n",argv[1]);
////		exit(1);
////	}
////	}
//    string _bamHeader = _bamReader.GetHeaderText();
//    BamTools::RefVector _references = _bamReader.GetReferenceData();
//
//    if (_bamHeader.empty() || _references.empty()) {
//    	cout << "This is not a bam file." << endl;
//    } else {
//    	cout << "This is a BAM file." << endl;
//    }
//	return 0;
//
//	ifstream myFile(argv[1]);
//	if (!myFile.good()) {
//		cerr << "Error: Can't open genome file" << argv[1] << "Exiting..." << endl;
//		exit(1);
//	}
//	string sLine;
//	vector<string> fields;
//	string chrName;
//
//	vector<string> chroms;
//	chroms.push_back("1");
//	chroms.push_back("2");
//	chroms.push_back("10");
//	chroms.push_back("11");
//
//	vector<int> chromCounts(4, 0);
//	int chromIdx = 0;
//	while (!myFile.eof()) {
//		sLine.clear();
//		fields.clear();
//		getline(myFile, sLine);
//		if (sLine[0] == '@') {
//			cout << sLine << endl;
//			continue;
//		}
//		Tokenize(sLine.c_str(), fields);
//		const string &currChrom = fields[2];
//		if (currChrom == chroms[chromIdx]) {
//			cout << sLine << endl;
//			chromCounts[chromIdx]++;
//			if (chromCounts[chromIdx] >= 3000) {
//				chromIdx++;
//			}
//			if (chromIdx > 3) {
//				break;
//			}
//		}
//	}
//
//	return 0;
//
//	ContextIntersect context;
//	context.addInputFile(argv[1]);
//	context.setSortedInput(true);
////	context.setObeySplits(true);
//
//	FileRecordMgr frm(0, &context);
////	frm.getBlockMgr()->setBreakOnSkipOps(true);
//	if (!frm.open()) {
//		cerr << "Error: couldn't open file " << argv[1] << ". Exiting." << endl;
//		exit(1);
//	}
//	cout << "File Type is : " << frm.getFileType() << ", " << frm.getFileTypeName() << "." << endl;
//	cout << "RecordType is : " << frm.getRecordType() << ", " << frm.getRecordTypeName() << "."  << endl;
//
//	bool headerFound = false;
//	string outbuf;
//	while (!frm.eof()) {
//		Record *record = frm.getNextRecord();
//		if (!headerFound && frm.hasHeader()) {
//			cout << frm.getHeader() << endl;
//			headerFound = true;
//		}
//		if (record == NULL) {
//			break;
//		}
//
//		outbuf.clear();
//		record->print(outbuf);
//		printf("%s\n", outbuf.c_str());
//
//		RecordKeyList recList(record);
//		int blockCount = frm.getBlockMgr()->getBlocks(recList);
//		printf("The %d blocks are:\n", blockCount);
//		for (RecordKeyList::const_iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
//			iter->value()->print();
//			printf("\n");
//		}
//		printf("\n\n");
//		frm.getBlockMgr()->deleteBlocks(recList);

//		frm.deleteRecord(record);
//	}
//	cout << "Final header is: " << frm.getHeader() << endl;
//	frm.close();

	return 0;
}

//int nek_sandbox1_main2(int argc,char** argv)
//{
//
//
//	vector<FormatGuess *> formats;
//	formats.push_back(new VCFGuess());
//	formats.push_back(new XMLGuess());
//
//	PushBackStreamBuf pbs(std::cin.rdbuf());
//	std::istream in(&pbs);
//	const char* format=0;
//	for(size_t i=0;i< 2;++i)
//	{
//		if( formats[i]->guess(&pbs))
//		{
//			format=formats[i]->format();
//			break;
//		}
//	}
//	std::string line;
//	while(getline(in,line))
//	{
//		if(format!=0) cout << format << "\t";
//		cout << line << endl;
//	}
//	return 0;
//}
//
#ifdef false
int nek_sandbox1_main(int argc, char **argv) {



	for (int i=0; i < argc; i ++) {
		cout << "Arg " << i << " is: " << argv[i] << endl;
	}
	test_streams();
	return 0;

	/////////////////////////////////
	//
	//	BLOCK FOR DUEL QUEUE
	//
	if (argc < 2) {
		cerr << "Error: need at least one file name." << endl;
		return 1;
	}

	Context *context = new Context();
	context->addInputFile(argv[1]);
	testDualQueue(context);
	return 0;
	//
	////////////////////////////////





	/////////////////////////////////
	//
	//	BLOCK FOR SWEEP
	//

	if (argc < 3) {
		cerr << "Error: need at least two data file names." << endl;
		return 1;
	}

	Context *context = new Context();
	context->parseCmdArgs(argc, argv, 1);

	doSweep(context);
	return 0;
	//
	////////////////////////////////

}

void doSweep(const char *file1, const char *file2, const string &genomeFile)
{
   Context *context = new Context();
   context->addInputFile(file1);
   context->addInputFile(file2);
   context->openGenomeFile(genomeFile);

   ChromSweep sweep = ChromSweep(context);

	if (!sweep.init()) {
		cerr << "ERROR: Failure to open files in jaccard's getIntersection method." << endl;
		return;
	}

   RecordKeyList hit_set;
   while (sweep.next(hit_set)) {
//	   _intersectionVal += getTotalIntersection(&hit_set);
	   continue;
	}
   unsigned long unionVal = sweep.getQueryTotalRecordLength() + sweep.getDatabaseTotalRecordLength();
   cout << endl << endl << "Union value is: " << unionVal << endl << endl;
}

void testDualQueue(Context *context) {
	DualQueue<Record *, DualQueueAscending > dqAsc;
//	DualQueue<Record *, DualQueueDescending> dqDesc;

	FileRecordMgr frm(context->getInputFileName(0), context);
	frm.open();

	printf("Original record order is:\n");
	while (!frm.eof()) {
		Record *record = frm.getNextRecord();
		if (record == NULL) {
			continue;
		}
//		printf("\n\nNext Record is:\n");
		record->print();
		printf("\n");
		dqAsc.push(record);
//		dqDesc.push(record);
	}

	printf("\nSupposedly ascending order is:\n");
	while (!dqAsc.empty()) {
		const Record *record = dqAsc.top();
		dqAsc.pop();
		record->print();
		printf("\n");
	}

//	printf("\nSupposedly descending order is:\n");
//	while (!dqDesc.empty()) {
//		const Record *record = dqDesc.top();
//		dqDesc.pop();
//		record->print();
//		printf("\n");
//	}
	frm.close();
}

void test_streams()
{
	char myBuf[10];
	memset(myBuf, 0, 10);

	cin >> noskipws >> setw(9) >>  myBuf;

	stringstream newStream;
	newStream << "myBuf =:" << myBuf << endl;

	newStream << "Full stream was:" << endl;
	newStream << myBuf;
	newStream << cin.rdbuf();


	cout << newStream.str() << endl;

}

#endif //ifdef false

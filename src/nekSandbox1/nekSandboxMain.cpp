using namespace std;

#include "Context.h"
#include "FileRecordMgr.h"
#include <iostream>
#include <cstdio>
#include "RecordKeyList.h"
#include "NewChromsweep.h"
#include "DualQueue.h"
#include "ParseTools.h"
#include <sstream>
#include <iomanip>
//#include "FormatGuess.h"

#include "PushBackStreamBuf.h"
#include "InflateStreamBuf.h"
#include "InputStreamMgr.h"
#include "BufferedStreamMgr.h"
//void doSweep(const Context *context);
//void testDualQueue(Context *context);
//
//void test_streams();


using namespace std;

int nek_sandbox1_main2(int argc,char** argv);

int nek_sandbox1_main(int argc,char** argv)
{

	if (argc < 2) {
		cerr << "Error: Need one input file. Use \"-\" for stdin." << endl;
	}


	Context context;
	context.addInputFile(argv[1]);
	context.setSortedInput(true);
//	context.setObeySplits(true);

	FileRecordMgr frm(0, &context);
//	frm.getBlockMgr()->setBreakOnSkipOps(true);
	if (!frm.open()) {
		cerr << "Error: couldn't open file " << argv[1] << ". Exiting." << endl;
		exit(1);
	}
	cout << "File Type is : " << frm.getFileType() << ", " << frm.getFileTypeName() << "." << endl;
	cout << "RecordType is : " << frm.getRecordType() << ", " << frm.getRecordTypeName() << "."  << endl;

	bool headerFound = false;
	QuickString outbuf;
	while (!frm.eof()) {
		Record *record = frm.allocateAndGetNextRecord();
		if (!headerFound && frm.hasHeader()) {
			cout << frm.getHeader() << endl;
			headerFound = true;
		}
		if (record == NULL) {
			break;
		}

		if (record->getStartPos() == 90647945) {
			printf("Breakpoint here.\n");
		}
		outbuf.clear();
		record->print(outbuf);
		printf("%s\n", outbuf.c_str());

//		RecordKeyList recList(record);
//		int blockCount = frm.getBlockMgr()->getBlocks(recList);
//		printf("The %d blocks are:\n", blockCount);
//		for (RecordKeyList::const_iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
//			iter->value()->print();
//			printf("\n");
//		}
//		printf("\n\n");
//		frm.getBlockMgr()->deleteBlocks(recList);

		frm.deleteRecord(record);
	}
	cout << "Final header is: " << frm.getHeader() << endl;
	frm.close();

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
		Record *record = frm.allocateAndGetNextRecord();
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


#include "RegressTest.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "string.h"

void usage()  {
	printf("Usage: bedtools regressTest sub-prog targetVersion configFile [optionsToTest]\n");
}

void setOneWayOptions(RegressTest *regressTest) {
	//Use this for programs whose main arguments are two files, preceded by -a and -b, respectively.

	regressTest->setFilesPerRun(1);
	regressTest->addFilePrecessorOption("-i");
}

void setTwoWayOptions(RegressTest *regressTest) {
	//Use this for programs whose main arguments are two files, preceded by -a and -b, respectively.

	regressTest->setFilesPerRun(2);
	regressTest->addFilePrecessorOption("-a");
	regressTest->addFilePrecessorOption("-b");
}

int regress_test_main(int argc, char **argv) {

	//usage: bedtools regressTest sub-prog targetVersion [optionsToTest]
	if (argc < 5) {
		usage();
		exit(1);
	}
	string program(argv[2]);

	RegressTest *regressTest = new RegressTest();

	//set specific options for each sub-program
	if (program == "intersect") {
		setTwoWayOptions(regressTest);
	} else if (program ==  "jaccard") {
		setTwoWayOptions(regressTest);
	} else if (program == "merge") {
		setOneWayOptions(regressTest);
	} else {
		//TBD: Handle all other programs eventually
		fprintf(stderr, "Sorry, sub-program %s is not yet supported.\n", argv[2]);
		delete regressTest;
		exit(1);
	}

	if (!regressTest->init(argc, argv)) {
		fprintf(stderr, "Error: could not initialize tests for %s.\n", argv[2]);
		delete regressTest;
		exit(1);
	}

	if (!regressTest->runTests()) {
		fprintf(stderr, "Error: Failed to run tests for %s.\n", argv[2]);
		delete regressTest;
		exit(1);
	}

	delete regressTest;
	exit(0);

}

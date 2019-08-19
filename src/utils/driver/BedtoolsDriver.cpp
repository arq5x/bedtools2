#include "BedtoolsDriver.h"
#include <iostream>
#include "RecordOutputMgr.h"
//contexts
#include "ContextIntersect.h"
#include "ContextMerge.h"
#include "ContextMap.h"
#include "ContextSample.h"
#include "ContextJaccard.h"
#include "ContextClosest.h"
#include "ContextComplement.h"
#include "ContextSubtract.h"
#include "ContextSpacing.h"
#include "ContextCoverage.h"
#include "ContextSummary.h"

//tools
#include "intersectFile.h"
#include "mapFile.h"
#include "closestFile.h"
#include "mergeFile.h"
#include "jaccard.h"
#include "subtractFile.h"
#include "sampleFile.h"
#include "spacingFile.h"
#include "fisher.h"
#include "coverageFile.h"
#include "complementFile.h"
#include "groupBy.h"
#include "summaryFile.h"

BedtoolsDriver::BedtoolsDriver()
: _hadError(false) {
	_supported.insert("intersect");
	_supported.insert("map");
	_supported.insert("closest");
	_supported.insert("merge");
	_supported.insert("jaccard");
	_supported.insert("subtract");
	_supported.insert("sample");
	_supported.insert("spacing");
	_supported.insert("fisher");
	_supported.insert("coverage");
	_supported.insert("complement");
	_supported.insert("groupby");
	_supported.insert("summary");
}


bool BedtoolsDriver::supports(const string &tool) {
	supportType::iterator iter = _supported.find(tool);
	return (iter != _supported.end());
}

bool BedtoolsDriver::subMain(int argc, char **argv) 
{
	_subCmd = argv[1];
	ContextBase *context = getContext();

	//process all command line arguments, check for valid usage,
	//show help and error messages if needed.
	if (!context->testCmdArgs(argc - 1, argv + 1)) {
		_hadError = context->errorEncountered();
		_errors = context->getErrorMessages();
		delete context;
		return false;
	}

	//establish which tool we're using (intersect, map, closest, etc).
	//initialize it.
	ToolBase *tool = getTool(context);
	if (!tool->init()) {
		delete context;
		return false;
	}

	//set-up output manager.
	RecordOutputMgr *outputMgr = new RecordOutputMgr();
	outputMgr->init(context);

	//process input
	RecordKeyVector hits;
	while (tool->findNext(hits)) {
		tool->processHits(outputMgr, hits);
		tool->cleanupHits(hits);
	}
	tool->finalizeCalculations();
	tool->giveFinalReport(outputMgr);
	delete outputMgr;
	delete tool;
	delete context;
	return true;
}


ContextBase *BedtoolsDriver::getContext()
{
	ContextBase *context = NULL;
	if (_subCmd == "intersect") {
		context = new ContextIntersect();
	} else if (_subCmd == "map") {
		context = new ContextMap();
	} else if (_subCmd == "merge") {
		context = new ContextMerge();
	} else if (_subCmd == "jaccard") {
		context = new ContextJaccard();
	} else if (_subCmd == "closest") {
		context = new ContextClosest();
	} else if (_subCmd == "subtract") {
		context = new ContextSubtract();
	} else if (_subCmd == "sample") {
		context = new ContextSample();
	} else if (_subCmd == "spacing") {
		context = new ContextSpacing();
	} else if (_subCmd == "fisher") {
		context = new ContextFisher();
	} else if (_subCmd == "coverage") {
		context = new ContextCoverage();
	} else if (_subCmd == "complement") {
		context = new ContextComplement();
	} else if (_subCmd == "groupby") {
		context = new ContextGroupBy();
	} else if (_subCmd == "summary") {
		context = new ContextSummary();
	} else {
		cerr << "Error: Tool " << _subCmd << " is not supported. Exiting..." << endl;
		exit(1);
	}
	return context;
}

ToolBase *BedtoolsDriver::getTool(ContextBase *context)
{
	ToolBase *tool = NULL;
	if (_subCmd == "intersect") {
		tool = new IntersectFile(static_cast<ContextIntersect *>(context));
	} else if (_subCmd == "map") {
		tool = new MapFile(static_cast<ContextMap *>(context));
	} else if (_subCmd == "closest") {
		tool = new ClosestFile(static_cast<ContextClosest *>(context));
	} else if (_subCmd == "merge") {
		tool = new MergeFile(static_cast<ContextMerge *>(context));
	} else if (_subCmd == "jaccard") {
		tool = new Jaccard(static_cast<ContextJaccard *>(context));
	} else if (_subCmd == "subtract") {
		tool = new SubtractFile(static_cast<ContextSubtract *>(context));
	} else if (_subCmd == "sample") {
		tool = new SampleFile(static_cast<ContextSample *>(context));
	} else if (_subCmd == "spacing") {
		tool = new SpacingFile(static_cast<ContextSpacing *>(context));
	} else if (_subCmd == "fisher") {
		tool = new Fisher(static_cast<ContextFisher *>(context));
	} else if (_subCmd == "coverage") {
		tool = new CoverageFile(static_cast<ContextCoverage *>(context));
	} else if (_subCmd == "complement") {
		tool = new ComplementFile(static_cast<ContextComplement *>(context));
	} else if (_subCmd == "groupby") {
		tool = new GroupBy(static_cast<ContextGroupBy *>(context));
	} else if (_subCmd == "summary") {
		tool = new SummaryFile(static_cast<ContextSummary *>(context));
	}

	else {
		cerr << "Error: Tool " << _subCmd << " is not supported. Exiting..." << endl;
		exit(1);
	}
	return tool;
}

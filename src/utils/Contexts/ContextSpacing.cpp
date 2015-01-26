/*
 * ContextSpacing.cpp
 *
 *  Created on: Jan 15, 2015
 *      Author: Aaron Quinlan
 */
#include "ContextSpacing.h"

ContextSpacing::ContextSpacing()
{
	setSortedInput(true);
}

ContextSpacing::~ContextSpacing()
{

}

bool ContextSpacing::parseCmdArgs(int argc, char **argv, int skipFirstArgs) {
	_argc = argc;
	_argv = argv;
	_skipFirstArgs = skipFirstArgs;
	if (_argc < 2) {
		setShowHelp(true);
		return false;
	}

	setProgram(_programNames[argv[0]]);

	_argsProcessed.resize(_argc - _skipFirstArgs, false);

	for (_i=_skipFirstArgs; _i < argc; _i++) {
		if (isUsed(_i - _skipFirstArgs)) {
			continue;
		}
	}
	return ContextBase::parseCmdArgs(argc, argv, 1);
}



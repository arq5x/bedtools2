/*
 * ContextGroupBy.h
 *
 *  Created on: Mar 26, 2014
 *      Author: nek3d
 */

#ifndef CONTEXTGROUPBY_H_
#define CONTEXTGROUPBY_H_

#include "ContextBase.h"

class ContextGroupBy: public ContextBase {
public:
	ContextGroupBy();
	~ContextGroupBy();
	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);
	virtual bool isValidState();

	bool printFullCols() const { return _printFullCols; }
	const string &getGroupCols() const { return _groupStr; }
	bool ignoreCase() const { return _ignoreCase; }
	const string &getDefaultHeader();

protected:
	string _groupStr;
	string _defaultHeader;
	bool _printFullCols;
	bool _ignoreCase;

	bool handle_g();
	bool handle_inheader();
	bool handle_outheader();
	bool handle_header();
	bool handle_full();
	bool handle_ignorecase();
};


#endif /* CONTEXTGROUPBY_H_ */

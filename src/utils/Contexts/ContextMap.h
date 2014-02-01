/*
 * ContextMap.h
 *
 *  Created on: Jan 7, 2014
 *      Author: arq5x
 */

#ifndef CONTEXTMAP_H_
#define CONTEXTMAP_H_

#include "ContextIntersect.h"

class ContextMap : public ContextIntersect {
public:
	ContextMap();
	virtual ~ContextMap();
	virtual bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

	int getColumn() const { return _column; }
	void setColumn(int column) { _column = column; }

	const QuickString & getColumnOperation() const { return _columnOperation; }
	void setColumnOperation(const QuickString & operation) { _columnOperation = operation; }

	const QuickString & getNullValue() const { return _nullValue; }
	void setNullValue(const QuickString & nullValue) { _nullValue = nullValue; }

    virtual bool hasIntersectMethods() const { return true; }

private:
	virtual bool handle_c();
	virtual bool handle_o();
	virtual bool handle_null();

};

#endif /* CONTEXTMAP_H_ */

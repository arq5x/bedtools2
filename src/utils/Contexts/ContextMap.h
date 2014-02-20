/*
 * ContextMap.h
 *
 *  Created on: Jan 7, 2014
 *      Author: arq5x
 */

#ifndef CONTEXTMAP_H_
#define CONTEXTMAP_H_

#include "ContextIntersect.h"
#include "KeyListOps.h"

class ContextMap : public ContextIntersect {
public:
	ContextMap();
	virtual ~ContextMap();
	virtual bool isValidState();

	virtual bool parseCmdArgs(int argc, char **argv, int skipFirstArgs);

	const QuickString &getColumns() const { return _columns; }
	void setColumns(const QuickString &columns) { _columns = columns; }

	const QuickString & getOperations() const { return _operations; }
	void setOperations(const QuickString & operation) { _operations = operation; }

	const QuickString & getNullValue() const { return _nullValue; }
	void setNullValue(const QuickString & nullValue) { _nullValue = nullValue; }

	const QuickString &getDelim() const { return _delimStr; }
    virtual bool hasIntersectMethods() const { return true; }

    typedef vector<pair<int, KeyListOps::OP_TYPES> > colOpsType;
	const colOpsType &getColOps() const { return _colOps; }

private:
    QuickString _operations;
    QuickString _columns;
    QuickString _nullValue;
    KeyListOps _keyListOps;
    colOpsType _colOps;
    QuickString _delimStr;

	virtual bool handle_c();
	virtual bool handle_o();
	virtual bool handle_null();
	virtual bool handle_delim();

};

#endif /* CONTEXTMAP_H_ */

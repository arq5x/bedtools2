#include "ContextBase.h"
#include "ToolBase.h"

class BedtoolsDriver {
public:
	BedtoolsDriver();
	bool subMain(int argc, char **argv);
	bool supports(const string &tool);
	ContextBase *getContext();
	ToolBase *getTool(ContextBase *context);
	bool hadError() const { return _hadError; }
	string getErrors() const { return _errors; }
protected:
	string _subCmd;
	typedef set<string> supportType;
	supportType _supported;
	bool _hadError;
	string _errors;
};

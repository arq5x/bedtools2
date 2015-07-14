#include "ContextBase.h"
#include "ToolBase.h"

class BedtoolsDriver {
public:
	BedtoolsDriver();
	bool subMain(int argc, char **argv);
	bool supports(const QuickString &tool);
	ContextBase *getContext();
	ToolBase *getTool(ContextBase *context);
	bool hadError() const { return _hadError; }
protected:
	QuickString _subCmd;
	typedef set<QuickString> supportType;
	supportType _supported;
	bool _hadError;

};

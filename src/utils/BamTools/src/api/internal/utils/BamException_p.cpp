// ***************************************************************************
// BamException_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a basic exception class for BamTools internals
// ***************************************************************************

#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;
using namespace std;

const string BamException::SEPARATOR = ": ";

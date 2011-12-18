// ***************************************************************************
// api_global.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides macros for exporting & importing BamTools API library symbols
// ***************************************************************************

#ifndef API_GLOBAL_H
#define API_GLOBAL_H

#include "shared/bamtools_global.h"

#ifdef BAMTOOLS_API_LIBRARY
#  define API_EXPORT BAMTOOLS_LIBRARY_EXPORT
#else
#  define API_EXPORT BAMTOOLS_LIBRARY_IMPORT
#endif

#endif // API_GLOBAL_H

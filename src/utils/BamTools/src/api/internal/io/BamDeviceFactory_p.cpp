// ***************************************************************************
// BamDeviceFactory_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 September 2011 (DB)
// ---------------------------------------------------------------------------
// Creates built-in concrete implementations of IBamIODevices
// ***************************************************************************

#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/io/BamFile_p.h"
#include "api/internal/io/BamFtp_p.h"
#include "api/internal/io/BamHttp_p.h"
#include "api/internal/io/BamPipe_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
using namespace std;

IBamIODevice* BamDeviceFactory::CreateDevice(const string& source) {

    // check for requested pipe
    if ( source == "-" || source == "stdin" || source == "stdout" )
        return new BamPipe;

    // check for HTTP prefix
    if ( source.find("http://") == 0 )
        return new BamHttp(source);

    // check for FTP prefix
    if ( source.find("ftp://") == 0 )
        return new BamFtp(source);

    // otherwise assume a "normal" file
    return new BamFile(source);
}

// ***************************************************************************
// BamDeviceFactory_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Creates built-in concrete implementations of IBamIODevices
// ***************************************************************************

#ifndef BAMDEVICEFACTORY_P_H
#define BAMDEVICEFACTORY_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/IBamIODevice.h"
#include <string>

namespace BamTools {
namespace Internal {

class BamDeviceFactory {
    public:
        static IBamIODevice* CreateDevice(const std::string& source);
};

} // namespace Internal
} // namespace BamTools

#endif // BAMDEVICEFACTORY_P_H

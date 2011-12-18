// ***************************************************************************
// NetUnix_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides common networking-related includes, etc. for all UNIX-like systems
// ***************************************************************************

#ifndef NETUNIX_P_H
#define NETUNIX_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#ifndef _WIN32 // <-- source files only include the proper Net*_p.h, but this is a double-check

#include <arpa/inet.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <netdb.h>
#include <unistd.h>

#ifndef   BT_SOCKLEN_T
#  define BT_SOCKLEN_T socklen_t
#endif

#endif // _WIN32
#endif // NETUNIX_P_H

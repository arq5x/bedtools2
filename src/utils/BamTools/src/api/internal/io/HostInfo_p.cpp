// ***************************************************************************
// HostInfo_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides DNS lookup functionality for hostname & its discovered addresses
// ***************************************************************************

#include "api/internal/io/HostInfo_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

// platorm-specifics
#ifdef _WIN32
#  include "api/internal/io/NetWin_p.h"
#else
#  include "api/internal/io/NetUnix_p.h"
#endif

// standard C++ includes
#include <cstdlib>
#include <cstring>
#include <set>
using namespace std;

// -------------------------
// HostInfo implementation
// -------------------------

HostInfo::HostInfo(void)
    : m_error(HostInfo::NoError)
{ }

HostInfo::HostInfo(const HostInfo& other)
    : m_hostName(other.m_hostName)
    , m_addresses(other.m_addresses)
    , m_error(other.m_error)
    , m_errorString(other.m_errorString)
{ }

HostInfo::~HostInfo(void) { }

vector<HostAddress> HostInfo::Addresses(void) const {
    return m_addresses;
}

HostInfo::ErrorType HostInfo::GetError(void) const {
    return m_error;
}

string HostInfo::GetErrorString(void) const {
    return m_errorString;
}

string HostInfo::HostName(void) const {
    return m_hostName;
}

void HostInfo::SetAddresses(const std::vector<HostAddress>& addresses) {
    m_addresses = addresses;
}

void HostInfo::SetError(const HostInfo::ErrorType error) {
    m_error = error;
}

void HostInfo::SetErrorString(const std::string& errorString) {
    m_errorString = errorString;
}

void HostInfo::SetHostName(const string& name) {
    m_hostName = name;
}

// ---------------------------------
// HostInfo::Lookup(host, port)
//  - the real "heavy-lifter" here
// ---------------------------------

HostInfo HostInfo::Lookup(const string& hostname, const string& port) {

    HostInfo result;
    result.SetHostName(hostname);
    set<HostAddress> uniqueAddresses;

#ifdef _WIN32
    WindowsSockInit init;
#endif

    HostAddress address;
    address.SetAddress(hostname);

    // if hostname is an IP string ('0.0.0.0' or IPv6 format)
    // do reverse lookup for host domain name
    //
    // TODO: might just remove this... not sure if proper 'hostname' from IP string is needed
    //
    //       so far, haven't been able to successfully fetch a domain name with reverse DNS
    //       getnameinfo() on test sites just returns original IP string. BUT this is likely a rare
    //       case that client code tries to use an IP string and the connection should work fine
    //       anyway. GetHostName() just won't quite show what I was hoping for. :(
    if ( address.HasIPAddress() ) {

        const uint16_t portNum = static_cast<uint16_t>( atoi(port.c_str()) );

        sockaddr_in  sa4;
        sockaddr_in6 sa6;
        sockaddr* sa = 0;
        BT_SOCKLEN_T saSize = 0;

        // IPv4
        if ( address.GetProtocol() == HostAddress::IPv4Protocol ) {
            sa = (sockaddr*)&sa4;
            saSize = sizeof(sa4);
            memset(&sa4, 0, sizeof(sa4));
            sa4.sin_family = AF_INET;
            sa4.sin_addr.s_addr = htonl(address.GetIPv4Address());
            sa4.sin_port = htons(portNum);
        }

        // IPv6
        else if ( address.GetProtocol() == HostAddress::IPv4Protocol ){
            sa = (sockaddr*)&sa6;
            saSize = sizeof(sa6);
            memset(&sa6, 0, sizeof(sa6));
            sa6.sin6_family = AF_INET6;
            memcpy(sa6.sin6_addr.s6_addr, address.GetIPv6Address().data, sizeof(sa6.sin6_addr.s6_addr));
            sa6.sin6_port = htons(portNum);
        }

        // unknown (should be unreachable)
        else BT_ASSERT_X(false, "HostInfo::Lookup: unknown network protocol");

        // lookup name for IP
        char hbuf[NI_MAXHOST];
        char serv[NI_MAXSERV];
        if ( sa && (getnameinfo(sa, saSize, hbuf, sizeof(hbuf), serv, sizeof(serv), 0) == 0) )
            result.SetHostName(string(hbuf));

        // if no domain name found, just use the original address's IP string
        if ( result.HostName().empty() )
            result.SetHostName(address.GetIPString());

        // store address in HostInfo
        uniqueAddresses.insert(address);
    }

    // otherwise, hostname is a domain name ('www.foo.bar')
    // do 'normal' lookup
    else {

        // setup address lookup 'hints'
        addrinfo hints;
        memset(&hints, 0, sizeof(hints));
        hints.ai_family   = AF_UNSPEC;   // allow either IPv4 or IPv6
        hints.ai_socktype = SOCK_STREAM; // for TCP
        hints.ai_protocol = IPPROTO_TCP;

        // fetch addresses for requested hostname/port
        addrinfo* res;
        int status = getaddrinfo(hostname.c_str(), port.c_str(), &hints, &res );

        // if everything OK
        if ( status == 0 ) {

            // iterate over all IP addresses found
            addrinfo* p = res;
            for ( ; p != NULL; p = p->ai_next ) {

                // IPv4
                if ( p->ai_family == AF_INET ) {
                    sockaddr_in* ipv4 = (sockaddr_in*)p->ai_addr;
                    HostAddress a( ntohl(ipv4->sin_addr.s_addr) );
                    uniqueAddresses.insert(a);
                }

                // IPv6
                else if ( p->ai_family == AF_INET6 ) {
                    sockaddr_in6* ipv6 = (sockaddr_in6*)p->ai_addr;
                    HostAddress a(ipv6->sin6_addr.s6_addr);
                    uniqueAddresses.insert(a);
                }
            }

            // if we iterated, but no addresses were stored
            if ( uniqueAddresses.empty() && (p == NULL) ) {
                result.SetError(HostInfo::UnknownError);
                result.SetErrorString("HostInfo: unknown address types found");
            }
        }

        // handle error cases
        else if (
#ifndef _WIN32
                     status == EAI_NONAME
                  || status == EAI_FAIL
#  ifdef EAI_NODATA
                  || status == EAI_NODATA  // officially deprecated, but just in case we happen to hit it
#  endif // EAI_NODATA

#else  // _WIN32
                     WSAGetLastError() == WSAHOST_NOT_FOUND
                  || WSAGetLastError() == WSANO_DATA
                  || WSAGetLastError() == WSANO_RECOVERY
#endif // _WIN32
                )
        {
            result.SetError(HostInfo::HostNotFound);
            result.SetErrorString("HostInfo: host not found");
        }
        else {
            result.SetError(HostInfo::UnknownError);
            result.SetErrorString("HostInfo: unknown error encountered");
        }

        // cleanup
        freeaddrinfo(res);
    }

    // store fetched addresses (converting set -> vector) in result & return
    result.SetAddresses( vector<HostAddress>(uniqueAddresses.begin(), uniqueAddresses.end()) );
    return result;
}

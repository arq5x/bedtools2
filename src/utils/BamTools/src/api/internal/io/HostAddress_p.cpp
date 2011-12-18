// ***************************************************************************
// HostAddress_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a generic IP address container
// ***************************************************************************

#include "api/internal/io/HostAddress_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <vector>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

namespace BamTools {
namespace Internal {

// split a string into fields, on delimiter character
static inline
vector<string> Split(const string& source, char delim) {
    stringstream ss(source);
    string field;
    vector<string> fields;
    while ( getline(ss, field, delim) )
        fields.push_back(field);
    return fields;
}

// return number of occurrences of @pattern in @source
static inline
uint8_t CountHits(const string& source, const string& pattern) {

    uint8_t count(0);
    size_t found = source.find(pattern);
    while ( found != string::npos ) {
        ++count;
        found = source.find(pattern, found+1);
    }
    return count;
}

static
bool ParseIp4(const string& address, uint32_t& maybeIp4 ) {

    // split IP address into string fields
    vector<string> addressFields = Split(address, '.');
    if ( addressFields.size() != 4 )
        return false;

    // convert each field to integer value
    uint32_t ipv4(0);
    for ( uint8_t i = 0; i < 4; ++i ) {

        const string& field = addressFields.at(i);
        const size_t fieldSize = field.size();
        for ( size_t j = 0; j < fieldSize; ++j ) {
            if ( !isdigit(field[j]) )
                return false;
        }

        int value = atoi( addressFields.at(i).c_str() );
        if ( value < 0 || value > 255 )
            return false;

        // append byte value
        ipv4 <<= 8;
        ipv4 += value;
    }

    // store 32-bit IP address & return success
    maybeIp4 = ipv4;
    return true;
}

static
bool ParseIp6(const string& address, uint8_t* maybeIp6 ) {

    string tmp = address;

    // look for '%' char (if found, lop off that part of address)
    // we're going to ignore any link-local zone index, for now at least
    const size_t percentFound = tmp.rfind('%');
    if ( percentFound != string::npos )
        tmp = tmp.substr(0, percentFound);

    // split IP address into string fields
    vector<string> fields = Split(tmp, ':');
    const uint8_t numFields = fields.size();
    if ( numFields < 3 || numFields > 8 )
        return false;

    // get number of '::' separators
    const uint8_t numColonColons = CountHits(tmp, "::");
    if ( numFields == 8 && numColonColons > 1 )
        return false;

    // check valid IPv6 'compression'
    // must be valid 'pure' IPv6 or mixed IPv4/6 notation
    const size_t dotFound = tmp.find('.');
    const bool isMixed = ( dotFound != string::npos );
    if ( numColonColons != 1 && (numFields < (isMixed ? 7 : 8)) )
        return false;

    // iterate over provided fields
    size_t index = 16;
    size_t fillCount = 9 - numFields;
    for ( int8_t i = numFields - 1; i >= 0; --i ) {
        if ( index == 0 )
            return false;
        const string& field = fields.at(i);

        // if field empty
        if ( field.empty() ) {

            // if last field empty
            if ( i == numFields - 1 ) {
                const string& previousField = fields.at(i-1);
                if ( previousField.empty() )
                    return false;
                maybeIp6[--index] = 0;
                maybeIp6[--index] = 0;
            }

            // if first field empty
            else if ( i == 0 ) {
                // make sure ':' isn't first character
                const string& nextField = fields.at(i+1);
                if ( nextField.empty() ) return false;
                maybeIp6[--index] = 0;
                maybeIp6[--index] = 0;
            }

            // fill in 'compressed' 0s
            else {
                for ( uint8_t j = 0; j < fillCount; ++j ) {
                    if ( index == 0 ) return false;
                    maybeIp6[--index] = 0;
                    maybeIp6[--index] = 0;
                }
            }
        }

        // field has data
        else {
            uint32_t value = static_cast<uint32_t>( strtoul(field.c_str(), 0, 16) );

            if ( value <= 0xffff ) {
                maybeIp6[--index] =  value       & 0xff;
                maybeIp6[--index] = (value >> 8) & 0xff;
            }

            // possible mixed IPv4/6 notation
            else {

                // mixed field must be last
                if ( i != numFields - 1 )
                    return false;

                // parse the IPv4 section
                uint32_t maybeIp4;
                if ( !ParseIp4(field, maybeIp4) )
                    return false;

                // store IPv4 fields in IPv6 container
                maybeIp6[--index] =  maybeIp4        & 0xff;
                maybeIp6[--index] = (maybeIp4 >> 8)  & 0xff;
                maybeIp6[--index] = (maybeIp4 >> 16) & 0xff;
                maybeIp6[--index] = (maybeIp4 >> 24) & 0xff;
                --fillCount;
            }
        }
    }

    // should have parsed OK, return success
    return true;
}

} // namespace Internal
} // namespace BamTools

// ----------------------------
// HostAddress implementation
// ----------------------------

HostAddress::HostAddress(void)
    : m_protocol(HostAddress::UnknownNetworkProtocol)
    , m_ip4Address(0)
    , m_hasIpAddress(true)
{ }

HostAddress::HostAddress(const uint32_t ip4Address)
    : m_protocol(HostAddress::UnknownNetworkProtocol)
    , m_ip4Address(0)
    , m_hasIpAddress(true)
{
    SetAddress(ip4Address);
}

HostAddress::HostAddress(const uint8_t* ip6Address)
    : m_protocol(HostAddress::UnknownNetworkProtocol)
    , m_ip4Address(0)
    , m_hasIpAddress(true)
{
    SetAddress(ip6Address);
}

HostAddress::HostAddress(const IPv6Address& ip6Address)
    : m_protocol(HostAddress::UnknownNetworkProtocol)
    , m_ip4Address(0)
    , m_hasIpAddress(true)
{
    SetAddress(ip6Address);
}

HostAddress::HostAddress(const std::string& address)
    : m_protocol(HostAddress::UnknownNetworkProtocol)
    , m_ip4Address(0)
{
    SetAddress(address);
}

HostAddress::HostAddress(const HostAddress& other)
    : m_protocol(other.m_protocol)
    , m_ip4Address(other.m_ip4Address)
    , m_ip6Address(other.m_ip6Address)
    , m_ipString(other.m_ipString)
    , m_hasIpAddress(other.m_hasIpAddress)
{ }

HostAddress::~HostAddress(void) { }

bool HostAddress::operator==(const HostAddress& other) const {

    // if self is IPv4
    if ( m_protocol == HostAddress::IPv4Protocol ) {
        return ( other.m_protocol == HostAddress::IPv4Protocol &&
                 m_ip4Address == other.m_ip4Address
               );
    }

    // if self is IPv6
    else if ( m_protocol == HostAddress::IPv6Protocol ) {
        return ( other.m_protocol == HostAddress::IPv6Protocol &&
                 memcmp(&m_ip6Address, &other.m_ip6Address, sizeof(IPv6Address)) == 0
               );
    }

    // otherwise compare protocols
    else return m_protocol == other.m_protocol;
}

bool HostAddress::operator<(const HostAddress& other) const {

    // if self is IPv4
    if ( m_protocol == HostAddress::IPv4Protocol ) {
        if ( other.m_protocol == HostAddress::IPv4Protocol )
            return m_ip4Address < m_ip4Address;
    }

    // if self is IPv6
    else if ( m_protocol == HostAddress::IPv6Protocol ) {
        if ( other.m_protocol == HostAddress::IPv6Protocol )
            return (memcmp(&m_ip6Address, &other.m_ip6Address, sizeof(IPv6Address)) < 0);
    }

    // otherwise compare protocol types
    return m_protocol < other.m_protocol;
}

void HostAddress::Clear(void) {

    m_protocol = HostAddress::UnknownNetworkProtocol;
    m_ip4Address = 0;
    memset(&m_ip6Address, 0, sizeof(IPv6Address));
    m_ipString.clear();

    // this may feel funny, but cleared IP (equivalent to '0.0.0.0') is technically valid
    // and that's not really what this flag is checking anyway
    //
    // this flag is false *iff* the string passed in is a 'plain-text' hostname (www.foo.bar)
    m_hasIpAddress = true;
}

bool HostAddress::HasIPAddress(void) const {
    return m_hasIpAddress;
}

bool HostAddress::IsNull(void) const {
    return m_protocol == HostAddress::UnknownNetworkProtocol;
}

uint32_t HostAddress::GetIPv4Address(void) const {
    return m_ip4Address;
}

IPv6Address HostAddress::GetIPv6Address(void) const {
    return m_ip6Address;
}

std::string HostAddress::GetIPString(void) const {

    stringstream ss("");

    // IPv4 format
    if ( m_protocol == HostAddress::IPv4Protocol ) {
        ss << ( (m_ip4Address>>24) & 0xff ) << '.'
           << ( (m_ip4Address>>16) & 0xff ) << '.'
           << ( (m_ip4Address>> 8) & 0xff ) << '.'
           << (  m_ip4Address      & 0xff );

    }

    // IPv6 format
    else if ( m_protocol == HostAddress::IPv6Protocol ) {
        for ( uint8_t i = 0; i < 8; ++i ) {
            if ( i != 0 )
                ss << ':';
                ss << hex << ( (uint16_t(m_ip6Address[2*i]) << 8) |
                               (uint16_t(m_ip6Address[2*i+1]))
                             );
        }
    }

    // return result (empty string if unknown protocol)
    return ss.str();
}

HostAddress::NetworkProtocol HostAddress::GetProtocol(void) const {
    return m_protocol;
}

bool HostAddress::ParseAddress(void) {

    // all IPv6 addresses should have a ':'
    string s = m_ipString;
    size_t found = s.find(':');
    if ( found != string::npos ) {
        // try parse IP6 address
        uint8_t maybeIp6[16];
        if ( ParseIp6(s, maybeIp6) ) {
            SetAddress(maybeIp6);
            m_protocol = HostAddress::IPv6Protocol;
            return true;
        }
    }

    // all IPv4 addresses should have a '.'
    found = s.find('.');
    if ( found != string::npos ) {
        uint32_t maybeIp4(0);
        if ( ParseIp4(s, maybeIp4) ) {
            SetAddress(maybeIp4);
            m_protocol = HostAddress::IPv4Protocol;
            return true;
        }
    }

    // else likely just a plain-text host name "www.foo.bar"
    // will need to look up IP address info later
    m_protocol = HostAddress::UnknownNetworkProtocol;
    return false;
}

void HostAddress::SetAddress(const uint32_t ip4Address) {
    m_ip4Address = ip4Address;
    m_protocol = HostAddress::IPv4Protocol;
    m_hasIpAddress = true;
}

void HostAddress::SetAddress(const uint8_t* ip6Address) {
    for ( uint8_t i = 0; i < 16; ++i )
        m_ip6Address[i] = ip6Address[i];
    m_protocol = HostAddress::IPv6Protocol;
    m_hasIpAddress = true;
}

void HostAddress::SetAddress(const IPv6Address& ip6Address) {
    m_ip6Address = ip6Address;
    m_ip4Address = 0;
    m_protocol = HostAddress::IPv6Protocol;
    m_hasIpAddress = true;
}

void HostAddress::SetAddress(const std::string& address) {
    m_ipString = address;
    m_hasIpAddress = ParseAddress();
}

// ***************************************************************************
// ByteArray_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a dynamic, variable-length byte buffer
// ***************************************************************************

#include "api/internal/io/ByteArray_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdlib>
#include <cstring>
using namespace std;

// --------------------------
// ByteArray implementation
// --------------------------

ByteArray::ByteArray(void)
    : m_data()
{ }

ByteArray::ByteArray(const string& value)
    : m_data(value.begin(), value.end())
{ }

ByteArray::ByteArray(const vector<char>& value)
    : m_data(value)
{ }

ByteArray::ByteArray(const char* value, size_t n) {
    const string s(value, n);
    m_data.assign(s.begin(), s.end());
}

ByteArray::ByteArray(const ByteArray& other)
    : m_data(other.m_data)
{ }

ByteArray::~ByteArray(void) { }

ByteArray& ByteArray::operator=(const ByteArray& other) {
    m_data = other.m_data;
    return *this;
}

void ByteArray::Clear(void) {
    m_data.clear();
}

const char* ByteArray::ConstData(void) const {
    return &m_data[0];
}

char* ByteArray::Data(void) {
    return &m_data[0];
}

const char& ByteArray::operator[](size_t i) const {
    return m_data[i];
}

char& ByteArray::operator[](size_t i) {
    return m_data[i];
}

size_t ByteArray::IndexOf(const char c, const size_t from, const size_t to) const {
    const size_t size = ( (to == 0 ) ? m_data.size() : to );
    for ( size_t i = from; i < size; ++i ) {
        if ( m_data.at(i) == c ) 
            return i;
    }
    return m_data.size();
}

ByteArray& ByteArray::Remove(size_t from, size_t n) {

    // if 'from' outside range, just return
    const size_t originalSize = m_data.size();
    if ( from >= originalSize )
        return *this;

    // if asked to clip from 'from' to end (or beyond), simply resize
    if ( from + n >= originalSize )
        Resize(from);

    // otherwise, shift data & resize
    else {
        memmove( &m_data[from], &m_data[from+n], (originalSize-from-n) );
        Resize(originalSize - n);
    }

    // return reference to modified byte array
    return *this;
}

void ByteArray::Resize(size_t n) {
    m_data.resize(n, 0);
}

size_t ByteArray::Size(void) const {
    return m_data.size();
}

void ByteArray::Squeeze(void) {
    vector<char> t(m_data);
    t.swap(m_data);
}

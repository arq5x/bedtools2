// ***************************************************************************
// ByteArray_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a dynamic, variable-length byte buffer
// ***************************************************************************

#ifndef BYTEARRAY_P_H
#define BYTEARRAY_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/api_global.h"
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

// provides a wrapper around a byte vector
class ByteArray {

    // ctors & dtor
    public:
        ByteArray(void);
        ByteArray(const std::string& value);
        ByteArray(const std::vector<char>& value);
        ByteArray(const char* value, size_t n);
        ByteArray(const ByteArray& other);
        ~ByteArray(void);

        ByteArray& operator=(const ByteArray& other);

    // ByteArray interface
    public:

        // data access
        const char* ConstData(void) const;
        char* Data(void);
        const char& operator[](size_t i) const;
        char& operator[](size_t i);

        // byte array manipulation
        void Clear(void);
        size_t IndexOf(const char c, const size_t from = 0, const size_t to = 0) const;
        ByteArray& Remove(size_t from, size_t n);
        void Resize(size_t n);
        size_t Size(void) const;
        void Squeeze(void);

    // data members
    private:
        std::vector<char> m_data;
};

} // namespace Internal
} // namespace BamTools

#endif // BYTEARRAY_P_H

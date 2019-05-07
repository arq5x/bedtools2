/*
 * BedtoolsTypes.h
 *
 *  Created on: Feb 13, 2013
 *      Author: nek3d
 */

#ifndef BEDTOOLSTYPES_H_
#define BEDTOOLSTYPES_H_

// Some old environments require this even though standard C++ doesn't
#define __STDC_FORMAT_MACROS

#include <stdint.h>
#include <inttypes.h>
#include <climits>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include "string.h"

//STL headers
#include <utility> //include pair and make_pair
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef int64_t CHRPOS;

#define PRId_CHRPOS PRId64


#endif /* BEDTOOLSTYPES_H_ */

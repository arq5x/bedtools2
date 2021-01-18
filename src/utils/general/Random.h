#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdint.h>

// Seeds the random number generator
void rand_set_seed(int seed);

// Returns a random integer in [0, limit)
uint64_t rand_range(uint64_t limit);

// Returns a random proportion in [0.0, 1.0]
double rand_proportion();

#endif

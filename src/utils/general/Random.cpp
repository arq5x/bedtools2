#include "Random.h"

#ifdef USE_RAND
#include <stdlib.h>

void rand_set_seed(int seed)
{
    srand(seed);
}

uint64_t rand_range(uint64_t limit)
{
    if (limit <= RAND_MAX) {
        int n, max = RAND_MAX - (RAND_MAX % limit);
        do {
            n = rand();
        } while (n >= max);

        return n % limit;
    }
    else {
        // We need to combine two consective calls to rand() because
        // RAND_MAX is likely only 2^31 (2147483648). Note that this
        // will not cover the range if RAND_MAX is smaller.
        uint64_t n, max = UINT64_MAX - (UINT64_MAX % limit);
        do {
            n = (((uint64_t) rand()) << 31) | rand();
        } while (n >= max);

        return n % limit;
    }
}

double rand_proportion()
{
    return (double) rand() / (double) RAND_MAX;
}

#else
#include <random>

static std::mt19937_64 mt_rand;

void rand_set_seed(int seed)
{
    mt_rand.seed(seed);
}

// Ideally we would just use std::uniform_int_distribution and
// std::generate_canonical here, but unlike std::mt19937_64 they
// do not have defined implementations so would produce different
// results on different platforms.

uint64_t rand_range(uint64_t limit)
{
    uint64_t n, max = mt_rand.max() - (mt_rand.max() % limit);
    do {
        n = mt_rand();
    } while (n >= max);

    return n % limit;
}

double rand_proportion()
{
    return (double) mt_rand() / (double) mt_rand.max();
}

#endif

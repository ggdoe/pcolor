#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef double real_t;

#define MALLOC(size) aligned_alloc(64, size)
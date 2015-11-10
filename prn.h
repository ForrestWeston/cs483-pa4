#ifndef __PNR_H__
#define __PNR_H__

#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

struct context {
	uint64_t A;
	uint64_t B;
	uint64_t P;
	uint64_t seed;
	uint64_t itt;
	uint64_t size;
	int rank;
	int numprocs;
};

uint64_t*
ParallelGen(
	struct context* ctx
);

#endif

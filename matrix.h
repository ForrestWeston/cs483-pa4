#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

struct matrix {
	uint64_t m[4];
};

struct vector {
	uint64_t v[2];
};

void
RecvMatrix(
	struct matrix* mrecv
);

void
SendMatrix(
	int sendto,
	struct matrix* m
);

struct
matrix* NewMatrix();

struct
vector* NewVector();

void
MultMatrices(
	struct matrix* x,
	struct matrix* y,
	struct matrix* nMatri
);

void
MultModMatrix(
	struct matrix* x,
	struct matrix* y,
	struct matrix* nMatrix,
	uint64_t mod
);

void MultVectorMatrix(
	struct vector *vec, 
	struct matrix* mat,
	struct vector* nvector
);

#endif

#include "matrix.h"

void RecvMatrix(struct matrix* mrecv)
{
	MPI_Status status;
	MPI_Recv(mrecv->m, 4, MPI_INT, MPI_ANY_SOURCE,
			MPI_ANY_TAG, MPI_COMM_WORLD, &status);
}

void SendMatrix(int sendto, struct matrix* m)
{
	MPI_Send(m->m, 4, MPI_INT, sendto, 0, MPI_COMM_WORLD);
}

struct matrix* NewMatrix()
{
	return malloc(sizeof(struct matrix));
}

struct vector* NewVector()
{
	return malloc(sizeof(struct vector));
}

void MultMatrices(struct matrix* x, struct matrix* y, struct matrix* nMatrix)
{
	nMatrix->m[0] = (x->m[0]*y->m[0]) + (x->m[1]*y->m[2]);
	nMatrix->m[1] = (x->m[0]*y->m[1]) + (x->m[1]*y->m[3]);
	nMatrix->m[2] = (x->m[2]*y->m[0]) + (x->m[3]*y->m[2]);
	nMatrix->m[3] = (x->m[2]*y->m[1]) + (x->m[3]*y->m[3]);
}

void MultModMatrix(struct matrix* x, struct matrix* y, struct matrix* nMatrix,
	uint64_t mod)
{
	nMatrix->m[0] = ((x->m[0]*y->m[0]) + (x->m[1]*y->m[2])) % mod;
	nMatrix->m[1] = ((x->m[0]*y->m[1]) + (x->m[1]*y->m[3])) % mod;
	nMatrix->m[2] = ((x->m[2]*y->m[0]) + (x->m[3]*y->m[2])) % mod;
	nMatrix->m[3] = ((x->m[2]*y->m[1]) + (x->m[3]*y->m[3])) % mod;
}

void MultVectorMatrix(struct vector *vec, struct matrix* mat,
	struct vector* nvector)
{
	nvector->v[0] = (vec->v[0] * mat->m[0]) + (vec->v[1] * mat->m[2]);
	nvector->v[1] = (vec->v[0] * mat->m[1]) + (vec->v[1] * mat->m[3]);
}

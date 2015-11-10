#include "matrix.h"
#include "pnr.h"

uint64_t* ParallelGen(struct context* ctx)
{
	uint64_t count = 1;
	struct matrix *mbase = NewMatrix();
	struct matrix **mlist = NULL;
	struct vector *vbase = NewVector();
	uint64_t* result =  NULL;

	// init base matrix for parallel prefix sum
	mbase->m[0] = ctx->A;
	mbase->m[1] = 0;
	mbase->m[2] = ctx->B;
	mbase->m[3] = 1;

	//init base vector for parallel prefix sum
	vbase->v[0] = ctx->seed;
	vbase->v[1] = 1;

	// create list to hold matrix entries
	mlist = malloc(sizeof(struct matrix*) * ctx->size);
	mlist[0] = NewMatrix();
	*mlist[0] = *mbase;

	// fill list with intial entries
	for (; count < ctx->size; count++) {
		mlist[count] = NewMatrix();
		MultModMatrix(mlist[count-1], mbase, mlist[count], ctx->P);
	}

	// matricies to share among processors
	struct matrix* recv = NewMatrix();
	struct matrix* send = NewMatrix();
	struct matrix* temp = NewMatrix();
	//we share the list element in our local list
	*send = *mlist[ctx->size - 1];

	//Communication step exchnage matricies hypercubic
	int till = (int)log2(ctx->numprocs);
	for (count = 1; count <= till; count++) {
		if (ctx->rank % (1 << count) < (1 << (count - 1))) {
			SendMatrix(ctx->rank + (1 << (count-1)), send);
			RecvMatrix(recv);
			MultModMatrix(send, recv, temp, ctx->P);
			*send = *temp;
		} else {
			RecvMatrix(temp);
			int to = ctx->rank - (1<<(count-1));
			SendMatrix(to, send);
			uint64_t i;
			for (i = 0; i < ctx->size; i++) {
				MultModMatrix(mlist[i], temp, recv, ctx->P);
				*mlist[i] = *recv;
			}
			MultModMatrix(send, temp, recv, ctx->P);
			*send = *recv;
		}
	}

	uint64_t *output = malloc(sizeof(uint64_t) * ctx->size);
	struct vector *vres = NewVector();
	vres->v[0] = 0;
	vres->v[1] = 0;
	for (count = 0; count < ctx->size; count++) {
		MultVectorMatrix(vbase, mlist[count], vres);
		output[count] = vres->v[0] % ctx->P;
	}

	if (ctx->rank == 0) {
		result = malloc(sizeof(uint64_t) * ctx->itt);
	}

	//MPI_Gather(output, ctx->size, MPI_INT, result, ctx->size,
	//MPI_INT, 0,MPI_COMM_WORLD);

	return result;
}

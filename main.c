#include "prn.h"
#include "matrix.h"

#define CONST_A 7
#define CONST_B 3
#define PRIME_P 1223
#define SEED	18
#define NUM_ITT 100

int mpiErr;

struct sort_context{
	int numNodes;
	int id;
	uint64_t* localSplitters;
	uint64_t* globalSplitters;
	uint64_t* myArray;
	uint64_t myArraySize;
};

struct rand_context{
	uint64_t A;
	uint64_t B;
	uint64_t P;
	uint64_t seed;
	uint64_t itt;
	uint64_t size;
	int rank;
	int numprocs;
};

void init_mpi(int argv, char ***argv, struct context *s_ctx);

void init_mpi(int argc, char ***argv, struct context *s_ctx)
{
	mpiErr = MPI_Init(argc, argv);
	mpiErr = MPI_Comm_size(MPI_COMM_WORLD, s_ctx->numNodes);
	mpiErr = MPI_Comm_rank(MPI_COMM_WORLD, &s_ctx->id);
}

int comp(const void *elem1, const void *elem2)
{
	int f = *((int*)elem1);
	int s = *((int*)elem2);
	if (f > s) return 1;
	if (f < s) return -1;
	return 0;
}

void localSort(struct sort_context *s_ctx)
{
	qsort(s_ctx->myArray, s_ctx->myArraySize, sizeof(*s_ctx->myArray), comp);
}

void chooseLocalPivots(struct sort_context *s_ctx)
{
	int numChoose = s_ctx->numNodes - 1;
	int i = s_ctx->myArraySize/numChoose;
	while (i < s_ctx->myArraySize){
		printf("i'm super smelly\n");
	}

}

int main(int argc, char *argv[])
{
	struct sort_context *s_ctx = malloc(sizeof(struct sort_context));
	init_mpi(argc, argv, s_ctx);

	int n = (int)log2(s_ctx->numNodes);
	if (1 << n != s_ctx->numNodes){
		printf("Number of procs must be a power of 2\n");
		exit(1);
	}

	struct rand_context r_ctx = {
		.A = CONST_A,
		.B = CONST_B,
		.P = PRIME_P,
		.seed = SEED,
		.itt = NUM_ITT,
		.size = (NUM_ITT/s_ctx->numNodes),
		.rank = s_ctx->id,
		.numprocs = s_ctx->numNodes
	};

	s_ctx->myArray = ParallelGen(&r_ctx);
	s_ctx->myArraySize = r_ctx->size;

	//sort local array
	localSort(s_ctx);

	//pick p-1 evenly spaced pivots from local sorted array
	chooseLocalPivots(s_ctx);
}




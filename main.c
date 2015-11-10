#include "prn.h"
#include "matrix.h"

#define CONST_A 7
#define CONST_B 3
#define PRIME_P 1223
#define SEED	18
#define NUM_ITT 10000000 //ten million

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
	uint64_t f = *((uint64_t*)elem1);
	uint64_t s = *((uint64_t*)elem2);
	if (f > s) return 1;
	if (f < s) return -1;
	return 0;
}

// sorts array inplace
void localSort(uint64_t *array, uint64_t size)
{
	qsort(array, size, sizeof(array), comp);
}

//chooses numPivots evenvly spaced pivots from array to out array
void choosePivots(uint64_t *array, uint64_t size, uint64_t numPivots, uint64_t *out)
{
	uint64_t i;
	uint64_t loc = size/numPivots;
	for(i=0;i<numPivots;i++) {
		out[i] = array[(i+1)*loc];
	}
}

//sorts s_ctx->localpivot arrays (locally) chooses pivots into s_ctx->globalPivots
void globalChoosePivots(struct sort_context *s_ctx)
{
	uint64_t *rbuf;
	uint64_t size = s_ctx->numNodes*(s_ctx->numNodes - 1);
	MPI_Comm comm;
	MPI_Comm_size(comm, s_ctx->numNodes);

	rbuf = malloc(sizeof(uint64_t) * size);
	MPI_Allgather(s_ctx->localSplitters, (s_ctx->numNodes-1), MPI_UINT64_T,
				  rbuf, (s_ctx->numNodes-1), MPI_UINT64_T, comm);

	localSort(rbuf,size);
	choosePivots(rbuf, size, s_ctx->myArraySize, s_ctx->globalSplitters);

}
int main(int argc, char *argv[])
{
	struct sort_context *s_ctx = malloc(sizeof(struct sort_context));
	s_ctx->localSplitters = malloc(sizeof(uint64_t)*s_ctx->numNodes-1);
	s_ctx->globalSplitters = malloc(sizeof(uint64_t)*s_ctx->numNodes-1);
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
	localSort(s_ctx->myArray, s_ctx->myArraySize);

	//pick p-1 evenly spaced pivots from local sorted array
	choosePivots(s_ctx->myArray, s_ctx->myArraySize,
					s_ctx->myArraySize-1, s_ctx->localSplitters);

	//Sort the p(p-1) pivots and choose p-1 evenly spaced
	//pivots for global pivots
	globalChoosePivots(s_ctx);
}




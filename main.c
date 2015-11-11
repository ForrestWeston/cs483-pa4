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
	uint64_t* localPivots;
	uint64_t* globalPivots;
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
	MPI_Allgather(s_ctx->localPivots, (s_ctx->numNodes-1), MPI_UINT64_T,
			rbuf, (s_ctx->numNodes-1), MPI_UINT64_T, comm);

	localSort(rbuf,size);
	choosePivots(rbuf, size, s_ctx->myArraySize, s_ctx->globalPivots);

}

void distributePartitions(struct sort_context *s_ctx)
{
	/* MPI_Alltoallv()
	 * sendBuf == s_ctx->myArray;
	 * sendcount[numNodes] [count to rank 0], [count to rank 1] ...
	 * senddispl[numNode] [0 to first pivt], [1st pivt to 2nd] ...
	 * sendtype == MPI_UINT64_T
	 *
	 * recvBuf == s_ctx->myArray
	 * recvcount[numNodes] [max possible count from rank 0] ....
	 * recvdispl[numNodes] [0 to recvcount[0]] [recvcount[0] to recvcount[1]]...
	 * recvtype == MPI_UINT64_T
	 */

	//calculate count and displacement for alltoallv
	uint64_t *sendcount = malloc(sizeof(uint64_t)*s_ctx->numNodes);
	uint64_t *senddispl = malloc(sizeof(uint64_t)*s_ctx->numNodes);
	uint64_t *recvcount = malloc(sizeof(uint64_t)*s_ctx->numNodes);
	uint64_t *recvdispl = malloc(sizeof(uint64_t)*s_ctx->numNodes);

	int count = 0, disp = 0, itt = 0, pivcount = 0;

	for (; itt < s_ctx->myArraySize; itt++) {
		if (s_ctx->myArray[itt] == s_ctx->globalPivots[pivcount]) {
			count++;
			sendcount[pivcount] = count;
			senddispl[pivcount] = disp;
			pivcount++;
			disp = count;
			count = 0;
			continue;

		}
		if (s_ctx->myArray[itt] > s_ctx->globalPivots[pivcount]) {
			sendcount[pivcount] = count = 0;
			senddispl[pivcount] = 0;
			pivcount++;
			continue;
		}
		count++;
	}
	//brodcast the count to all proccesses
	MPI_Alltoall(sendcount, s_ctx->numNodes, MPI_UINT64_T,
				 recvcount, s_ctx->numNodes, MPI_UINT64_T, MPI_COMM_WORLD);

	//derive recvdispl from recvcount
	recvdispl[0] = 0;
	for (itt=1; itt < s_ctx->numNodes; itt++) {
		recvdispl[itt] = recvcount[itt-1] + recvdispl[itt-1]; //this is neat
	}

	//redistribute arrays with alltoallv
	MPI_Alltoallv(MPI_IN_PLACE, sendcount, senddispl, MPI_UINT64_T,
				  s_ctx->myArray, recvcount, recvdispl, MPI_UINT64_T,
				  MPI_COMM_WORLD);

}


uint64_t* main(int argc, char *argv[])
{
	struct sort_context *s_ctx = malloc(sizeof(struct sort_context));
	s_ctx->localPivots = malloc(sizeof(uint64_t)*s_ctx->numNodes-1);
	s_ctx->globalPivots = malloc(sizeof(uint64_t)*s_ctx->numNodes-1);
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
					s_ctx->myArraySize-1, s_ctx->localPivots);

	//Sort the p(p-1) pivots and choose p-1 evenly spaced
	//pivots for global pivots
	globalChoosePivots(s_ctx);

	//use global pivots to divide up the data among all procs
	distributePartitions(s_ctx);

	//after recieveing the correct elements, sort them locally
	localSort(s_ctx->myArray, s_ctx->myArraySize);

	//every processor now has a sorted list s.t. rank0[arraysize] < rank1[arraysize].....
	//now all we need to do is combine them in order from lowest rank to highest

	if (s_ctx->id == 0) {
		uint64_t *result = malloc(sizeof(uint64_t)
							* s_ctx->myArraySize * s_ctx->numNodes);
	}

}




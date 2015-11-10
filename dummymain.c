#include "prn.h"
#include "matrix.h"
#include <sys/time.h>
#include <time.h>
#define BILLION 1000000000L

struct timer {
	struct timespec start;
	struct timespec end;
};

void start(struct timer *t) {
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t->start);
}

void stop(struct timer *t) {
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t->end);
}

struct timespec difftimer(struct timer *t) {
	struct timespec diff;

	if ( t->end.tv_nsec - t->start.tv_nsec < 0 ) {
		t->end.tv_sec  -= 1;
		t->end.tv_nsec += BILLION;
	}

	diff.tv_sec  = t->end.tv_sec - t->start.tv_sec;
	diff.tv_nsec = t->end.tv_nsec - t->start.tv_nsec;

	return diff;
}

int get_ms(struct timer *t) {
	struct timespec d;
	d = difftimer(t);
	return (d.tv_sec * 1000) + (d.tv_nsec / 1000000);
}

uint64_t* SerialGen(struct context* ctx)
{
	uint64_t index = 0;
	uint64_t* result = malloc(sizeof(uint64_t) *ctx->itt);
	result[index] = ((ctx->A * ctx->seed) + ctx->B) % ctx->P;
	index++;
	while (ctx->itt != index) {
		result[index] = ((ctx->A * result[index-1]) + ctx->B) % ctx->P;
		index++;
	}
	return result;
}

int main(int argc, char* argv[])
{
	int rank, p;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (argc != 6) {
		printf("ERROR: Incorrect number of args\n");
		printf("A = Constant\n");
		printf("B = Constant\n");
		printf("P = Prime Constant\n");
		printf("seed = seed for generation (x)\n");
		printf("itt = number of itterations\n");
		return -1;
	}

	struct context* ctx = malloc(sizeof(struct context));
	struct timer tParallel;
	struct timer tSerial;

	ctx->rank = rank;
	ctx->numprocs = p;

	ctx->A = atoi(argv[1]);
	ctx->B = atoi(argv[2]);
	ctx->P = atoi(argv[3]);
	ctx->seed = atoi(argv[4]);
	ctx->itt = atoi(argv[5]);
	ctx->size = ctx->itt/ctx->numprocs;
	//printf("ctx->size: %" PRIu64 "\n", ctx->size);

	if(rank == 0) {
		printf("Number procs: %d\n", ctx->numprocs);
	}

	int n = (int)log2(ctx->numprocs);
	if (1 << n != ctx->numprocs) {
		printf("number of procs must be power of 2\n");
		exit(1);
	}

	uint64_t* parallelresults = malloc(sizeof(uint64_t*)* ctx->itt);
	uint64_t* serialresults = malloc(sizeof(uint64_t*)* ctx->itt);


	if (rank == 0) {
		start(&tSerial);
		serialresults = SerialGen(ctx);
		stop(&tSerial);
	}

	start(&tParallel);
	parallelresults = ParallelGen(ctx);
	stop(&tParallel);

	if (rank == 0) {
		int pass = memcmp(parallelresults, serialresults, sizeof(parallelresults));
		printf("PASS: %d\n", pass);
		printf("Parallel Time: %dms\n", get_ms(&tParallel));
		printf("Serial Time: %dms\n",get_ms(&tSerial));
	}


	MPI_Finalize();
	return 0;
}

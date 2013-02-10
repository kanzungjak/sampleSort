#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "sort-util.c"
#include "sort.h"
#include "math.h"
#include "random.h"
#include "my_alltoallv.c"

int nSamples; //number of samples
int n; //number of values to be sorted
int nP, iP, c;

// randInit(key, n) generates random numbers
// intSort(item, size) sequential sort
//MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)

void sampleSort(int *vals, int *lSamples) {
	int *pivots, *samples, *sendcounts, *recv, *recvcounts, *rdispls;
	int **start;
	int i, j, err, size;
	nSamples = c * log(nP); //number of local samples
	sendcounts = malloc(nP * sizeof(int));
	recv = malloc(2 * nSamples * sizeof(int));  // This is very big
	//printf("nsamples %d\n", nSamples);

	initParallelRandomLEcuyer( time(NULL) * n, iP, nP); //I'mportant!
	randInit(vals, n); //generate random values
	
	if(!iP) {
		samples = malloc(nSamples * nP * sizeof(int));
	}
	pivots = malloc(nP * sizeof(int));

	for (i = 0; i < nSamples; i++) { //randomly select nSamples samples 
		lSamples[i] = vals[(int) (nextRandomLEcuyer() * n)]; 
		//printf("%d lSamples %d \n", iP, lSamples[i]);
	}

	//collect all samples to root
	MPI_Gather( lSamples, nSamples, MPI_INT, 
    		    samples, nSamples, MPI_INT, 
		    	0, MPI_COMM_WORLD);
	if(!iP) {
		for (i = 0; i < nSamples*nP; i++) { 
		//	printf("%d samples %d \n", iP, samples[i]);
		}
	}

	
	if(!iP) { //select pivots
		intSort(samples, nSamples * nP);
		for (i = nSamples, j = 1; i < nSamples * (nP - 1); i += nSamples) {
			pivots[j++] = samples[i]; 
		}
		pivots[0] = INT_MIN; //-infinity
		pivots[nP - 1] = INT_MAX; //+infinity
	}
	
	MPI_Bcast(pivots, nP , MPI_INT, 0, MPI_COMM_WORLD);
	for (i = 0; i < nP; i++) { 
			printf("%d pivots %d \n", iP, pivots[i]);
	}

	
	start = malloc(n * sizeof(int));
	partition(vals, n, pivots, start, nP);

	for (i = 0; i < nP; i++) { 
			printf("%d start %p | %d\n", iP, start[i], start[i][0]);
	}

	/*
	for (i = 0; i < nP; i++) {
		sendcounts[i] = start[i + 1] - start[i]; //displacements
	}

	for(i = 0; i < sizeof(sendcounts)/sizeof(int); i++)
		printf("%d: %d\n", iP,sendcounts[i]);

	recv = malloc(2 * nSamples * sizeof(int));
	recvcounts = malloc(nP * sizeof(int));
	rdispls = malloc(nP * sizeof(int));
	err = my_Alltoallv(vals, sendcounts, *start, MPI_INT, 
                     recv, 2 * nSamples, &size,
                     recvcounts, rdispls, MPI_INT,
                     MPI_COMM_WORLD);

	if(err) {
		puts("Error: Bad Sample!");
		printf("nSamples: %d \n", nSamples);
		MPI_Abort(MPI_COMM_WORLD, 99);
	}

	intSort(recv, recvcounts[iP]);

	free(samples);*/
}


int main(int argc, char** argv) {
	int *lSamples, *samples; //local samples, all samples (root)
	int *vals;

	int debug = 0;
	c = atoi(argv[2]); //constant value for nSamples computation
	n = atoi(argv[1]); //number of values to be sorted

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nP);
	MPI_Comm_rank(MPI_COMM_WORLD, &iP);

	if(argc != 3) {
		if(!iP) printf("Usage: sampleSort <number of values> <constant for samples>\n");
		MPI_Finalize();
		exit(0);
	}

	lSamples = malloc( nSamples * sizeof(int) ); 
	vals = malloc(n * sizeof(int));

	sampleSort(vals, lSamples);

	/*Debug*/
	if(debug && !iP) {
		printf("nSamples %d\n", nSamples);
	}
	
	/*cleanup*/
	free(lSamples);	
	free(vals);
	MPI_Finalize();
	return 0;
}

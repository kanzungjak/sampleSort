#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include "mpi.h"
#include "sort.h"
#include "math.h"
#include "random.h"


int nSamples; //number of samples
int n; //number of values to be sorted
int nP, iP, c;

// randInit(key, n) generates random numbers
// intSort(item, size) sequential sort
//MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)

void sampleSort(int *vals, int *lSamples, int **start) {
	int i, j, *pivots, *samples, *sendcounts, *recv, size, *recvcounts, *rdispls;
	nSamples = c * log(nP); //number of local samples
	sendcounts = malloc(nP * sizeof(int));
	recv = malloc(2 * nSamples * sizeof(int));  // This is very big
	
	randInit(vals, n); //generate random values
	
	if(!iP) {
		samples = malloc(nSamples * nP * sizeof(int));
	}
	pivots = malloc(nP * sizeof(int));

	for (i = 0; i < nSamples; i++) {
		lSamples[i] = vals[(int) (nextRandomLEcuyer() * nSamples)];
	}

	//collect all samples to root
	MPI_Gather( lSamples, nSamples, MPI_INT, 
    		    samples, nSamples, MPI_INT, 
		    0, MPI_COMM_WORLD);

	if(!iP) {
		intSort(samples, nSamples * nP);
		for (i = nSamples, j = 0; i < nSamples * nP; i += nSamples) {
			pivots[j++] = samples[i];
		}
	}
	MPI_Bcast(pivots, nP, MPI_INT, 0, MPI_COMM_WORLD);

	partition(vals, n, pivots, start, nP);

	for (i = 0; i < nP; i++) {
		sendcounts[i] = start[i + 1] - start[i];
	}
	my_Alltoallv(vals, sendcounts, start, MPI_INT, 
                     recv, 2 * nSamples, &size,
                     recvcounts, rdispls, MPI_INT,
                     MPI_COMM_WORLD);

	intSort(vals, n / nP);

	free(samples);
}


int main(int argc, char** argv) {
	int *lSamples, *samples; //local samples, all samples (root)
	int *vals, **res;

	int debug = 1;
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
	res = malloc(n * n * sizeof(int));

	sampleSort(vals, lSamples, res);

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

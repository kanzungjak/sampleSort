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

void sampleSort(int *vals) {
	int *pivots, *samples, *sdispls, *recv, *sendcounts, *recvcounts, *rdispls;
	int *lSamples;
	int **start;
	int i, j, err, size;
	nSamples = c * log(nP); //number of local samples
	lSamples = malloc(nSamples * sizeof(int)); 
	sdispls = malloc((nP - 1) * sizeof(int));
	recv = malloc(2 * nSamples * sizeof(int));  // This is very big
	sendcounts = malloc(nP * sizeof(int));
	
	//printf("nsamples %d\n", nSamples);

	initParallelRandomLEcuyer( time(NULL) * n, iP, nP); //I'mportant!
	randInit(vals, n); //generate random values
	
	if(!iP) {
		samples = malloc(nSamples * nP * sizeof(int));
	}
	pivots = malloc((nP + 1) * sizeof(int));

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
		for (i = nSamples, j = 1; i < nSamples * nP; i += nSamples) {
			pivots[j++] = samples[i]; 
		}
		pivots[0] = INT_MIN; //-infinity
		pivots[nP] = INT_MAX; //+infinity
	}
	
	MPI_Bcast(pivots, nP + 1, MPI_INT, 0, MPI_COMM_WORLD);

	start = malloc((nP + 1) * sizeof(int *)); //enough memory?!
	partition(vals, n, pivots, start, nP);

	for (i = 0; i < nP; i++) { 
	//	printf("%d start %p | %d\n", iP, start[i], start[i][0]);
	}

	
	for (i = 0; i < nP - 1; i++) {
		sdispls[i] = abs(start[i + 1] - start[i]); //displacements
	}
	
	for (i = 0; i < nP - 1; i++) { 
	//	printf("%d sdispls %d\n", iP, sdispls[i]);
	}
	for (i = 0; i < nP ; i++) {
		sendcounts[i] = (start[i + 1] - start[i]) / sizeof(int); //displacements
		printf("%d sCounts: %d.\n", iP, sendcounts[i]);
	}
	/*int piv = 1;
	int ctr = 0;
	for(i = 0; i < n; i++) {
		//printf("%d Vals: %d.\n", iP, vals[i]);
		if(vals[i] < pivots[piv]) {
			ctr++;
			if(i == (n - 1)) {
				printf("IF# I am %d for value %d (%d) with pivot %d.\n", iP, i, vals[i], piv);
				sendcounts[piv-1] = ctr;
			}
		} else if (piv <= nP){
			printf("ELSE# I am %d for value %d (%d) with pivot %d.\n", iP, i, vals[i], piv);
			sendcounts[piv-1] = ctr;
			piv++;
			ctr = 0;
		}
	}*/
	/*for (i = 0; i < nP; i++) {
		printf("%d Start: %d\n", iP, start[i][0]);
	}*/
	
	//printf("piv: %d\n", piv);
	recv = malloc(2 * nSamples * sizeof(int));
	recvcounts = malloc(nP * sizeof(int));
	rdispls = malloc(nP * sizeof(int));
	/*err = my_Alltoallv(vals, sendcounts, sdispls, MPI_INT, 
                     recv, 2 * nSamples, &size,
                     recvcounts, rdispls, MPI_INT,
                     MPI_COMM_WORLD);*/

	//intSort(recv, recvcounts[iP]);
	//free(lSamples);	
	//if(!iP) free(samples);
}


int main(int argc, char** argv) {
	int *samples; //local samples, all samples (root)
	int *vals;

	int debug = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nP);
	MPI_Comm_rank(MPI_COMM_WORLD, &iP);
	if(argc != 3) {
		if(!iP) printf("Usage: sampleSort <number of values> <constant for samples>\n");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	c = atoi(argv[2]); //constant value for nSamples computation
	n = atoi(argv[1]); //number of values to be sorted

	vals = malloc(n * sizeof(int));

	sampleSort(vals);

	/*Debug*/
	/*if(debug && !iP) {
		printf("nSamples %d\n", nSamples);
	}*/

	/*cleanup*/
	//free(vals);
	MPI_Barrier(MPI_COMM_WORLD);
		printf("myiP %d\n", iP);

	MPI_Finalize();
	return 0;
}
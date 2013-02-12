#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "sort-util.c"
#include "sort.h"
#include "math.h"
#include "random.h"
#include "my_alltoallv.c"

int n, c; // Number of values to be sorted (per PE), parameter for amount of samples
int nP, iP; // Amount of PEs, number of own PE
MPI_Comm comm; // Communicator for Message Passing

void printSorted(int *res, int size) {
	int lSorted, sorted;
	sorted = 0;
	lSorted = isGloballySorted(res, size);
	MPI_Reduce(&lSorted, &sorted, 1, MPI_INT, MPI_LAND, 0, comm);
	if(!iP && sorted) {
		printf("Sorted!\n");
	} else if(!iP) {
		printf("Unsorted!\n");
	}
}

int sampleSort(int *vals, int *res, int maxrcount) {
	/* Declaration */
	int size, nSamples, i, j, ctr, pos;
	int *samples, *lSamples, *pivots, *sendcounts, *recvcounts, *sdispls, *rdispls;
	int **start;

	/* Initialisation */
	size = maxrcount; // Variable for local array size, currently set to maxrcount
	nSamples = ceil(c * log(nP)); // Number of local samples
	if(!iP) {
		samples =
		     malloc(nSamples * nP * sizeof(int	));
	}
	lSamples =   malloc(nSamples 	  * sizeof(int	));
	pivots =     malloc((nP + 1) 	  * sizeof(int	));
	start =      malloc((nP + 1)  	  * sizeof(int *)); //Pointer!
	sendcounts = malloc( nP	     	  * sizeof(int	));
	recvcounts = malloc( nP	     	  * sizeof(int	));
	sdispls =    malloc((nP + 1) 	  * sizeof(int	));
	rdispls =    malloc((nP + 1) 	  * sizeof(int	));

	/* Error Handling */
	if(!iP && !samples) {
		perror("Failure: ");
	}
	if(!lSamples || !pivots || !sendcounts || !recvcounts || !sdispls || !rdispls || !start) {
		perror("Failure: ");
	}

	/* Local Sample Generation */
	initParallelRandomLEcuyer(time(NULL) * n, iP, nP);
	for (i = 0; i < nSamples; i++) { // Random Selection
		lSamples[i] = vals[(int) (nextRandomLEcuyer() * n)]; 
	}

	// Collect all samples in root
	MPI_Gather(lSamples, nSamples, MPI_INT, 
    		   samples,  nSamples, MPI_INT, 
		   0, MPI_COMM_WORLD);

	/* Pivot Element Generation */
	if(!iP) {
		intSort(samples, nSamples * nP);
		for (i = nSamples, j = 1; i < nSamples * nP; i += nSamples) {
			pivots[j++] = samples[i];
		}
		pivots[0] = INT_MIN; // - Infinity
		pivots[nP] = INT_MAX; // + Infinity
	}

	MPI_Bcast(pivots, nP + 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Partitioning */
	partition(vals, n, pivots, start, nP);

	/* Send Block Count & Displacement Calculation */
	ctr = 0;
	pos = 0;
	sdispls[0] = 0;
	for(i = 1; i <= nP; i++) {
		for(j = pos; j < n; j++) {
			if(vals[j] < pivots[i]) {
				ctr++;
				pos++;
				if(pos == n) {
					sendcounts[i-1] = ctr;
					sdispls[i] = &vals[j] - vals;
				}
			} else {
				sendcounts[i-1] = ctr;
				sdispls[i] = &vals[j] - vals;
				ctr = 0;
				break;
			}
		}
	}

	/* Block Reallocation to PEs */
	my_Alltoallv(vals, sendcounts, sdispls, MPI_INT, 
                     res, maxrcount, &size,
                     recvcounts, rdispls, MPI_INT,
                     comm);

	intSort(res, size);
	
	/*for(i = 0; i < nP; i++) {
		if(iP == i) {
			for(j = 0; j < n; j++){
				printf("%d ", vals[j]);
			}
			puts("");
			printf("pivots: ");
			for(j = 1; j < nP; j++){
				printf("%d ", pivots[j]);
			}
			puts("");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/

	/* Cleanup */
	if(!iP && !samples) {
		free(samples);
	}
	if(!lSamples) {
		free(lSamples);
	}
	if(!pivots) {
		free(pivots);
	}
	if(!start) {
		free(start);
	}
	if(!sendcounts) {
		free(sendcounts);
	}
	if(!recvcounts) {
		free(recvcounts);
	}
	if(!sdispls) {
		free(sdispls);
	}
	if(!rdispls) {
		free(rdispls);
	}

	return size;
}


int main(int argc, char** argv) {
	/* Declaration */
	int *vals, *res, size, limit;

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &nP);
	MPI_Comm_rank(comm, &iP);

	if(argc != 3) {
		if(!iP) puts("Usage: sampleSort <number of values> <constant for computation of sample amount>");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	c = atoi(argv[2]); // Constant value for computation of nSamples
	n = atoi(argv[1]); // Number of values to be sorted (per PE)
	limit = 2 * n;

	vals = malloc(n     * sizeof(int));
	res =  malloc(limit * sizeof(int));
	if(!vals || !res) {
		perror("Failure: ");
	}

	/* Initialisation With Random Values */
	initParallelRandomLEcuyer(time(NULL) * n, iP, nP);
	randInitBound(vals, n, 2 * n);
	// randInit(vals, n); // Generate random values

	/* Method Calls */
	size = sampleSort(vals, res, limit);
	printSorted(res, size);
	printItemsGlobally(res, size);

	/* Cleanup */
	if(!vals) {
		free(vals);
	}
	if(!res) {
		free(res);
	}

	MPI_Finalize();
	return 0;
}

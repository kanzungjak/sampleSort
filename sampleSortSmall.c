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
	int size, nSamples, i, j, ctr, pos;
	int *pivots, *samples, /* *lSamples, */ *sdispls, *sendcounts, *recvcounts, *rdispls;
	int **start;

	/* Initialisation */
	size = maxrcount; // Variable for local array size, currently set to maxrcount
	nSamples = 3; // Number of local samples
	if(!iP) {
		samples =
		     malloc(nSamples * nP * sizeof(int	));
	}
	//lSamples =   malloc(nSamples 	  * sizeof(int	));
	pivots =     malloc((nP + 1) 	  * sizeof(int	));
	start =      malloc((nP + 1)  	  * sizeof(int *));
	sendcounts = malloc(nP 	     	  * sizeof(int	));
	recvcounts = malloc(nP	     	  * sizeof(int	));
	sdispls =    malloc((nP + 1) 	  * sizeof(int	));
	rdispls =    malloc((nP + 1) 	  * sizeof(int	));
	int lSamples[3];

	/* Error Handling */
	if(!iP && !samples) {
		perror("Failure: ");
	}
	if(!lSamples || !pivots || !sendcounts || !recvcounts || !sdispls || !rdispls || !start) {
		perror("Failure: ");
	}

	/* Local Sample Generation */
	if(iP == 0) {
		lSamples[0] = 7;
		lSamples[1] = 13,
		lSamples[2] = 25;
	}
	if(iP == 1) {
		lSamples[0] = 6;
		lSamples[1] = 17,
		lSamples[2] = 10;
	}
	if(iP == 2) {
		lSamples[0] = 20;
		lSamples[1] = 18,
		lSamples[2] = 21;
	}

	/*initParallelRandomLEcuyer(time(NULL) * n, iP, nP);
	for (i = 0; i < nSamples; i++) { // Random Selection
		lSamples[i] = vals[(int) (nextRandomLEcuyer() * n)]; 
	}*/

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
					// sdispls[i] = &vals[j] - vals;
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
	int /* *vals,*/ *res, size, limit;
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &nP);
	MPI_Comm_rank(comm, &iP);

	if(argc != 1) {
		if(!iP) puts("Usage: sampleSort");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	if(nP != 3){
		if(!iP) puts("Wrong PE number");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	c = 1; // Constant value for computation of nSamples
	n = 9; // Number of values to be sorted (per PE)
	limit = 2 * n;

	//vals = malloc(n     * sizeof(int));
	res =  malloc(limit * sizeof(int));
	if(/*!vals ||*/ !res) {
		perror("Failure: ");
	}

	/* Initialisation With Book Values */
	int vals[9];
	if(iP == 0) {
		vals[0] = 19;
		vals[1] = 7;
		vals[2] = 12;
		vals[3] = 1;
		vals[4] = 9;
		vals[5] = 13;
		vals[6] = 25;
		vals[7] = 4;
		vals[8] = 2;
	}
	if(iP == 1) {
		vals[0] = 6;
		vals[1] = 30;
		vals[2] = 17;
		vals[3] = 13;
		vals[4] = 10;
		vals[5] = 11;
		vals[6] = 16;
		vals[7] = 27;
		vals[8] = 22;
	}
	if(iP == 2) {
		vals[0] = 3;
		vals[1] = 20;
		vals[2] = 14;
		vals[3] = 18;
		vals[4] = 5;
		vals[5] = 16;
		vals[6] = 15;
		vals[7] = 21;
		vals[8] = 8;
	}

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

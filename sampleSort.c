#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "sort-util.c"
#include "sort.h"
#include "math.h"
#include "random.h"
#include "my_alltoallv.c"

int n, c; // Number of values to be sorted (per PE), parameter for amount of samples
int nP, iP; // Amount of PEs, rank of own PE
MPI_Comm comm; // Communicator for Message Passing

int printGloballySorted(int *res, int size) {
	int lSorted, sorted;
	sorted = 0;
	lSorted = isGloballySorted(res, size);
	if(nP > 1) {
		MPI_Reduce(&lSorted, &sorted, 1, MPI_INT, MPI_LAND, 0, comm);
	} else {
		sorted = lSorted;
	}
	if(!iP && sorted) {
		printf("Sorted!\n");
	} else if(!iP) {
		printf("Unsorted!\n");
	}

	return sorted;
}

int sampleSort(int *vals, int *res, int maxrcount) {
	/* Declaration */
	int size, m, i, j;
	int *samples, *lSamples, *pivots, *sendcounts, *recvcounts, *sdispls, *rdispls;
	int **start;

	/* Initialisation */
	if(nP > 1) {
		size = maxrcount; // Variable for local array size, currently set to maxrcount
	} else {
		size = n;
	}

	if(nP > 1) {
		m = c * log(nP); // Number of local samples
		if(!iP) {
			samples = 
			     malloc(m * nP	  * sizeof(int	));
		}
		lSamples =   malloc(m	 	  * sizeof(int	));
		pivots =     malloc((nP + 1) 	  * sizeof(int	));
		start =      malloc((nP + 1)  	  * sizeof(int *)); // Pointers!
		sendcounts = malloc( nP	     	  * sizeof(int	));
		recvcounts = malloc( nP	     	  * sizeof(int	));
		sdispls =    malloc( nP	 	  * sizeof(int	));
		rdispls =    malloc( nP	 	  * sizeof(int	));	

		/* Error Handling */
		if(!iP && !samples) {
			perror("Failure: ");
		}
		if(!lSamples || !pivots || !sendcounts || !recvcounts || !sdispls || !rdispls || !start) {
			perror("Failure: ");
		}

		/* Local Sample selection */
		initParallelRandomLEcuyer(time(NULL) * n, iP, nP);
		for (i = 0; i < m; i++) { // Random Selection
			lSamples[i] = vals[(int) (nextRandomLEcuyer() * n)];
		}

		// Collect all samples in root
		MPI_Gather(lSamples, m, MPI_INT,
				samples,  m, MPI_INT,
				0, MPI_COMM_WORLD);

		/* Pivot Element selection */
		if(!iP) {
			intSort(samples, m * nP);
			for (j = 1, i = m; i < m * nP; i += m) {
				pivots[j++] = samples[i];
			}
			pivots[0] =  INT_MIN; // - Infinity
			pivots[nP] = INT_MAX; // + Infinity
		}

		MPI_Bcast(pivots, nP + 1, MPI_INT, 0, MPI_COMM_WORLD);

		/* Partitioning */
		partition(vals, n, pivots, start, nP);

		/* Send Block Count Calculation */
		for(i = 0; i < nP; i++) {
			sendcounts[i] = start[i+1] - start[i];
		}

		/* Send Block Displacement Calculation */
		sdispls[0] = 0;
		for(i = 1; i < nP; i++) {
			sdispls[i] = sendcounts[i-1] + sdispls[i-1];
		}

		/* Block Reallocation to PEs */
		my_Alltoallv(vals, sendcounts, sdispls, MPI_INT,
				res, maxrcount, &size,
				recvcounts, rdispls, MPI_INT,
				comm);
	} else {
		for(i = 0; i < n; i++) {
			res[i] = vals[i];
		}
	}

	/* Local Sorting */
	intSort(res, size);

	/* Cleanup */
	if(nP > 1) {
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
	}

	return size;
}


int main(int argc, char** argv) {
	/* Declaration */
	int *vals, *res, size, limit, sorted, it, i, draw;	
	double startTime, diffTime, avgTime;
	char pFile[20], npFile[20];
	FILE *pData, *npData;

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &nP);
	MPI_Comm_rank(comm, &iP);
	
	if(argc != 3) {
		if(!iP) {
			puts("Usage: sampleSort <number of values> <constant for computation of sample amount>");
		}
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	c = atoi(argv[2]); // Constant value for computation of m
	// For n = 2 ^ 18, 700 is a good value for c
	n = atoi(argv[1]); // Number of values to be sorted (per PE)


	it = 1; // Number of iterations to get "proper" measurement results
	draw = 0; // Should we draw measurement curves?

	limit = 1.1 * n; // Maximum amount of values for one PE

	vals = malloc(n     * sizeof(int));
	res =  malloc(limit * sizeof(int));
	if(!vals || !res) {
		perror("Failure: ");
	}

	/* Create data files for curve drawing */
	if(draw) {		
		sprintf(pFile, "p-%d.log", nP);
		sprintf(npFile, "np-%d.log", n * nP);
		pData = fopen(pFile, "a");
		npData = fopen(npFile, "a");
	}
	


	if(nP > 1) {
		/* Initialisation With Random Values */
		initParallelRandomLEcuyer(time(NULL), iP, nP);
	} else {
		initRandomLEcuyer(time(NULL));
	}

	avgTime = 0.0;
	for(i = 0; i < it; i++) {
		if(nP > 1) { //parallel case
			/* Initialisation With Random Values */
			randInit(vals, n);
			// randInitBound(vals, n, nP * n);

			/* Sorting And Measuring */
			startTime = MPI_Wtime();
			size = sampleSort(vals, res, limit);
			diffTime = MPI_Wtime() - startTime;
			// printItemsGlobally(res, size);
			sorted = printGloballySorted(res, size); // Sorted or Unsorted?
		} else { //sequential case
			//random init
			randInit(vals, n);

			startTime = MPI_Wtime();
			intSort(vals, n);
			diffTime = MPI_Wtime() - startTime;
			sorted = isSorted(vals, n);
			if(sorted)
				puts("Sorted!");
			else
				puts("Unsorted!");

		}

		if (!iP && it < 2) {
			printf("Time: %f seconds\n", diffTime);
		}
		avgTime += (double) (diffTime / it);
	}

	/* Create data for curve drawing */
	if (!iP && draw) {
		fprintf(pData,"%d %f\n", n * nP, avgTime);
		fprintf(npData,"%d %f\n", nP, avgTime);
		if(!sorted) {
			fprintf(pData,"%f for %d was incorrect.\n", avgTime, n * nP);
			fprintf(npData,"%f for %d was incorrect.\n", avgTime, nP);
		}
	}

	if (!iP && it > 1) {
		printf("Average Time: %f seconds\n", avgTime);
	}

	/* Cleanup */
	if(!vals) {
		free(vals);
	}
	if(!res) {
		free(res);
	}
/*
	if(!pData) {
		fclose(pData);
	}
	if(!npData) {
		fclose(npData);
	}
*/

	MPI_Finalize();
	return 0;
}

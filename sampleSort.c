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
	int i, j, size, sorted, lSorted, ctr, pos;
	int *pivots, *samples, *lSamples, *sdispls, *res, *sendcounts, *recvcounts, *rdispls;
	int **start;
	nSamples = ceil(c * log(nP)); //number of local samples
	
	lSamples = malloc(nSamples * sizeof(int));
	
	res = malloc(2 * n * sizeof(int));
	
	sendcounts = malloc(nP * sizeof(int));
	recvcounts = malloc(nP * sizeof(int));
	
	sdispls = malloc((nP + 1) * sizeof(int));
	rdispls = malloc((nP + 1) * sizeof(int));
	
	pivots =  malloc((nP + 1) * sizeof(int));
	start =   malloc((nP + 1) * sizeof(int *)); 
	
	if(!iP) {
		samples = malloc(nSamples * nP * sizeof(int));
	}

	initParallelRandomLEcuyer( time(NULL) * n, iP, nP); //I'mportant!
	randInitBound(vals, n, 2 * n);

//	randInit(vals, n); //generate random values
	
	for (i = 0; i < nSamples; i++) { //randomly select nSamples samples 
		lSamples[i] = vals[(int) (nextRandomLEcuyer() * n)]; 
	}

	//collect all samples to root
	MPI_Gather( lSamples, nSamples, MPI_INT, 
    		    samples, nSamples, MPI_INT, 
		    0, MPI_COMM_WORLD);
	
	if(!iP) { //select pivots
		intSort(samples, nSamples * nP);
		for (i = nSamples, j = 1; i < nSamples * nP; i += nSamples) {
			pivots[j++] = samples[i]; 
		}
		pivots[0] = INT_MIN; //-infinity
		pivots[nP] = INT_MAX; //+infinity
	}

	MPI_Bcast(pivots, nP + 1, MPI_INT, 0, MPI_COMM_WORLD);

	partition(vals, n, pivots, start, nP);

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
	
	my_Alltoallv(vals, sendcounts, sdispls, MPI_INT, 
                     res, 2 * n, &size,
                     recvcounts, rdispls, MPI_INT,
                     MPI_COMM_WORLD);

	intSort(res, size);
	
	for(i = 0; i < nP; i++) {
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
	}
	
	printItemsGlobally(res, size);
	lSorted = isGloballySorted(res, size);
	//int numSorted = 0;
	//lSorted = isSorted(res, size);
	//MPI_Reduce(&lSorted, &numSorted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&lSorted, &sorted, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	if(!iP) printf("Success? %d\n", sorted);
	//if(!iP) printf("From %d PEs %d are sorted.\n", nP, numSorted);
	if(!lSamples) free(lSamples);
	if(!iP && !samples) free(samples);
}


int main(int argc, char** argv) {
	int *vals;

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

	/*cleanup*/
	if(!vals) free(vals);

	MPI_Finalize();
	return 0;
}
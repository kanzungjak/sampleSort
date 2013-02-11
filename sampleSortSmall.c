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

void sampleSort() {
	int *pivots, *samples, *sdispls, *recv, *sendcounts, *recvcounts, *rdispls;
	int **start;
	int i, j, err, size;
	nSamples = 3; //number of local samples
	
	sdispls = malloc((nP) * sizeof(int));
	sendcounts = malloc(nP * sizeof(int));
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

	int lSamples[3];
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
	
	if(!iP) {
		samples = malloc(nSamples * nP * sizeof(int));
	}
	pivots = malloc((nP + 1) * sizeof(int));

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
	
	start = malloc((nP + 1) * sizeof(int *)); //enough memory?!
	partition(vals, n, pivots, start, nP);


	int ctr = 0;
	int pos = 0;
	sdispls = malloc(nP * sizeof(int));
	sdispls[0] = 0;
	for(i = 1; i <= nP; i++) {
		for(j = pos; j < n; j++) {
			if(vals[j] < pivots[i]) {
				ctr++;
				pos++;
				if(pos == n) {
					sendcounts[i-1] = ctr;
				}
			} else {
				sendcounts[i-1] = ctr;
				sdispls[i] = &vals[j] - vals;
				ctr = 0;
				break;
			}
		}
	}
	
	recv = malloc(2 * n * sizeof(int));
	recvcounts = malloc(nP * sizeof(int));
	rdispls = malloc(nP * sizeof(int));
	
	err = my_Alltoallv(vals, sendcounts, sdispls, MPI_INT, 
                     recv, 2 * n, &size,
                     recvcounts, rdispls, MPI_INT,
                     MPI_COMM_WORLD);

	intSort(recv, size);
	printItemsGlobally(recv, size);
	int lSorted = isGloballySorted(recv, size);
	int sorted;
	MPI_Reduce(&lSorted, &sorted, 1,	MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	if(!iP) printf("Success? %d\n", sorted );
	if(!lSamples) free(lSamples);	
	if(!iP && !samples) free(samples);
}


int main(int argc, char** argv) {
	int *samples; //local samples, all samples (root)
	int *vals;

	int debug = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nP);
	MPI_Comm_rank(MPI_COMM_WORLD, &iP);
	if(argc != 1) {
		if(!iP) printf("Usage: sampleSort \n");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	
	if(nP != 3){
		if(!iP) printf("Wrong PE number\n");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
		
	c = 1; //constant value for nSamples computation
	n = 9; //number of values to be sorted

	sampleSort();

	MPI_Finalize();
	return 0;
}
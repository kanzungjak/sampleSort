#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include "mpi.h"
#include "sort.h"
#include "math.h"


int nSamples; //number of samples
int n; //number of values to be sorted

// randInit(key, n) generates random numbers
// intSort(item, size) sequential sort
//MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)




int main(int argc, char** argv) {
	int iP, nP;
	int *lSamples, *samples; //local samples, all samples (root)
	int *data;

	int debug = 1;
	int c = atoi(argv[2]); //constant value for nSamples computation
	n = atoi(argv[1]); //number of values to be sorted

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nP);
	MPI_Comm_rank(MPI_COMM_WORLD, &iP);

	if(argc != 3) {
		if(!iP) printf("Usage: sampleSort <number of values> <constant for samples>\n");
		MPI_Finalize();
		exit(0);
	}

	nSamples = c * log(nP); //number of local samples
	lSamples = malloc( nSamples * sizeof(int) ); 
	data = malloc( n * sizeof(int));
	
//	randInit(data, n); //generate random values
	
	if(!iP)
		samples = malloc( nSamples * nP * sizeof(int) );
	
	//collect all samples to root
	MPI_Gather( lSamples, nSamples, MPI_INT, 
    		    samples, nSamples, MPI_INT, 
		    0, MPI_COMM_WORLD);


	/*Debug*/
	if(debug && !iP) {
		printf("nSamples %d\n", nSamples);
		
	}
	
	/*cleanup*/
	free(lSamples);
	free(samples);
	free(data);
	MPI_Finalize();
	return 0;
}

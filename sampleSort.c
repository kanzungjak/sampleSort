#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "sort.h"
#include "math.h"


int nSamples; //number of samples

// randInit(key, n) generates random numbers
// intSort(item, size) sequential sort
//MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)




int main(int argc, char** argv) {
	int iP, nP;
	int *lSamples, *samples; //local samples, all samples (root)
	int *data;
	int c = atoi(argv[1]); //constant value for nSamples computation
	
	MPI_Comm_size(MPI_COMM_WORLD, &nP);
	MPI_Comm_rank(MPI_COMM_WORLD, &iP);

	nSamples = c * log(nP); //number of local samples
	lSamples = malloc( nSamples * sizeof(int) ); 
	data = malloc(  * sizeof(int));
	//randInit(samples, nSamples);
	
	if(!iP)
		samples = malloc( nSamples * nP * sizeof(int) );
	
	//collect all samples to root
	MPI_Gather( lSamples, nSamples, MPI_INT, 
				samples, nSamples, MPI_INT, 
				0, MPI_COMM_WORLD)

	
	
	/*cleanup*/
	free(lSamples);
	free(samples);
	MPI_Finalize();
	return 1;
}
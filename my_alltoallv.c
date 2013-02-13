/* Function similar to MPI_Alltoallv except that 
 * the receiver does not need to know the exact
 * displacements and lengths beforehand.
 * Instead, only a maximum total length needs to
 * be given. The received data are concatenated.
 * Parameters are used as in MPI_ALLTOALLV except:
  IN   maxRCount   specifies the maximum total 
                   number of items (of type recvtype)
                   that would fit into recvbuf
  OUT  size        total number of received items (not bytes)
  OUT  recvcounts  
  OUT  rdispls     number of elements and displacements
                   for each process - as in MPI_Alltoallv - 
                   but they are an output now 
*/

int my_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, 
                 void *recvbuf, int maxrcount, int *size,
                 int *recvcounts, int *rdispls, MPI_Datatype recvtype,
                 MPI_Comm comm) {
	/* Declaration & Initialisation */
	int iP, nP, i;
	MPI_Comm_rank(comm, &iP);
	MPI_Comm_size(comm, &nP);

	/* Receive Block Count Calculation */
 	MPI_Alltoall(sendcounts, 1, sendtype,
                     recvcounts, 1, recvtype,
                     comm);

	/* Receive Block Displacement Calculation */
	rdispls[0] = 0;
	for(i = 1; i < nP; i++) {
		rdispls[i] = recvcounts[i-1] + rdispls[i-1];
	}

	/* Actual Block Size Calculation */
	*size = 0;
	for (i = 0; i < nP; i++) {
		*size += recvcounts[i];
	}

	/* Error Handling For Bad Sample */
	if(*size > maxrcount) {
		puts("Error: Bad Sample!");
		MPI_Abort(MPI_COMM_WORLD, 42);
	}

	/* Block Reallocation to PEs */
	MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, 
		      recvbuf, recvcounts, rdispls, recvtype,
		      comm);
	return 0;
}

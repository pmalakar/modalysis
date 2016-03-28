#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <error.h>
#include <spi/include/kernel/memory.h>
#include "modalysis.h"

Modalysis::Modalysis() {
	nlocal = 0;
	x = v = f = NULL;

	atevery = NULL;
	atsteps = NULL;
	afname = NULL;
}

Modalysis::~Modalysis() {}

/*
 * initialize nlocal, nglobal, ntimesteps 
 */
void Modalysis::init(int me, int nprocs, int postprocess, long long int totalAtoms, int timesteps) {
	
	comm = MPI_COMM_WORLD;

	myrank = me;
	nprocs = nprocs;
	postprocess = postprocess;
	nglobal = totalAtoms;
	ntimesteps = timesteps;

	nlocal = nglobal/nprocs;
	//manage number of atoms in last rank
	if (myrank == nprocs-1)
		nlocal = nglobal - (nlocal*(nprocs-1));

//#ifdef DEBUG
	printf ("%d: Init %d %d %d %d\n", myrank, nprocs, nglobal, nlocal, ntimesteps);

	return;

}

void Modalysis::setupPostprocess() {

	uint64_t heapavail;

	//positions for all timesteps
	x = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) 
	{
		x[n] = NULL;
		x[n] = (double *) malloc (3 * nlocal * sizeof(double));
		if (x[n] == NULL)
			printf("%d: Failed to allocate %d bytes to x[%d]\n", myrank, 3*nlocal*sizeof(double), n);
	}

	//Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heapavail);
	//printf("%d: after allocating x free %lld\n", myrank, heapavail);

	//velocities for all timesteps
	v = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) {
		v[n] = NULL;
		v[n] = (double *) malloc (3 * nlocal * sizeof(double));
		if (v[n] == NULL)
			printf("%d: Failed to allocate %d bytes to v[%d]\n", myrank, 3*nlocal*sizeof(double), n);
		//else
		//	printf("%d: Allocated %d bytes to v[%d] free %lld\n", myrank, 3*nlocal*sizeof(double), n, heapavail);
	}

//	f = (double **) malloc (ntimesteps * sizeof(double *));

	xoriginal = (double *) malloc (3 * nlocal * sizeof(double));
	voriginal = (double *) malloc (3 * nlocal * sizeof(double));

	vacf = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) 
		vacf[n] = (double *) malloc (4 * sizeof(double));
	msd = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) 
		msd[n] = (double *) malloc (4 * sizeof(double));

	return;

}

/*
 * readFile
 *
 * reads atom positions
 * reads atom velocities
 *
 */
 
void Modalysis::readFile() {

	char *posfile = "positions.txt";
	char *velfile = "velocities.txt";

	MPI_Offset offset, mpifo;
	MPI_Status status;
	MPI_File posfh, velfh;
	int rcount;

	MPI_Barrier(comm);

	MPI_Scan(&nlocal, &nPartialSum, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

  MPI_File_open(comm, posfile, MPI_MODE_RDONLY, MPI_INFO_NULL, &posfh);
	if (posfh == NULL) 
		perror("positions file open failed");
  MPI_File_open(comm, velfile, MPI_MODE_RDONLY, MPI_INFO_NULL, &velfh);
	if (velfh == NULL)
		perror("velocities file open failed");

	if (posfh == NULL || velfh == NULL)
		return;

	MPI_Barrier(comm);

	int numelem;
	for (int n = 0; n<ntimesteps ; n++) {

		offset = n*nglobal + nPartialSum - nlocal;
		mpifo = offset * 3 * sizeof(double);

		numelem = 3*nlocal;

		if (MPI_File_read_at_all(posfh, mpifo, x[n], numelem, MPI_DOUBLE, &status) != MPI_SUCCESS)
			perror("positions file read error");
		if(MPI_File_read_at_all(velfh, mpifo, v[n], numelem, MPI_DOUBLE, &status) != MPI_SUCCESS)
			perror("velocities file read error");
		if (MPI_Get_count (&status, MPI_DOUBLE, &rcount) != MPI_SUCCESS)
			perror("MPI get count error");

#ifdef DEBUG
		//for (int i=0; i<nlocal; i++) {
		for (int i=0; i<2; i++) {
			printf("%d x[%d][%d] = %lf\n", myrank, n, i*PAD+0, x[n][i*PAD+0]);
			printf("%d v[%d][%d] = %lf\n", myrank, n, i*PAD+0, v[n][i*PAD+0]);
		}
#endif
	}

  MPI_File_close(&posfh);
  MPI_File_close(&velfh);

}

void Modalysis::initAnalysis() {

	for (int i = 0; i<nlocal*3 ; i++) {
		xoriginal[i] = x[0][i];
		voriginal[i] = v[0][i];
	}

}

long long int Modalysis::getnlocal() {
	return nlocal;
}

long long int Modalysis::getnglobal() {
	return nglobal;
}



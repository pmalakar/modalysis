#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#ifndef NERSC_HOST
#include <spi/include/kernel/memory.h>
#endif
#include "modalysis.h"

void Modalysis::setupPostprocess() {

	uint64_t heapavail;

	//positions for all timesteps
	x = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) 
	{
		x[n] = NULL;
		x[n] = (double *) malloc (3 * nlocal * sizeof(double));
		if (x[n] == NULL)
			printf("%d: Failed to allocate %d bytes for x[%d]\n", myrank, 3*nlocal*sizeof(double), n);
	}

	//Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heapavail);
	//printf("%d: after allocating x free %lld\n", myrank, heapavail);

	//velocities for all timesteps
	v = (double **) malloc (ntimesteps * sizeof(double *));
	for (int n = 0; n<ntimesteps ; n++) {
		v[n] = NULL;
		v[n] = (double *) malloc (3 * nlocal * sizeof(double));
		if (v[n] == NULL)
			printf("%d: Failed to allocate %d bytes for v[%d]\n", myrank, 3*nlocal*sizeof(double), n);
		//else
		//	printf("%d: Allocated %d bytes to v[%d] free %lld\n", myrank, 3*nlocal*sizeof(double), n, heapavail);
	}

//	f = (double **) malloc (ntimesteps * sizeof(double *));

	return;

}

/*
 * readFile for postprocess
 *
 * reads atom positions
 * reads atom velocities
 *
 */
 
void Modalysis::readFile() {

	char *posfile = "pos.txt";
	char *velfile = "vel.txt";

	MPI_Offset offset, mpifo;
	MPI_Status status;
	MPI_File posfh, velfh;
	int rcount;

	if (comm == 0) {
		printf("\n%d: Comm null %d\n", myrank, comm);
		exit(1);		
	}

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
	//		printf("%d v[%d][%d] = %lf\n", myrank, n, i*PAD+0, v[n][i*PAD+0]);
		}
#endif
	}

  MPI_File_close(&posfh);
  MPI_File_close(&velfh);

	for (int i = 0; i<nlocal*3 ; i++) {
		xoriginal[i] = x[0][i];
		voriginal[i] = v[0][i];
	}
}

void Modalysis::vacf_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<ntimesteps ; n++) { 
		compute_vacf(n, v[n]);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
	if (myrank == 0) printf("%lld Time to compute vacf: %lf\n", nglobal, time);

}

void Modalysis::msd_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<ntimesteps ; n++) { 
		compute_msd(n, x[n]);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
	if (myrank == 0) printf("%lld Time to compute msd: %lf\n", nglobal, time);

}

void Modalysis::histo_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<ntimesteps ; n=n+1) { 
		compute_histo(n, x[n]);	
		compute_histo(n, v[n]);	
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
	if (myrank == 0) printf("%lld Time to compute histo: %lf\n", nglobal, time);

}

void Modalysis::fft_() {
	
	int tend;	
	if (ntimesteps>=100)
		tend = 100;
	else
		tend = ntimesteps;

	double time;
	double stime = MPI_Wtime();

	for (int i=0; i<nlocal ; i++) 
		compute_fft_1d(0, tend, i, x);

	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
	if (myrank == 0) printf("%lld Time to compute fft: %lf\n", nglobal, time);

}

void Modalysis::postprocessdata() {

	setupPostprocess();
	allocate_();
	readFile();

	vacf_();
	msd_();
	histo_();
	fft_();

}



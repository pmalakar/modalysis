#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <spi/include/kernel/memory.h>
#include <mpi.h>

#include "modalysis.h"

#define PAD 3

int me, nprocs;
long long int totalAtoms = 4000;
int timesteps = 10;

Modalysis modalysis;

void vacf_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<timesteps ; n++) { 
		modalysis.compute_vacf(n);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("Time taken to compute vacf: %lf\n", time);

}

void msd_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<timesteps ; n++) { 
		modalysis.compute_msd(n);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("Time taken to compute msd: %lf\n", time);

}

void histo_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<timesteps ; n=n+1) { 
		modalysis.compute_histo(n, POSITION);	
		modalysis.compute_histo(n, VELOCITY);	
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("Time taken to compute histo: %lf\n", time);

}

void execute() {

	modalysis.myrank = me;
	modalysis.nprocs = nprocs;

	//post-processing mode setup
	modalysis.init(totalAtoms, timesteps);
	modalysis.setup();
	modalysis.readFile();

	modalysis.initAnalysis();

	//vacf_();
	//msd_();
	histo_();

}

int main (int argc, char** argv) {

		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		if (argc == 2) 
			totalAtoms = atoi(argv[1]);

		execute();

		MPI_Finalize();
		return 0;

}


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <spi/include/kernel/memory.h>
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

#include "modalysis.h"

#define PAD 3

int me, nprocs;
long long int totalAtoms;
int timesteps = 100;

Modalysis modalysis;

void vacf_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<timesteps ; n++) { 
		modalysis.compute_vacf(n);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("%lld Time taken to compute vacf: %lf\n", modalysis.getnglobal(), time);

}

void msd_() {

	double time;
	double stime = MPI_Wtime();
	for (int n=0; n<timesteps ; n++) { 
		modalysis.compute_msd(n);
	}
	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("%lld Time taken to compute msd: %lf\n", modalysis.getnglobal(), time);

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
	if (me == 0) printf("%lld Time taken to compute histo: %lf\n", modalysis.getnglobal(), time);

}

void fft_() {
	
	//for (int n=0; n<1 ; n++) { 
	//	modalysis.compute_fft(n);
	//}
	int tend;	
	if (timesteps>=100)
		tend = 100;
	else
		tend = timesteps;

	double time;
	double stime = MPI_Wtime();

	for (int i=0; i<modalysis.getnlocal() ; i++) 
		modalysis.compute_fft_1d(0, tend);

	stime = MPI_Wtime() - stime;
	MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, modalysis.comm);
	if (me == 0) printf("%lld Time taken to compute fft: %lf\n", modalysis.getnglobal(), time);

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
	//histo_();
	fft_();

}

int main (int argc, char** argv) {

		MPI_Init(&argc, &argv);
		//fftw_mpi_init();
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		if (argc == 2) 
			totalAtoms = atoi(argv[1]);

		execute();

		//fftw_mpi_cleanup();
		MPI_Finalize();
		return 0;

}


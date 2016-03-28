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
int timesteps = 20;
int postprocess = 0;

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

	if (postprocess) {

		modalysis.setupPostprocess();
		modalysis.readFile();
		modalysis.initAnalysis();

		vacf_();
		msd_();
		histo_();
		fft_();

	}

	else {		//coanalysis

	}

}


int main (int argc, char** argv) {

		char *analysiscfg = NULL;
		postprocess = 0;

		MPI_Init(&argc, &argv);
		//fftw_mpi_init();
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		for(int i = 0; i < argc; i++) {

     if((strcmp(argv[i], "-a") == 0) || (strcmp(argv[i], "--num_atoms") == 0)) {
       totalAtoms = atoi(argv[++i]);
       continue;
     }   

     if((strcmp(argv[i], "-n") == 0) || (strcmp(argv[i], "--num_steps") == 0)) {
       timesteps = atoi(argv[++i]);
       continue;
     }   

     if((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--postprocess") == 0)) {
       postprocess = atoi(argv[++i]);
       continue;
     }

		 if((strcmp(argv[i], "-acfg") == 0) || (strcmp(argv[i], "--analysis_config_file") == 0)) {
      if (analysiscfg == NULL) analysiscfg = new char[256];
      if (analysiscfg != NULL) strcpy(analysiscfg, argv[++i]);
      else printf("whats wrong!\n");
      continue;
     }

		}

		modalysis.init(me, nprocs, postprocess, totalAtoms, timesteps);

		printf("%s postprocess %d\n", analysiscfg, postprocess);

		if (postprocess == 0) 
			modalysis.readConfig(analysiscfg);

		execute();

		//fftw_mpi_cleanup();
		MPI_Finalize();
		return 0;

}



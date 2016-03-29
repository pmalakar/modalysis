/*
 * Developed at Argonne National Laboratory
 *
 * Contact: pmalakar@anl.gov, malakar.preeti@gmail.com
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <spi/include/kernel/memory.h>
#include <mpi.h>
//#include <fftw3.h>
//#include <fftw3-mpi.h>

#include "modalysis.h"

#define PAD 3

int me, nprocs;
long long int totalAtoms;
int timesteps = 10;
int postprocess = 0;

//Modalysis modalysis;

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

		Modalysis modalysis;

		modalysis.init(me, nprocs, postprocess, totalAtoms, timesteps);

		if (postprocess == 0) 
			modalysis.coanalyze(analysiscfg);
		else 
			modalysis.postprocessdata();

		//fftw_mpi_cleanup();
		MPI_Finalize();
		return 0;

}



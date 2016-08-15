#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include "modalysis.h"

#define numbins 10

//void Modalysis::compute_histo(int tstart, int tend, int attribute, double **arr) {
void Modalysis::compute_histo(int ts, double *arr) {

	double *vector;
	double value, minval, maxval;

	double MIN= 999;
	double MAX=-999;
	
/*
	switch(attribute) {
		case 1:
			vector = x[aindex][ts]; break;
		case 2:
			vector = v[aindex][ts]; break;
		case 3:
			vector = f[aindex][ts]; break;
		default:
			break;
	}
*/

#ifdef DEBUG
	printf("%d %d %d %d\n", myrank, ts, attribute, nlocal); 
#endif

	//MIN = MAX = vector[0];
	for (int i = 0; i < nlocal; i++) {
		value = arr[i*PAD+0];	//x-direction
		if (value <= MIN) MIN = value;
		else MAX = value;
  }

	MPI_Allreduce(&MIN, &minval, 1, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(&MAX, &maxval, 1, MPI_DOUBLE, MPI_MAX, comm);

	float binsize = (MAX-MIN)/numbins;
	double **bin;	//[numbins][2];
	bin = (double **) malloc (sizeof(double *) * numbins);	
	for (int i = 0; i < numbins; i++) {
		bin[i] = (double *) malloc (sizeof(double) * 3);	
		bin[i][2] = 0.0;
	}

	value = MIN;
	for (int i = 0; i < numbins; i++) {
		bin[i][0] = value;
		value += binsize;
		bin[i][1] = value;
	}

	for (int i = 0; i < nlocal; i++) {
		value = arr[i*PAD+0];	//x-direction
				
		for (int i = 0; i < numbins; i++) {
			if (value >= bin[i][0] && value <= bin[i][1]) bin[i][2] += 1;
		}
	
	}

#ifdef DEBUG
	if (myrank == 0) {
		for (int i = 0; i < numbins; i++) 
			printf("bin-%d Min %lf Max %lf count %lf\n", i, bin[i][0], bin[i][1], bin[i][2]);

		printf("%d gMin %lf gMax %lf\n", ts, minval, maxval);
	}
#endif

	return;

}


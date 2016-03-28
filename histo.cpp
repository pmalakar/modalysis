#include <mpi.h>
#include "modalysis.h"

void Modalysis::compute_histo(int ts, int attribute) {

	double *vector;
	double value, minval, maxval;

	double MIN= 999;
	double MAX=-999;
	
	switch(attribute) {
		case 1:
			vector = x[ts]; break;
		case 2:
			vector = v[ts]; break;
		case 3:
			vector = f[ts]; break;
		default:
			break;
	}

#ifdef DEBUG
	printf("%d %d %d %d\n", myrank, ts, attribute, nlocal); 
#endif

	//MIN = MAX = vector[0];
	for (int i = 0; i < nlocal; i++) {
		value = vector[i*PAD+0];	//x-direction
		if (value <= MIN) MIN = value;
		else MAX = value;
  }

	MPI_Allreduce(&MIN, &minval, 1, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(&MAX, &maxval, 1, MPI_DOUBLE, MPI_MAX, comm);

	if (myrank == 0) 
		printf("%d %d: Min %lf Max %lf\n", ts, attribute, minval, maxval);
	
	return;

}


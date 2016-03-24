#include "modalysis.h"

void Modalysis::compute_histo(int ts, int attribute) {

	double *vector;
	double value, minval, maxval;

	double MIN;
	double MAX;
	
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

	printf("%d %d %d %d\n", myrank, ts, attribute, nlocal);
 
	MIN = MAX = vector[0];
	for (int i = 0; i < nlocal; i++) {
		value = vector[i*PAD+0];
		if (value < MIN)
			MIN = value;
		else if (value > MAX)
			MAX = value;
  }

	MPI_Allreduce(&value, &minval, 1, MPI_DOUBLE, MPI_MIN, comm);
	MPI_Allreduce(&value, &maxval, 1, MPI_DOUBLE, MPI_MAX, comm);

	if (myrank == 0) printf("%d %d: Min %lf Max %lf\n", ts, attribute, minval, maxval);
	
	return;

}


#include <mpi.h>
#include <stdio.h>
#include "modalysis.h"

void Modalysis::compute_vacf(int ts, double *arr) {

	double vxsq, vysq, vzsq;
	double vector[4] = {0};

	for (int i = 0; i < nlocal; i++) {

		vxsq = arr[i*PAD+0] * voriginal[i*PAD+0];
		vysq = arr[i*PAD+1] * voriginal[i*PAD+1];
		vzsq = arr[i*PAD+2] * voriginal[i*PAD+2];
		vector[0] += vxsq;
		vector[1] += vysq;
		vector[2] += vzsq;
		vector[3] += vxsq + vysq + vzsq;

	}

	MPI_Allreduce(vector, vacf[ts], 4, MPI_DOUBLE, MPI_SUM, comm);

	int nvacf = nglobal;
  if (nvacf) {
    vacf[ts][0] /= nvacf;
    vacf[ts][1] /= nvacf;
    vacf[ts][2] /= nvacf;
    vacf[ts][3] /= nvacf;
  }

#ifdef DEBUG
	if (myrank == 0)
		printf("vacf %d: %lf %lf %lf %lf | %lf %lf %lf %lf\n", ts, vector[0], vector[1], vector[2], vector[3], vacf[0], vacf[1], vacf[2], vacf[3]);
#endif

	return;

}


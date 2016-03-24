#include "modalysis.h"

void Modalysis::compute_msd(int ts) {

	double dx, dy, dz;

	double vector[4] = {0};

	for (int i = 0; i < nlocal; i++) {

	    //dx = x[ts][i*PAD+0] - cm[0] - xoriginal[i*PAD+0];
	    dx = x[ts][i*PAD+0] - xoriginal[i*PAD+0];
      dy = x[ts][i*PAD+1] - xoriginal[i*PAD+1];
      dz = x[ts][i*PAD+2] - xoriginal[i*PAD+2];

      vector[0] += dx*dx;
      vector[1] += dy*dy;
      vector[2] += dz*dz;
      vector[3] += dx*dx + dy*dy + dz*dz;

			MPI_Allreduce(vector, msd[ts], 4, MPI_DOUBLE, MPI_SUM, comm);

			//nmsd - #atoms in a group 
			int nmsd = nlocal;

		  if (nmsd) {
				msd[ts][0] /= nmsd;
				msd[ts][1] /= nmsd;
				msd[ts][2] /= nmsd;
				msd[ts][3] /= nmsd;
      }

  }

	return;

}


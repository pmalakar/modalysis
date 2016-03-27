#include <fftw3.h>
#include <fftw3-mpi.h>
#include "modalysis.h"

void Modalysis::compute_fft(int ts) {

	fftw_plan plan;
	double *data = new double[3*nlocal];
	fftw_complex *out; 
	ptrdiff_t alloc_local, local_nx, local_x_start;
	ptrdiff_t i, j, k;

	ptrdiff_t dim = nglobal ; //nlocal;	//dimensions must be of ptrdiff_t in fftw_mpi* calls

	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_3d(dim, dim, dim/2+1, comm, &local_nx, &local_x_start);
	if (myrank < 10) printf("%d: dim=%i local_nx=%i local_x_start=%i alloc_local=%i\n", myrank, dim, local_nx, local_x_start, alloc_local);

	data = fftw_alloc_real(2 * alloc_local);
  out = fftw_alloc_complex(alloc_local);
	//out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);

	/* create plan for forward DFT */
	plan = fftw_mpi_plan_dft_r2c_3d(dim, dim, dim, data, out, comm, FFTW_MEASURE);

	j = -1;
#ifdef DEBUG
	double t = MPI_Wtime();
#endif
/*
	for (int i = 0; i < nlocal; i++) { 
		data[++j] = x[ts][i*PAD+0];
		data[j+nlocal] = x[ts][i*PAD+1];
		data[j+2*nlocal] = x[ts][i*PAD+2];
	}
*/
 	 if(myrank == 0)
   for (i = 0; i < local_nx; ++i)
   	for (j = 0; j < dim; ++j)
   		for (k = 0; k < dim; ++k)
   			printf("dim %d %d %d %d %d\n", i, j, k, (i*dim + j) * (2*(dim/2+1)) + k, (i*dim + j) * dim + k);
   			//data[index] = x[ts][];
   			////rin[(i*M + j) * (2*(N/2+1)) + k] = my_func(local_0_start+i, j, k);

#ifdef DEBUG
	t = MPI_Wtime() - t;
	if (myrank < 3) printf("%d: FFT: array copy time %lf\n", myrank, t);
#endif

	/* compute transforms */
	fftw_execute(plan);


	fftw_destroy_plan(plan);
	fftw_free(data);


	return;

}


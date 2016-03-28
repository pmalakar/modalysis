#ifndef MODALYSIS_H_
#define MODALYSIS_H_

#include <stdio.h>
#include <iostream>
#include <mpi.h>

#define PAD 3
#define FILENAMELEN 64

enum attributes {POSITION=1, VELOCITY, FORCE};

class Modalysis
{
	private:
		
		double **x;
		double **v;
		double **f;

		double *xoriginal;
		double *voriginal;

		long long int nlocal;
		long long int nglobal;
		long long int nPartialSum;

		int ntimesteps;

		int postprocess;

		double **vacf;
		double **msd;

	public:

		int myrank;
		int nprocs;

		MPI_Comm comm;

		int anum;
		int *atevery;
		int *atsteps;
		char *afname;

		Modalysis();
		~Modalysis();

		void init(int, int, int, long long int, int);

		void setupPostprocess();
		void readFile();
		void initAnalysis();

		void aalloc(int);
		
		void readConfig(char *);

		long long int getnlocal();
		long long int getnglobal();

		void compute_vacf(int);
		void compute_msd(int);
		void compute_histo(int, int);
		void compute_fft(int);
		void compute_fft_1d(int, int);

};

#endif


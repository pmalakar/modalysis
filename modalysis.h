#ifndef MODALYSIS_H_
#define MODALYSIS_H_

#include <stdio.h>
#include <iostream>
#include <mpi.h>

#define PAD 3

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

		int preprocess;

		double **vacf;
		double **msd;

	public:

		int myrank;
		int nprocs;

		MPI_Comm comm;

		Modalysis();
		~Modalysis();

		void init(long long int, int);
		void setup();
		void readFile();
		void initAnalysis();

		void compute_vacf(int);
		void compute_msd(int);
		void compute_histo(int, int);

};

#endif

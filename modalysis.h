/*
 * Developed at Argonne National Laboratory
 *
 * Contact: pmalakar@anl.gov, malakar.preeti@gmail.com
 *
 */

#ifndef MODALYSIS_H_
#define MODALYSIS_H_

#include <stdio.h>
#include <iostream>
#include <mpi.h>

#define PAD 3
#define ANAMELEN 16
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

		int myrank;
		int nprocs;
		int ntimesteps;

		long long int nlocal;
		long long int nglobal;
		long long int nPartialSum;

		int postprocess;

		double **vacf;
		double **msd;

		char *configFile;

		int anum;
		int *adim;
		int *atevery;
		int *atsteps;
		int *acurrstep;
		char *afname;
		char *aname;

		MPI_File *afh;

		int *current_ts;
		int *newts;
		double ***array;

		MPI_Comm comm;

	public:

		Modalysis();
		~Modalysis();

		void init(int, int, int, long long int, int);

		void postprocessdata(); 
		void setupPostprocess();
		void readFile();
		void vacf_();
		void msd_();
		void histo_();
		void fft_();

		void coanalyze(char *);	
		void readConfig(char *);
		void initAnalyses();
		void finiAnalyses();
		void processTimeStep(int, int);
		void process();
		void aalloc(int);
		bool check_new_timestep(int);

		void allocate_();

		long long int getnlocal();
		long long int getnglobal();

		void compute_vacf(int, double *);
		void compute_msd(int, double *);
		void compute_histo(int, double *);
		void compute_fft(int);
		void compute_fft_1d(int, int, int, double **);

};

#endif


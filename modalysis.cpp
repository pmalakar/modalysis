/*
 * Developed at Argonne National Laboratory
 *
 * Contact: pmalakar@anl.gov, malakar.preeti@gmail.com
 *
 */

#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <error.h>
#include <spi/include/kernel/memory.h>
#include "modalysis.h"

Modalysis::Modalysis() {

	comm = 0;

	nlocal = 0;
	x = v = f = NULL;

	atevery = NULL;
	atsteps = NULL;
	afname = NULL;
	newts = NULL;
	current_ts = NULL;

}

Modalysis::~Modalysis() {}

/*
 * initialize nlocal, nglobal, ntimesteps 
 */
void Modalysis::init(int me, int nprocs, int postprocess, long long int totalAtoms, int timesteps) {
	
	comm = MPI_COMM_WORLD;

	myrank = me;
	nprocs = nprocs;
	postprocess = postprocess;
	nglobal = totalAtoms;
	ntimesteps = timesteps;

	nlocal = nglobal/nprocs;
	//manage number of atoms in last rank
	if (myrank == nprocs-1)
		nlocal = nglobal - (nlocal*(nprocs-1));

#ifdef DEBUG
	printf ("%d: Init %d %d %d %d\n", myrank, nprocs, nglobal, nlocal, ntimesteps);
#endif

	if (comm == 0) {
		printf("\n%d: init Comm null %d\n", myrank, comm);
		exit(1);		
	}
	return;

}


long long int Modalysis::getnlocal() {
	return nlocal;
}

long long int Modalysis::getnglobal() {
	return nglobal;
}



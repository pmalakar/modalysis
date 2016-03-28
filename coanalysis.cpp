#include <mpi.h>
#include <errno.h>
#include <stdlib.h>
#include "modalysis.h"

/*
 * read config file for analysis (simulation output) file names and number of timesteps
 * process 0 reads, broadcasts
 */
void Modalysis::readConfig(char *analysiscfg)
{
	int i;

	if (analysiscfg == NULL) { 
		printf("Config file NULL error %d\n", myrank);
		return;
	}

	if (myrank == 0) {

    FILE *fp = fopen (analysiscfg, "r");
    if (fp == NULL) {
      printf("Config file open error %d %s\n", errno, strerror(errno));
      exit(1);
    }   
    
    fscanf(fp, "%d", &anum);
    aalloc(anum);

    i = 0;
    while(i<anum) {
      fscanf(fp, "%d %d %d %s", &adim[i], &atevery[i], &atsteps[i], afname+i*FILENAMELEN);
      ++i;
    }   
  
    fclose(fp);
  }

	if (comm == 0) {
		printf("\n%d: Comm null %d\n", myrank, comm);
		exit(1);		
	}

	MPI_Bcast(&anum, 1, MPI_INT, 0, comm);

	//allocate anum elements in all processes but 0 (reader)
	if(atsteps == NULL) aalloc(anum);

	if (MPI_Bcast(adim, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config adim bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atevery, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config atevery bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atsteps, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config atsteps bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(afname, anum*FILENAMELEN, MPI_CHAR, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis file name bcast error %d %s\n", errno, strerror(errno));

#ifdef DEBUG
	for (i=0; i<anum; i++) 
		printf("%d %d | %d %d %d %s\n", myrank, i, adim[i], atevery[i], atsteps[i], afname+i*FILENAMELEN);
#endif

}

void Modalysis::initAnalyses() {

	int i, j, retval;
	long long int numbytes;
	
	//allocate 
	for (i=0; i<anum; i++) {

		current_ts[i] = 0;
		newts[i]  = false;
		array[i] = new double*[atsteps[i]];

		for (j=0; j<atsteps[i]; j++) {
			numbytes = adim[i] * nlocal * sizeof(double);
			array[i][j] = (double *) malloc (numbytes);
			if(array[i][j] == NULL) 
				printf("%d: Failed to allocate %d bytes for array[%d][%d]\n", myrank, numbytes, i, j);
		}
			
		retval = MPI_File_open (comm, afname+i*FILENAMELEN, MPI_MODE_RDONLY, MPI_INFO_NULL, &afh[i]);
		if (retval != MPI_SUCCESS) 
			printf("\nAnalysis file open error %d %s\n", errno, strerror(errno));
	}


}

void Modalysis::finiAnalyses() {

	for (int i=0; i<anum; i++) {
		MPI_File_close (&afh[i]);
	}

	delete atsteps;
	delete atevery;
	delete afname;
	delete newts;
	delete current_ts;

	free(afh);

}

void Modalysis::processTimeStep(int aindex) {

	MPI_Offset offset, mpifo;
	MPI_Status status;
	long long int numelem;
	int n = current_ts[aindex];

	MPI_Scan(&nlocal, &nPartialSum, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

	offset = n * nglobal + nPartialSum - nlocal;
	mpifo = offset * adim[aindex] * sizeof(double);
	numelem = adim[aindex] * nlocal;

	if (MPI_File_read_at_all(afh[aindex], mpifo, array[aindex][n], numelem, MPI_DOUBLE, &status) != MPI_SUCCESS)
		perror("file read error");

	for (int j=0; j<3; j++) 
		printf("%d array[%d][%d][%d] = %lf\n", myrank, aindex, n, j, array[aindex][n][j]);

}

void Modalysis::process() {

	int i;
		
	for (i=0; i<anum; i++) {
		newts[i] = -1;
		if (atevery[i] == 1) {
		//	while (newts[i] != current_ts[i]) {
				//check if new timestep is written
		//	}
			processTimeStep(i);				
			current_ts[i] ++;
		}
	}

}

void Modalysis::aalloc(int anum)
{
	adim = new int[anum];
	atevery = new int[anum];
  atsteps = new int[anum];
  afname = new char[anum*FILENAMELEN];

	afh = (MPI_File *) malloc(anum * sizeof(MPI_File)); 
  current_ts = new int[anum];
  newts = new int[anum];

	array = new double**[anum];
	for (int n = 0; n<anum ; n++) 
		array[n] = NULL;
}

void Modalysis::coanalyze(char *cfg) {

	readConfig(cfg);
	initAnalyses();
	process();
	finiAnalyses();

}


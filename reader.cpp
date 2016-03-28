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
	else
		printf("%d: %s postprocess %d\n", myrank, analysiscfg, postprocess);

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
      fscanf(fp, "%d %d %s", &atevery[i], &atsteps[i], afname+i*FILENAMELEN);
      ++i;
    }   
  
    fclose(fp);
  }

	MPI_Bcast(&anum, 1, MPI_INT, 0, comm);

	//allocate anum elements in all processes but 0 (reader)
	if(atsteps == NULL) aalloc(anum);

	if (MPI_Bcast(atevery, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config atevery bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atsteps, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config atsteps bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(afname, anum*FILENAMELEN, MPI_CHAR, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis file name bcast error %d %s\n", errno, strerror(errno));

//#ifdef DEBUG
	for (i=0; i<anum; i++) 
		printf("%d %d | %d %d filename=%s\n", myrank, i, atevery[i], atsteps[i], afname+i*FILENAMELEN);
//#endif


}



void Modalysis::aalloc(int anum)
{
    atevery = new int[anum];
    atsteps = new int[anum];
    afname = new char[anum*FILENAMELEN];
}



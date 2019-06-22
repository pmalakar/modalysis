/*
 * Developed at Argonne National Laboratory
 *
 * Contact: pmalakar@anl.gov, malakar.preeti@gmail.com
 *
 */

#include <mpi.h>
#include <time.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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
		exit(1);
	}

	configFile = new char[strlen(analysiscfg)];
	strcpy(configFile, analysiscfg);

	if (myrank == 0) {

    FILE *fp = fopen (analysiscfg, "r");
    if (fp == NULL) {
      printf("Config file open error %d %s\n", errno, strerror(errno));
      exit(1);
    }   
    
    fscanf(fp, "%d", &anum);
    aalloc(anum);

		//Read number of dimensions | compute after every n steps | total N steps | file name | analysis name 
    i = 0;
    while(i<anum) {
      fscanf(fp, "%d %d %d %d %s %s", &adim[i], &atevery[i], &atsteps[i], &acurrstep[i], afname+i*FILENAMELEN, aname+i*ANAMELEN);
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
	if (MPI_Bcast(acurrstep, anum, MPI_INT, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis config acurrstep bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(afname, anum*FILENAMELEN, MPI_CHAR, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis file name bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(aname, anum*ANAMELEN, MPI_CHAR, 0, comm) != MPI_SUCCESS)
		printf("\nAnalysis name bcast error %d %s\n", errno, strerror(errno));

#ifdef DEBUG
	if (myrank < 2)
	for (i=0; i<anum; i++) 
		printf("%d %d | %d %d %d %d %s %s\n", myrank, i, adim[i], atevery[i], atsteps[i], acurrstep[i], afname+i*FILENAMELEN, aname+i*ANAMELEN);
#endif

}

void Modalysis::initAnalyses() {

	int i, j, retval;
	long long int numbytes;
	
	//allocate 
	for (i=0; i<anum; i++) {

		current_ts[i] = -1;

		array[i] = new double*[atsteps[i]];

		for (j=0; j<atsteps[i]; j++) {
			numbytes = adim[i] * nlocal * sizeof(double);
			array[i][j] = (double *) malloc (numbytes);
			if(array[i][j] == NULL) 
				if (myrank == 0) printf("%d: Failed to allocate %lld bytes for array[%d][%d]\n", myrank, numbytes, i, j);
		}
			
		retval = MPI_File_open (comm, afname+i*FILENAMELEN, MPI_MODE_RDONLY, MPI_INFO_NULL, &afh[i]);
		if (retval != MPI_SUCCESS) 
			if (myrank == 0) printf("\nAnalysis file open error %d %s\n", errno, strerror(errno));

	}

}


void Modalysis::processTimeStep(int aindex, int n) {

	MPI_Offset offset, mpifo;
	MPI_Status status;
	long long int numelem;
	int nend = n + atevery[aindex];
	double time, stime;

	if (myrank ==0)
		printf("Process %d to %d step for %d analysis [%s]\n", n, nend, aindex, aname+aindex*ANAMELEN);

	MPI_Scan(&nlocal, &nPartialSum, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

	offset = n * nglobal + nPartialSum - nlocal;
	mpifo = offset * adim[aindex] * sizeof(double);
	numelem = adim[aindex] * nlocal;

#ifdef DEBUG
	if (myrank ==0) {
		printf("test %lld %lld %lld: %lld\n", nlocal, nPartialSum, nglobal, numelem);
		printf("Read from %lld for %d analysis [%s]\n", mpifo, aindex, aname+aindex*ANAMELEN);
	}
#endif

	double read_time = MPI_Wtime();
	for (int nidx=n; nidx<nend ; nidx++) 
		if (MPI_File_read_at_all(afh[aindex], mpifo, array[aindex][nidx], numelem, MPI_DOUBLE, &status) != MPI_SUCCESS)
			perror("File read error");
	read_time = MPI_Wtime() - read_time;
	if (myrank == 0 && n<4) printf("%d Time to read for %s: %lf\n", n, aname+aindex*ANAMELEN, read_time);

#ifdef DEBUG
	for (int j=0; j<3; j++) 
		printf("%d %s array[%d][%d][%d] = %lf\n", myrank, aname+aindex*ANAMELEN, aindex, n, j, array[aindex][n][j]);
#endif

	//Call analysis functions
	
	if (strcmp(aname+aindex*ANAMELEN, "vacf") == 0) {	

		if (n == 0)
		 for (int k = 0; k<nlocal*adim[aindex] ; k++) 
			voriginal[k] = array[aindex][0][k];
		
		stime = MPI_Wtime();
		compute_vacf(n, array[aindex][n]);
		stime = MPI_Wtime() - stime;
		MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
		if (myrank == 0) printf("%lld Time to compute vacf[%d]: %lf\n", nglobal, n, time);
	
	}

	else if (strcmp(aname+aindex*ANAMELEN, "msd") == 0)	{

		if (n == 0)
		 for (int k = 0; k<nlocal*adim[aindex] ; k++) 
			xoriginal[k] = array[aindex][0][k];
		
		stime = MPI_Wtime();
		compute_msd(n, array[aindex][n]);
		stime = MPI_Wtime() - stime;
		MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
		if (myrank == 0) printf("%lld Time to compute msd[%d]: %lf\n", nglobal, n, time);
	}

	else if (strncmp(aname+aindex*ANAMELEN, "histo", 5) == 0)	{

		stime = MPI_Wtime();
		compute_histo(n, array[aindex][n]);
			//compute_histo(n, nend, k, array[aindex]);
		stime = MPI_Wtime() - stime;
		MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
		if (myrank == 0) printf("%lld Time to compute histo[%d]: %lf\n", nglobal, n, time);
	}

	else if (strncmp(aname+aindex*ANAMELEN, "fft", 3) == 0)	{

		stime = MPI_Wtime();
		for (int k=0; k<nlocal ; k++) 
			compute_fft_1d(n, nend, k, array[aindex]);
		stime = MPI_Wtime() - stime;
		MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, comm);
		if (myrank == 0) printf("%lld Time to compute fft[%d, %d): %lf\n", nglobal, n, nend, time);
	}

}

void Modalysis::process() {

	int aindex;
	bool alldone = false;
	bool done[anum];
		
	while (!alldone) {
		
	 alldone = true;
	 for (aindex=0; aindex<anum; aindex++) {

	 	done[aindex] = false;
		//if (acurrstep[aindex]+1 == atsteps[aindex])
		if (current_ts[aindex]+1 == atsteps[aindex])
			done[aindex] = true;
		else
			done[aindex] = false, alldone = false;

#ifdef DEBUG
		if(myrank < 1) 
			if (done[aindex]) printf ("done for %d: %d %d\n", aindex, acurrstep[aindex], atsteps[aindex]);
#endif
		
	 }

	 if (alldone) break;	 

	 for (aindex=0; aindex<anum; aindex++) {

			//Check if all timesteps processed, round robin for all
			if (current_ts[aindex]+1 == atsteps[aindex]) continue;

			//while (check_new_timestep(aindex) != true) sleep(10);  
			if (check_new_timestep(aindex) != true) {
			 	sleep(5); 
				continue;	//check for the other analysis timesteps
			} 
			processTimeStep(aindex, current_ts[aindex]+1);				
			current_ts[aindex] += atevery[aindex]; 

	 }

	}
}

bool Modalysis::check_new_timestep(int aidx) {

	bool value = false;
	
	if (myrank == 0) {
#ifdef DEBUG
		printf("0 check %s for new timestep: current_ts[%d]=%d acurrstep[%d]=%d\n", configFile, aidx, current_ts[aidx], aidx, acurrstep[aidx]);
#endif
		if (configFile == NULL)
			printf("configFile NULL\n");

		int num;
		int a, b, c, currstep;
		char str1[FILENAMELEN], str2[FILENAMELEN];

    FILE *fp = fopen (configFile, "r");
    if (fp == NULL) {
      printf("Config file open error %d %s\n", errno, strerror(errno));
      exit(1);
    }   
    
    fscanf(fp, "%d", &num);

		//Read number of dimensions | compute after every n steps | total N steps | file name | analysis name 
    int i = 0;
    while(i<num) {
      fscanf(fp, "%d %d %d %d %s %s", &a, &b, &c, &currstep, str1, str2);
			if (i == aidx) break;
      ++i;
    }   
		fclose(fp);

		acurrstep[aidx] = currstep;
		//printf ("%d current step read %d\n", aidx, acurrstep[aidx]);

	}

	MPI_Bcast(&acurrstep[aidx], 1, MPI_INT, 0, comm);

	if (current_ts[aidx] + atevery[aidx] <= acurrstep[aidx])
		value = true;
	else
		value = false;

	return value;

}


void Modalysis::aalloc(int anum)
{
	adim = new int[anum];
	atevery = new int[anum];
  atsteps = new int[anum];
  acurrstep = new int[anum];
  afname = new char[anum*FILENAMELEN];
  aname = new char[anum*ANAMELEN];

	afh = (MPI_File *) malloc(anum * sizeof(MPI_File)); 

  current_ts = new int[anum];

	array = new double**[anum];
	for (int n = 0; n<anum ; n++) 
		array[n] = NULL;
}

void Modalysis::finiAnalyses() {

	for (int i=0; i<anum; i++) {
		MPI_File_close (&afh[i]);
	}

	delete adim;
	delete atevery;
	delete atsteps;
	delete acurrstep;
	delete afname;
	delete aname;

	delete current_ts;

	free(afh);

}

void Modalysis::coanalyze(char *cfg) {

	readConfig(cfg);

	allocate_();
	initAnalyses();
	process();
	finiAnalyses();

}


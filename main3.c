/*
   stress test program
   version 03
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include <pthread.h>

#include <omp.h>
#include <mpi.h>

#include <math.h>

#include <time.h>
#include <sys/time.h>

struct sysinfo *info;
int i;
int j;
int k;
char name[MPI_MAX_PROCESSOR_NAME];
int resultlen;
int myrank;
int nprocs;
int ret;
unsigned long memorysize;
unsigned long memorysize_5;
unsigned long memorysize_dimension;
unsigned long MAXI;
unsigned long MAXJ;
unsigned long MAXK;
int cpunumber_node;
int dimension = 3;
int thread_number = 2;
pthread_t *thread2;
int arg1 = 0;
int arg2 = 1;
int ret1;
size_t filemaxsize;
int ret2;
double **a;
char filename[20];
typedef struct {
	int id;
	unsigned long loop;
	unsigned long runtime;
	unsigned long maxi;
	int MAX_THREADS;
	FILE *fp;
	unsigned int WriteLen;
} THREADS_STRUCT;

double cbrt_newton(double a, double x) {
	double e;
	do {
		e = (x * x * x - a) / (3.0 * x * x);
		x = x - e;
	} while ( fabs(e) > 1.0e-16 );
	return x;
}

double cbrt_simple(double x) {
	return pow(x, 1.0 / 3.0);
}

int ogafile (size_t gmin, 
             size_t gmax, 
             size_t lmin, 
             size_t lmax, 
             size_t mod_nmax, 
             size_t mod_minsize, 
             size_t mod_maxsize,
             size_t mod_dummynum);
/*
const size_t	gmin	     = 20;
const size_t	gmax	     = 30;
const size_t	lmin	     = 5;
const size_t	lmax	     = 20;
const size_t	mod_nmax     = 10000;
const size_t	mod_minsize  = 512;
const size_t	mod_maxsize  = 1536;
const size_t	mod_dummynum = 100000;
*/

void *cal_thread(void *p)
{
	/*
	int *ips = (int *)p;
	sscanf(p->id, "%d", ip);
	*/
	int ic;
	int jc;
	int kc;
	int ip;
	int ii;
	int im;
	unsigned int WriteLen;
	unsigned long looptmp;
	unsigned long loop;
	unsigned long maxi;
	int counter=0;
	THREADS_STRUCT *pt = (THREADS_STRUCT *)p;
	double sum=0.0;
	double Totalrealsec=0.0;
	double Totaliorealsec=0.0;
	double Totalsec=0.0;
	double Totalmemsize=0.0;
	double *b;
	struct timeval startTime, endTime, ioTime;
	clock_t startClock, endClock, ioClock;

	ip      = pt->id;
	looptmp = pt->loop;
	maxi    = pt->maxi;
#ifdef DEBUG
	printf(" looptmp = %ld , ip = %d \n",looptmp,ip);
	printf(" maxi    = %ld , ip = %d \n",maxi,ip);
	printf(" maxi pt = %ld , ip = %d \n",pt->maxi,ip);
#endif
	Totalrealsec   = 0.0;
	Totaliorealsec = 0.0;
	Totalsec       = 0.0;
	Totalmemsize   = 0.0;

	if(looptmp == 0){
		loop = 10000000;
	} else {
		loop = looptmp;
	}
	for(ii = 0 ;ii < loop ; ii++){

#ifdef DEBUG
	printf(" ii = %d ,ip = %d \n",ii,ip);
	printf(" loop = %ld ,ip = %d \n",loop,ip);
	printf(" maxi = %ld ,ip = %d \n",maxi,ip);
#endif
		sum = 0.0;

		gettimeofday(&startTime, NULL);
		startClock = clock();

		if(maxi < 3){
			maxi = 4;
		}
		for(im = 3;im<maxi;im++){
			b = (double *)malloc((sizeof(double))*im*im*im);
			if( b == NULL) {
				err(EXIT_FAILURE, "can not create thread 2" );
			}
				
			for(ic=0;ic<im;ic++){
				for(jc=0;jc<im;jc++){
					for(kc=0;kc<im;kc++){
						b[ic*jc*kc] = 0.0; 
					}
				}
			}
			for(ic=0;ic<im;ic++){
				b[ic] = 1.0; 
			}
			for(ic=1;ic<(im-1);ic++){
				for(jc=1;jc<(im-1);jc++){
					for(kc=1;kc<(im-1);kc++){
						sum = sum + (b[ic*jc*kc] + b[(ic-1)*jc*kc] + b[(ic+1)*jc*kc] 
						                     	+ b[ic*(jc-1)*kc] + b[ic*(jc+1)*kc] 
                                                                     	+ b[ic*jc*(kc-1)] + b[ic*jc*(kc+1)] 
                                                    	)/7.0; 
					}
				}
			}
			for(ic=0;ic<im;ic++){
				b[ic] = sum; 
			}
	
			gettimeofday(&ioTime, NULL);
			ioClock = clock();
	
			WriteLen = (int)((im*im*im)-1);
			fwrite((char *)b,sizeof(char),WriteLen,pt->fp);
			fseek(pt->fp,  0L, SEEK_SET);

			free(b);
		}
		gettimeofday(&endTime, NULL);
		endClock = clock();
	
		time_t diffsec = difftime(endTime.tv_sec, startTime.tv_sec);
		suseconds_t diffsub = endTime.tv_usec - startTime.tv_usec;
		if (diffsub < 0) {
			diffsec -= 1;
			diffsub = 1000000 + diffsub;
		}
		double realsec = diffsec+diffsub*1e-6;
		double cpusec = (endClock - startClock)/(double)CLOCKS_PER_SEC;
		double percent = 100.0*cpusec/realsec;
	
		time_t iodiffsec = difftime(endTime.tv_sec, ioTime.tv_sec);
		suseconds_t iodiffsub = endTime.tv_usec - ioTime.tv_usec;
		if (iodiffsub < 0) {
			iodiffsec -= 1;
			iodiffsub = 1000000 + iodiffsub;
		}
		double iorealsec = iodiffsec+iodiffsub*1e-6;
		double iocpusec = (endClock - ioClock)/(double)CLOCKS_PER_SEC;
		double iopercent = 100.0*iocpusec/iorealsec;

		Totalrealsec   += realsec;
		Totaliorealsec += iorealsec;
		Totalsec       += realsec+ iorealsec;
		Totalmemsize   += im;

		counter++;
	}
	printf("******************************** \n");
	printf("**** Real Time  %f s\n", Totalrealsec);
	printf("**** I/O  Time  %f s\n", Totaliorealsec/Totalmemsize);
	printf("**** Ave.       %f %\n", Totalrealsec/counter);
	printf("******************************** \n");

#ifdef DEBUG
	printf(" write start WriteLen = %d ,ip = %d \n",pt->WriteLen,ip);
#endif
	return 0;
}
	
void error_arg()
{
	printf(" argment is less \n");
	printf(" ./stress sigle memsize loop time \n");
	printf(" Example ... \n");
	printf(" ./stress 2 500 30 0 2000 \n");
	printf(" \n");
	printf(" 1st argument ... IO or Calculation or Mix ? \n");
	printf(" if IO only, argument is 0 \n");
	printf(" if Calculation only, argument is 1 \n");
	printf(" if mix mode, argument is 2\n");
	printf(" \n");
	printf(" 2nd argument ... Array size (memory size). \n");
	printf(" if memsize not need, argument is 0 \n");
	printf(" \n");
	printf(" 3rd argument ... Loop number. \n");
	printf(" if loop not need, argument is 0 \n");
	printf(" \n");
	printf(" 4th argument ... Time. \n");
	printf(" if time not need, argument is 0 \n");
	printf(" \n");
	printf(" 5th argument ... IO size. \n");
	printf(" unit is Byte. \n");
	printf(" if size not need, argument is 0 \n");
	printf(" if 1st argument is 1, argument is 0 \n");
}

int main(int argc,char *argv[])
{
	THREADS_STRUCT *ts;

	if(argc != 6){
		error_arg();
		exit(-1);
	}
	info = malloc(sizeof(struct sysinfo));
	sysinfo(info);

#ifdef DEBUG
	printf(" %ld kb\n ",info->freeram);
#endif
	
#ifdef MULTI
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Get_processor_name(name, &resultlen);
	printf(" My processor is %s Total Free Memory is %ld \n",name,info->freeram);
	printf(" Total processor is %d myrank is %d \n",nprocs,myrank);
#endif

	printf(" Memory is %ld \n",info->totalram);
	printf(" Free Memory is %ld \n",info->freeram);

	cpunumber_node = (sysconf(_SC_NPROCESSORS_CONF) * 1 );
	memorysize = info->totalram * 0.7;
	memorysize_5 = (memorysize / cpunumber_node) * 0.7;
	/*
	memorysize_dimension = memorysize_5^(1/dimension);
	memorysize_dimension = (int)(memorysize_5)^(1/4);
	memorysize_dimension = (int)(cbrt_newton(memorysize_5, 10.0));
	memorysize_dimension = (int)((memorysize_5)/4);
	MAXI = memorysize_dimension/cpunumber_node/sizeof(double);
	MAXJ = memorysize_dimension/cpunumber_node/sizeof(double);
	MAXK = memorysize_dimension/cpunumber_node/sizeof(double);
	printf(" memorysize_dimension  is %ld \n",memorysize_dimension);
	*/

	memorysize_dimension = (cbrt_simple(memorysize_5))*0.7;

	if(atoi(argv[2]) == 0){
		MAXI = memorysize_dimension;
		MAXJ = memorysize_dimension;
		MAXK = memorysize_dimension;
	} else if( atoi(argv[2]) > memorysize_dimension ){
		MAXI = memorysize_dimension;
		MAXJ = memorysize_dimension;
		MAXK = memorysize_dimension;
		printf(" ****** Attension! ******* \n");
		printf(" Array size may be over ... \n");
	} else {
		MAXI = atoi(argv[2]);
		MAXJ = atoi(argv[2]);
		MAXK = atoi(argv[2]);
	}

	printf(" MAXI is %ld \n",MAXI);
	printf(" MAXJ is %ld \n",MAXJ);
	printf(" MAXK is %ld \n",MAXK);

	filemaxsize = atoi(argv[5]);
	printf(" File Size is %ld \n",filemaxsize);

#ifdef DEBUG
	printf(" MAXI is %ld \n",MAXI);
	printf(" MAXJ is %ld \n",MAXJ);
	printf(" MAXK is %ld \n",MAXK);
#endif
	printf(" CPU number is %d \n",cpunumber_node);

	thread2 = (pthread_t *)malloc((sizeof(pthread_t)*(cpunumber_node+1)));
	ts = (THREADS_STRUCT *)malloc( (sizeof(THREADS_STRUCT)*(cpunumber_node+1)) );

	a = (double **)malloc( (sizeof(double *)*cpunumber_node) );

	if(atoi(argv[1]) == 1){
		for(i=0;i<cpunumber_node;i++){
			if(nprocs == 0){
				sprintf(filename, "tmpfile_%03d.txt", i);
			} else {
				sprintf(filename, "tmpfile_%03d_%04d.txt", i,myrank);
			}
			if( (ts[i].fp= fopen(filename,"wb")) == NULL){
				printf("Open Error \n");
				exit(-1);
			}
			ts[i].id = i;
			ts[i].maxi = MAXI;
			ts[i].loop = atoi(argv[3]);
			ts[i].runtime = atoi(argv[4]);
			ts[i].MAX_THREADS = cpunumber_node;
#ifdef DEBUG
			printf(" ID local = %d \n",ts[i].id);
			printf(" ID local = %d \n",ts[i].MAX_THREADS);
#endif
			ret2 = pthread_create(&thread2[i],NULL,(void *)cal_thread,(void *) &ts[i]);
			if (ret2 != 0) {
				err(EXIT_FAILURE, "can not create thread 2: %s", strerror(ret2) );
			}
		}
	}else if(atoi(argv[1]) == 0){
	
		ret = ogafile(20, 30, 5, 20, 10000, filemaxsize*0.8, filemaxsize, 100000);
/*
const size_t	gmin	     = 20;
const size_t	gmax	     = 30;
const size_t	lmin	     = 5;
const size_t	lmax	     = 20;
const size_t	mod_nmax     = 10000;
const size_t	mod_minsize  = 512;
const size_t	mod_maxsize  = 1536;
const size_t	mod_dummynum = 100000;
*/
	}else if(atoi(argv[1]) == 2 ){
		for(i=0;i<cpunumber_node;i++){
			if(nprocs == 0){
				sprintf(filename, "tmpfile_%03d.txt", i);
			} else {
				sprintf(filename, "tmpfile_%03d_%04d.txt", i,myrank);
			}
			if( (ts[i].fp= fopen(filename,"wb")) == NULL){
				printf("Open Error \n");
				exit(-1);
			}
			ts[i].id = i;
			ts[i].maxi = MAXI;
			ts[i].loop = atoi(argv[3]);
			ts[i].runtime = atoi(argv[4]);
			ts[i].MAX_THREADS = cpunumber_node;
#ifdef DEBUG
			printf(" ID local = %d \n",ts[i].id);
			printf(" ID local = %d \n",ts[i].MAX_THREADS);
#endif
			ret2 = pthread_create(&thread2[i],NULL,(void *)cal_thread,(void *) &ts[i]);
			if (ret2 != 0) {
				err(EXIT_FAILURE, "can not create thread 2: %s", strerror(ret2) );
			}
		}
		ret = ogafile(20, 30, 5, 20, 10000, filemaxsize*0.8, filemaxsize, 100000);
	}

	if(atoi(argv[1]) == 1 || atoi(argv[1]) == 2){
		for(i=0;i<cpunumber_node;i++){
			ret2 = pthread_join(thread2[i],NULL);
			if (ret2 != 0) {
				err(EXIT_FAILURE, "can not join thread 2: %s", strerror(ret2) );
			}
		}
	}

#ifdef MULTI
	MPI_Finalize();
#endif
	free(info);
	return 0;
}


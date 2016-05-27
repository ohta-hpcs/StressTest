/*
   stress test program
   Need MPI Library
   version 02
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

#define VERSION 02

struct sysinfo *info;
int i;
int j;
int k;
char name[MPI_MAX_PROCESSOR_NAME];
int resultlen;
int myrank;
int nprocs;
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
int ret2;
double **a;
char filename[20];
typedef struct {
	int id;
	unsigned long loop;
	unsigned long runtime;
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
	unsigned long loop;
	int counter=0;
	THREADS_STRUCT *pt = (THREADS_STRUCT *)p;
	double sum=0.0;
	double Totalrealsec=0.0;
	double Totaliorealsec=0.0;
	double Totalsec=0.0;
	struct timeval startTime, endTime, ioTime;
	clock_t startClock, endClock, ioClock;

	ip   = pt->id;
	loop = pt->loop;
#ifdef DEBUG
	printf(" ip = %d \n",ip);
#endif

	a[ip] = (double *)malloc((sizeof(double))*MAXI*MAXJ*MAXK);
	if(loop == 0){
		loop = 10000000;
	} 
	for(ii = 0 ;ii < loop ; ii++){

		gettimeofday(&startTime, NULL);
		startClock = clock();

		for(ic=0;ic<MAXI;ic++){
			for(jc=0;jc<MAXJ;jc++){
				for(kc=0;kc<MAXK;kc++){
					a[ip][ic*jc*kc] = 0.0; 
				}
			}
		}
		for(ic=0;ic<MAXI;ic++){
			a[ip][ic] = 1.0; 
		}
		for(ic=0;ic<MAXI;ic++){
			for(jc=0;jc<MAXJ;jc++){
				for(kc=0;kc<MAXK;kc++){
					sum = sum + (a[ip][ic*jc*kc] + a[ip][(ic-1)*jc*kc] + a[ip][(ic+1)*jc*kc] 
						                     + a[ip][ic*(jc-1)*kc] + a[ip][ic*(jc+1)*kc] 
                                                                     + a[ip][ic*jc*(kc-1)] + a[ip][ic*jc*(kc+1)] 
                                                    )/7.0; 
				}
			}
		}
		for(ic=0;ic<MAXI;ic++){
			a[ip][ic] = sum; 
		}
	
		gettimeofday(&ioTime, NULL);
		ioClock = clock();

		pt->WriteLen = (int)((MAXI*MAXJ*MAXK)-1);
		fwrite((char *)a[ip],sizeof(double),pt->WriteLen,pt->fp);
		fseek(pt->fp,  0L, SEEK_SET);
		
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

		Totalrealsec += realsec;
		Totaliorealsec += iorealsec;
		Totalsec  += realsec+ iorealsec;

		counter++;
		/*
		if(Totalsec > pt->runtime) break;
		*/
	}
	printf("******************************** \n");
	printf("**** Real Time  %f s\n", Totalrealsec);
	printf("**** I/O  Time  %f s\n", Totaliorealsec);
	printf("**** Ave.       %f %\n", Totalrealsec/counter);
	printf("******************************** \n");

#ifdef DEBUG
			printf(" write start WriteLen = %d ,ip = %d \n",pt->WriteLen,ip);
#endif
	free(a[ip]);
	return 0;
}
	
void help_arg()
{
	printf(" help \n");
	printf(" Example ... \n");
	printf(" ./stress 0 500 300 0 \n");
}
void error_arg()
{
	printf(" argment is less \n");
	printf(" ./stress sigle memsize loop time \n");
	printf(" Example ... \n");
	printf(" ./stress 0 500 300 0 \n");
	printf(" \n");
	printf(" 1st argument ... Single or Parallel? \n");
	printf(" if single process, argument is 0 \n");
	printf(" if parallel process, argument is 1 \n");
	printf(" \n");
	printf(" 2nd argument ... Array size (memory size). \n");
	printf(" if memsize not need, argument is 0 \n");
	printf(" \n");
	printf(" 3rd argument ... Loop number. \n");
	printf(" if loop not need, argument is 0 \n");
	printf(" \n");
	printf(" 4th argument ... Time. \n");
	printf(" if time not need, argument is 0 \n");
}

int main(int argc,char *argv[])
{
	THREADS_STRUCT *ts;

	if(argc != 5){
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

	if(argv[2] == 0){
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

#ifdef DEBUG
	printf(" MAXI is %ld \n",MAXI);
	printf(" MAXJ is %ld \n",MAXJ);
	printf(" MAXK is %ld \n",MAXK);
#endif
	printf(" CPU number is %d \n",cpunumber_node);

	while((result=getopt(argc,argv,"s:m:l:t:hv"))!=-1){
		switch(result){

			/* 値をとらないオプション */
			case 'h':
				/* getoptの返り値は見付けたオプションである. */
				help_arg();
				fprintf(stdout,"%c\n",result);
				break;
			case 'v':
				fprintf(stdout,"%c\n",result);
				fprintf(stdout," Version is %d \n",VERSION);
				break;
	
			case 's':
				sa = atoi(optarg);
				break;
			case 'l':
				la = atoi(optarg);
				break;
			case 't':
				ta = atoi(optarg);
				break;
			case 'm':
				/* 値を取る引数の場合は外部変数optargにその値を格納する. */
				fprintf(stdout,"%c %s\n",result,optarg);
				ma = atoi(optarg);
				break;

				/* 以下二つのcaseは意味がないようだ.
				getoptが直接エラーを出力してくれるから.
				プログラムを終了するなら意味があるかも知れない */
			case ':':
				/* 値を取る引数に値がなかった場合:を返す. */
				fprintf(stdout,"%c needs value\n",result);
				break;

			/* getoptの引数で指定されなかったオプションを受け取ると?を返す. */
			case '?':
				fprintf(stdout,"unknown\n");
				break;
		}
	}

	thread2 = (pthread_t *)malloc((sizeof(pthread_t)*(cpunumber_node+1)));
	ts = (THREADS_STRUCT *)malloc( (sizeof(THREADS_STRUCT)*(cpunumber_node+1)) );

	a = (double **)malloc( (sizeof(double *)*cpunumber_node) );

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
		ts[i].loop = atoi(argv[3]);
		ts[i].runtime = atoi(argv[4]);
		ts[i].MAX_THREADS = cpunumber_node;
		/*
		printf(" ID local = %d \n",ts[i].id);
		printf(" ID local = %d \n",ts[i].MAX_THREADS);
		*/
		ret2 = pthread_create(&thread2[i],NULL,(void *)cal_thread,(void *) &ts[i]);
		if (ret2 != 0) {
			err(EXIT_FAILURE, "can not create thread 2: %s", strerror(ret2) );
		}
	}

	for(i=1;i<cpunumber_node;i++){
		ret2 = pthread_join(thread2[i],NULL);
		if (ret2 != 0) {
			err(EXIT_FAILURE, ret2, "can not join thread 2");
		}
	}

#ifdef MULTI
	MPI_Finalize();
#endif
	free(info);
	return 0;
}


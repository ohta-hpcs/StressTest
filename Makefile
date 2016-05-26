CC = icc
MPIHOME = /opt/intel/impi/4.1.3.048/intel64
MPIINCS = ${MPIHOME}/include
MPILIBS = ${MPIHOME}/lib
CFLAGS = -O3 -openmp -pthread -openmp-report=2 -I${MPIINCS} -parallel 
LFLAGS = -L${MPILIBS} -lmpi -pthread -openmp -openmp-report=2 -parallel

TARGET=stress

TARGET_SOURCES=main3.c
TARGETOBJ=${TARGET_SOURCES:.c=.o}

${TARGET}:${TARGETOBJ}
	${CC} ${CFLAGS} ${LFLAGS} -o $@ $^ -lm
clean: 
	rm -rf *.o ${TARGET} tmpfile_*.txt


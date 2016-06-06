CC = icc
CXX = icpc
MPIHOME = /opt/intel/impi/4.1.3.048/intel64
MPIINCS = ${MPIHOME}/include
MPILIBS = ${MPIHOME}/lib
CFLAGS = -O3 -openmp -pthread -openmp-report=2 -I${MPIINCS}
CXXFLAGS = -lstdc++
LFLAGS = -L${MPILIBS} -lmpi -pthread -openmp -openmp-report=2 -lstdc++


TARGET=stress

TARGET_SOURCES=main3.c
TARGETOBJ=${TARGET_SOURCES:.c=.o}

TARGETCXX_SOURCES=ogafile.cpp
TARGETCXXOBJ=${TARGETCXX_SOURCES:.cpp=.o}

${TARGET}:${TARGETOBJ} ${TARGETCXXOBJ}
	${CC} ${CFLAGS} ${LFLAGS} -o $@ $^ -lm
clean: 
	rm -rf *.o ${TARGET} tmpfile_*.txt bi_*.txt
ioclean: 
	rm -rf tmpfile_*.txt bi_*.txt


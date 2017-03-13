EXECS = dotnet
OBJS = MPI_Version_2.o Network.o Crack.o

MPICC = mpic++
MPIFLAGS = -std=c++11 
LD = g++
CFLAGS = -std=c++11 -c 
OFLAGS = -std=c++11 -o

all : ${EXECS}

$(EXECS) : $(OBJS)
	$(MPICC) $(MPIFLAGS) $(OBJS) -o $(EXECS)

MPI_Version_2.o : MPI_Version_2.cpp Network.cpp Network.h
	${MPICC} $(CFLAGS) MPI_Version_2.cpp

Network.o : Network.cpp Network.h crack.h
	${LD} $(CFLAGS) Network.cpp

Crack.o : Crack.cpp crack.h
	${LD} $(CFLAGS) Crack.cpp

clean :
	rm -f *.o ${EXECS}

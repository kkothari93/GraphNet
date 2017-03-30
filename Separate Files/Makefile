EXECS = dotnet
OBJS = main.o network.o crack.o helper_funcs.o mpi_network.o sac_network.o

MPICC = mpic++
MPIFLAGS = -std=c++11
LD = g++
CFLAGS = -std=c++11 -c
OFLAGS = -std=c++11 -o

all : ${EXECS}

$(EXECS) : $(OBJS)
	$(MPICC) $(MPIFLAGS) $(OBJS) -o $(EXECS)

main.o : main.cpp network.cpp mpi_network.cpp
	$(MPICC) $(CFLAGS) main.cpp

mpi_network.o : mpi_network.cpp mpi_network.h network.h sac_network.h
	$(LD) $(CFLAGS) mpi_network.cpp

network.o : network.cpp network.h helper_funcs.h crack.h
	$(LD) $(CFLAGS) network.cpp

crack.o : crack.cpp crack.h
	$(LD) $(CFLAGS) crack.cpp

helper_funcs.o : helper_funcs.cpp helper_funcs.h
	$(LD) $(CFLAGS) helper_funcs.cpp

sac_network.o : sac_network.cpp sac_network.h network.h
	$(LD) $(CFLAGS) sac_network.cpp


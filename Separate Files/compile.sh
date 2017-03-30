#!/bin/bash
make clean
rm forcesMPI.txt
make

if [ $2 == 'y' ];
then
	mpirun -np $3 ./dotnet $1 1
else
	./dotnet $1 0
fi 

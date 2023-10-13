#!/bin/bash

sed -e "s|output.nc|${USER}/output.nc|" namelist.in > namelist.tmp

mkdir /tmp/${USER}

if [ -z "$1" ]
then
	mpiexec -n 1 ./main.exe namelist.tmp
else
	mpiexec -n $1 ./main.exe namelist.tmp
fi
rm namelist.tmp 




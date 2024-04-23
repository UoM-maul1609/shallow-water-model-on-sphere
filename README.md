[![Pylint](https://github.com/petrafisher-code/ShallowWaterModelCode/actions/workflows/pylint.yml/badge.svg)](https://github.com/petrafisher-code/ShallowWaterModelCode/actions/workflows/pylint.yml)

## 1. Getting started
Clone the repository
```bash
git clone https://github.com/petrafisher-code/ShallowWaterModelCode.git
```

Move into the repository
```bash
cd ShallowWaterModelCode
```

<br>

## 2. Compilation and running the code

You can either compile and run using [Docker](https://www.docker.com/get-started/) or directly on your machine.
Docker is recommended. 

### Running using Docker
Install docker on your machine. Open docker, go back to the code/your IDE and run
the following command in terminal

```bash
docker build -t fortran_code ./ --progress=plain
```

Once complete, run the container
```bash
docker run fortran_code:latest
```

In the docker application under "Containers", copy the latest container ID to clipboard.

In another terminal (or after stopping the container) copy the output file to your local folder
```bash
docker cp [CONTAINER_ID]:/tmp/output.nc /path/to/your/tests/directory
```

### Running without Docker

You will need to install NETCDF on your system and note the directories where the 
library and module reside. 

In the Makefile change the variables NETCDF_FOR, NETCDF_C and NETCDF_LIB to suit 
your system. Sometimes when compiling NETCDF you may have to install both the 
fortran and c versions. If you do not need to do this, set NETCDF_FOR and NETCDF_C
equal to the same value. 

If you install both c and fortran versions, set NETCDF_LIB to -lnetcdff
If you only need to install one version, set NETCDF_LIB to lnetcdf
Finally, once netcdf is installed and the make file is edited, type "make" to compile
the code

Note, these variables are commented out in the Makefile. You may set them as environment
variables in your .bashrc, or .bash_profile and the Makefile will pick them up.

Type the following to run the code on a single processor
```bash
./main.exe ../config/namelist.in
```
To run the code with multiple processors as an MPI job
```bash
mpiexec -n 32 ./main.exe ../config/namelist.in
```

<br>

## 3. Animating results
There a python files provided to animate the resulting output.nc netCDF file. 
Move to the python directory,
```bash
cd src/python
```

Run the python animate_output.py file
```bash
python animate_output.py
```

This will produce frames in the ShallowWaterModelCode/output/frames directory
and animations in the ShallowWaterModelCode/output/animations directory.
(If you have already produced animations, the terminal will ask you whether you
want to overwrite the previous animations.)

<br>


## 4. Doxygen pages
Self generating documentation (from the source code) can be created by typing
"doxygen fortran.dxg" at the command line.

You can then view the html pages by opening the file at "doxygen/html/index.html"

<br>

## 5. Contributing development
At the moment the code is under private repo with no development branches. Contact 
the code owner to contribute ideas at 
[https://github.com/UoM-maul1609/shallow-water-model-on-sphere](https://github.com/UoM-maul1609/shallow-water-model-on-sphere).

<br>

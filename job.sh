#!/bin/bash
####################################
#
# Executing a scenario on the cluster 
#
# with one argument: name of the scenario (without file extension)
#
####################################

# request resources
#SBATCH --job-name=submission
#SBATCH --output=results_${1}.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00

# load modules
module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

# run code
rm -rf build/
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make install
./finalproject ../ini/${1}.txt > logs_${1}.txt &

# show job number
jobs
# fg <job_number>
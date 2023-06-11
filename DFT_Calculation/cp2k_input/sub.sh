#!/bin/bash -l
#PBS -q normal3
#PBS -N sp
#PBS -l nodes=1:ppn=24
#PBS -l walltime=1000:00:00

cd $PBS_O_WORKDIR

source /home/pengchao/App/gcc/9.3.0/env.sh
source /home/pengchao/App/cp2k/cp2k-9.1_plumed-grifoni/tools/toolchain/install/setup
export PATH=$PATH:/home/pengchao/App/cp2k/cp2k-9.1_plumed-grifoni/exe/local

workpath=$PBS_O_WORKDIR
index="001"

for (( i=10001; i<10002; i=i+1 )); do
	cd $workpath/$index/$i/pbe
	test $? -ne 0 && exit 1
	{ if [ ! -f tag_pbe_finished ] ;then
	mpirun -np 24 cp2k.popt pbe.inp |tee pbe.out 
    { if test $? -ne 0; then touch tag_pbe_failure
	else
	  touch tag_pbe_finished; fi }
	fi }
	
	cd $workpath/$index/$i/m062x
	test $? -ne 0 && exit 1
	{ if [ ! -f tag_m062x_finished ] ;then
	mpirun -np 24 cp2k.popt m062x.inp |tee m062x.out 
    { if test $? -ne 0; then touch tag_m062x_failure 
	else
	  touch tag_m062x_finished; fi } 
	fi }
done

cd $workpath
test $? -ne 0 && exit 1
wait
touch $PBS_JOBNAME

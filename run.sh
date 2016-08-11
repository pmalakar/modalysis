#!/bin/bash -l
#SBATCH -p regular
#SBATCH -t 00:05:00
#SBATCH -A m1248
##SBATCH -N 1
#SBATCH -J my_job
#SBATCH -o my_job.o%j
DIR=/var/opt/cray/dws/mounts/batch/2154246/ss//output
#DW stage_out source=$DIR destination=/global/cscratch1/sd/preeti/minimd/dout type=directory

#Edison has 24 cores per compute node

nodes=$SLURM_JOB_NUM_NODES
#echo $nodes

if [ $nodes -lt 1 ]; then
	echo $nodes 
	exit
fi

ppn=4 #24

#for nodes in 1 #2 3 4 5 # 6 7 8 9 10 11 12 13 14 #24 48 #72 96
#for dim in 128 256 512 1024
#do
		echo 
		echo "* * * * *"
		procs=`echo "$nodes*$ppn"|bc`
    OUTPUT=${procs}_${dim}
		echo $OUTPUT
		ARGS= #${dim}
###DW jobdw capacity=20GB access_mode=striped type=scratch
##mkdir $DW_JOB_STRIPED/outputdir
    srun -n $procs -N $nodes ./modalysis -p 1 #$ARGS #> $OUTPUT
		echo "* * * * *"
#done


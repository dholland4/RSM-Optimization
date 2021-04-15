#!/bin/bash
# Created by Darren Holland
# Modified by Darren Holland 2020-11-02

# File to overwrites submission variables and run Dakota/project

# Input to ./batchDakota.sh $1 $2 $3 $4 $5
# $1 is the Dakota/project input name
# $2 is the number of nodes
# $3 is the number of processors per node
# $4 is the number of total processors
# $5 is the folder name

# Overwrite submission settings
sed -i -- "s?nodes=1:ppn=4?nodes=$2:ppn=$3?g" DakotarunTEMP.pbs
sed -i -- "s?OMP_NUM_THREADS=4?OMP_NUM_THREADS=$4?g" DakotarunTEMP.pbs
sed -i -- "s?filename -t 4?dakota -i $1 -w $5".rst" > $5".log"?g" DakotarunTEMP.pbs

# Submit job
qsub DakotarunTEMP.pbs

#!/bin/bash
#SBATCH -Jtest                    # Job name
#SBATCH --account=gts-creinhard3-fy20phase2               # account to which job is charged
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=1G                # memory per core
#SBATCH -t00:05:00          # duration of the job (ex: 15 min)
#SBATCH -qinferno                               # QOS Name
#SBATCH -otest${1}.out                               # QOS Name

cd $SLURM_SUBMIT_DIR  
# ./${INPUT1} ${INPUT2} ${INPUT3} ${INPUT4} ${INPUT5}
# ${INPUT1} ${INPUT2} ${INPUT3} ${INPUT4} ${INPUT5}
echo $1 $2 $3 $4 $5 $6 $7 $8
$1 $2 $3 $4 $5 $6 $7 $8


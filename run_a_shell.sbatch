#!/bin/bash
#SBATCH --account=gts-creinhard3-fy20phase2               # account to which job is charged
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=1G                # memory per core
#SBATCH -t00:05:00          # duration of the job (ex: 15 min)
#SBATCH -qinferno                               # QOS Name

cd $SLURM_SUBMIT_DIR 
file=$0
jobidinfo="${file//[!0-9]/}"
echo $jobidinfo
# ./${INPUT1} ${INPUT2} ${INPUT3} ${INPUT4} ${INPUT5}
# ${INPUT1} ${INPUT2} ${INPUT3} ${INPUT4} ${INPUT5}
echo ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21} ${22} ${23} ${24} ${25}
${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21} ${22} ${23} ${24} ${25}


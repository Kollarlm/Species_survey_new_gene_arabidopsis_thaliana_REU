#!/bin/bash --login
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name exonerate_L1
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}

#Export paths to conda
export PATH="${conda}/envs/exonerate/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/exonerate/lib:$LD_LIBRARY_PATH"

#contig=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./assembly_pseudomolecules/S1_inversion_contigs.txt)
#contig=`echo $contig`

exonerate Laurens_gene ./S1/S1-v1.fa --bestn 1 --showtargetgff yes 

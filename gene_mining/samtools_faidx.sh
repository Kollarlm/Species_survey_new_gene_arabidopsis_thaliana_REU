#!/bin/bash --login
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name samtools_index
#SBATCH --output=%x-%j.SLURMout
##SBATCH --partition=niederhu
##SBATCH --account=niederhu
##SBATCH --array=1-17%17

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}

#Export paths to conda
export PATH="${conda}/envs/exonerate/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/exonerate/lib:$LD_LIBRARY_PATH"

#contig=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./assembly_pseudomolecules/L1_inversion_contigs.txt)
#contig=`echo $contig`

#Running for chromosomes
#for i in {1..14}
#do
#samtools faidx ./assembly_pseudomolecules/S1-v1_allmaps6.fa chr${i} -o chr${i}.S1
#done

#Running for contigs
#samtools faidx ./assembly_pseudomolecules/L1_ragtag_IM62.fa "$contig"_RagTag -o "$contig".L1

# Removing breakpoints

#breakpoint=$(sed -n "$SLURM_ARRAY_TASK_ID"p breakpoints_merged_L1.txt)
#breakpoint=`echo $breakpoint`

samtools faidx L1-v1_pseudochromosome.fa  chr5:13629641-13645049 -o L1_inversion5_start.fa
samtools faidx L1-v1_pseudochromosome.fa  chr5:17847181-17857181 -o L1_inversion5_end.fa  

samtools faidx S1-v1_pseudochromosome.fa  chr5:14094885-14125402 -o S1_inversion5_start.fa
samtools faidx S1-v1_pseudochromosome.fa  chr5:18136144-18146144 -o S1_inversion5_end.fa


samtools faidx S1-v1_pseudochromosome.fa  chr8:846808-858736 -o S1_inversion8_start.fa
samtools faidx S1-v1_pseudochromosome.fa  chr8:6465310-6527136 -o S1_inversion8_end.fa

samtools faidx L1-v1_pseudochromosome.fa  chr8:849409-851381 -o L1_inversion8_start.fa
samtools faidx L1-v1_pseudochromosome.fa  chr8:7604769-7703000 -o L1_inversion8_end.fa

samtools faidx L1-v1_pseudochromosome.fa  chr14:5308079-5329939 -o L1_inversion14_start.fa
samtools faidx L1-v1_pseudochromosome.fa  chr14:7791347-7820218 -o L1_inversion14_end.fa

samtools faidx S1-v1_pseudochromosome.fa  chr14:5868251-5892556 -o S1_inversion14_start.fa
samtools faidx S1-v1_pseudochromosome.fa  chr14:8671180-8677481 -o S1_inversion14_end.fa


samtools faidx S1-v1_pseudochromosome.fa  chr8:5996241-5999201 -o S1_inversionsmall8_start.fa
samtools faidx S1-v1_pseudochromosome.fa  chr8:6260288-6263731 -o S1_inversionsmall8_end.fa

samtools faidx L1-v1_pseudochromosome.fa  chr8:1026660-1032696 -o L1_inversionsmall8_start.fa
samtools faidx L1-v1_pseudochromosome.fa  chr8:1242892-1249416 -o L1_inversionsmall8_end.fa

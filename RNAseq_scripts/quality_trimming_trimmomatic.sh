#!/bin/bash --login
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --job-name rna_cleaning_L1_field1
#SBATCH --output=rna_cleaning_L1_field1.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:${LD_LIBRARY_PATH}"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/polishing/share/trimmomatic-0.39-2/adapters"

#################
### VARIABLES ###
#################

genotype=L1_tissues

sample=Field1

#path2="/mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018"

adapters=NexteraPE-PE.fa

## DONT TOUCH Fastq files, these should not have to be changed, but should set automatically
path3="/mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/${genotype}/${sample}/"
r1="${path3}/combined.1.fastq.gz"
r2="${path3}/combined.2.fastq.gz"
t1="${path3}/trimmed.1.fastq.gz"
t2="${path3}/trimmed.2.fastq.gz"
t3="${path3}/trimmed.1.single.fastq.gz"
t4="${path3}/trimmed.2.single.fastq.gz"

if ls ${path3}/*R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*R1_001.fastq.gz > $r1
			cat ${path3}/*R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*R1_001.fastq.gz > $r1
		fi
	fi
elif ls ${path3}/*R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*R2_001.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*R1_001.fastq.gz > $r1
			cat ${path3}/*R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*R1_001.fastq.gz > $r1
		fi
	fi
else
	echo "Data Missing"
fi

#Trim & QC reads
if [ -f ${t1} ]
then
	if [ ${PE} = "TRUE" ]
	then
		echo "To rerun this step, please delete ${t1} & ${t2} and resubmit"
	else
		echo "To rerun this step, please delete ${t1} and resubmit"
	fi
else
	if [ ${PE} = "TRUE" ]
	then
		echo "Running trimmomatic PE"
		trimmomatic PE \
			-threads 20 \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
			ILLUMINACLIP:${adapter_path}/${adapters}:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${t2} ${r1} ${r2}
	elif [ ${PE} = "FALSE" ]
	then
		echo "Running trimmomatic SE"
		trimmomatic SE \
			-threads 20 \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${t1} \
			ILLUMINACLIP:${adapter_path}/${adapters}:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:50 
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${r1}
	fi
fi
rm ${r1} ${r2}

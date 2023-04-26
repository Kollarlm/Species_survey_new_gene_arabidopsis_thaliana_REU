#!/bin/bash --login
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name star_align_S1
#SBATCH --output=star_align_S1.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/rnaseq/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/rnaseq/lib:${LD_LIBRARY_PATH}"

################
## Variables ###
################

path1="/mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/"
index="L1_index"
genotype="L1"
tissue="Bud Leaf Root Field1 Field2"
t1="${path1}/trimmed.1.fastq.gz"
t2="${path1}/trimmed.2.fastq.gz"
output1="rnaseq1"
output2="rnaseq2"

#Submit in ${genotype}_tissue

################     
### Comand  ####
################

#Star 1st pass
echo "Running star 1st pass"
for i in ${tissue}
do
	echo "Running sample: $i"
	mkdir $i/$output1 $i/$output2
	cd $i/$output1
	
STAR \
--genomeDir ${path1}/${genotype}_index_star/ \
--runThreadN 6 \
--runMode alignReads \
--readFilesIn /mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/${genotype}_tissues/${i}/trimmed.1.fastq.gz /mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/${genotype}_tissues/${i}/trimmed.2.fastq.gz \
--readFilesCommand zcat \
--sjdbGTFfile $path1/${genotype}_tissues/${genotype}_index/${genotype}-v1.2.gtf \
--outFilterMultimapNmax 20 
--outFileNamePrefix out_$i_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
--outSAMstrandField intronMotif \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04
junctions="../../"$i"/"$output1"/SJ.out.tab $junctions"
	cd ../../
done
echo $junctions

#Star 2nd pass
echo "Running star 2nd pass"
for i in $tissue
do
	echo "Running sample: $i"
	cd $i/$output2
	STAR \
		--runThreadN 6 \
		--runMode alignReads \
		--genomeDir ${path1}/${genotype}_index_star/ \
		--readFilesIn /mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/${genotype}_tissues/${i}/trimmed.1.fastq.gz /mnt/gs21/scratch/kollarle/mimulus-assembly/data/Mguttatus/final/final_pseudo_assemblies/rnaseq_alignments/gould_et_al_2018/${genotype}_tissues/${i}/trimmed.2.fastq.gz \
		--sjdbFileChrStartEnd $junctions \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandField intronMotif \
		--outFilterType BySJout \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--quantMode GeneCounts
	cut -f1,4 ReadsPerGene.out.tab | sed '1,4d' > "$i"_counts.tsv 
	cd ../../
done

echo "Done"

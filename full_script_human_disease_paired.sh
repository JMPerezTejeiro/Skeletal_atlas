#!/usr/bin/env bash
# Leave only one comment symbol on selected options
# Those with two commets will be ignored:
# The name to show in queue lists for this job:
##SBATCH -J fastqc

# Number of desired cpus (can be in any node):
#SBATCH --ntasks=1

# Number of desired cpus (all in same node):
#SBATCH --cpus-per-task=8

# Amount of RAM needed for this job:
#SBATCH --mem=64gb

# The time the job will be running:
#######SBATCH --time=72:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features uncomment the desired constraint line:
# * to request only the machines with 80 cores and 2TB of RAM
##SBATCH --constraint=bigmem
# * to request only machines with 16 cores and 64GB with InfiniBand network
#SBATCH --constraint=cal
# * to request only machines with 24 cores and Gigabit network
##SBATCH --constraint=dx
##SBATCH --constraint=ssd

# Set output and error files
#SBATCH --error=qc.%J.err
#SBATCH --output=qc.%J.out

# Leave one comment in following line to make an array job. Then N jobs will be launched. In each one SLURM_ARRAY_TASK_ID will take one value from 1 to 100
##SBATCH --array=1-100

# To load some software (you can show the list with 'module avail'):

###################################  HUMAN DISEASES PAIRED  ################################################################################################

############################# FastQC of RAW DATA  ##############################################################################################################

# module load software
module purge
module load fastqc/0.11.9

fastqc ../00_raw_data/human_disease/hd_paired/*.fastq.gz -o ../01_raw_fastQC/human_disease/hd_paired/ -t 16
#multiqc .


###################################################### TRIMMING THE DATA WITH TRIMMOMATIC ########################################################################

# module load software
module purge
module load trimmomatic/0.39

# the program to execute with its parameters:
for seq1 in ../00_raw_data/human_disease/hd_paired/*_1.fastq.gz
do
  echo $seq1

  seq2=${seq1/_1\.fastq\.gz/_2\.fastq\.gz}

  echo $seq2

  seq=${seq1/00_raw_data/02_trimmed_data}
  seq_r1=${seq/_1\.fastq\.gz/_1_trimmed\.fastq\.gz}
  seq_r1u=${seq/_1\.fastq\.gz/_1_trim_unpaired\.fastq\.gz}
  seq_r2=${seq/_1\.fastq\.gz/_2_trimmed\.fastq\.gz}
  seq_r2u=${seq/_1\.fastq\.gz/_2_trim_unpaired\.fastq\.gz}

  java -jar /mnt/home/soft/trimmomatic/programs/x86_64/0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 $seq1 $seq2 $seq_r1 $seq_r1u $seq_r2 $seq_r2u ILLUMINACLIP:/mnt/home/users/bio_369_uma/jmperez/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

  done

############################################################################## FastQC OF TRIMMED DATA #################################################################################################

# module load software
module purge
module load fastqc/0.11.9

fastqc ../02_trimmed_data/human_disease/hd_paired/*.fastq.gz -o ../03_trimmed_fastQC/human_disease/hd_paired/ -t 16


################################################################ MAPPING THE READS WITH HISAT2 AND OBTAINING TPM MATRIX WITH FEATURE COUNTS FROM SUBREAD PACKAGE #########################################################

# module load software
module unload gcc
module load hisat/2.2.1
module unload gcc
module load samtools/1.16
module unload gcc
module load subread/2.0.3

#Variable definition for Hisat2
for seq1 in ../02_trimmed_data/human_disease/hd_paired/*_1_trimmed.fastq.gz
do

echo $seq1
seq2=${seq1/_1_trimmed\.fastq\.gz/_2_trimmed\.fastq\.gz}
echo $seq2

hisat_output1=${seq1/_1_trimmed\.fastq\.gz/_s.bam}
hisat_report1=${seq1/_1_trimmed\.fastq\.gz/.txt}
hisat_output=${hisat_output1/02_trimmed_data/04_hisat2_output/}
hisat_report=${hisat_report1/02_trimmed_data/04_hisat2_output/}

#Allignment

hisat2 -p 16 --summary $hisat_report --dta -x /mnt/home/users/bio_369_uma/jmperez/ref_genome/human_ENS_index/human_ENS_index -1 $seq1 -2 $seq2 | samtools view -F 0x4 -b - | samtools sort -o $hisat_output

done
###################### featureCOUNTS #############################

#Ejecutar featureCounts:

for sample in ../04_hisat2_output/human_disease/hd_paired/*_s.bam;

do

echo $sample

sample1=${sample/_s.bam/.txt}
sample_output=${sample1/04_hisat2_output/05_TPMs_matrix}

featureCounts -p --countReadPairs -C -T 16 -a /mnt/home/users/bio_369_uma/jmperez/ref_genome/Homo_sapiens.GRCh38.109_ENSEMBL.gtf -t exon -g gene_id -o $sample_output $sample
#rm $sample

done

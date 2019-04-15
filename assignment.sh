###########################################
#  - NGS1 Course - Assignment             #
#  - Bash Script                          #
#  - April 15,2019                        #
#  - Copyright: Mohamed AboelEla          #
#  - Nile University                      #
###########################################

#!/bin/bash

### 01.Downloading data file ####
#Creatinging base working directory and downloading data
mkdir ngs_assignment && cd ngs_assignment
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

### 02.Preparing the data ###
#Converting SRA file to FastQ file and spliting R1 and R2 using SRA Toolkit
fastq-dump -I --split-files SRR8797509.sra

#Creating specific directores for shuffled and unshuffled fastq files
mkdir -p unshuffled/ shuffled/

#Moving and renaming R1 and R2 files
mv SRR8797509_1.fastq unshuffled/SRR8797509_1_unshuffled.fastq
mv SRR8797509_2.fastq unshuffled/SRR8797509_2_unshuffled.fastq

#Shuffling R1 and R2 files
seqkit shuffle --threads 6 unshuffled/SRR8797509_1_unshuffled.fastq > shuffled/SRR8797509_1_shuffled.fastq
seqkit shuffle --threads 6 unshuffled/SRR8797509_2_unshuffled.fastq > shuffled/SRR8797509_2_shuffled.fastq

#Splitting 1M reads from each file
seqkit head -n 5000000 unshuffled/SRR8797509_1_unshuffled.fastq | seqkit split2 -s 1000000 -O unshuffled/r1/
seqkit head -n 5000000 unshuffled/SRR8797509_2_unshuffled.fastq | seqkit split2 -s 1000000 -O unshuffled/r2/
seqkit head -n 5000000 shuffled/SRR8797509_1_shuffled.fastq | seqkit split2 -s 1000000 -O shuffled/r1/
seqkit head -n 5000000 shuffled/SRR8797509_2_shuffled.fastq | seqkit split2 -s 1000000 -O shuffled/r2/

#Renaming sample files
for dir in shuffled unshuffled;
do
	for dirc in r1 r2;
	do
		for i in {1..5};
		do
			mv $dir/$dirc/stdin.part_00$i.fastq  $dir/$dirc/s$i-$dirc-$dir.fastq
		done
	done
done

### 03.FASTQ Quality Control ###
#Creating new output directory for FastQC
mkdir fastqc/
for dir in shuffled unshuffled;
do
	for dirc in r1 r2;
	do
		fastqc -o fastqc/ -t 1 -f fastq -noextract $dir/$dirc/s1*.fastq;
	done
done
#Merging output reports into one report
multiqc -z fastqc/

### 04.Trimming ###
#Creating directories for trimmming
for dir in unshuffled shuffled;
do
	for dirc in r1 r2;
	do
		mkdir $dir/$dirc/trimmed/
	done
done

#Applying mild trimming for unshuffled samples
for i in {1..5};
do
	f1="unshuffled/r1/s$i*.fastq"
	f2="unshuffled/r2/s$i*.fastq"
	nf1="unshuffled/r1/s$i-r1-unshuffled-pe-trim.fastq"
	nf2="unshuffled/r2/s$i-r2-unshuffled-pe-trim.fastq"
	nf1u="unshuffled/r1/s$i-r1-unshuffled-se-trim.fastq"
	nf2u="unshuffled/r2/s$i-r2-unshuffled-se-trim.fastq"
	trimmomatic PE -threads 12 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $nf1 $nf1u $nf2 $nf2u  SLIDINGWINDOW:4:20 MINLEN:36
done

#Applying aggressive trimming for shuffled samples
for i in {1..5};
do
	f1="shuffled/r1/s$i*.fastq"
	f2="shuffled/r2/s$i*.fastq"
        nf1="shuffled/r1/s$i-r1-shuffled-pe-trim.fastq"
        nf2="shuffled/r2/s$i-r2-shuffled-pe-trim.fastq"
        nf1u="shuffled/r1/s$i-r1-shuffled-se-trim.fastq"
        nf2u="shuffled/r2/s$i-r2-shuffled-se-trim.fastq"
        trimmomatic PE -threads 12 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $nf1 $nf1u $nf2 $nf2u  SLIDINGWINDOW:4:30 MINLEN:36
done

### 05.Alignment ###
#Downloading reference genome
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

#Aligning unshuffled data with BWA
#Generating Reference Genome Index
bwa index -a bwtsw Homo_sapiens.GRCh38.dna.chromosome.22.fa
#Aligning reads to reference
for i in {1..5};
do
	r1="unshuffled/r1/s$i-r1-unshuffled-pe-trim.fastq"
	r2="unshuffled/r2/s$i-r2-unshuffled-pe-trim.fastq"
	bwa mem Homo_sapiens.GRCh38.dna.chromosome.22.fa $r1 $r2 > s$i-unshuffled-pe-trim.sam
done

#Reporting results stat
for i in {1..5};
do
	samtools flagstat s$i-unshuffled-pe-trim.sam > s$i-unshuffled-pe-trim.stat
done

#Aligning shuffled data with HISAT
#Generating Reference Genome Index
hisat2-build  Homo_sapiens.GRCh38.dna.chromosome.22.fa  Homo_sapiens.GRCh38.dna.chromosome.22.ht2

#Aligning reads to reference
for i in {1..5};
do
	r1="shuffled/r1/s$i-r1-shuffled-pe-trim.fastq"
	r2="shuffled/r2/s$i-r2-shuffled-pe-trim.fastq"
	hisat2 -q -x Homo_sapiens.GRCh38.dna.chromosome.22.ht2 -1$r1 -2$r2 -S s$i-shuffled-pe-trim.sam
done

#Reporting results stat
for i in {1..5};
do
	samtools flagstat s$i-shuffled-pe-trim.sam > s$i-shuffled-pe-trim.stat
done

### 06.Assembly ###
#Converting SAM to BAM
for i in unshuffled shuffled;
do
	for x in {1..5};
	do
		samtools view -bS s$x-$i-pe-trim.sam > s$x-$i-pe-trim.bam
		samtools sort s$x-$i-pe-trim.bam -o s$x-$i-pe-trim.sorted.bam
	done
done

#Assembly without known annotation
for i in unshuffled shuffled;
do
	for x in {1..5};
	do
		stringtie s$x-$i-pe-trim.sorted.bam --rf -l ref_free -o s$x-$i-pe-trim.sorted.gtf
	done
done

#Assembly with known annotation
#Downloading GFF3 annotation file
wget ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.chromosome.22.gff3.gz
gunzip Homo_sapiens.GRCh38.96.chromosome.22.gff3.gz

for i in unshuffled shuffled;
do
       	for x in {1..5};
	do
		stringtie s$x-$i-pe-trim.sorted.bam --rf -l ref_sup -G Homo_sapiens.GRCh38.96.chromosome.22.gff3 -o s$x-$i-pe-trim.sorted_r.gtf
	done
done



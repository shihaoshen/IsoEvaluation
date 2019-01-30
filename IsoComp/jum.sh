#!/bin/bash

# qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-1:1 jum.sh

# -- our name ---
#$ -N jum.sh
#$ -S /bin/bash
#$ -V

conda activate turbo

od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output/
bin=/mnt/isilon/xing_lab/shens/bin/github/JUM_2.0.2/

#g1=PC3E-1Aligned.sort.bam.fasta.Aligned.out_sorted.bam,PC3E-2Aligned.sort.bam.fasta.Aligned.out_sorted.bam,PC3E-3Aligned.sort.bam.fasta.Aligned.out_sorted.bam
#g2=GS689.LI-1Aligned.sort.bam.fasta.Aligned.out_sorted.bam,GS689.LI-2Aligned.sort.bam.fasta.Aligned.out_sorted.bam,GS689.LI-3Aligned.sort.bam.fasta.Aligned.out_sorted.bam

g1=PC3E-1Aligned.sort.bam.fasta.,PC3E-2Aligned.sort.bam.fasta.,PC3E-3Aligned.sort.bam.fasta.
g2=GS689.LI-1Aligned.sort.bam.fasta.,GS689.LI-2Aligned.sort.bam.fasta.,GS689.LI-3Aligned.sort.bam.fasta.


#JUM Example
#/user/home/JUM_2.0.2/JUM_A.sh --Folder "directory" --JuncThreshold "junction_read_threshold" --Condition1_fileNum_threshold "thre_file_num_1" --Condition2_fileNum_threshold "thre_file_num_2" --IRthreshold "IR_read_threshold" --Readlength "read_length" --Thread "thread_num" --Condition1SampleName "condition_1sample" --Condition2SampleName "condition_2sample"
cp $bin/* $od/
cd $od/
./JUM_A.sh --Folder /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output/ --JuncThreshold 5 --Condition1_fileNum_threshold 2 --Condition2_fileNum_threshold 2 --IRthreshold 5 --Readlength 100 --Thread 3 --Condition1SampleName $g1 --Condition2SampleName $g2

cd $od/JUM_diff/
Rscript $bin/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout

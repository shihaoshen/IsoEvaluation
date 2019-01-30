#!/bin/bash

# qsub -pe smp 3 -l h_vmem=20G -l m_mem_free=20G -t 1-6:1 jum_star_2ndpass.sh

# -- our name ---
#$ -N jum_star_2ndpass.sh
#$ -S /bin/bash
#$ -V

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation
list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list.txt
list2=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list2.txt
od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output/
bin=/home/shens/xing_lab/bin/github/STAR/bin/Linux_x86_64_static/

#JUM Example
#STAR --runThreadN 3 --genomeDir genome_index_STAR_r1 --outFileNamePrefix ctrl_1 --readFilesIn ctrl_1_R01.fastq ctrl_1_R02.fastq --outSJfilterReads Unique --outSAMstrandField intronMotif --outFilterMultimapNmax 1 -sjdbFileChrStartEnd 1st_SJ/ctrl_1SJ.out.tab 1st_SJ/ctrl_2SJ.out.tab 1st_SJ/ctrl_3SJ.out.tab 1st_SJ/treat_1SJ.out.tab 1st_SJ/treat_2SJ.out.tab 1st_SJ/treat_3SJ.out.tab

echo 'Star run. Output:'
echo $od
$bin/STAR --runThreadN 3 --genomeDir /mnt/isilon/xing_lab/shens/Annotation/hg19_genome/ --outFileNamePrefix $od/`sed -n ${SGE_TASK_ID}p $list2`. --readFilesIn `sed -n ${SGE_TASK_ID}p $list` --outSJfilterReads Unique --outSAMstrandField intronMotif --outFilterMultimapNmax 1 -sjdbFileChrStartEnd $od/1st_SJ/GS689.LI-1.fasta.SJ.out.tab $od/1st_SJ/GS689.LI-1.fasta.SJ.out.tab $od/1st_SJ/GS689.LI-2Aligned.sort.bam.fasta.SJ.out.tab $od/1st_SJ/GS689.LI-3Aligned.sort.bam.fasta.SJ.out.tab $od/1st_SJ/PC3E-1.fasta.SJ.out.tab $od/1st_SJ/PC3E-2Aligned.sort.bam.fasta.SJ.out.tab $od/1st_SJ/PC3E-3Aligned.sort.bam.fasta.SJ.out.tab

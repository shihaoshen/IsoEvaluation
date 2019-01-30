#!/bin/bash

# qsub -pe smp 3 -l h_vmem=20G -l m_mem_free=20G -t 1-6:1 jum_star_2ndpass_clean.sh

# -- our name ---
#$ -N jum_star_2ndpass_clean.sh
#$ -S /bin/bash
#$ -V

/home/shens/.bash_profile

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation
list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list.txt
list2=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list2.txt
od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output/
bin=/home/shens/xing_lab/bin/github/STAR/bin/Linux_x86_64_static/

samtools view -bS $od/`sed -n ${SGE_TASK_ID}p $list2`.Aligned.out.sam > $od/`sed -n ${SGE_TASK_ID}p $list2`.Aligned.out.bam
samtools sort -o $od/`sed -n ${SGE_TASK_ID}p $list2`.Aligned.out_sorted.bam -T $od/`sed -n ${SGE_TASK_ID}p $list2`_temp $od/`sed -n ${SGE_TASK_ID}p $list2`.Aligned.out.bam 
samtools index $od/`sed -n ${SGE_TASK_ID}p $list2`.Aligned.out_sorted.bam

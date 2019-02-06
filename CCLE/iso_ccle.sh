#!/bin/bash

# qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-2:1 iso_ccle.sh
# qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 5-7:1 iso_ccle.sh

# -- our name ---
#$ -N iso_ccle.sh
#$ -S /bin/bash
#$ -V

conda activate iso

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/em_stat/

list=/mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/em_stat/tissue.list
gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
od=/mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/em_stat

python /mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/rMATS-ISO.py stat --bam $od/`sed -n ${SGE_TASK_ID}p $list`/input.list -o $od/`sed -n ${SGE_TASK_ID}p $list` >iso_stat_${SGE_TASK_ID}_log.txt 2>iso_stat_${SGE_TASK_ID}_log_error.txt

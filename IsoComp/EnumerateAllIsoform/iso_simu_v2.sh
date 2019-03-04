#!/bin/bash

# qsub -pe smp 1 -l h_vmem=80G -l m_mem_free=80G -t 1-1:1 iso_simu_v2.sh

# -- our name ---
#$ -N iso_simu_v2.sh
#$ -S /bin/bash
#$ -V

conda activate iso

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation

list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_input.list
gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
bin=/mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/
log=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/log
od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output_r2/

thisdate="$(date +"%M-%H-%m-%d-%Y")"
python $bin/rMATS-ISO.py module --gtf $gtf --bam $list -o $od >$log/iso_log_$thisdate.txt 2>$log/iso_log_error_$thisdate.txt
python $bin/rMATS-ISO.py stat --psi-low 0.01 --iso-thresh 10 --bam $list -o $od >$log/iso_stat_log_$thisdate.txt 2>$log/iso_stat_log_error_$thisdate.txt

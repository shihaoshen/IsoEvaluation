#!/bin/bash

# qsub -l h_vmem=20G -l m_mem_free=20G -t 1-1:1 iso.sh

# -- our name ---
#$ -N iso.sh
#$ -S /bin/bash
#$ -V

conda activate iso

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation
bin=/mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO
od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation

echo 'Start iso run. Output:'
echo $od/iso_simu_output/
python $bin/rMATS-ISO.py stat --bam $od/simu_input.list -o $od/iso_simu_output/


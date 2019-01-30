#!/bin/bash

# qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 2-3:1 GSpsiFromPro.sh

# -- our name ---
#$ -N GSpsiFromPro.sh
#$ -S /bin/bash
#$ -V

conda activate turbo

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/truepsi_iso/
list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list.txt
list2=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_fasta_list2.txt
od=/mnt/isilon/xing_lab/shens/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/od/
bin=/home/shens/xing_lab/bin/github/STAR/bin/Linux_x86_64_static/
gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
pro=/mnt/isilon/xing_lab/shens/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/pro/
psi=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/truepsi_turbo/
iso=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output/ISO_module/

#JUM Example
#python calc_psi_keycoord.py SE rMATS-ISO_Flux_Simulation/od/SE.MATS.JCEC.txt $GTF rMATS-ISO_Flux_Simulation/pro/PC3E-1Aligned.sort.bam.pro rMATS-ISO_Flux_Simulation/rMATS_psi/PC3E-1Aligned.sort.bam.psi
echo 'Calculate true psi for turbo SE:'
python -O GSpsiFromPro.py $pro/PC3E-${SGE_TASK_ID}Aligned.sort.bam.pro $gtf $iso/GS689.LI-1_STARAligned.out.sort.bam.IsoExon > PC3E-${SGE_TASK_ID}Aligned.sort.bam.psi
python -O GSpsiFromPro.py $pro/GS689.LI-${SGE_TASK_ID}Aligned.sort.bam.pro $gtf $iso/GS689.LI-1_STARAligned.out.sort.bam.IsoExon > GS689.LI-${SGE_TASK_ID}Aligned.sort.bam.psi

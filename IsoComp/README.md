# IsoEvaluation
## Comparison of ISO with JUM, majiq, leafcutter in Flux simulated data

#CHOP Folder location
    
    /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation

#Run rMATS-ISO 
        
    ls /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_reads/*.bam > simu_input.list
    list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_input.list
    gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
    python /mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/rMATS-ISO.py module --gtf $gtf --bam $list -o iso_simu_output/ >iso_log.txt 2>iso_log_error.txt
    python /mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/rMATS-ISO.py stat --bam $list -o iso_simu_output/ >iso_stat_log.txt 2>iso_stat_log_error.txt
    
#Run rMATS with splicing difference c cutoff = 1% (gold standard null)
  
    GTF=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
    TMP=/mnt/isilon/xing_lab/shens/ISO_ReRun/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/tmp
    OD=/mnt/isilon/xing_lab/shens/ISO_ReRun/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/od
    BAM=/mnt/isilon/xing_lab/shens/ISO_ReRun/rMATS-iso-bak/Simulation/
    python /mnt/isilon/xing_lab/shens/bin/github/rmats_pipeline/rmats.py --b1 $BAM/bam1.list --b2 $BAM/bam2.list --gtf $GTF --od $OD --tmp $TMP --readLength 101 --cstat 0.01 --task both

#Run leafcutter (V0.2)

    conda activate turbo
    #Convert bam to junction files
    export PATH=$PATH:/mnt/isilon/xing_lab/shens/bin/github/leafcutter/scripts/
    if [ -e test_juncfiles.txt ]; then rm test_juncfiles.txt; fi
    for bamfile in `ls simu_reads/*.bam`; do
        echo Converting $bamfile to $bamfile.junc
        sh /mnt/isilon/xing_lab/shens/bin/github/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc
        echo $bamfile.junc >> test_juncfiles.txt
    done
    # Finds intron clusters and quantifies junction usage within them.
    python /mnt/isilon/xing_lab/shens/bin/github/leafcutter/clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o leafcutter_output -l 500000
    # Differential splicing analysis.
    /mnt/isilon/xing_lab/shens/bin/github/leafcutter/scripts/leafcutter_ds.R --num_threads 4 leafcutter_output/leafcutter_output_perind_numers.counts.gz simu_reads/groups_file.txt -i 3
    
#Run JUM (V2.0.2)

    qsub -pe smp 3 -l h_vmem=20G -l m_mem_free=20G -t 1-6:1 jum_star.sh
    ./jum_star_clean.sh
    qsub -pe smp 3 -l h_vmem=20G -l m_mem_free=20G -t 1-6:1 jum_star_2ndpass.sh
    ./jum_star_2ndpass_clean.sh
    qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-1:1 jum.sh

#Run majiq (V1.1.7a) and voila (V1.1.9)
#1. By default majiq only outputs the events likely to have splicing changes, I add 'show-all' option to output all events.
#2. majiq default posterior probability is for P(|delta psi|>0.2), I change it to P(|delta psi|>0.1).

    majiq build ensembl.hg19.gff3 -c settings.ini -j 1 -o output/
    majiq psi output/PC3E-1_STARAligned.out.sort.majiq output/PC3E-2_STARAligned.out.sort.majiq output/PC3E-3_STARAligned.out.sort.majiq output/GS689.LI-1_STARAligned.out.sort.majiq output/GS689.LI-2_STARAligned.out.sort.majiq output/GS689.LI-3_STARAligned.out.sort.majiq -j 2 -o output/ -n g
    majiq deltapsi -grp1 output/PC3E-1_STARAligned.out.sort.majiq output/PC3E-2_STARAligned.out.sort.majiq output/PC3E-3_STARAligned.out.sort.majiq -grp2 output/GS689.LI-1_STARAligned.out.sort.majiq output/GS689.LI-2_STARAligned.out.sort.majiq output/GS689.LI-3_STARAligned.out.sort.majiq -j 2 -o output/ -n pc3e gs689 --default-prior
    voila deltapsi output/pc3e_gs689.deltapsi.voila -s output/splicegraph.sql -o output_showall_p01/ --show_all --threshold 0.1
    
#Calculate true psi for ISO modules from Flux simulation parameters
#Revision 02/07/19: Require the total read counts of all isoforms that are averaged among the replicates to be greater than 10

    cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/truepsi_iso
    qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-3:1 GSpsiFromPro.sh
    #get glmm logistic regression p value for the isoform with largest difference
    python isopsidiff.py isopsidiff10glmm.txt 1
    #get t-test p value for Flux psi estimates
    python isopsidiff_ttest.py isopsidiff10_ttest.txt 1

#Match turbo, ISO, majiq, leafcutter and jum

    #match the ISO and rMATS results in the simulation data
    #this has the min rMATS p value of the module
    python findASM.py events.txt eventsp.txt /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output/ISO_module/GS689.LI-1_STARAligned.out.sort.bam.IsoExon ASM2SE.112319.txt

    #match the ISO and majiq results in the simulation data
    python findASM_majiq.py events_majiq_p001.txt /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output/ISO_module/GS689.LI-1_STARAligned.out.sort.bam.IsoExon ASM2majiq_p001.txt

    #match the ISO and leafcutter results in the simulation data
    python findASM_leafcutter.py events_leafcutter.txt eventsp_leafcutter.txt /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output/ISO_module/GS689.LI-1_STARAligned.out.sort.bam.IsoExon ASM2leafcutter.102218.txt

    #match the ISO and leafcutter results in the simulation data
    python findASM_jum.py events_jum.txt /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output/ISO_module/GS689.LI-1_STARAligned.out.sort.bam.IsoExon ASM2jum.102218.txt

Roc of all modules, modules of SE events, and modules of ASS events. 

1. Gold standard H1: |delta psi| > 10%. Gold standard H0: |delta psi|<1%. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

    Rscript roc.R


[ROC with majiq delta psi 10%, rMATS-turbo JC, Flux delta psi gold standard](https://drive.google.com/file/d/1CEr6e8QwDqw-ue2nIbk0xBRx4-W_zl-b/view?usp=sharing)

[ROC with majiq delta psi 1%, rMATS-turbo JC, Flux delta psi gold standard](https://drive.google.com/open?id=1hzNcLh1bWaiJ_2FTH27-lp1JsRrLJ3Nw)

2. Gold standard H1: GLMM logistic of Flux counts P < 0.005. Gold standard H0: logistic P > 0.5, |delta psi|<1%. GLMM calculated on the isoform with largest Flux psi difference. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

    Rscript roc_logit.R

[ROC with majiq delta psi 10%, rMATS-turbo JC, logistic gold standard](https://drive.google.com/open?id=1zDU43xx0j9QB8hLYCxRIz0fA-pUf9cJM)

[ROC with majiq delta psi 1%, rMATS-turbo JC, logistic gold standard](https://drive.google.com/open?id=1O_dsLrKfOh4ncP6wwREz3PdJvcGlArWR)

3. Gold standard H1: T-test of Flux psi P < 0.05. Gold standard H0: T-test P > 0.5, |delta psi|<1%. T-test calculated on the isoform with largest Flux psi difference. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

    Rscript roc_ttest.R

[ROC with majiq delta psi 10%, rMATS-turbo JC, ttest gold standard](https://drive.google.com/file/d/17pIHZw17XgQV1Irxp9L-LxN_uz_eIzBS/view?usp=sharing)

[ROC with majiq delta psi 1%, rMATS-turbo JC, ttest gold standard](https://drive.google.com/open?id=1dZixAxbma8ivQaZdDUp_W3Ply78V6VGr)


![](https://github.com/shihaoshen/IsoEvaluation/blob/master/IsoComp/rMATS_ISO_Flux_roc.png)

![](https://github.com/shihaoshen/IsoEvaluation/blob/master/IsoComp/rMATS_ISO_Flux_roc_SE.png)

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
    
#Run rMATS 
  
    GTF=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
    TMP=/mnt/isilon/xing_lab/shens/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/tmp
    OD=/mnt/isilon/xing_lab/shens/rMATS-iso-bak/Simulation/rMATS-ISO_Flux_Simulation/od
    BAM=/mnt/isilon/xing_lab/shens/rMATS-iso-bak/Simulation/
    python /mnt/isilon/xing_lab/shens/bin/github/rmats_pipeline/rmats.py --b1 $BAM/bam1.list --b2 $BAM/bam2.list --gtf $GTF --od $OD -t paired --readLength 76 --cstat 0.0001 --libType fr-unstranded

#Run leafcutter

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

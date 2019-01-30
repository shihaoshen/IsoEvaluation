# IsoEvaluation

#CHOP Folder location
    
    /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation

#Run rMATS-ISO 
        
    ls /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_reads/*.bam > simu_input.list
    list=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/simu_input.list
    gtf=/mnt/isilon/xing_lab/shens/Annotation/Homo_sapiens.GRCh37.87.new.gtf
    python /mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/rMATS-ISO.py module --gtf $gtf --bam $list -o iso_simu_output/ >iso_log.txt 2>iso_log_error.txt
    python /mnt/isilon/xing_lab/shens/bin/github/rMATS-ISO/rMATS-ISO.py stat --bam $list -o iso_simu_output/ >iso_stat_log.txt 2>iso_stat_log_error.txt

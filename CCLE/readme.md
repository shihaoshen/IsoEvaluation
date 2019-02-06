# ISO runs on 377 CCLE samples

    /mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE

#ISO module detection runs on CCLE samples

    qsub -l h_vmem=60G -l m_mem_free=60G iso.sh
    
#Prepare the samples for statistical test between EM samples

    cd /mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/em_stat
    
    #EM sample information is available at /mnt/isilon/xing_lab/shens/CCLE
    #get the list of tissues with at least 5 E and 5 M samples
    Rscript prepEMT.R
    
    [1] "breast"
     E  M
    36  8
    [1] "endometrium"
     E  M
    18  8
    [1] "lung"
     E  M
    74 89
    [1] "ovary"
     E  M
    26 16
    [1] "pancreas"
     E  M
    31  9
    [1] "stomach"
     E  M
    27 10
    [1] "urinary_tract"
     E  M
    16  8

# Copy the ISO samples into each folder
    python sampleEMT.py
    
# ISO statistical analysis between E and M samples in 7 tissues
    qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-7:1 iso_ccle.sh

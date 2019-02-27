# rMATS-ISO run on real PC3E data

/mnt/isilon/xing_lab/shens/ISO_ReRun/PC3E

#Run ISO on the real PC3E data, using stock GTF
```
qsub -pe smp 1 -l h_vmem=80G -l m_mem_free=80G -t 1-1:1 iso.sh
```

#Identify module type (SE or ASS) based on GTF

    qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-3:1 GSpsiFromPro2.sh

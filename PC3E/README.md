# rMATS-ISO run on real PC3E data

/mnt/isilon/xing_lab/shens/ISO_ReRun/PC3E

```
qsub -pe smp 1 -l h_vmem=80G -l m_mem_free=80G -t 1-1:1 iso.sh
```

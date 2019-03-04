#update ISO results when enumerating all isoforms (having isoforms cut from IsoExon output leads to incomplete isoforms);

qsub -pe smp 1 -l h_vmem=80G -l m_mem_free=80G -t 1-1:1 iso_simu.sh

#Update the true psi values based on ISO

cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/truepsi_iso_allisoforms
qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-3:1 GSpsiFromPro2.sh
python isopsidiff.py isopsidiff.txt
python isopsidiff_ttest.py isopsidiff_ttest.txt 1
python isopsidiff_logit.py isopsidiff_logit.txt 1
cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/

#Match rMATS-turbo, majiq, leafcutter and jum results to ISO modules
cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_roc/
isoexon=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/iso_simu_output_noveljuc_allisoforms/ISO_module/GS689.LI-1_STARAligned.out.sort.bam.IsoExon
python findASM.py events.txt eventsp.JC.txt $isoexon ASM2SE.JC.txt
python findASM_majiq.py events_majiq_p001.txt $isoexon ASM2majiq_p001.txt
python findASM_leafcutter.py events_leafcutter.txt eventsp_leafcutter.txt $isoexon ASM2leafcutter.txt
python findASM_jum.py events_jum.txt $isoexon ASM2jum.txt

Roc of all modules, modules of SE events, and modules of ASS events. 

1. Gold standard H1: |delta psi| > 10%. Gold standard H0: |delta psi|<1%. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

        Rscript roc.R


[ROC with majiq delta psi 1%, rMATS-turbo JC delta 1%, Flux delta psi gold standard](https://drive.google.com/open?id=1uRmRe-XVNnQcmSBlaBiBJZFYptsq3smg)

2. Gold standard H1: GLMM logistic of Flux counts P < 0.005. Gold standard H0: logistic P > 0.5, |delta psi|<1%. GLMM calculated on the isoform with largest Flux psi difference. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

        Rscript roc_logit.R

[ROC with majiq delta psi 1%, rMATS-turbo JC delta 1%, logistic gold standard](https://drive.google.com/open?id=1UgVUlBScd97POvSeBbXBNKC6tI__Jhia)

3. Gold standard H1: T-test of Flux psi P < 0.05. Gold standard H0: T-test P > 0.5, |delta psi|<1%. T-test calculated on the isoform with largest Flux psi difference. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

        Rscript roc_ttest.R

[ROC with majiq delta psi 1%, rMATS-turbo JC delta 1%, t-test gold standard](https://drive.google.com/open?id=16jIRMUjSKBwQbXO0ruO_oe1krx_0UkOo)


4. Gold standard H1: T-test of Flux psi P < 0.05, |delta psi|>5% or 10% (two sets). Gold standard H0: T-test P > 0.5, |delta psi|<1%. T-test calculated on the isoform with largest Flux psi difference. Gold standard required sum of module counts of all transcript averaged over all replicates to be >=10.

        Rscript roc_ttest_deltapsi5.R
        Rscript roc_ttest_deltapsi10.R

[ROC with majiq delta psi 1%, rMATS-turbo JC delta 1%, DeltaPsi H1 5% + t-test gold standard](https://drive.google.com/open?id=1YUcBgKqCL9dpe1J2ITmJru_DGMH_yYqw)

[ROC with majiq delta psi 1%, rMATS-turbo JC delta 1%, DeltaPsi H1 10% + t-test gold standard](https://drive.google.com/open?id=1Fiv1olgv4CUxV2LeH-E0dqf7eaxzZ4Gw)


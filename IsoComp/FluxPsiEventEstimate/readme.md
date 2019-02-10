#This script outputs the psi estimate from Flux simulation and check whether the module is SE or ASS by GTF annotation

    qsub -pe smp 1 -l h_vmem=20G -l m_mem_free=20G -t 1-3:1 GSpsiFromPro2.sh

#Details:

#For each module, the script first finds all the gtf transcripts overlapped with the coordinate range of module

#1. 

#For psi value estimate, the script requires the gtf transcript to be 100% match with all the coordinates of at least of module isoform inside the range of the module. 

#Since one isoform can have multiple transcripts to be 100% match inside the module range, all matched transcripts' Flux molecule numbers are summed up as the molecule numbers of the isoform.

#The gold standard isoform psi value is estimated by dividing one isoform's molecule numbers by numbers from all isoforms.

#2.

#For identification of simple event (SE or ASS), the script finds all the gtf transcripts with at least one exon overlapped inside the module range.

#If multiple gtf transcript annotations have 100% identical exon cooridnates inside the module, these identical transcripts are merged.

#After merging identical transcript inside the module, if there are 2 and only 2 unique transcripts inside the module, the scripts continue to check whether it is SE or ASS.

#SE events have only one exon different between the two transcripts. The exon is a middle exon. Both of its splice sites are absent in one transcript and present in another.

#ASS events also has only one exon different between the two transcripts. The exon has one splice site in both transcript and another splice site different between the two transcripts.

#Note that the transcript is NOT required to be 100% matched with any of the rMATS-ISO module isoforms, since rMATS-ISO can miss detection of certain isoforms due to lack of reads.

#The identification of SE or ASS is completely determined by GTF transcript annotation over the module region. 

#If the GTF indicates more than 2 transcripts with different exon coordinates inside the module, this will not be identified as either SE or ASS, despite the number of rMATS-ISO module isoforms.

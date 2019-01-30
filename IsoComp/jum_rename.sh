cd /mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output
for file in PC3E-1.fasta*; do
    mv "$file" "$(basename "$file" PC3E-1.fasta)PC3E-2Aligned.sort.bam.fasta"
done

# Align
for filename in *.fastq; do
/bigdata/messaoudilab/arivera/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --genomeLoad NoSharedMemory --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --readFilesCommand cat --genomeDir /big$
done

samtools view -h -q 255 file.fastqAligned.sortedByCoord.out.bam > test.bam

# Align
for i in *.read1_val_1.fq.gz; do
        /bigdata/messaoudilab/arivera/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --genomeLoad NoSharedMemory --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 255 --readFilesCommand zcat --genomeDir /bigdata/messaoudilab/arivera/viral_genomes/Rhesus_SVV/STAR_genome --runThreadN 10 --readFilesIn $i ${i%.read1_val_1.fq.gz}.read2_val_2.fq.gz --outFileNamePrefix ${i%.read1_val_1.fq.gz}
done

samtools view -h -q 255 file.fastqAligned.sortedByCoord.out.bam > test.bam

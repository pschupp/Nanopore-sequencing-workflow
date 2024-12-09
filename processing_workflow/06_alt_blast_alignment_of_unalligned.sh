./makeblastdb -in ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa -dbtype nucl 

/home/patrick/bin/blast/ncbi-blast-2.7.1+/bin/blastn -query ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.fa -task blastn -db ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa  -out ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.blast.sam  -word_size 6 -gapopen 0 -gapextend 2 -penalty -1 -reward 1 -outfmt 5 -num_threads 6 -evalue 0.001



perl /home/patrick/bin/samtools/misc/blast2sam.pl -sd ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.blast.sam  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/test.sam



samtools view -hT ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/transcriptome/gencode.v27.lncRNA_and_transcripts.fa  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/test.sam > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/test2.sam



blast2sam.pl is a script for parsing output of NCBI's blastn output (default format) into sam format


blast2sam.pl -sd out.blast > out.blast.sam
samtools view -hT your_ref.fasta your_file.sam > your_file_with_header.sam'
$ blast2bam [options] blast.xml ref.fasta FastQ_1 [FastQ_2] > out.sam
~/bin/Blast2Bam/bin/blast2bam  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.blast.sam  ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.fa > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/test3.sam
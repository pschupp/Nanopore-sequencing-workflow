for file in *.fa; do
/home/patrick/bin/blast/ncbi-blast-2.7.1+/bin/blastn -query $file -task blastn -db ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa  -out $file.sam  -word_size 16 -gapopen 0 -gapextend 2 -penalty -1 -reward 1 -outfmt 5 -num_threads 6 -evalue 0.00001 -max_target_seqs 10
done
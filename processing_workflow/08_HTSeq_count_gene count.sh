for f in *alignment.sam; do var="$f.counts.txt";echo $f;/usr/local/bin/htseq-count -m union --nonunique all -t exon -a 0 -i gene_id --additional-attr=gene_name -s no ~/Documents/24_sec_seq/root/08_read_assigned_to_gene/$f ~/Documents/hg_align_db/GRCh37.p13_Gencode/gencode.v19.annotation.gff > $var; done

for f in *alignment.sam; do var="$f.counts.txt";echo $f;/usr/local/bin/htseq-count -m intersection-nonempty -t exon -a 0 -i gene_id --additional-attr=gene_name -s no ~/Documents/24_sec_seq/root/07_alignment_split_by_BC/GRCh38/$f ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/gencode.v27.chr_patch_hapl_scaff.annotation.gff > $var; done



/usr/local/bin/htseq-count -m intersection-nonempty -t exon -a 0 -i gene_id --additional-attr=gene_name -s no ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh37.p7_alignment/genome_alignment/BC_all_reads_output_k14.sam ~/Documents/hg_align_db/GRCh37.p13_Gencode_Genome/gencode.v19.annotation.gff > unique_GRCh37_combo.txt
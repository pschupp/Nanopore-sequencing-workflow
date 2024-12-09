/home/patrick/bin/minimap2/minimap2 -ax splice -k14 /home/patrick/Documents/hg_align_db/GRCh37.p13_Gencode_Genome/GRCh37.p13.genome.fa /home/patrick/Documents/24_sec_seq/root/04_alignment_to_barcode/BC_all_reads.fasta  > BC_all_reads_output_k14.sam 


./samtools stats /home/patrick/Documents/24_sec_seq/root/04_alignment_to_barcode/BC_all_reads_output_k14.sam >  /home/patrick/Documents/24_sec_seq/root/04_alignment_to_barcode/BC_all_reads_output_k14.sam.summary_ercc.txt

d
/home/patrick/bin/minimap2/minimap2 -ax splice -k14 /home/patrick/Documents/hg_align_db/ercc.fa /home/patrick/Documents/24_sec_seq/root/05_read_fasta_renamed_with_barcode/BC_all_reads.fasta  > BC_all_reads_output_k14_ercc.sam 

/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_ercc.sam >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14.sam.summary_ercc.txt

/home/patrick/bin/samtools/samtools sort ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_ercc.sam >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14.ercc.sort.sam

/home/patrick/bin/samtools/samtools depth  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14.ercc.sort.sam > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14.ercc.sort.sam.depth

/home/patrick/bin/minimap2/minimap2 -ax splice -k14 /home/patrick/Documents/hg_align_db/GRCh37.p13_Gencode_Genome/Transcriptome/genecode.v19.lnc_and_pc_transcripts.custom.fa /home/patrick/Documents/24_sec_seq/root/05_read_fasta_renamed_with_barcode/BC_all_reads.fasta  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome.sam 

/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome.sam >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome.sam.summary.txt

cat BC_all_reads_output_k14_pc_lnc_txome.sam.summary.txt | grep ^SN
SN	raw total sequences:	1088145
SN	filtered sequences:	0
SN	sequences:	1088145
SN	is sorted:	0
SN	1st fragments:	1088145
SN	last fragments:	0
SN	reads mapped:	636651
SN	reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	451494
SN	reads properly paired:	0	# proper-pair bit set
SN	reads paired:	0	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	157852	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	527164
SN	total length:	1979057432	# ignores clipping
SN	bases mapped:	1334451331	# ignores clipping
SN	bases mapped (cigar):	1033879710	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	164270871	# from NM fields
SN	error rate:	1.588878e-01	# mismatches / bases mapped (cigar)
SN	average length:	1818
SN	maximum length:	20905
SN	average quality:	255.0
SN	insert size average:	0.0
SN	insert size standard deviation:	0.0
SN	inward oriented pairs:	0
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0

/home/patrick/bin/minimap2/minimap2 -ax map-ont -k14 /home/patrick/Documents/hg_align_db/GRCh37.p13_Gencode_Genome/Transcriptome/genecode.v19.lnc_and_pc_transcripts.custom.fa /home/patrick/Documents/24_sec_seq/root/05_read_fasta_renamed_with_barcode/BC_all_reads.fasta  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam 

/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam  >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam.summary.txt

cat BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam.summary.txt | grep ^SN
SN	raw total sequences:	1088145
SN	filtered sequences:	0
SN	sequences:	1088145
SN	is sorted:	0
SN	1st fragments:	1088145
SN	last fragments:	0
SN	reads mapped:	624786
SN	reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	463359
SN	reads properly paired:	0	# proper-pair bit set
SN	reads paired:	0	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	155151	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	533703
SN	total length:	1979057432	# ignores clipping
SN	bases mapped:	1312538593	# ignores clipping
SN	bases mapped (cigar):	998303335	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	145436897	# from NM fields
SN	error rate:	1.456841e-01	# mismatches / bases mapped (cigar)
SN	average length:	1818
SN	maximum length:	20905
SN	average quality:	255.0
SN	insert size average:	0.0
SN	insert size standard deviation:	0.0
SN	inward oriented pairs:	0
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0


/home/patrick/bin/samtools/samtools sort ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/transcriptome_alignment/BC_all_reads_output_k14_pc_lnc_txome.sam >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/transcriptome_alignment/BC_all_reads_output_k14_pc_lnc_txome.sort.sam 

/home/patrick/bin/samtools/samtools index  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/transcriptome_alignment/BC_all_reads_output_k14_pc_lnc_txome.sort.sam 
> ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/transcriptome_alignment/BC_all_reads_output_k14_pc_lnc_txome.sort.index




/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam  >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_all_reads_output_k14_pc_lnc_txome_no_splice.sam.summary.txt



/home/patrick/bin/minimap2/minimap2 -ax splice -k14 /home/patrick/Documents/hg_align_db/GRCh37.p13_Gencode_Genome/Transcriptome/genecode.v19.lnc_and_pc_transcripts.custom.fa /home/patrick/Documents/24_sec_seq/root/05_read_fasta_renamed_with_barcode/unmapped_reads.fa  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_unmapped_reads_output_k14_pc_lnc_txome.sam 


/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_unmapped_reads_output_k14_pc_lnc_txome.sam   >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/BC_unmapped_reads_output_k14_pc_lnc_txome.sam.summary.txt

cat BC_unmapped_reads_output_k14_pc_lnc_txome.sam.summary.txt | grep '^SN'
SN	raw total sequences:	267380
SN	filtered sequences:	0
SN	sequences:	267380
SN	is sorted:	0
SN	1st fragments:	267380
SN	last fragments:	0
SN	reads mapped:	1420
SN	reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	265960
SN	reads properly paired:	0	# proper-pair bit set
SN	reads paired:	0	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	485	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	1862
SN	total length:	269506612	# ignores clipping
SN	bases mapped:	1923313	# ignores clipping
SN	bases mapped (cigar):	1285001	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	496792	# from NM fields
SN	error rate:	3.866083e-01	# mismatches / bases mapped (cigar)
SN	average length:	1007
SN	maximum length:	20905
SN	average quality:	255.0
SN	insert size average:	0.0
SN	insert size standard deviation:	0.0
SN	inward oriented pairs:	0
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0


/home/patrick/bin/minimap2/minimap2 -ax splice -k14 ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa /home/patrick/Documents/24_sec_seq/root/05_read_fasta_renamed_with_barcode/BC_all_reads.fasta  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/BC_all_reads_output_k14_38.p10.sam 

/home/patrick/bin/samtools/samtools stats ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/BC_all_reads_output_k14_38.p10.sam    >  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/BC_all_reads_output_k14_38.p10.sam.summary.txt

cat BC_all_reads_output_k14_38.p10.sam.summary.txt | grep '^SN'
SN	raw total sequences:	1088145
SN	filtered sequences:	0
SN	sequences:	1088145
SN	is sorted:	0
SN	1st fragments:	1088145
SN	last fragments:	0
SN	reads mapped:	821516
SN	reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	266629
SN	reads properly paired:	0	# proper-pair bit set
SN	reads paired:	0	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	29383	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	174605
SN	total length:	1979057432	# ignores clipping
SN	bases mapped:	1710450217	# ignores clipping
SN	bases mapped (cigar):	1550988201	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	228925626	# from NM fields
SN	error rate:	1.475999e-01	# mismatches / bases mapped (cigar)
SN	average length:	1818
SN	maximum length:	20905
SN	average quality:	255.0
SN	insert size average:	0.0
SN	insert size standard deviation:	0.0
SN	inward oriented pairs:	0
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0

/home/patrick/bin/minimap2/minimap2 -ax splice -k10 ~/Documents/hg_align_db/GRCh38.p10_Gencode_Genome/genomic/GRCh38.p10.genome.fa  ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.fa  > ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh38.p10_alignment/unmapped_reads_GRCh38.p10.k_10.minmap2.sam 
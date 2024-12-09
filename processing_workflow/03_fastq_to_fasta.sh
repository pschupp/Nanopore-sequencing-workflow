#!/bin/bash
#Part of FASTX Toolkit 0.0.13 by A. Gordon (gordon@cshl.edu)
/home/patrick/Downloads/fast_tools/bin/fastq_to_fasta -n -v -Q33 -i /home/patrick/Single_Cell_Run/read_combo/basecalled/workspace/pass/all.fastq -o /home/patrick/Single_Cell_Run/read_combo/basecalled/workspace/pass/all.fasta

#It assumes Illumina (ASCII offset 64), but you have Sanger (offset 33), hence -Q33 See: http://seqanswers.com/forums/showthread.php?t=7399


/home/patrick/Downloads/fast_tools/bin/fastq_to_fasta -n -v -Q33 -i ~/Documents/24_sec_seq/root/01_raw_reads/all_reads.fastq -o ~/Documents/24_sec_seq/root/01_raw_reads/all_reads.fasta



for d in *.fastq; do
    mv $d ${d}_fail_basecld.fastq
done

for d in *.fastq; do
    echo $d
    fastq_to_fasta -Q33 -n -i $d -o ${d}.fasta
done

dest="/nas100/backups/servers/z/zebra/mysql.tgz"
## get file name i.e. basename such as mysql.tgz
tempfile="${dest##*/}"
 
## display filename 
echo "${tempfile%.*}"

dest=/home/patrick/Documents/all_combo/fasta/*fasta
## get file name i.e. basename such as mysql.tgz
tempfile="${dest##*/}"
 
## display filename 
echo "${tempfile%.*}"

for f in test_suite/*.args;do 
	fn="${f%.*}" 
	printf "filename without path and extension: %s, extension only: %s\n" "${fn##*/}" "${f##*.}"  
done

for f in *fasta;do 
	fn="${f%.*}" 
	echo "${fn##*/}" >> output_fasta_files.txt
done

for f in *fastq;do 
	fn="${f%.*}" 
	echo "${fn##*/}" >> output_fastq_files.txt
done
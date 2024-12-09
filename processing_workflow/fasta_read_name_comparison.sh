for barc in {20..49}; do
echo "BC$barc"
cat unmapped_reads_GRCh38.p10.fa | grep BC$barc | wc -l
done

diff unaligned_names.txt aligned_names.txt > diff.txt

comm -12 unaligned_names.txt aligned_names.txt > share.txt
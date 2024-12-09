let "start_s=1"
let "end_s=20000"
for mult in {1..25}; do
sed -n "$start_s,$end_s"p unmapped_reads_GRCh38.p10.fa > $mult.fa
let "start_s=(20000*$mult)+1"
let "end_s=20000*($mult+1)"
done

let "start_s=520001"
let "end_s = 533258"
sed -n "$start_s,$end_s"p unmapped_reads_GRCh38.p10.fa > 26.fa

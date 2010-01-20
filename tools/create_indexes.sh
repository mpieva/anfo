#! /bin/bash

# creates ANFO genome files and indices, split into manageable chunks;
# currently for human and chimp

humdesc="homo sapiens sapiens, rev. 18"
chimpdesc="pan troglodytes, rev. 2"

declare -a chroms
set -e

i=0
while read -a chroms ; do
    for c in "${chroms[@]}" ; do
        twoBitToFa "../hg18.2bit:$c" stdout
    done | fa2dna -v -o hg18_${i}.dna -m 2 -g hg18 - -d "$humdesc"
    dnaindex -v -o hg18_${i}_12_4.idx -g hg18_${i}.dna -s 12 -S 4 -r 6
    i=$((i+1)) ;
done << EOF
chr17_random chrX chr7 chr6 chr5 chr4 chr3 chr2 chr1
chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrY chr22 chr21 chr6_random chrX_random chr21_random chr1_random chr9_random chr8_random chr4_random chr15_random chr3_random chr7_random chr19_random chr22_random chr11_random chr13_random chr2_random chr5_random chr10_random chr16_random chrM chr18_random
EOF

i=0
while read -a chroms ; do
    for c in "${chroms[@]}" ; do
        twoBitToFa "../pt2.2bit:$c" stdout
    done | fa2dna -v -o pt2_${i}.dna -m 2 -g pt2 - -d "$chimpdesc"
    dnaindex -v -o pt2_${i}_12_4.idx -g pt2_${i}.dna -s 12 -S 4 -r 6
    i=$((i+1)) ;
done << EOF
chr8 chrX chr7 chr6 chr5 chr4 chr3 chr1 chr2b
chr9 chr12 chr10 chr11 chr13 chr2a chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrUn chr22 chr21 chrY chr1_random chr6_random chr13_random chr11_random chr10_random chr9_random chr8_random chr7_random chr4_random chr16_random chr17_random chrX_random chr3_random chr5_random chr15_random chr2a_random chr19_random chr12_random chr2b_random chr14_random chr18_random chr20_random chr22_random chrY_random chrM
EOF




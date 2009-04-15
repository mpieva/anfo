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
        twoBitToFa "hg18.2bit:$c" stdout
    done | fa2dna -v -o hg18_${i}.dna -m 2 -g hg18 - -d "$humdesc"
    dnaindex -v -o hg18_${i}_10.idx -g hg18_${i}.dna -s 10
    i=$((i+1)) ;
done << EOF
  chr7 chr1
  chrX chr2
  chr4 chr3
  chr18_random chrM chr16_random chr10_random chr5_random chr2_random \
    chr13_random chr11_random chr22_random chr19_random chr7_random \
    chr3_random chr15_random chr4_random chr8_random chr9_random chr6 \
    chr1_random chr21_random chrX_random chr6_random chr17_random chr5
  chr13 chr9 chr8
  chr12 chr11 chr10
  chr17 chr16 chr15 chr14
  chr18 chr19 chr20 chrY chr22 chr21
EOF

i=0
while read -a chroms ; do
    for c in "${chroms[@]}" ; do
        twoBitToFa "pt2.2bit:$c" stdout
    done | fa2dna -v -o pt2_${i}.dna -m 2 -g pt2 - -d "$chimpdesc"
    dnaindex -v -o pt2_${i}_10.idx -g pt2_${i}.dna -s 10
    i=$((i+1)) ;
done << EOF
  chr5 chr2b
  chr3 chr1
  chr22 chr6 chr4
  chr13 chrX chr7
  chr12 chr9 chr8
  chrY chr2a chr11 chr10
  chr17_random chr6_random chr1_random chr17 chr16 chr15 chr14
  chr18 chr19 chr20 chrUn chr21 chr13_random chr11_random chr10_random \
    chr9_random chr8_random chr7_random chr4_random chr16_random \
    chrX_random chr3_random chr5_random chr15_random chr2a_random \
    chr19_random chr12_random chr2b_random chr14_random chr18_random \
    chr20_random chr22_random chrY_random chrM
EOF




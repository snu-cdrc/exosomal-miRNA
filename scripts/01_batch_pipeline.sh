#!/bin/bash
# ==============================================================================
# Small RNA-seq Analysis Pipeline
# ==============================================================================
# Identifies differentially expressed miRNAs between primary tumor-derived EVs
# (pTDEs) and metastatic tumor-derived EVs (mTDEs) from canine mammary
# carcinoma cell lines (CHMp vs CHMm).
#
# Reference genome: CanFam3.1
# Sequencing: Illumina small RNA-seq (NEXTflex Small RNA-Seq Kit v3)
# Samples: CHMp (primary, n=2), CHMm (metastatic, n=2)
#
# Adapter sequences (NEXTflex Small RNA-Seq Kit):
#   3' adapter: TGGAATTCTCGGGTGCCAAGG
#   5' adapter: GATCGTCGGACTGTAGAACTCTGAAC
#   Note: NEXTflex uses 4bp random sequences (NNNN) at both ends of the insert
# ==============================================================================

THREADS=23
GENOME_INDEX="canFam3v01"
GENOME_FA="canFam3v01.fa"
MATURE_CFA="mature_cfa.fa"
HAIRPIN_CFA="hairpin_cfa.fa"

# ==============================================================================
# Step 1: Quality Control (before trimming)
# ==============================================================================

for fq in *.fq.gz; do
    pre=${fq%%.*}
    mkdir -p qc/before_trim/${pre}
    fastqc -o qc/before_trim/${pre} -t ${THREADS} ${fq}
done

# ==============================================================================
# Step 2: Adapter Trimming (cutadapt)
# ==============================================================================
# -a: 3' adapter sequence
# -u 4: remove 4bp from 5' end (NEXTflex random adapter)
# -m 22 / -M 30: keep reads between 22-30bp (mature miRNA size range)
# -q 20: quality threshold (Phred33)
# -O 6: minimum overlap for adapter matching

mkdir -p trim

# Forward reads
for fq1 in *_1.fq.gz; do
    pre=${fq1%%.*}
    cutadapt --quality-base 33 -u 4 -m 22 -M 30 -f fastq -q 20 -O 6 \
        -j ${THREADS} -a TGGAATTCTCGGGTGCCAAGG \
        ${fq1} -o trim/trim_${fq1}
done

# Reverse reads
for fq2 in *_2.fq.gz; do
    pre=${fq2%%.*}
    cutadapt --quality-base 33 -u 4 -m 22 -M 30 -f fastq -q 20 -O 6 \
        -j ${THREADS} -a GATCGTCGGACTGTAGAACTCTGAAC \
        ${fq2} -o trim/trim_${fq2}
done

# ==============================================================================
# Step 3: Remove 3' Random Adapter Sequences (NEXTflex 4bp)
# ==============================================================================
# NEXTflex adds 4bp random sequences at both ends of the insert.
# 5' random bases: removed by cutadapt (-u 4) in Step 2.
# 3' random bases: removed by this custom script from sequence/quality tails.

python scripts/02_trim_read_tails.py trim/ --trim-length 4

# ==============================================================================
# Step 4: Quality Control (after trimming)
# ==============================================================================

for fq in trim/cut4bp_*.fq.gz; do
    pre=$(basename ${fq%%.*})
    mkdir -p qc/after_trim/${pre}
    fastqc -o qc/after_trim/${pre} -t ${THREADS} ${fq}
done

# ==============================================================================
# Step 5: Merge Paired-End Reads to Single-End
# ==============================================================================

mkdir -p merged

for sample in CHMm_rep1 CHMm_rep2 CHMp_rep1 CHMp_rep2; do
    cat trim/cut4bp_trim_*_${sample}_1.fq trim/cut4bp_trim_*_${sample}_2.fq \
        > merged/${sample}.fq
done

# ==============================================================================
# Step 6: miRDeep2 - Read Collapsing & Known miRNA Quantification
# ==============================================================================

mkdir -p collapsed

# Collapse reads
for fq in merged/*.fq; do
    pre=$(basename ${fq%%.fq})
    mapper.pl ${fq} -e -h -j -m -s collapsed/${pre}_collapsed.fa
done

# Quantify known miRNAs
for fa in collapsed/*_collapsed.fa; do
    quantifier.pl \
        -p ${HAIRPIN_CFA} -m ${MATURE_CFA} \
        -r ${fa} -t cfa -g 2 -e 2 -f 5
done

# ==============================================================================
# Step 7: Novel miRNA Detection (miRDeep2)
# ==============================================================================
# Merge all samples for higher sensitivity in novel miRNA discovery.

cat merged/*.fq > merged/all_samples_merged.fq

mapper.pl merged/all_samples_merged.fq \
    -e -h -j -m \
    -p ${GENOME_INDEX} \
    -s collapsed/all_merged_collapsed.fa \
    -t collapsed/all_merged_vs_genome.arf -v

miRDeep2.pl \
    collapsed/all_merged_collapsed.fa \
    ${GENOME_FA} \
    collapsed/all_merged_vs_genome.arf \
    ${MATURE_CFA} mature_hsa_mmu.fa ${HAIRPIN_CFA} \
    -t cfa 2> logs/mirdeep2_novel_report.log

# ==============================================================================
# Step 8: Extract miRNA Primary Transcript Regions
# ==============================================================================

awk '$3 == "miRNA_primary_transcript" {print $1"\t"$4"\t"$5}' \
    reference/cfa_miRBase.gff3 > reference/cfa_miRNA_primary_transcript.bed

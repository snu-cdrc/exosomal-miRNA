# Exosomal miRNA Analysis Pipeline

Analysis code for:

**"Extracellular vesicle-mediated transfer of miRNA-1 from primary tumors represses the growth of distant metastases"**

Chae-Yi Kim<sup>1,2,3†</sup>, Kang-Hoon Lee<sup>1,2,3</sup>, [Keun Hong Son](https://keun-hong.github.io/)<sup>1,2,3</sup>, Tae-Jin Shin<sup>1,2,3</sup>, and [**Je-Yoel Cho**](https://vetbio.snu.ac.kr/)<sup>1,2,3*</sup>

<sup>1</sup> Department of Biochemistry, College of Veterinary Medicine, Seoul National University, Seoul, Korea<br>
<sup>2</sup> Comparative Medicine and Disease Research Center (CDRC), Science Research Center (SRC), Seoul National University, Seoul, Korea<br>
<sup>3</sup> BK21 PLUS Program for Creative Veterinary Science Research and Research Institute for Veterinary Science, Seoul National University, Seoul, Korea

<sup>†</sup> First author · <sup>*</sup> Corresponding author

Published in *Experimental & Molecular Medicine* (2024)<br>
DOI: [10.1038/s12276-024-01181-7](https://doi.org/10.1038/s12276-024-01181-7)

## Abstract

Metastases originate from primary tumors and reach distant organs. Growing evidence suggests that metastases are under the control of primary tumors even outside the primary site; however, the mechanisms by which primary tumors remotely control metastases remain unclear. Here, we discovered a molecular mechanism by which primary tumors suppress metastatic growth. We found that extracellular vesicles (EVs) derived from the primary tumor can inhibit the growth of metastases both *in vitro* and *in vivo*. miR-1 was particularly enriched in primary tumor-derived EVs (pTDEs) and was found to be responsible for the suppression of metastatic growth. Mechanistically, intracellular reactive oxygen species (ROS) production and DNA damage were induced, which led to cell cycle arrest. Collectively, our data demonstrate that primary tumors restrict the growth of distant metastases via miR-1 in pTDEs and that miR-1 could potentially be used as an antimetastatic agent.

## Overview

Small RNA sequencing analysis pipeline for identifying differentially expressed miRNAs between primary tumor-derived extracellular vesicles (pTDEs) and metastatic tumor-derived extracellular vesicles (mTDEs) from canine mammary carcinoma cell lines.

## Pipeline

```
Raw small RNA reads (Illumina, NEXTflex Small RNA-Seq Kit v3)
  → FastQC (quality control)
  → cutadapt (adapter trimming, quality filtering)
  → trim_read_tails.py (NEXTflex 4bp random adapter removal)
  → Bowtie (alignment to CanFam3.1)
  → miRDeep2 (miRNA identification & quantification)
  → edgeR (differential expression analysis)
```

## Samples

| Sample | Description | Replicates |
|--------|-------------|------------|
| CHMp | Primary tumor-derived EVs (CHMp cell line) | 2 |
| CHMm | Metastatic tumor-derived EVs (CHMm cell line) | 2 |

Raw sequencing data: [GEO GSE213969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213969)

## Directory Structure

```
├── scripts/
│   ├── 01_batch_pipeline.sh       # Main analysis pipeline
│   └── 02_trim_read_tails.py      # Trim N bp from FASTQ read tails
├── reference/
│   ├── adapter_list.txt           # Adapter sequences for FastQC
│   ├── cfa_miRBase.gff3           # Canine miRNA annotations (miRBase)
│   ├── cfa_miRNA_primary_transcript.bed  # Canine miRNA genomic regions
│   ├── cfa_hairpin_ref.fa         # Canine miRNA hairpin sequences (miRBase)
│   └── cfa_mature_ref.fa          # Canine mature miRNA sequences (miRBase)
├── results/
│   ├── Differentially_expressed_miRNAs.txt       # edgeR DE results
│   └── miRNAs_expressed_all_samples_quantifier.csv  # miRDeep2 quantification
└── README.md
```

## Key Parameters

### cutadapt
```bash
cutadapt \
    --quality-base 33 \
    -u 4           `# remove 4bp from 5' end (NEXTflex random adapter)` \
    -m 22 -M 30    `# keep reads between 22-30bp` \
    -q 20          `# quality threshold` \
    -O 6           `# minimum adapter overlap` \
    -a TGGAATTCTCGGGTGCCAAGG \
    input.fq.gz -o output.fq.gz
```

### miRDeep2
```bash
# Read mapping & collapsing
mapper.pl reads.fq -e -h -j -m -s collapsed.fa

# Known miRNA quantification
quantifier.pl -p hairpin.fa -m mature.fa -r collapsed.fa -t cfa -g 2 -e 2 -f 5

# Novel miRNA detection
miRDeep2.pl collapsed.fa canFam3v01.fa collapsed_vs_genome.arf \
    mature_cfa.fa mature_hsa_mmu.fa hairpin_cfa.fa -t cfa
```

## Citation

```bibtex
@article{kim2024extracellular,
  title     = {Extracellular vesicle-mediated transfer of miRNA-1 from primary
               tumors represses the growth of distant metastases},
  author    = {Kim, Chae-Yi and Lee, Kang-Hoon and Son, Keun Hong and
               Shin, Tae-Jin and Cho, Je-Yoel},
  journal   = {Experimental \& Molecular Medicine},
  year      = {2024},
  doi       = {10.1038/s12276-024-01181-7}
}
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

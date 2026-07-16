# Human RNome Project

UCSC Genome Browser assembly hub visualizing RNA modifications across the human
RNome: poly(A) RNA (hg38), ribosomal RNA (rRNA), and transfer RNA (tRNA).
Modifications are detected by three modalities — short-read sequencing (SRS,
Illumina), long-read sequencing (LRS, Oxford Nanopore/ONT), and mass
spectrometry (MS). Each modality is rendered as a single flat track per RNA
type, colored by modification type with saturation proportional to frequency.
Each RNA type also has a tiered cross-platform consensus track (tier1/tier2
confidence). The poly(A) RNA assembly additionally includes a native RNA
sequence track (bigPSL) with consensus modification sites overlaid as
decorators. Track descriptions are sourced from `RNA_modifications_manifest.tsv`.

---

## Repo structure

```
RNome_dataviz_repo/
├── RNA_modifications_manifest.tsv          # Source of truth: Filename → Description, per modality
├── final_bedRmods/                         # Pipeline input: pre-merged, genome-mapped BED files
│   ├── Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed
│   ├── ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed
│   └── MS_rRNA_tRNA.bed
│
├── concensus/                              # Pipeline input: combined tiered consensus bed
│   └── consensus_draft_sequence.bed           # header: #chr/start/end/name/tier/strand
│
├── draft_sequence/                         # Native RNA sequence alignment files
│   ├── HRP.psl                                # PSL alignment of reads to hg38
│   ├── HRP_reads.fa                           # Read sequences (used to embed in bigPSL)
│   ├── HRP.bigPsl.bb                          # Indexed bigPSL (gitignored; hosted on GCS)
│   ├── bigPsl.as                              # AutoSQL schema for bigPSL
│   ├── HRP_sequence_draft.bam(.bai)           # Source BAM
│   └── draft_native_reference.*              # Reference FASTA, FAI, 2bit
│
├── reference/                              # Reference sequences & annotations
│   ├── hs_rRNAs_NR_046235.fa(.fai)            # Human rRNA reference (18S/28S/5.8S/5S)
│   ├── Step3_hg38_Homo_sapiens_47seq...fa     # tRNA reference sequences
│   └── gencode.v49...annotation.gtf.gz        # GENCODE hg38 annotation
│
├── final_data/                             # Pipeline output, organized by RNA type
│   ├── RNA_modifications_manifest.xlsx
│   ├── polyA_RNA_hg38/
│   │   └── polyA_RNA_hg38_{SRS,LRS,consensus_tiered}.bed
│   ├── rRNA/
│   │   └── rRNA_{SRS,LRS,MS,consensus_tiered}.bed
│   └── tRNA/
│       └── tRNA_{SRS,MS,consensus_tiered}.bed
│
├── scripts/
│   ├── pipeline.py                         # End-to-end hub generation (see below)
│   ├── make_decorator.py                   # Builds HRP_decorator.bigBed for the sequence track
│   └── BedPyLift.py                        # Transcriptome → genome coordinate lifting
│
├── ucsc_hub/                               # Public UCSC assembly hub (on GitHub)
│   ├── hub.txt / genomes.txt / hubDescription.html
│   ├── hg38_polyA-RNA/                        # hg38: bigPSL track + modifications + NA12878 variants
│   ├── rrna/                                  # hs_rRNA custom assembly
│   └── trna/                                  # hs_tRNA custom assembly
│
├── ucsc_hub_internal/                      # Internal hub — no HTML files (gitignored)
│   ├── hg38_polyA-RNA/                        # Standalone hub (own hub.txt + genomes.txt)
│   ├── rrna/                                  # Standalone hub
│   └── trna/                                  # Standalone hub
│
└── utils/                                  # UCSC command-line tools
    ├── bedToBigBed
    ├── bigBedToBed
    ├── faToTwoBit
    ├── pslToBigPsl
    └── hg38.chrom.sizes
```

---

## Pipeline

### Requirements

- Python 3.7+
- `openpyxl` (optional — used for the Excel manifest; falls back to TSV if missing)
- UCSC command-line utilities (included in `utils/` for macOS Apple Silicon):

```bash
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/pslToBigPsl
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigBedToBed
chmod +x bedToBigBed faToTwoBit pslToBigPsl bigBedToBed
```

### Run the pipeline

```bash
python scripts/pipeline.py
```

This will:

1. Split each combined SRS/LRS/MS BED file by RNA type into
   `final_data/{RNA_type}/{RNA_type}_{Modality}.bed`
2. Normalize and convert each split BED to bigBed (`bed9+4`); also split and
   convert `concensus/consensus_draft_sequence.bed` to tiered consensus bigBed
   (`bed9+1`) per RNA type
3. Generate one HTML description per track, using text from
   `RNA_modifications_manifest.tsv`
4. Generate `final_data/RNA_modifications_manifest.xlsx`
5. Generate the UCSC hub config (`hub.txt`, `genomes.txt`, per-assembly
   `trackDb.txt`)
6. Run `scripts/make_decorator.py` to build `HRP_decorator.bigBed` — the
   modification decorator overlay for the RNA sequence track

### Building the RNA sequence track (one-time)

The bigPSL track is built from the source BAM and is too large for GitHub
(166 MB); it is hosted on GCS. To regenerate:

```bash
# 1. Convert BAM → PSL
utils/bamToPsl draft_sequence/HRP_sequence_draft.bam draft_sequence/HRP.psl

# 2. Convert PSL → bigPSL (embed read sequences for base-level rendering)
utils/pslToBigPsl -fa=draft_sequence/HRP_reads.fa draft_sequence/HRP.psl \
  draft_sequence/HRP.bigPsl.as stdout \
  | sort -k1,1 -k2,2n \
  | awk 'NF==24{print $0"\t1"}' \
  > draft_sequence/HRP.filtered.bigPsl

utils/bedToBigBed -type=bed12+13 -tab -as=draft_sequence/bigPsl.as \
  draft_sequence/HRP.filtered.bigPsl utils/hg38.chrom.sizes \
  draft_sequence/HRP.bigPsl.bb

# 3. Upload to GCS
gsutil cp draft_sequence/HRP.bigPsl.bb gs://hrp_proj/HRP.bigPsl.bb
```

Then run `make_decorator.py` to rebuild the decorator:
```bash
python scripts/make_decorator.py
```

### One-time setup — custom assembly sequences (.2bit)

```bash
utils/faToTwoBit reference/hs_rRNAs_NR_046235.fa ucsc_hub/rrna/hs_rRNAs_NR_046235.2bit
utils/faToTwoBit reference/Step3_hg38_Homo_sapiens_47seq20251209.fa ucsc_hub/trna/hs_tRNA.2bit
```

### Push to GitHub

```bash
git add .
git commit -m "Update hub"
git push
```

---

## Loading the hub in UCSC

1. Go to [genome.ucsc.edu](https://genome.ucsc.edu)
2. Click **My Data** → **Track Hubs** → **My Hubs**
3. Paste the hub URL and click **Add Hub**

| Hub | URL |
|---|---|
| Public (with HTML descriptions) | `https://raw.githubusercontent.com/baihesun/RNome_UCSC_hub/main/ucsc_hub/hub.txt` |
| Internal poly(A) RNA | `ucsc_hub_internal/hg38_polyA-RNA/hub.txt` (local) |
| Internal rRNA | `ucsc_hub_internal/rrna/hub.txt` (local) |
| Internal tRNA | `ucsc_hub_internal/trna/hub.txt` (local) |

---

## Assemblies & tracks

| Assembly | UCSC name | Tracks | Default position |
|---|---|---|---|
| Poly(A) RNA | `hg38` | RNA Sequence, SRS, LRS, Consensus | `chr6:32,659,538-32,659,638` |
| rRNA | `hs_rRNA` | SRS, LRS, MS, Consensus | `hs_rRNA_18S:1-1869` |
| tRNA | `hs_tRNA` | SRS, MS, Consensus | `hs_tRNAAla_AGC:1-76` |

### Track types

**RNA Sequence (poly(A) RNA only)** — native RNA reads aligned to hg38 in
bigPSL format. Read sequences are embedded for base-level rendering. Consensus
modification sites are overlaid as decorators colored by modification type.
The bigPSL file is hosted on GCS (`gs://hrp_proj/HRP.bigPsl.bb`).

**Per-platform tracks (SRS / LRS / MS)** — all modification sites for that
RNA type and modality, colored by modification type with saturation
proportional to frequency. Filterable by modification type, score, coverage,
and frequency.

**Consensus track** — tiered cross-platform consensus merging evidence across
platforms. Saturation encodes confidence tier: vivid = tier1 (strongest
evidence), pale = tier2 (≥2 datasets). Filterable by modification type and
tier.

**NA12878 variants (poly(A) RNA only)** — WGS variant calls for NA12878,
intersected from DeepVariant v1.9.0 and GATK GenotypeGVCFs. Hosted on GCS.

---

## BED columns (bed9+4, `rna_mods.as`)

| Column | Field | Description |
|---|---|---|
| 1 | chrom | Chromosome or RNA molecule |
| 2–3 | chromStart / chromEnd | Position of the modified site (0-based, half-open) |
| 4 | name | Modification type (e.g. `m6A`, `Y`) |
| 5 | score | Detection score, clamped 0–1000 |
| 6 | strand | `+` or `-` |
| 7–8 | thickStart / thickEnd | Same as start/end |
| 9 | itemRgb | Color — hue per modification type, saturation ∝ frequency |
| 10 | coverage | Reads/spectra covering this position |
| 11 | frequency | Modification frequency (%), 0–100 |
| 12 | single_letter_code | Single-letter modification code (bedRmod) |
| 13 | mod_id | Modification identifier |

## Tiered consensus BED columns (bed9+1, `tiered_mods.as`)

| Column | Field | Description |
|---|---|---|
| 1 | chrom | Chromosome or RNA molecule |
| 2–3 | chromStart / chromEnd | Position of the modified site |
| 4 | name | Modification type |
| 5 | score | `1000` = tier1, `500` = tier2 |
| 6 | strand | `+` or `-` |
| 7–8 | thickStart / thickEnd | Same as start/end |
| 9 | itemRgb | Color — hue per modification type, vivid = tier1, pale = tier2 |
| 10 | tier | `tier1` or `tier2` |

---

## Color scheme

| Modification | Full name | Hex |
|---|---|---|
| m6A | N6-methyladenosine | `#003f5c` |
| Am | 2'-O-methyladenosine | `#1464a0` |
| Ino / I | Inosine | `#0e8fd4` |
| m1A | N1-methyladenosine | `#00b4d8` |
| m66A | N6,N6-dimethyladenosine | `#00d4b4` |
| m1acp3Y | 1-methyl-3-(3-amino-3-carboxypropyl)pseudouridine | `#00c875` |
| Y | Pseudouridine | `#90be6d` |
| Ym | 2'-O-methylpseudouridine | `#c9b84c` |
| Um | 2'-O-methyluridine | `#ffa600` |
| m3U | 3-methyluridine | `#ff7c43` |
| Gm | 2'-O-methylguanosine | `#ff6361` |
| m7G | N7-methylguanosine | `#dd7371` |
| Cm | 2'-O-methylcytidine | `#ef476f` |
| m5C | 5-methylcytidine | `#d45087` |
| ac4C | N4-acetylcytidine | `#7b2d8b` |

Modifications not in this table are rendered in gray (`#808080`). On
per-platform tracks, saturation encodes stoichiometry (vivid = high frequency).
On the consensus track, saturation encodes tier (vivid = tier1, pale = tier2).

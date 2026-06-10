# RNome RNA Modification Atlas — UCSC Track Hub

UCSC Genome Browser assembly hub visualizing RNA modifications across the human
"RNome": poly-A RNA (hg38), ribosomal RNA (rRNA, 18S/28S/5.8S/5S), and transfer
RNA (tRNA). Modifications are detected by three modalities — Illumina
sequencing, Oxford Nanopore (ONT), and mass spectrometry (MS) — and combined
into per-RNA-type, per-modality consensus tracks plus per-experiment subtracks.

---

## Repo structure

```
RNome_dataviz_repo/
├── reference/                              # Reference sequences & annotations
│   ├── hs_rRNAs_NR_046235.fa(.fai)            # Human rRNA reference (18S/28S/5.8S/5S)
│   ├── Step3_hg38_Homo_sapiens_47seq...fa     # tRNA reference sequences
│   └── gencode.v49...annotation.gtf.gz        # GENCODE hg38 annotation (for BedPyLift)
│
├── processed_data/                         # Per-modality intermediate files
│   ├── illumina/ ont/ ms/
│   │   ├── *_header.txt                       # bedRmod header (mod-code legend, provenance)
│   │   ├── *_rna.bed                          # Modifications in transcriptomic coords
│   │   └── *_genome.bed                       # Lifted to genomic coords via BedPyLift.py
│   └── consensus_genome.bed                   # Modalities merged, duplicate sites collapsed
│
├── final_data/                             # Pipeline inputs/outputs, organized by RNA type
│   ├── Illumina_combined_polyARNA_tRNA_rRNA.bed     # Combined Illumina input
│   ├── ONT_polyARNA_rRNA_combined.filtered_*.bed    # Combined ONT input
│   ├── MS_rRNA_tRNA.bed                             # Combined MS input
│   ├── RNA_modifications_manifest.xlsx              # Generated manifest of all tracks
│   ├── polyA-RNA_hg38/
│   │   ├── polyA-RNA_hg38_{Illumina,ONT}.bed          # Consensus tracks
│   │   └── experiments_{illumina,ont}/                # Individual raw experiments
│   ├── rRNA/
│   │   ├── rRNA_{Illumina,ONT,MS}.bed
│   │   └── experiments_{illumina,ont,ms}/
│   └── tRNA/
│       ├── tRNA_{Illumina,MS}.bed
│       └── experiments_{illumina,ms}/
│
├── scripts/
│   ├── pipeline.py                      # End-to-end hub generation (see below)
│   └── BedPyLift.py                        # Transcriptome → genome coordinate lifting
│
├── ucsc_hub/                                # Generated UCSC assembly hub
│   ├── hub.txt / genomes.txt / hubDescription.html
│   ├── hg38_polyA-RNA/                        # hg38 assembly: trackDb, bigBeds, HTML
│   ├── rrna/                                  # hs_rRNA custom assembly: 2bit, trackDb, bigBeds
│   └── trna/                                  # hs_tRNA custom assembly: 2bit, trackDb, bigBeds
│
├── utils/                                    # UCSC command-line tools
│   ├── bedToBigBed
│   └── faToTwoBit
│
└── README.md
```

---

## Pipeline

### Requirements

- Python 3.7+
- [`polars`](https://pola.rs/) (used by `BedPyLift.py`)
- `openpyxl` (optional — used for the Excel manifest; falls back to a TSV if missing)
- UCSC command line utilities (already included in `utils/` for macOS Apple Silicon):

```bash
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit
chmod +x bedToBigBed faToTwoBit
```

### 1. Lift transcriptomic coordinates to genomic coordinates (`BedPyLift.py`)

For each modality, convert per-transcript (ENST) bedRmod positions in
`processed_data/{modality}/{modality}_rna.bed` to genomic coordinates using a
GTF transcript model:

```bash
python scripts/BedPyLift.py \
  --gtf reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
  --bed processed_data/illumina/illumina_rna.bed \
  --output processed_data/illumina/illumina_genome.bed
```

This builds an exon-based transcript index from the GTF, maps each
transcriptomic position to `(chrom, genomic_pos, strand)`, merges duplicate
genomic sites using a coverage-weighted average frequency, and logs mapping
statistics (conversion rate, unmappable transcripts, duplicate sites).

### 2. Run the main pipeline (`pipeline.py`)

Update `SOURCE_DIR`, `COMBINED_FILES`, and `EXPERIMENT_FILES` at the top of
`scripts/pipeline.py` to point at your source bedRmod files, then from the
repo root:

```bash
python scripts/pipeline.py
```

This will:

1. Split the combined Illumina/ONT/MS BED files by RNA type (`rRNA`, `tRNA`,
   `polyA-RNA_hg38`) into `final_data/{RNA_type}/{RNA_type}_{Modality}.bed`
2. Copy individual raw experiment files into
   `final_data/{RNA_type}/experiments_{modality}/`
3. Normalize each BED file for `bedToBigBed`: clamp scores to 0–1000, recolor
   by modification type (hue) × frequency (saturation), remap Ensembl
   scaffold names to UCSC names where needed, and convert to bigBed
   (`bedToBigBed -type=bed9+4 -as=rna_mods.as`)
4. Generate the UCSC hub config (`hub.txt`, `genomes.txt`, per-assembly
   `trackDb.txt`) with consensus tracks, per-experiment composite subtracks,
   and frequency / modification-type filters
5. Generate placeholder per-track HTML descriptions
6. Generate `final_data/RNA_modifications_manifest.xlsx` (or `.tsv` if
   `openpyxl` is unavailable) listing every track and its hub bigBed URL

### One-time setup — convert reference sequences to .2bit

```bash
./utils/faToTwoBit reference/hs_rRNAs_NR_046235.fa ucsc_hub/rrna/hs_rRNAs_NR_046235.2bit
./utils/faToTwoBit reference/Step3_hg38_Homo_sapiens_47seq20251209.fa ucsc_hub/trna/hs_tRNA.2bit
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
3. Paste the following URL and click **Add Hub**:

```
https://raw.githubusercontent.com/baihesun/RNome_UCSC_hub/main/ucsc_hub/hub.txt
```

---

## Assemblies & tracks

| Assembly | UCSC name | Modalities | Notes |
|---|---|---|---|
| Poly-A RNA | `hg38` | Illumina, ONT | Standard hg38; default position `chr21:8400001-8500000` |
| rRNA | `hs_rRNA` | Illumina, ONT, MS | Custom 2bit assembly (18S/28S/5.8S/5S) |
| tRNA | `hs_tRNA` | Illumina, MS | Custom 2bit assembly |

For each modality present on an assembly, the hub provides:

- a **consensus track** (`{RNA_type}_{Modality}.bigBed`) — all sites for that
  RNA type/modality, colored and filterable by modification type and frequency
- an **experiments composite** — one subtrack per individual raw experiment
  file, under `experiments_{modality}/`

The hg38 assembly additionally carries an externally hosted track,
**SRS000090 variants** (WGS merged variant calls from two callers), streamed
from `https://storage.googleapis.com/srs000090-wgs/`.

---

## BED columns (bed9+4, `rna_mods.as`)

| Column | Field | Description |
|---|---|---|
| 1 | chrom | Chromosome or RNA molecule (e.g. `chr1`, `hs_rRNA_18S`, `hs_tRNAAla_AGC`) |
| 2–3 | chromStart / chromEnd | Position of the modified site (0-based, half-open) |
| 4 | name | Modification type (e.g. `Am`, `Y`, `m6A`, `I`) |
| 5 | score | Detection score, clamped to 0–1000 |
| 6 | strand | `+` or `-` |
| 7–8 | thickStart / thickEnd | Same as start/end (display only) |
| 9 | itemRgb | Display color — hue per modification type, saturation ∝ frequency |
| 10 | coverage | Number of reads/spectra covering this position |
| 11 | frequency | Modification frequency (%), 0–100 |
| 12 | single_letter_code | Single-letter modification code (bedRmod convention) |
| 13 | mod_id | Modification identifier |

Clicking a site in UCSC shows all 13 fields in a popup.

---

## Color scheme

Each modification type has a fixed hue. Saturation encodes stoichiometry:
vivid = high frequency, pale/gray = low frequency. Modification types not in
this table are rendered in gray (`#808080`).

| Modification | Type | Hex |
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

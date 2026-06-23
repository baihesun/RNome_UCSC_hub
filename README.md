# RNome RNA Modification Atlas — UCSC Track Hub

UCSC Genome Browser assembly hub visualizing RNA modifications across the human
"RNome": poly-A RNA (hg38), ribosomal RNA (rRNA, 18S/28S/5.8S/5S), and transfer
RNA (tRNA). Modifications are detected by three modalities — short-read
sequencing (SRS, Illumina), long-read sequencing (LRS, Oxford Nanopore/ONT),
and mass spectrometry (MS) — each provided as a single pre-merged,
genome-mapped consensus BED file. Each modality is
rendered as one flat track per RNA type (no per-experiment subtracks); track
descriptions are sourced directly from `RNA_modifications_manifest.tsv`.

---

## Repo structure

```
RNome_dataviz_repo/
├── RNA_modifications_manifest.tsv          # Source of truth: Filename → Description, per modality
├── final_bedRmods/                         # Pipeline input: 3 pre-merged, genome-mapped BED files
│   ├── Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed
│   ├── ONT_polyARNA_rRNA_combined.filtered_rmchrY.bed
│   └── MS_rRNA_tRNA.bed
│
├── reference/                              # Reference sequences & annotations
│   ├── hs_rRNAs_NR_046235.fa(.fai)            # Human rRNA reference (18S/28S/5.8S/5S)
│   ├── Step3_hg38_Homo_sapiens_47seq...fa     # tRNA reference sequences
│   └── gencode.v49...annotation.gtf.gz        # GENCODE hg38 annotation (only needed by BedPyLift.py)
│
├── final_data/                             # Pipeline output, organized by RNA type
│   ├── RNA_modifications_manifest.xlsx        # Generated manifest of all tracks + hub URLs
│   ├── polyA-RNA_hg38/
│   │   └── polyA-RNA_hg38_{SRS,LRS}.bed
│   ├── rRNA/
│   │   └── rRNA_{SRS,LRS,MS}.bed
│   └── tRNA/
│       └── tRNA_{SRS,MS}.bed
│
├── scripts/
│   ├── pipeline.py                         # End-to-end hub generation (see below)
│   └── BedPyLift.py                        # Transcriptome → genome coordinate lifting (not
│                                            # needed for this hub — final_bedRmods/ is already
│                                            # genome-mapped; kept for future transcript-coord inputs)
│
├── ucsc_hub/                                # Generated UCSC assembly hub
│   ├── hub.txt / genomes.txt / hubDescription.html
│   ├── hg38_polyA-RNA/                        # hg38 assembly: trackDb, bigBeds, HTML
│   ├── rrna/                                  # hs_rRNA custom assembly: 2bit, trackDb, bigBeds
│   └── trna/                                  # hs_tRNA custom assembly: 2bit, trackDb, bigBeds
│
├── utils/                                    # UCSC command-line tools
│   ├── bedToBigBed
│   ├── faToTwoBit
│   └── hg38.chrom.sizes
│
├── *_old/                                    # Previous hub generation (experiment-level subtracks);
│                                              # kept for reference, not used by the current pipeline
└── README.md
```

---

## Pipeline

### Requirements

- Python 3.7+
- `openpyxl` (optional — used for the Excel manifest; falls back to a TSV if missing)
- UCSC command line utilities (already included in `utils/` for macOS Apple Silicon):

```bash
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit
chmod +x bedToBigBed faToTwoBit
```

### Run the pipeline (`pipeline.py`)

The 3 files in `final_bedRmods/` (SRS/Illumina, LRS/ONT, MS) are already in
genomic coordinates. From the repo root:

```bash
python scripts/pipeline.py
```

This will:

1. Split each combined SRS/LRS/MS BED file by RNA type (`rRNA`, `tRNA`,
   `polyA-RNA_hg38`) into `final_data/{RNA_type}/{RNA_type}_{Modality}.bed`
2. Normalize each split BED file for `bedToBigBed`: clamp scores to 0–1000,
   recolor by modification type (hue) × frequency (saturation), remap
   Ensembl scaffold names to UCSC names where needed, and convert to bigBed
   (`bedToBigBed -type=bed9+4 -as=rna_mods.as`)
3. Generate one HTML description per `{RNA_type}_{Modality}` track, using the
   matching modality's description text from `RNA_modifications_manifest.tsv`
4. Generate `final_data/RNA_modifications_manifest.xlsx` (or `.tsv` if
   `openpyxl` is unavailable) listing every track and its hub bigBed URL
5. Generate the UCSC hub config (`hub.txt`, `genomes.txt`, per-assembly
   `trackDb.txt`) — one flat track per RNA type/modality, each with
   score / coverage / frequency / modification-type filters

To change a track's description, edit `RNA_modifications_manifest.tsv`
(matched to `COMBINED_FILES` by filename) and re-run the pipeline.

### One-time setup — custom assembly sequences (.2bit)

Needed once for the custom `hs_rRNA` / `hs_tRNA` assemblies (already present
under `ucsc_hub/rrna/` and `ucsc_hub/trna/`; regenerate only if the reference
sequences change):

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
| Poly-A RNA | `hg38` | SRS, LRS | Standard hg38; default position `chr21:8400001-8500000` |
| rRNA | `hs_rRNA` | SRS, LRS, MS | Custom 2bit assembly (18S/28S/5.8S/5S) |
| tRNA | `hs_tRNA` | SRS, MS | Custom 2bit assembly |

Each present modality is a single flat consensus track
(`{RNA_type}_{Modality}.bigBed`) — all sites for that RNA type/modality,
colored and filterable by modification type, score, coverage, and frequency.
There are no per-experiment subtracks in this hub.

The hg38 assembly additionally carries an externally hosted track,
**NA12878 variants** (WGS merged variant calls from two callers), streamed
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

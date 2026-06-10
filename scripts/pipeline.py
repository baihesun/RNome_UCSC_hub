"""
process_files.py — v2 data organization and UCSC hub generation

Steps:
  1. Split combined BED files by RNA type → final_data/{RNA_type}/{RNA_type}_{Modality}.bed
  2. Copy raw experiment files → final_data/{RNA_type}/experiments_{modality}/
  3. Convert split BED files to bigBed
  4. Generate UCSC hub v2 (ucsc_hub/) with hub.txt, genomes.txt, per-assembly trackDb
  5. Generate per-file HTML readmes (placeholder descriptions)
  6. Generate Excel manifest

RNA type classification:
  chrom contains "rRNA"  →  rRNA
  chrom contains "tRNA"  →  tRNA
  else                   →  polyA-RNA_hg38

Run from repo root:
  python scripts/process_files.py
"""

import colorsys
import os
import re
import shutil
import subprocess

# ── Source data ────────────────────────────────────────────────────────────────

SOURCE_DIR = "/Users/baihesun/Desktop/RNome_final_bedRmods"

COMBINED_FILES = {
    "Illumina": os.path.join(SOURCE_DIR, "illumina",  "Illumina_combined_polyARNA_tRNA_rRNA.bed"),
    "ONT":      os.path.join(SOURCE_DIR, "ont",       "ONT_polyARNA_rRNA_combined.filtered.bed"),
    "MS":       os.path.join(SOURCE_DIR, "mass_spec", "MS_rRNA_tRNA.bed"),
}

# {RNA_type: {Modality: [src_path, ...]}}
EXPERIMENT_FILES = {
    "rRNA": {
        "Illumina": [
            os.path.join(SOURCE_DIR, "illumina",  "raw", "rRNA_Illumina.bed"),
        ],
        "ONT": [
            os.path.join(SOURCE_DIR, "ont",       "raw", "rRNA.native_vs_ivt_fisher.filtered.bed"),
        ],
        "MS": [
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "HRP_C_001_rRNA_001.bed"),
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "HRP_C_002_rRNA_001.bed"),
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "HRP_C_005_rRNA_003+Kaiser_20260213_28S_rRNA_3.5uL.bed"),
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "Limbach_group_Protiowizard convert_rRNA_T1_G2.bed"),
        ],
    },
    "tRNA": {
        "Illumina": [
            os.path.join(SOURCE_DIR, "illumina",  "raw", "combined_YM+TP_bedRmod_RNOME_tRNA_noScoreAndFrequency.bed"),
        ],
        "MS": [
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "HRP_C_001_tRNA_003.bed"),
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "Kaiser_20260206_tRNA_1000ng_1.bed"),
            os.path.join(SOURCE_DIR, "mass_spec", "raw", "Limbach_group_Proteowizard convert tRNA_HRP.bed"),
        ],
    },
    "polyA-RNA_hg38": {
        "Illumina": [
            os.path.join(SOURCE_DIR, "illumina",  "raw", "polyA_m6A.bed"),
            os.path.join(SOURCE_DIR, "illumina",  "raw", "polyA_m5C.bed"),
            os.path.join(SOURCE_DIR, "illumina",  "raw", "polyA_psiU.bed"),
            os.path.join(SOURCE_DIR, "illumina",  "raw", "Inosine_JACUSA_REDItools.bed"),
        ],
        "ONT": [
            os.path.join(SOURCE_DIR, "ont",       "raw", "polyA.native_vs_ivt_fisher.filtered.bed"),
        ],
    },
}

# ── Output paths ───────────────────────────────────────────────────────────────

FINAL_DATA_DIR = "final_data"
HUB_DIR    = "ucsc_hub"
UTILS_DIR      = "utils"

# ── Reference files ────────────────────────────────────────────────────────────

RRNA_FAI    = os.path.join("reference", "hs_rRNAs_NR_046235.fa.fai")
RRNA_2BIT   = os.path.join("ucsc_hub",  "rrna", "hs_rRNAs_NR_046235.2bit")
HG38_SIZES  = os.path.join(UTILS_DIR,   "hg38.chrom.sizes")
BEDTOBIGBED = os.path.join(UTILS_DIR,   "bedToBigBed")

TRNA_2BIT   = os.path.join(HUB_DIR, "trna", "hs_tRNA.2bit")
TRNA_SIZES  = os.path.join(HUB_DIR, "trna", "hs_tRNA.chrom.sizes")

# ── Hub metadata ───────────────────────────────────────────────────────────────

GITHUB_USER = "baihesun"
GITHUB_REPO = "RNome_UCSC_hub"
HUB_BASE_URL = f"https://raw.githubusercontent.com/{GITHUB_USER}/{GITHUB_REPO}/main"
PAPER_URL    = "PLACEHOLDER_DOI"  # replace with actual DOI/URL

# ── Modification colors (hue per mod type, saturation ∝ frequency) ─────────────

COLOR_CODE = {
    "m6A":     "#003f5c",
    "Am":      "#1464a0",
    "Ino":     "#0e8fd4",
    "I":       "#0e8fd4",
    "m1A":     "#00b4d8",
    "m66A":    "#00d4b4",
    "m1acp3Y": "#00c875",
    "Y":       "#90be6d",
    "Ym":      "#c9b84c",
    "Um":      "#ffa600",
    "m3U":     "#ff7c43",
    "Gm":      "#ff6361",
    "m7G":     "#dd7371",
    "Cm":      "#ef476f",
    "m5C":     "#d45087",
    "ac4C":    "#7b2d8b",
}
DEFAULT_COLOR = "#808080"

RNA_TYPES  = ["polyA-RNA_hg38", "tRNA", "rRNA"]
MODALITIES = ["Illumina", "ONT", "MS"]

# ── Assembly config ────────────────────────────────────────────────────────────

ASSEMBLY_CFG = {
    "rRNA": {
        "ucsc_name":  "hs_rRNA",
        "hub_dir":    "rrna",
        "2bit_src":   RRNA_2BIT,
        "2bit_dest":  "hs_rRNAs_NR_046235.2bit",
        "sizes":      None,          # derived from FAI at runtime
        "label":      "rRNA",
        "organism":   "Homo sapiens rRNA",
    },
    "polyA-RNA_hg38": {
        "ucsc_name":  "hg38",
        "hub_dir":    "hg38_polyA-RNA",
        "2bit_src":   None,          # standard UCSC assembly
        "2bit_dest":  None,
        "sizes":      HG38_SIZES,
        "label":      "poly-A RNA",
        "organism":   "Homo sapiens",
    },
    "tRNA": {
        "ucsc_name":  "hs_tRNA",
        "hub_dir":    "trna",
        "2bit_src":   TRNA_2BIT,
        "2bit_dest":  "hs_tRNA.2bit",
        "sizes":      TRNA_SIZES,
        "label":      "tRNA",
        "organism":   "Homo sapiens tRNA",
    },
}

# ── AutoSQL schema (bed9 + 4 extra fields) ─────────────────────────────────────

AUTOSQL_SCHEMA = """\
table rna_mods
"RNA modification sites"
(
string  chrom;              "Chromosome or RNA name"
uint    chromStart;         "Start position (0-based)"
uint    chromEnd;           "End position"
string  name;               "Modification short name"
uint    score;              "Detection score (0-1000)"
char[1] strand;             "Strand (+ or -)"
uint    thickStart;         "Same as chromStart"
uint    thickEnd;           "Same as chromEnd"
uint    itemRgb;            "Display color (R,G,B)"
uint    coverage;           "Read coverage"
float   frequency;          "Modification frequency (%)"
string  single_letter_code; "Single-letter modification code"
lstring mod_id;             "Modification identifier"
)
"""

# ── Placeholder descriptions (update with actual text for publication) ─────────

FILE_DESCRIPTIONS = {
    "polyA-RNA_hg38_Illumina.bed": (
        "Illumina sequencing-based RNA modifications in human polyadenylated RNA, "
        "mapped to the hg38 reference genome."
    ),
    "polyA-RNA_hg38_ONT.bed": (
        "Oxford Nanopore sequencing-based RNA modifications in human polyadenylated RNA, "
        "mapped to the hg38 reference genome."
    ),
    "rRNA_Illumina.bed": (
        "Illumina sequencing-based RNA modifications in human ribosomal RNA "
        "(18S, 28S, 5.8S, 5S)."
    ),
    "rRNA_ONT.bed": (
        "Oxford Nanopore sequencing-based RNA modifications in human ribosomal RNA "
        "(18S, 28S, 5.8S, 5S)."
    ),
    "rRNA_MS.bed": (
        "Mass spectrometry-based RNA modifications in human ribosomal RNA "
        "(18S, 28S, 5.8S, 5S)."
    ),
    "tRNA_Illumina.bed": (
        "Illumina sequencing-based RNA modifications in human transfer RNA."
    ),
    "tRNA_MS.bed": (
        "Mass spectrometry-based RNA modifications in human transfer RNA."
    ),
}


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def mod_to_rgb(mod_type, frequency):
    """Return 'R,G,B' string: hue fixed per mod type, saturation ∝ frequency."""
    hex_color = COLOR_CODE.get(mod_type.strip(), DEFAULT_COLOR).lstrip("#")
    r = int(hex_color[0:2], 16) / 255
    g = int(hex_color[2:4], 16) / 255
    b = int(hex_color[4:6], 16) / 255
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s_scaled = s * max(0.0, min(100.0, float(frequency))) / 100.0
    r2, g2, b2 = colorsys.hsv_to_rgb(h, s_scaled, v)
    return f"{int(r2 * 255)},{int(g2 * 255)},{int(b2 * 255)}"


def clamp_score(val):
    try:
        return str(int(max(0, min(1000, float(val)))))
    except (ValueError, TypeError):
        return "0"


def build_ensembl_to_ucsc_map(sizes_path):
    """
    Parse hg38.chrom.sizes and return {ensembl_accession: ucsc_name}.

    UCSC embeds the accession in names like chrUn_GL000195v1 and
    chr6_GL000251v2_alt.  The regex extracts e.g. GL000195 + v1 → GL000195.1.
    Standard chr* entries have no version suffix and are not added to the map.
    """
    mapping = {}
    pat = re.compile(r'([A-Z]{2}\d+)v(\d+)')
    with open(sizes_path) as fh:
        for line in fh:
            ucsc = line.split('\t')[0].strip()
            m = pat.search(ucsc)
            if m:
                ensembl = f"{m.group(1)}.{m.group(2)}"
                mapping[ensembl] = ucsc
    return mapping


def fai_to_sizes(fai_path):
    """Parse .fai → {chrom: size}."""
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2:
                sizes[parts[0].strip()] = int(parts[1])
    return sizes


def write_sizes_file(sizes_dict, out_path):
    with open(out_path, "w") as fh:
        for chrom, size in sorted(sizes_dict.items()):
            fh.write(f"{chrom}\t{size}\n")


def classify_rna_type(chrom):
    if "rRNA" in chrom:
        return "rRNA"
    if "tRNA" in chrom:
        return "tRNA"
    return "polyA-RNA_hg38"


def safe_id(s):
    """Return a UCSC-safe track identifier (alphanumeric + underscore, ≤60 chars)."""
    s = re.sub(r'[^A-Za-z0-9]', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    if s and s[0].isdigit():
        s = 'x' + s
    return s[:60]


def get_mod_names(bed_path):
    """Collect sorted unique modification names (col 3) from a BED file."""
    names = set()
    try:
        with open(bed_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) >= 4:
                    names.add(cols[3].strip())
    except OSError:
        pass
    return ",".join(sorted(names))


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 — Split combined BED files
# ══════════════════════════════════════════════════════════════════════════════

def split_combined_files():
    print("\n── Step 1: Splitting combined files ──")
    for modality, src_path in COMBINED_FILES.items():
        if not os.path.exists(src_path):
            print(f"  [SKIP] not found: {src_path}")
            continue
        print(f"  Reading {os.path.basename(src_path)} ...")
        buckets = {rna: [] for rna in RNA_TYPES}
        with open(src_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                chrom = line.split("\t")[0]
                buckets[classify_rna_type(chrom)].append(line)
        for rna_type, lines in buckets.items():
            if not lines:
                continue
            out_dir = os.path.join(FINAL_DATA_DIR, rna_type)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"{rna_type}_{modality}.bed")
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")
            print(f"    {len(lines):>7,} lines → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — Copy raw experiment files
# ══════════════════════════════════════════════════════════════════════════════

def copy_experiment_files():
    print("\n── Step 2: Copying experiment files ──")
    for rna_type, modality_map in EXPERIMENT_FILES.items():
        for modality, src_paths in modality_map.items():
            dest_dir = os.path.join(FINAL_DATA_DIR, rna_type, f"experiments_{modality.lower()}")
            os.makedirs(dest_dir, exist_ok=True)
            for src in src_paths:
                if not os.path.exists(src):
                    print(f"  [SKIP] not found: {src}")
                    continue
                dest_name = os.path.basename(src).replace(" ", "_")
                dest = os.path.join(dest_dir, dest_name)
                shutil.copy2(src, dest)
                print(f"  Copied → {dest}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — Convert split BED → bigBed
# ══════════════════════════════════════════════════════════════════════════════

def process_bed_for_bigbed(src_path, dst_path, chrom_remap=None, valid_chroms=None):
    """
    Normalise a 13-col bedRmod file for bedToBigBed:
      - Score clamped to 0-1000
      - itemRgb recomputed from mod name + frequency
      - Exactly 13 columns, sorted by chrom then start
      - chrom_remap:   dict to rename chromosomes (Ensembl → UCSC); rows
                       whose non-chr chrom is absent are dropped
      - valid_chroms:  set of allowed chrom names (checked after remap);
                       rows with other chroms are dropped
    """
    rows = []
    skipped = 0
    with open(src_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            if chrom_remap is not None:
                chrom = cols[0]
                if not chrom.startswith("chr"):
                    ucsc = chrom_remap.get(chrom)
                    if ucsc is None:
                        skipped += 1
                        continue  # unmapped scaffold — skip
                    cols[0] = ucsc
            if valid_chroms is not None and cols[0] not in valid_chroms:
                skipped += 1
                continue
            while len(cols) < 13:
                cols.append(".")
            mod_name = cols[3].strip()
            try:
                freq = float(cols[10])
            except ValueError:
                freq = 0.0
            cols[4] = clamp_score(cols[4])
            cols[5] = cols[5] if cols[5] in ("+", "-") else "+"
            try:
                cols[6] = str(int(float(cols[6])))
                cols[7] = str(int(float(cols[7])))
            except ValueError:
                pass
            cols[8] = mod_to_rgb(mod_name, freq)
            try:
                cols[9] = str(int(float(cols[9])))
            except ValueError:
                cols[9] = "0"
            rows.append(cols[:13])
    rows.sort(key=lambda c: (c[0], int(c[1])))
    with open(dst_path, "w", encoding="utf-8") as fh:
        for cols in rows:
            fh.write("\t".join(cols) + "\n")
    return len(rows), skipped


def write_autosql(out_dir):
    as_path = os.path.join(out_dir, "rna_mods.as")
    with open(as_path, "w") as fh:
        fh.write(AUTOSQL_SCHEMA)
    return as_path


def run_bedtobigbed(bed_path, sizes_path, as_path, bigbed_path):
    cmd = [
        BEDTOBIGBED,
        "-type=bed9+4",
        "-tab",
        f"-as={as_path}",
        bed_path,
        sizes_path,
        bigbed_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] bedToBigBed: {result.stderr.strip()}")
        return False
    return True


def convert_split_to_bigbed():
    print("\n── Step 3: Converting split BED → bigBed ──")

    # Build rRNA chrom sizes from FAI
    rrna_sizes_path = os.path.join(HUB_DIR, "rrna", "hs_rRNA.chrom.sizes")
    if os.path.exists(RRNA_FAI):
        rrna_sizes = fai_to_sizes(RRNA_FAI)
        os.makedirs(os.path.dirname(rrna_sizes_path), exist_ok=True)
        write_sizes_file(rrna_sizes, rrna_sizes_path)
    else:
        print(f"  [WARN] rRNA FAI not found ({RRNA_FAI}); rRNA bigBed skipped")
        rrna_sizes_path = None

    # Build Ensembl→UCSC name map for hg38 scaffolds
    hg38_remap = None
    if os.path.exists(HG38_SIZES):
        hg38_remap = build_ensembl_to_ucsc_map(HG38_SIZES)
        print(f"  Built hg38 chrom remap: {len(hg38_remap)} scaffold entries")

    sizes_for = {
        "rRNA":            rrna_sizes_path,
        "polyA-RNA_hg38":  HG38_SIZES if os.path.exists(HG38_SIZES) else None,
        "tRNA":            TRNA_SIZES,  # None until tRNA reference is provided
    }

    for rna_type in RNA_TYPES:
        cfg = ASSEMBLY_CFG[rna_type]
        hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
        os.makedirs(hub_dir, exist_ok=True)
        sizes_path = sizes_for[rna_type]

        if sizes_path is None:
            print(f"  [TODO] No chrom sizes for {rna_type} — bigBed skipped")
            continue

        as_path = write_autosql(hub_dir)

        for modality in MODALITIES:
            src_bed    = os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_{modality}.bed")
            fixed_bed  = src_bed.replace(".bed", ".fixed.bed")
            bigbed_out = os.path.join(hub_dir, f"{rna_type}_{modality}.bigBed")

            if not os.path.exists(src_bed):
                continue

            remap = hg38_remap if rna_type == "polyA-RNA_hg38" else None
            valid  = set(fai_to_sizes(sizes_path).keys()) if sizes_path else None
            n, skipped = process_bed_for_bigbed(src_bed, fixed_bed,
                                                chrom_remap=remap,
                                                valid_chroms=valid)
            skip_str = f" ({skipped} skipped)" if skipped else ""
            print(f"  Processed {n:>6,} rows{skip_str}: {os.path.basename(src_bed)}")

            if run_bedtobigbed(fixed_bed, sizes_path, as_path, bigbed_out):
                print(f"  bigBed → {bigbed_out}")

            if os.path.exists(fixed_bed):
                os.remove(fixed_bed)


# ══════════════════════════════════════════════════════════════════════════════
# Step 3b — Convert experiment BED files to bigBed
# ══════════════════════════════════════════════════════════════════════════════

def _exp_bigbeds(rna_type, modality, cfg):
    """Return sorted list of (stem, filename) for experiment bigBed files."""
    exp_dir = os.path.join(HUB_DIR, cfg["hub_dir"], f"experiments_{modality.lower()}")
    if not os.path.isdir(exp_dir):
        return []
    return sorted(
        (os.path.splitext(f)[0], f)
        for f in os.listdir(exp_dir)
        if f.endswith(".bigBed")
    )


def convert_experiments_to_bigbed():
    print("\n── Step 3b: Converting experiment BED → bigBed ──")

    hg38_remap = build_ensembl_to_ucsc_map(HG38_SIZES) if os.path.exists(HG38_SIZES) else {}

    rrna_sizes_path = os.path.join(HUB_DIR, "rrna", "hs_rRNA.chrom.sizes")
    rrna_valid  = set(fai_to_sizes(RRNA_FAI).keys()) if os.path.exists(RRNA_FAI) else None
    trna_valid  = set(fai_to_sizes(TRNA_SIZES).keys()) if TRNA_SIZES and os.path.exists(TRNA_SIZES) else None

    sizes_for = {
        "rRNA":            rrna_sizes_path if os.path.exists(rrna_sizes_path) else None,
        "polyA-RNA_hg38":  HG38_SIZES     if os.path.exists(HG38_SIZES)      else None,
        "tRNA":            TRNA_SIZES     if TRNA_SIZES and os.path.exists(TRNA_SIZES) else None,
    }
    valid_for = {"rRNA": rrna_valid, "polyA-RNA_hg38": None, "tRNA": trna_valid}
    remap_for = {"rRNA": None, "polyA-RNA_hg38": hg38_remap, "tRNA": None}

    for rna_type, modality_map in EXPERIMENT_FILES.items():
        cfg     = ASSEMBLY_CFG[rna_type]
        sizes   = sizes_for[rna_type]
        as_path = os.path.join(HUB_DIR, cfg["hub_dir"], "rna_mods.as")

        if not sizes:
            print(f"  [SKIP] {rna_type} — no chrom sizes for experiment bigBed")
            continue

        for modality, src_paths in modality_map.items():
            exp_hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"], f"experiments_{modality.lower()}")
            os.makedirs(exp_hub_dir, exist_ok=True)

            for src in src_paths:
                fname   = os.path.basename(src).replace(" ", "_")
                src_bed = os.path.join(FINAL_DATA_DIR, rna_type, f"experiments_{modality.lower()}", fname)
                if not os.path.exists(src_bed):
                    continue
                stem       = os.path.splitext(fname)[0]
                fixed_bed  = src_bed.replace(".bed", ".fixed.bed")
                bigbed_out = os.path.join(exp_hub_dir, f"{stem}.bigBed")

                n, skipped = process_bed_for_bigbed(
                    src_bed, fixed_bed,
                    chrom_remap=remap_for[rna_type],
                    valid_chroms=valid_for[rna_type],
                )
                skip_str = f" ({skipped} skipped)" if skipped else ""
                if run_bedtobigbed(fixed_bed, sizes, as_path, bigbed_out):
                    print(f"  bigBed ({n} rows{skip_str}) → {bigbed_out}")
                if os.path.exists(fixed_bed):
                    os.remove(fixed_bed)


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — HTML readmes
# ══════════════════════════════════════════════════════════════════════════════

def write_html_readme(filename, description, out_dir):
    stem      = os.path.splitext(filename)[0]
    html_path = os.path.join(out_dir, f"{stem}.html")
    content = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{stem}</title>
</head>
<body>
<h2>{stem}</h2>
<p>{description}</p>

<p>
  <strong>Reference:</strong>
  <a href="{PAPER_URL}">{PAPER_URL}</a>
</p>

</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return html_path


def write_experiments_html(rna_type, modality, exp_filenames, hub_dir):
    stem      = f"{rna_type}_{modality}_experiments"
    html_path = os.path.join(hub_dir, f"{stem}.html")
    desc = (
        f"Individual {modality} sequencing experiment tracks for {rna_type} RNA modifications. "
        f"Each subtrack represents one experiment or sample."
    )
    items = "".join(f"    <li>{f}</li>\n" for f in exp_filenames)
    content = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{rna_type} {modality} Experiments</title>
</head>
<body>
<h2>{rna_type} — {modality} Experiments</h2>
<p>{desc}</p>
<p>
  <strong>Reference:</strong>
  <a href="{PAPER_URL}">{PAPER_URL}</a>
</p>
<h3>Included experiments</h3>
<ul>
{items}</ul>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return html_path


def generate_html_readmes():
    print("\n── Step 4: Generating HTML readmes ──")
    for rna_type in RNA_TYPES:
        cfg     = ASSEMBLY_CFG[rna_type]
        hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
        os.makedirs(hub_dir, exist_ok=True)
        for modality in MODALITIES:
            filename = f"{rna_type}_{modality}.bed"
            bed_path = os.path.join(FINAL_DATA_DIR, rna_type, filename)
            if not os.path.exists(bed_path):
                continue
            desc = FILE_DESCRIPTIONS.get(
                filename,
                f"RNA modification sites — {rna_type}, {modality}. (PLACEHOLDER description)"
            )
            html_path = write_html_readme(filename, desc, hub_dir)
            print(f"  HTML → {html_path}")

            # Experiments composite HTML
            exp_src_paths = EXPERIMENT_FILES.get(rna_type, {}).get(modality, [])
            if exp_src_paths:
                exp_names = [os.path.basename(p).replace(" ", "_") for p in exp_src_paths]
                exp_html  = write_experiments_html(rna_type, modality, exp_names, hub_dir)
                print(f"  HTML → {exp_html}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 — Excel manifest
# ══════════════════════════════════════════════════════════════════════════════

def _manifest_rows():
    """Yield (filename, rna_type, modality, track_type, description, hub_url) for every track."""
    for rna_type in RNA_TYPES:
        cfg = ASSEMBLY_CFG[rna_type]
        # Consensus rows
        for modality in MODALITIES:
            filename = f"{rna_type}_{modality}.bed"
            if not os.path.exists(os.path.join(FINAL_DATA_DIR, rna_type, filename)):
                continue
            hub_url = f"{HUB_BASE_URL}/{HUB_DIR}/{cfg['hub_dir']}/{rna_type}_{modality}.bigBed"
            desc    = FILE_DESCRIPTIONS.get(filename, f"RNA modification sites — {rna_type}, {modality}.")
            yield filename, rna_type, modality, "consensus", desc, hub_url
        # Experiment rows
        for modality, src_paths in EXPERIMENT_FILES.get(rna_type, {}).items():
            for src in src_paths:
                fname   = os.path.basename(src).replace(" ", "_")
                stem    = os.path.splitext(fname)[0]
                src_bed = os.path.join(FINAL_DATA_DIR, rna_type, f"experiments_{modality.lower()}", fname)
                if not os.path.exists(src_bed):
                    continue
                hub_url = (f"{HUB_BASE_URL}/{HUB_DIR}/{cfg['hub_dir']}/"
                           f"experiments_{modality.lower()}/{stem}.bigBed")
                desc    = f"Individual {modality} experiment — {rna_type} RNA modifications. (PLACEHOLDER)"
                yield fname, rna_type, modality, "experiment", desc, hub_url


def generate_excel_manifest():
    print("\n── Step 5: Generating Excel manifest ──")
    try:
        from openpyxl import Workbook
    except ImportError:
        print("  [SKIP] openpyxl not installed — run: pip install openpyxl")
        _write_tsv_manifest()
        return

    wb = Workbook()
    ws = wb.active
    ws.title = "Manifest"
    ws.append(["Filename", "RNA Type", "Modality", "Track Type", "Description", "Hub BigBed URL", "Paper"])
    for row in _manifest_rows():
        ws.append([*row, PAPER_URL])

    out_path = os.path.join(FINAL_DATA_DIR, "RNA_modifications_manifest.xlsx")
    wb.save(out_path)
    print(f"  Manifest → {out_path}")


def _write_tsv_manifest():
    """Fallback TSV manifest when openpyxl is unavailable."""
    out_path = os.path.join(FINAL_DATA_DIR, "RNA_modifications_manifest.tsv")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Filename\tRNA Type\tModality\tTrack Type\tDescription\tHub BigBed URL\tPaper\n")
        for row in _manifest_rows():
            fh.write("\t".join(row) + f"\t{PAPER_URL}\n")
    print(f"  Manifest (TSV) → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 6 — UCSC hub v2 config
# ══════════════════════════════════════════════════════════════════════════════

MODALITY_COLOR = {
    "Illumina": "0,120,180",
    "ONT":      "180,80,0",
    "MS":       "60,160,60",
}


def write_hub_txt():
    desc_html = os.path.join(HUB_DIR, "hubDescription.html")
    if not os.path.exists(desc_html):
        with open(desc_html, "w") as fh:
            fh.write("""\
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>Human RNome RNA Modification Atlas</title></head>
<body>
<h2>Human RNome RNA Modification Atlas (v2)</h2>
<p>
  A comprehensive map of RNA modifications across the human transcriptome,
  integrating Illumina sequencing, Oxford Nanopore (ONT), and mass spectrometry (MS) data
  for poly-A RNA, ribosomal RNA (rRNA), and transfer RNA (tRNA).
</p>
<p><strong>Reference:</strong> PLACEHOLDER_DOI</p>
</body>
</html>
""")
    content = f"""\
hub RNome
shortLabel RNome Modifications v2
longLabel Human RNome RNA Modification Atlas (v2)
genomesFile genomes.txt
email baihe_sun@mail.dfci.harvard.edu
descriptionUrl hubDescription.html
"""
    path = os.path.join(HUB_DIR, "hub.txt")
    with open(path, "w") as fh:
        fh.write(content)
    print(f"  hub.txt → {path}")


def _write_assembly_html(hub_dir, label):
    """Write a minimal assembly description HTML for custom assembly hubs."""
    html_path = os.path.join(HUB_DIR, hub_dir, "assemblyDescription.html")
    if not os.path.exists(html_path):
        with open(html_path, "w") as fh:
            fh.write(f"""\
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>{label} Assembly</title></head>
<body>
<h2>{label} — Custom Assembly</h2>
<p>
  Custom reference sequences for human {label} used in the RNome RNA Modification Atlas.
</p>
<p><strong>Reference:</strong> PLACEHOLDER_DOI</p>
</body>
</html>
""")


def write_genomes_txt():
    lines = []
    for rna_type in RNA_TYPES:
        cfg       = ASSEMBLY_CFG[rna_type]
        ucsc_name = cfg["ucsc_name"]
        hub_dir   = cfg["hub_dir"]
        label     = cfg["label"]
        organism  = cfg["organism"]

        if rna_type == "polyA-RNA_hg38":
            # Standard UCSC assembly — no custom 2bit needed
            block = (
                f"genome {ucsc_name}\n"
                f"trackDb {hub_dir}/trackDb.txt\n"
                f"defaultPos chr21:8400001-8500000\n"
                f"organism \"Homo sapiens\"\n"
                f"description \"Human RNA modifications (poly-A RNA, hg38)\"\n"
            )
        elif rna_type == "tRNA":
            _write_assembly_html(hub_dir, label)
            block = (
                f"genome {ucsc_name}\n"
                f"trackDb {hub_dir}/trackDb.txt\n"
                f"twoBitPath {hub_dir}/hs_tRNA.2bit\n"
                f"defaultPos hs_tRNAAla_AGC:1-76\n"
                f"htmlPath {hub_dir}/assemblyDescription.html\n"
                f"organism \"{organism}\"\n"
                f"description \"{label}\"\n"
                f"scientificName \"Homo sapiens\"\n"
            )
        else:  # rRNA
            twobit_dest = cfg["2bit_dest"]
            _write_assembly_html(hub_dir, label)
            block = (
                f"genome {ucsc_name}\n"
                f"trackDb {hub_dir}/trackDb.txt\n"
                f"twoBitPath {hub_dir}/{twobit_dest}\n"
                f"defaultPos hs_rRNA_18S:1-1869\n"
                f"htmlPath {hub_dir}/assemblyDescription.html\n"
                f"organism \"{organism}\"\n"
                f"description \"{label}\"\n"
                f"scientificName \"Homo sapiens\"\n"
            )
        lines.append(block)

    path = os.path.join(HUB_DIR, "genomes.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"  genomes.txt → {path}")


def _bigbed_exists(rna_type, modality, cfg):
    return os.path.exists(
        os.path.join(HUB_DIR, cfg["hub_dir"], f"{rna_type}_{modality}.bigBed")
    )


def _bed_exists(rna_type, modality):
    return os.path.exists(
        os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_{modality}.bed")
    )


# Shared filter directives for consensus and experiment tracks
_CONSENSUS_FILTER = (
    "scoreMin 0\n"
    "scoreMax 1000\n"
    "filter.frequency 0\n"
    "filterByRange.frequency on\n"
    "filterLimits.frequency 0:100\n"
    "filterLabel.frequency Modification frequency (%)\n"
    "filterLabel.name Modification type\n"
)
_EXP_FILTER = "".join(f"\t{ln}\n" for ln in _CONSENSUS_FILTER.splitlines())


def write_trackdb(rna_type):
    cfg     = ASSEMBLY_CFG[rna_type]
    hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
    os.makedirs(hub_dir, exist_ok=True)
    path    = os.path.join(hub_dir, "trackDb.txt")

    present = [m for m in MODALITIES if _bed_exists(rna_type, m)]
    if not present:
        with open(path, "w") as fh:
            fh.write(f"# No data available for {rna_type}\n")
        print(f"  trackDb (empty) → {path}")
        return

    safe_rna = safe_id(rna_type)
    stanzas  = []

    for priority, modality in enumerate(present, 1):
        super_tid    = f"{safe_rna}_{modality}_super"
        consensus_tid = f"{safe_rna}_{modality}_consensus"
        bigbed_fname  = f"{rna_type}_{modality}.bigBed"
        html_ref      = f"{rna_type}_{modality}"
        bed_path      = os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_{modality}.bed")
        mod_names     = get_mod_names(bed_path)

        # ── SuperTrack ────────────────────────────────────────────────────────
        stanzas.append(
            f"track {super_tid}\n"
            f"superTrack on show\n"
            f"shortLabel {modality}\n"
            f"longLabel {rna_type} — {modality}\n"
            f"html {html_ref}\n"
            f"priority {priority}\n"
        )

        # ── Consensus track ───────────────────────────────────────────────────
        if _bigbed_exists(rna_type, modality, cfg):
            stanzas.append(
                f"track {consensus_tid}\n"
                f"parent {super_tid}\n"
                f"bigDataUrl {bigbed_fname}\n"
                f"shortLabel {modality} consensus\n"
                f"longLabel {rna_type} {modality} modifications\n"
                f"type bigBed 9 +\n"
                f"itemRgb on\n"
                f"visibility pack\n"
                + _CONSENSUS_FILTER
                + f"filterValues.name {mod_names}\n"
                  f"html {html_ref}\n"
                  f"priority 1\n"
            )
        else:
            stanzas.append(f"# TODO: {bigbed_fname} not yet generated\n")

        # ── Experiments composite ─────────────────────────────────────────────
        exp_tracks = _exp_bigbeds(rna_type, modality, cfg)
        if not exp_tracks:
            continue

        exp_comp_tid = f"{safe_rna}_{modality}_experiments_composite"
        exp_html_ref = f"{rna_type}_{modality}_experiments"

        sg_entries = " ".join(
            f"{safe_id(stem).lower()}={stem[:20].replace(' ', '_')}"
            for stem, _ in exp_tracks
        )

        stanzas.append(
            f"track {exp_comp_tid}\n"
            f"compositeTrack on\n"
            f"parent {super_tid}\n"
            f"shortLabel {modality} experiments\n"
            f"longLabel Individual {modality} Experiment Tracks\n"
            f"type bigBed 9 +\n"
            f"itemRgb on\n"
            f"visibility hide\n"
            f"subGroup1 sample Sample {sg_entries}\n"
            f"dimensions dimX=sample\n"
            f"filterComposite on\n"
            f"html {exp_html_ref}\n"
            f"priority 2\n"
        )

        for j, (stem, fname) in enumerate(exp_tracks, 1):
            exp_tid   = f"{safe_rna}_{modality}_exp_{safe_id(stem)}"
            sub_key   = safe_id(stem).lower()
            rel_url   = f"experiments_{modality.lower()}/{fname}"
            exp_bed   = os.path.join(FINAL_DATA_DIR, rna_type,
                                     f"experiments_{modality.lower()}", f"{stem}.bed")
            exp_mods  = get_mod_names(exp_bed)
            short_lbl = stem[:24]

            stanzas.append(
                f"\ttrack {exp_tid}\n"
                f"\tparent {exp_comp_tid}\n"
                f"\tbigDataUrl {rel_url}\n"
                f"\tsubGroups sample={sub_key}\n"
                f"\tshortLabel {short_lbl}\n"
                f"\tlongLabel {rna_type} {modality} — {stem}\n"
                f"\ttype bigBed 9 +\n"
                f"\titemRgb on\n"
                + _EXP_FILTER
                + f"\tfilterValues.name {exp_mods}\n"
                  f"\thtml {exp_html_ref}\n"
                  f"\tpriority {j}\n"
            )

    if rna_type == "polyA-RNA_hg38":
        vcf_html = os.path.join(HUB_DIR, cfg["hub_dir"], "SRS000090_variants.html")
        if not os.path.exists(vcf_html):
            with open(vcf_html, "w") as fh:
                fh.write("""\
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>SRS000090 WGS Variants</title></head>
<body>
<h2>SRS000090 — WGS Merged Variant Calls (two-caller)</h2>
<p>Whole-genome sequencing variant calls for sample SRS000090, merged from two independent variant callers.</p>
</body>
</html>
""")
        stanzas.append(
            "track SRS000090_variants\n"
            "bigDataUrl https://storage.googleapis.com/srs000090-wgs/SRS000090.two_caller_merged.vcf.gz\n"
            "shortLabel SRS000090 variants\n"
            "longLabel SRS000090 WGS merged variant calls (two-caller)\n"
            "type vcfTabix\n"
            "visibility squish\n"
            "maxWindowToDraw 200000\n"
            f"html SRS000090_variants\n"
        )

    with open(path, "w") as fh:
        fh.write("\n".join(stanzas) + "\n")
    print(f"  trackDb → {path}")


def generate_hub_config():
    print("\n── Step 6: Generating UCSC hub v2 config ──")
    os.makedirs(HUB_DIR, exist_ok=True)
    write_hub_txt()
    write_genomes_txt()

    for rna_type in RNA_TYPES:
        cfg = ASSEMBLY_CFG[rna_type]
        hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
        os.makedirs(hub_dir, exist_ok=True)
        write_trackdb(rna_type)

        # Copy 2bit into hub_dir if available
        src_2bit  = cfg.get("2bit_src")
        dest_name = cfg.get("2bit_dest")
        if src_2bit and dest_name and os.path.exists(src_2bit):
            dest = os.path.join(hub_dir, dest_name)
            if not os.path.exists(dest):
                shutil.copy2(src_2bit, dest)
                print(f"  Copied 2bit → {dest}")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    split_combined_files()
    copy_experiment_files()
    convert_split_to_bigbed()
    convert_experiments_to_bigbed()
    generate_html_readmes()
    generate_excel_manifest()
    generate_hub_config()
    print("\nDone.")
    print(f"\nNext steps:")
    print(f"  1. Update PAPER_URL with actual DOI")
    print(f"  2. Update FILE_DESCRIPTIONS with actual text")
    print(f"  3. Push ucsc_hub/ to GitHub and load hub.txt in UCSC")


if __name__ == "__main__":
    main()

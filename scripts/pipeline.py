"""
process_files.py -- data organization and UCSC hub generation

Steps:
  1. Split combined BED files by RNA type → final_data/{RNA_type}/{RNA_type}_{Modality}.bed
  2. Convert split BED files to bigBed (per-platform), plus tiered consensus
     TSVs (concensus_tsvs/) to bigBed (per RNA type, already split)
  3. Generate UCSC hub (ucsc_hub/) with hub.txt, genomes.txt, per-assembly trackDb
  4. Generate per-file HTML readmes (descriptions sourced from RNA_modifications_manifest.tsv)
  5. Generate Excel manifest

Each modality (SRS, LRS, MS) is a single pre-merged consensus BED file —
there are no per-experiment subtracks in this hub. Each RNA type also gets a
tiered cross-platform consensus track (concensus_tsvs/), colored by
modification type with vivid/pale saturation for tier1/tier2 confidence.

RNA type classification:
  chrom contains "rRNA"  →  rRNA
  chrom contains "tRNA"  →  tRNA
  else                   →  polyA_RNA_hg38

Run from repo root:
  python scripts/pipeline.py
"""

import colorsys
import math
import os
import re
import shutil
import subprocess

# ── Source data ────────────────────────────────────────────────────────────────

SOURCE_DIR = "final_bedRmods"

# Modality labels: SRS = short-read sequencing (Illumina), LRS = long-read
# sequencing (ONT). Source filenames keep their original Illumina/ONT names.
COMBINED_FILES = {
    "SRS": os.path.join(SOURCE_DIR, "Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed"),
    "LRS": os.path.join(SOURCE_DIR, "ONT_polyARNA_rRNA_combined.filtered_rmchrY.bed"),
    "MS":  os.path.join(SOURCE_DIR, "MS_rRNA_tRNA.bed"),
}

# Tiered cross-platform consensus tracks: one pre-merged file per RNA type,
# already split (no further splitting needed). Columns: chr/start/end/name/
# tier/strand (+ 2 extra rRNA-only columns that are ignored — see
# process_tiered_for_bigbed).
TIERED_DIR = "concensus_tsvs"

TIERED_FILES = {
    "rRNA":           os.path.join(TIERED_DIR, "tiered_rRNA_only.tsv"),
    "tRNA":           os.path.join(TIERED_DIR, "tiered_tRNA.tsv"),
    "polyA_RNA_hg38": os.path.join(TIERED_DIR, "tiered_polyA.tsv"),
}

# Saturation (%) used by mod_to_rgb for each tier — tier1 is vivid, tier2 is pale.
TIER_SATURATION = {
    "tier1": 100,
    "tier2": 30,
}

MANIFEST_PATH = "RNA_modifications_manifest.tsv"

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

MOD_FULL_NAMES = {
    "m6A":     "N6-methyladenosine",
    "Am":      "2'-O-methyladenosine",
    "Ino":     "Inosine",
    "m1A":     "N1-methyladenosine",
    "m66A":    "N6,N6-dimethyladenosine",
    "m1acp3Y": "1-methyl-3-(3-amino-3-carboxypropyl)pseudouridine",
    "Y":       "Pseudouridine",
    "Ym":      "2'-O-methylpseudouridine",
    "Um":      "2'-O-methyluridine",
    "m3U":     "3-methyluridine",
    "Gm":      "2'-O-methylguanosine",
    "m7G":     "N7-methylguanosine",
    "Cm":      "2'-O-methylcytidine",
    "m5C":     "5-methylcytidine",
    "ac4C":    "N4-acetylcytidine",
}

RNA_TYPES  = ["polyA_RNA_hg38", "tRNA", "rRNA"]
MODALITIES = ["SRS", "LRS", "MS"]

# Spelled-out modality names, used in track HTML page titles/headings.
MODALITY_FULL_NAMES = {
    "SRS": "short-read sequencing",
    "LRS": "long-read sequencing",
    "MS":  "mass spectrometry",
}

# Track shortLabel/longLabel text — fixed across all RNA types (no RNA class in the name).
TRACK_LABELS = {
    "SRS": "Human RNome Project — Short-Read Sequencing Modifications",
    "LRS": "Human RNome Project — Long-Read Sequencing Modifications",
    "MS":  "Human RNome Project — Mass Spectrometry Modifications",
}
TIERED_TRACK_LABEL = "Human RNome Project — Concensus Sequence"

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
    "polyA_RNA_hg38": {
        "ucsc_name":  "hg38",
        "hub_dir":    "hg38_polyA-RNA",
        "2bit_src":   None,          # standard UCSC assembly
        "2bit_dest":  None,
        "sizes":      HG38_SIZES,
        "label":      "poly(A) RNA",
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
uint    chromStart;         "Start position"
uint    chromEnd;           "End position"
string  name;               "Modification short name"
uint    score;              "Score (0-1000)"
char[1] strand;             "Strand (+ or -)"
uint    thickStart;         "Same as chromStart"
uint    thickEnd;           "Same as chromEnd"
uint    itemRgb;            "Display color (R,G,B)"
uint    coverage;           "Coverage"
float   frequency;          "Frequency (%)"
string  single_letter_code; "Single-letter modification code"
lstring mod_id;             "Modification identifier"
)
"""

# ── AutoSQL schema for tiered consensus tracks (bed9 + tier) ──────────────────
# itemRgb needs the full bed9 core (thickStart/thickEnd/itemRgb), not bed6.

TIERED_AUTOSQL_SCHEMA = """\
table tiered_mods
"Tiered cross-platform consensus RNA modification sites"
(
string  chrom;      "Chromosome or RNA name"
uint    chromStart; "Start position"
uint    chromEnd;   "End position"
string  name;       "Modification short name"
uint    score;      "Display score (1000=tier1, 500=tier2)"
char[1] strand;     "Strand (+ or -)"
uint    thickStart; "Same as chromStart"
uint    thickEnd;   "Same as chromEnd"
uint    itemRgb;    "Display color (R,G,B) — hue by modification type, vivid=tier1/pale=tier2"
string  tier;       "Confidence tier: tier1 = strongest cross-platform evidence, tier2 = supported by >=2 datasets"
)
"""

# ── Track descriptions, sourced from the manifest ───────────────────────────────

def load_descriptions_from_manifest(file_map):
    """
    Map dict key → description text, read from MANIFEST_PATH.

    The manifest lists one row per source file (Filename, Description); rows
    are matched back to file_map's values by basename.
    """
    desc_by_filename = {}
    with open(MANIFEST_PATH, encoding="utf-8") as fh:
        next(fh, None)  # header
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fname, _, desc = line.partition("\t")
            desc_by_filename[fname.strip()] = desc.strip()

    descriptions = {}
    for key, src_path in file_map.items():
        fname = os.path.basename(src_path)
        if fname in desc_by_filename:
            descriptions[key] = desc_by_filename[fname]
        else:
            print(f"  [WARN] No manifest description found for {fname}")
    return descriptions


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
    return "polyA_RNA_hg38"


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


def get_field_values(bed_path, col_idx):
    """Collect sorted unique values of column col_idx (0-based) from a BED file."""
    values = set()
    try:
        with open(bed_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) > col_idx:
                    values.add(cols[col_idx].strip())
    except OSError:
        pass
    return ",".join(sorted(values))


def get_field_max(bed_path, col_idx):
    """Return the max value of column col_idx (0-based) across a BED file, or 0 if unavailable."""
    max_val = 0.0
    try:
        with open(bed_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) > col_idx:
                    try:
                        val = float(cols[col_idx])
                    except ValueError:
                        continue
                    if val > max_val:
                        max_val = val
    except OSError:
        pass
    return max_val


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
# Step 2 — Convert split BED → bigBed
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


def process_tiered_for_bigbed(src_path, dst_path, chrom_remap=None, valid_chroms=None):
    """
    Normalise a tiered consensus TSV (header: chr/start/end/name/tier/strand,
    plus optional trailing columns that are ignored) into a 10-col bed9+1 for
    bedToBigBed:
      - score derived from tier (tier1=1000, tier2=500)
      - itemRgb: hue from mod name, saturation from TIER_SATURATION[tier]
        (vivid for tier1, pale for tier2) — same hue table as the other tracks
      - chrom_remap / valid_chroms: same semantics as process_bed_for_bigbed
    """
    rows = []
    skipped = 0
    with open(src_path, encoding="utf-8") as fh:
        next(fh, None)  # header row
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, name, tier, strand = cols[:6]
            if chrom_remap is not None and not chrom.startswith("chr"):
                ucsc = chrom_remap.get(chrom)
                if ucsc is None:
                    skipped += 1
                    continue  # unmapped scaffold — skip
                chrom = ucsc
            if valid_chroms is not None and chrom not in valid_chroms:
                skipped += 1
                continue
            name = name.strip()
            tier = tier.strip()
            score = "1000" if tier == "tier1" else "500"
            strand = strand if strand in ("+", "-") else "+"
            item_rgb = mod_to_rgb(name, TIER_SATURATION.get(tier, 0))
            rows.append([chrom, start, end, name, score, strand, start, end, item_rgb, tier])
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


def write_tiered_autosql(out_dir):
    as_path = os.path.join(out_dir, "tiered_mods.as")
    with open(as_path, "w") as fh:
        fh.write(TIERED_AUTOSQL_SCHEMA)
    return as_path


def run_bedtobigbed(bed_path, sizes_path, as_path, bigbed_path, bed_type="bed9+4"):
    cmd = [
        BEDTOBIGBED,
        f"-type={bed_type}",
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


def _build_sizes_and_remap():
    """Build per-RNA-type chrom-sizes paths and the hg38 Ensembl→UCSC remap, shared
    by both the per-platform and tiered-consensus bigBed conversion steps."""
    rrna_sizes_path = os.path.join(HUB_DIR, "rrna", "hs_rRNA.chrom.sizes")
    if os.path.exists(RRNA_FAI):
        rrna_sizes = fai_to_sizes(RRNA_FAI)
        os.makedirs(os.path.dirname(rrna_sizes_path), exist_ok=True)
        write_sizes_file(rrna_sizes, rrna_sizes_path)
    else:
        print(f"  [WARN] rRNA FAI not found ({RRNA_FAI}); rRNA bigBed skipped")
        rrna_sizes_path = None

    hg38_remap = None
    if os.path.exists(HG38_SIZES):
        hg38_remap = build_ensembl_to_ucsc_map(HG38_SIZES)
        print(f"  Built hg38 chrom remap: {len(hg38_remap)} scaffold entries")

    sizes_for = {
        "rRNA":            rrna_sizes_path,
        "polyA_RNA_hg38":  HG38_SIZES if os.path.exists(HG38_SIZES) else None,
        "tRNA":            TRNA_SIZES,  # None until tRNA reference is provided
    }
    return sizes_for, hg38_remap


def convert_split_to_bigbed():
    print("\n── Step 2: Converting split BED → bigBed ──")

    sizes_for, hg38_remap = _build_sizes_and_remap()

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

            remap = hg38_remap if rna_type == "polyA_RNA_hg38" else None
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


def convert_tiered_to_bigbed():
    print("\n── Step 2b: Converting tiered consensus TSV → bigBed ──")

    sizes_for, hg38_remap = _build_sizes_and_remap()

    for rna_type, src_tsv in TIERED_FILES.items():
        if not os.path.exists(src_tsv):
            print(f"  [SKIP] not found: {src_tsv}")
            continue

        cfg        = ASSEMBLY_CFG[rna_type]
        hub_dir    = os.path.join(HUB_DIR, cfg["hub_dir"])
        sizes_path = sizes_for[rna_type]
        if sizes_path is None:
            print(f"  [TODO] No chrom sizes for {rna_type} — tiered bigBed skipped")
            continue
        os.makedirs(hub_dir, exist_ok=True)
        as_path = write_tiered_autosql(hub_dir)

        out_dir   = os.path.join(FINAL_DATA_DIR, rna_type)
        os.makedirs(out_dir, exist_ok=True)
        dst_bed   = os.path.join(out_dir, f"{rna_type}_consensus_tiered.bed")
        bigbed_out = os.path.join(hub_dir, f"{rna_type}_consensus_tiered.bigBed")

        remap = hg38_remap if rna_type == "polyA_RNA_hg38" else None
        valid  = set(fai_to_sizes(sizes_path).keys()) if sizes_path else None
        n, skipped = process_tiered_for_bigbed(src_tsv, dst_bed,
                                                chrom_remap=remap,
                                                valid_chroms=valid)
        skip_str = f" ({skipped} skipped)" if skipped else ""
        print(f"  Processed {n:>6,} rows{skip_str}: {os.path.basename(src_tsv)}")

        if run_bedtobigbed(dst_bed, sizes_path, as_path, bigbed_out, bed_type="bed9+1"):
            print(f"  bigBed → {bigbed_out}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — HTML readmes
# ══════════════════════════════════════════════════════════════════════════════

# Filter controls shown on each track type's details page.
PLATFORM_FILTER_ROWS = [
    ("Modification type", "show only selected modification types"),
    ("Detection score", "filter by score, 0–1000"),
    ("Read coverage", "filter by number of reads/spectra covering the site"),
    ("Modification frequency (%)", "filter by stoichiometry, 0–100%"),
]

TIERED_FILTER_ROWS = [
    ("Modification type", "show only selected modification types"),
    (
        "Confidence tier",
        "show only tier1 (strongest cross-platform evidence) and/or tier2 "
        "(supported by ≥2 datasets) sites",
    ),
]

PLATFORM_SATURATION_NOTE = (
    "within each modification type, color saturation is proportional to "
    "modification frequency — vivid colors mark high-frequency "
    "(high-stoichiometry) sites, pale colors mark low-frequency sites."
)

TIERED_SATURATION_NOTE = (
    "within each modification type, color saturation encodes confidence tier "
    "— tier1 sites (strongest cross-platform evidence) are vivid, tier2 "
    "sites (supported by ≥2 datasets) are pale."
)


def build_color_key_html():
    """Render the modification-type color key (hue table) as an HTML legend."""
    seen_hex = set()
    rows = []
    for name, hexcode in COLOR_CODE.items():
        if hexcode in seen_hex:
            continue  # skip duplicate aliases (e.g. "I" alongside "Ino")
        seen_hex.add(hexcode)
        display_name = "Ino / I" if name == "Ino" else name
        full_name = MOD_FULL_NAMES.get(name, "")
        rows.append(
            f"<tr><td>{display_name}</td><td>{full_name}</td>"
            f'<td><span style="display:inline-block;width:14px;height:14px;'
            f'background-color:{hexcode};border:1px solid #999;"></span></td>'
            f"<td>{hexcode}</td></tr>"
        )
    return (
        "<table border=\"1\" cellpadding=\"4\" cellspacing=\"0\">\n"
        "<tr><th>Modification</th><th>Full name</th><th>Color</th><th>Hex</th></tr>\n"
        + "\n".join(rows) + "\n"
        "</table>"
    )


def build_display_conventions_html(saturation_note, filter_rows):
    """Build a 'Display Conventions and Configuration' section: how items are
    colored, plus the filter controls available on the track details page."""
    filter_items = "".join(
        f"<li><strong>{label}</strong> — {desc}</li>\n" for label, desc in filter_rows
    )
    return (
        "<h3>Display Conventions and Configuration</h3>\n"
        "<p>\n"
        "  Items are colored by modification type; "
        f"{saturation_note}\n"
        "  Modification types not shown in the color key below are rendered in gray (#808080).\n"
        "</p>\n"
        f"{build_color_key_html()}\n"
        "<p>\n"
        "  This track can be filtered using the controls on the track description page "
        "(click the track name in the Genome Browser, then scroll down):\n"
        "</p>\n"
        f"<ul>\n{filter_items}</ul>"
    )


def write_html_readme(filename, description, out_dir, title=None, extra_html=""):
    stem      = os.path.splitext(filename)[0]
    html_path = os.path.join(out_dir, f"{stem}.html")
    title     = title or stem
    content = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{title}</title>
</head>
<body>
<h2>{title}</h2>
<p>{description}</p>
{extra_html}
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


def generate_html_readmes(descriptions, tiered_descriptions):
    print("\n── Step 3: Generating HTML readmes ──")
    for rna_type in RNA_TYPES:
        cfg     = ASSEMBLY_CFG[rna_type]
        label   = cfg["label"]
        hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
        os.makedirs(hub_dir, exist_ok=True)
        codon_table = build_codon_table_html() if rna_type == "tRNA" else ""

        for modality in MODALITIES:
            filename = f"{rna_type}_{modality}.bed"
            bed_path = os.path.join(FINAL_DATA_DIR, rna_type, filename)
            if not os.path.exists(bed_path):
                continue
            desc = descriptions.get(
                modality,
                f"RNA modification sites — {label}, {modality}. (PLACEHOLDER description)"
            )
            extra_html = (
                codon_table + "\n"
                + build_display_conventions_html(PLATFORM_SATURATION_NOTE, PLATFORM_FILTER_ROWS)
            )
            html_path = write_html_readme(
                filename, desc, hub_dir,
                title=f"{label} {MODALITY_FULL_NAMES[modality]}", extra_html=extra_html
            )
            print(f"  HTML → {html_path}")

        tiered_bed = os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_consensus_tiered.bed")
        if os.path.exists(tiered_bed):
            desc = tiered_descriptions.get(
                rna_type,
                f"Tiered consensus RNA modification sites — {label}. (PLACEHOLDER description)"
            )
            extra_html = (
                codon_table + "\n"
                + build_display_conventions_html(TIERED_SATURATION_NOTE, TIERED_FILTER_ROWS)
            )
            html_path = write_html_readme(
                f"{rna_type}_consensus_tiered.bed", desc, hub_dir,
                title=f"{label} tiered consensus", extra_html=extra_html
            )
            print(f"  HTML → {html_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — Excel manifest
# ══════════════════════════════════════════════════════════════════════════════

def _manifest_rows(descriptions, tiered_descriptions):
    """Yield (filename, rna_type, modality, description, hub_url) for every consensus track."""
    for rna_type in RNA_TYPES:
        cfg = ASSEMBLY_CFG[rna_type]
        for modality in MODALITIES:
            filename = f"{rna_type}_{modality}.bed"
            if not os.path.exists(os.path.join(FINAL_DATA_DIR, rna_type, filename)):
                continue
            hub_url = f"{HUB_BASE_URL}/{HUB_DIR}/{cfg['hub_dir']}/{rna_type}_{modality}.bigBed"
            desc    = descriptions.get(modality, f"RNA modification sites — {cfg['label']}, {modality}.")
            yield filename, rna_type, modality, desc, hub_url

        tiered_filename = f"{rna_type}_consensus_tiered.bed"
        if os.path.exists(os.path.join(FINAL_DATA_DIR, rna_type, tiered_filename)):
            hub_url = f"{HUB_BASE_URL}/{HUB_DIR}/{cfg['hub_dir']}/{rna_type}_consensus_tiered.bigBed"
            desc    = tiered_descriptions.get(rna_type, f"Tiered consensus RNA modification sites — {cfg['label']}.")
            yield tiered_filename, rna_type, "Tiered", desc, hub_url


def generate_excel_manifest(descriptions, tiered_descriptions):
    print("\n── Step 4: Generating Excel manifest ──")
    try:
        from openpyxl import Workbook
    except ImportError:
        print("  [SKIP] openpyxl not installed — run: pip install openpyxl")
        _write_tsv_manifest(descriptions, tiered_descriptions)
        return

    wb = Workbook()
    ws = wb.active
    ws.title = "Manifest"
    ws.append(["Filename", "RNA Type", "Modality", "Description", "Hub BigBed URL", "Paper"])
    for row in _manifest_rows(descriptions, tiered_descriptions):
        ws.append([*row, PAPER_URL])

    out_path = os.path.join(FINAL_DATA_DIR, "RNA_modifications_manifest.xlsx")
    wb.save(out_path)
    print(f"  Manifest → {out_path}")


def _write_tsv_manifest(descriptions, tiered_descriptions):
    """Fallback TSV manifest when openpyxl is unavailable."""
    out_path = os.path.join(FINAL_DATA_DIR, "RNA_modifications_manifest.tsv")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Filename\tRNA Type\tModality\tDescription\tHub BigBed URL\tPaper\n")
        for row in _manifest_rows(descriptions, tiered_descriptions):
            fh.write("\t".join(row) + f"\t{PAPER_URL}\n")
    print(f"  Manifest (TSV) → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 — UCSC hub config
# ══════════════════════════════════════════════════════════════════════════════


def write_hub_txt():
    desc_html = os.path.join(HUB_DIR, "hubDescription.html")
    if not os.path.exists(desc_html):
        with open(desc_html, "w") as fh:
            fh.write("""\
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>Human RNome Project</title></head>
<body>
<h2>Human RNome RNA Project </h2>
<p>
  A comprehensive map of RNA modifications across the human transcriptome,
  integrating short-read sequencing (SRS), long-read sequencing (LRS), and mass spectrometry (MS) data
  for poly(A) RNA, ribosomal RNA (rRNA), and transfer RNA (tRNA).
</p>
<p><strong>Reference:</strong> PLACEHOLDER_DOI</p>
</body>
</html>
""")
    content = f"""\
hub RNome
shortLabel Human RNome Project
longLabel Human RNome Project
genomesFile genomes.txt
email baihe_sun@mail.dfci.harvard.edu
descriptionUrl hubDescription.html
"""
    path = os.path.join(HUB_DIR, "hub.txt")
    with open(path, "w") as fh:
        fh.write(content)
    print(f"  hub.txt → {path}")


# Standard genetic code (mRNA codons, 5'→3') — used for the tRNA documentation's codon table.
STANDARD_GENETIC_CODE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Stop": "*",
}

AA_3_TO_FULL = {
    "Ala": "Alanine", "Arg": "Arginine", "Asn": "Asparagine", "Asp": "Aspartic acid",
    "Cys": "Cysteine", "Gln": "Glutamine", "Glu": "Glutamic acid", "Gly": "Glycine",
    "His": "Histidine", "Ile": "Isoleucine", "Leu": "Leucine", "Lys": "Lysine",
    "Met": "Methionine", "Phe": "Phenylalanine", "Pro": "Proline", "Ser": "Serine",
    "Thr": "Threonine", "Trp": "Tryptophan", "Tyr": "Tyrosine", "Val": "Valine",
    "Stop": "Stop",
}


def build_codon_table_html():
    """Render the standard genetic code as an HTML table, grouped by amino acid."""
    codons_by_aa = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        codons_by_aa.setdefault(aa, []).append(codon)

    rows = []
    for aa in sorted(codons_by_aa, key=lambda a: (a == "Stop", a)):
        codons = ", ".join(sorted(codons_by_aa[aa]))
        rows.append(
            f"<tr><td>{AA_3_TO_FULL[aa]}</td><td>{aa}</td><td>{AA_3_TO_1[aa]}</td><td>{codons}</td></tr>"
        )

    return (
        "<h3>Standard genetic code</h3>\n"
        "<table border=\"1\" cellpadding=\"4\" cellspacing=\"0\">\n"
        "<tr><th>Amino acid</th><th>3 letter abbreviation</th>"
        "<th>1 letter abbreviation</th><th>Codons (mRNA, 5'→3')</th></tr>\n"
        + "\n".join(rows) + "\n"
        "</table>"
    )


def _write_assembly_html(hub_dir, label, extra_html=""):
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
  Custom reference sequences for human {label} used in the Human RNome Project.
</p>
{extra_html}
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

        if rna_type == "polyA_RNA_hg38":
            # Standard UCSC assembly — no custom 2bit needed
            block = (
                f"genome {ucsc_name}\n"
                f"trackDb {hub_dir}/trackDb.txt\n"
                f"defaultPos chr21:8400001-8500000\n"
                f"organism \"Homo sapiens\"\n"
                f"description \"Human RNome Project\"\n"
            )
        elif rna_type == "tRNA":
            _write_assembly_html(hub_dir, label, extra_html=build_codon_table_html())
            block = (
                f"genome {ucsc_name}\n"
                f"trackDb {hub_dir}/trackDb.txt\n"
                f"twoBitPath {hub_dir}/hs_tRNA.2bit\n"
                f"defaultPos hs_tRNAAla_AGC:1-76\n"
                f"htmlPath {hub_dir}/assemblyDescription.html\n"
                f"organism \"{organism}\"\n"
                f"description \"Human RNome Project\"\n"
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
                f"description \"Human RNome Project\"\n"
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


def _tiered_bed_exists(rna_type):
    return os.path.exists(
        os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_consensus_tiered.bed")
    )


def _tiered_bigbed_exists(rna_type, cfg):
    return os.path.exists(
        os.path.join(HUB_DIR, cfg["hub_dir"], f"{rna_type}_consensus_tiered.bigBed")
    )


# Continuous range filters (score, coverage, frequency) for consensus tracks
def build_filter_block(bed_path):
    """Build trackDb filter directives for score, coverage, and frequency (all continuous sliders)."""
    cov_max = max(1, int(math.ceil(get_field_max(bed_path, 9))))
    lines = [
        "scoreMin 0",
        "scoreMax 1000",
        "filter.score 0",
        "filterByRange.score on",
        "filterLimits.score 0:1000",
        "filterLabel.score Detection score",
        "filter.coverage 0",
        "filterByRange.coverage on",
        f"filterLimits.coverage 0:{cov_max}",
        "filterLabel.coverage Read coverage",
        "filter.frequency 0",
        "filterByRange.frequency on",
        "filterLimits.frequency 0:100",
        "filterLabel.frequency Modification frequency (%)",
        "filterLabel.name Modification type",
    ]
    return "".join(f"{ln}\n" for ln in lines)


def write_trackdb(rna_type):
    cfg     = ASSEMBLY_CFG[rna_type]
    hub_dir = os.path.join(HUB_DIR, cfg["hub_dir"])
    os.makedirs(hub_dir, exist_ok=True)
    path    = os.path.join(hub_dir, "trackDb.txt")

    present     = [m for m in MODALITIES if _bed_exists(rna_type, m)]
    has_tiered  = _tiered_bed_exists(rna_type)
    if not present and not has_tiered:
        with open(path, "w") as fh:
            fh.write(f"# No data available for {rna_type}\n")
        print(f"  trackDb (empty) → {path}")
        return

    safe_rna = safe_id(rna_type)
    stanzas  = []

    # ── Tiered cross-platform consensus track (vivid=tier1, pale=tier2) ─────────
    # Listed first / highest priority so it appears at the top of the track list.
    if has_tiered:
        tiered_tid  = f"{safe_rna}_consensus_tiered"
        tiered_bed  = os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_consensus_tiered.bed")
        tiered_html = f"{rna_type}_consensus_tiered"
        mod_names   = get_mod_names(tiered_bed)
        tier_values = get_field_values(tiered_bed, 9)

        if _tiered_bigbed_exists(rna_type, cfg):
            stanzas.append(
                f"track {tiered_tid}\n"
                f"bigDataUrl {rna_type}_consensus_tiered.bigBed\n"
                f"shortLabel {TIERED_TRACK_LABEL}\n"
                f"longLabel {TIERED_TRACK_LABEL}\n"
                f"type bigBed 9 +\n"
                f"itemRgb on\n"
                f"visibility pack\n"
                f"filterValues.name {mod_names}\n"
                f"filterLabel.name Modification type\n"
                f"filterValues.tier {tier_values}\n"
                f"filterLabel.tier Confidence tier\n"
                f"html {tiered_html}\n"
                f"priority 1\n"
            )
        else:
            stanzas.append(f"# TODO: {rna_type}_consensus_tiered.bigBed not yet generated\n")

    for priority, modality in enumerate(present, 2 if has_tiered else 1):
        consensus_tid = f"{safe_rna}_{modality}_consensus"
        bigbed_fname  = f"{rna_type}_{modality}.bigBed"
        html_ref      = f"{rna_type}_{modality}"
        bed_path      = os.path.join(FINAL_DATA_DIR, rna_type, f"{rna_type}_{modality}.bed")
        mod_names     = get_mod_names(bed_path)

        # ── Consensus track (flat, no sub-tracks) ───────────────────────────────
        if _bigbed_exists(rna_type, modality, cfg):
            stanzas.append(
                f"track {consensus_tid}\n"
                f"bigDataUrl {bigbed_fname}\n"
                f"shortLabel {TRACK_LABELS[modality]}\n"
                f"longLabel {TRACK_LABELS[modality]}\n"
                f"type bigBed 9 +\n"
                f"itemRgb on\n"
                f"visibility pack\n"
                + build_filter_block(bed_path)
                + f"filterValues.name {mod_names}\n"
                  f"html {html_ref}\n"
                  f"priority {priority}\n"
            )
        else:
            stanzas.append(f"# TODO: {bigbed_fname} not yet generated\n")

    if rna_type == "polyA_RNA_hg38":
        vcf_html = os.path.join(HUB_DIR, cfg["hub_dir"], "NA12878_variants.html")
        if not os.path.exists(vcf_html):
            with open(vcf_html, "w") as fh:
                fh.write("""\
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>NA12878 WGS Variants</title></head>
<body>
<h2>NA12878 — WGS Variants</h2>
<p>Whole-genome sequencing variant calls for sample NA12878, which was utilized throughout the Human RNome Project. Generated by merging calls from two independent variant callers: DeepVariant and GATK HaplotypeCaller.</p>
</body>
</html>
""")
        stanzas.append(
            "track NA12878_variants\n"
            "bigDataUrl https://storage.googleapis.com/srs000090-wgs/SRS000090.two_caller_merged.vcf.gz\n"
            "shortLabel NA12878 variants\n"
            "longLabel NA12878 WGS Variants\n"
            "type vcfTabix\n"
            "visibility squish\n"
            "maxWindowToDraw 200000\n"
            f"html NA12878_variants\n"
        )

    with open(path, "w") as fh:
        fh.write("\n".join(stanzas) + "\n")
    print(f"  trackDb → {path}")


def generate_hub_config():
    print("\n── Step 5: Generating UCSC hub config ──")
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
    descriptions = load_descriptions_from_manifest(COMBINED_FILES)
    tiered_descriptions = load_descriptions_from_manifest(TIERED_FILES)
    split_combined_files()
    convert_split_to_bigbed()
    convert_tiered_to_bigbed()
    generate_html_readmes(descriptions, tiered_descriptions)
    generate_excel_manifest(descriptions, tiered_descriptions)
    generate_hub_config()
    print("\nDone.")
    print(f"\nNext steps:")
    print(f"  1. Update PAPER_URL with actual DOI")
    print(f"  2. Push ucsc_hub/ to GitHub and load hub.txt in UCSC")


if __name__ == "__main__":
    main()

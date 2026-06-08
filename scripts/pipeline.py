"""
pipeline.py

End-to-end pipeline to generate a UCSC assembly hub from rRNA bedRmod files.

Steps:
    1. Fix BED files — remap chrom names, recolor, rescale scores, sort
    2. Write AutoSql (.as) file
    3. Convert to bigBed using UCSC's bedToBigBed
    4. Generate hub config files — hub.txt, genomes.txt, trackDb.txt

Three technology groups are supported, each producing a consensus track and a
composite of three simulated individual experiment tracks:
    - Illumina   (Illumina_combined_orthogonal.bed)
    - Mass Spec  (massspec_consensus_statistics_q0.05_2minsam_100percoverlap.bedrmod)
    - ONT        (ONT_combined_orthogonal.bed)

Input BED format (bedRmod):
    chrom chromStart chromEnd name strand score thickStart thickEnd itemRgb
    coverage frequency [extra cols...]
    Note: strand is at col 4 and score at col 5 (opposite to standard BED).

DATA_MODE:
    "real"  — use data/real/ (gitignored, actual experimental data)
    "demo"  — use data/demo/ (committed, scrambled placeholder data)
              Run with --scramble first to generate demo data from the real FAI.

Usage:
    python scripts/pipeline.py              # run with real data
    python scripts/pipeline.py --scramble   # generate demo data then run in demo mode

Requirements:
    - bedToBigBed and faToTwoBit in supplementary/ or on PATH
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
               https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit

After running:
    1. Convert FASTA to .2bit (if not done already):
           supplementary/faToTwoBit reference/hs_rRNAs_NR_046235.fa ucsc_hub/rrna/hs_rRNAs_NR_046235.2bit
    2. Push to GitHub:
           git add .
           git commit -m "Updated patch labels"
           git push -u origin main
    3. In UCSC: My Data > Track Hubs > My Hubs
           Paste: https://baihesun.github.io/RNome_UCSC_hub/ucsc_hub/hub.txt
"""

import csv
import os
import random
import colorsys
import subprocess
import sys
from collections import defaultdict

# ── Configuration ─────────────────────────────────────────────────────────────
DATA_MODE   = "real"

REAL_DIR    = "data/real"
DEMO_DIR    = "data/demo"
SUPP_DIR    = "supplementary"   # intermediate files + binary tools (gitignored)
REF_DIR   = "reference"     # reference FASTA + FAI (gitignored)
OUTPUT_DIR  = "ucsc_hub"
RRNA_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "rrna")
HG38_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "hg38")
GITHUB_USER = "baihesun"
GITHUB_REPO = "RNome_UCSC_hub"

# Remap short chromosome names in input BEDs to the hs_rRNA assembly names
# used in the .2bit / FAI file.
CHROM_MAP = {
    "18S":   "hs_rRNA_18S",
    "28S":   "hs_rRNA_28S",
    "5S":    "hs_rRNA_5S",
    "5.8S":  "hs_rRNA_5.8S",
}

# One entry per sequencing technology.
# experiment_ids names three simulated experiment subtracks; all three point to
# the same consensus_file.  Replace with real per-experiment file paths when
# individual-experiment BEDs become available.
TECHNOLOGY_GROUPS = {
    "illumina": {
        "consensus_file": "Illumina_combined_orthogonal.bed",
        "short_label":    "Illumina",
        "long_label":     "Illumina",
        "experiment_ids": ["ILL01", "ILL02", "ILL03"],
        # col indices (0-based) for std columns in the input file
        "std_cols": {"score": 12, "coverage": 13, "frequency": 14},
    },
    "massspec": {
        "consensus_file": "massspec_consensus_statistics_q0.05_2minsam_100percoverlap.bedrmod",
        "short_label":    "Mass Spec",
        "long_label":     "Mass Spectrometry",
        "experiment_ids": ["MS01", "MS02", "MS03"],
        # massspec has frequency_std=11, scores_std=12; no coverage_std
        "std_cols": {"score": 12, "coverage": None, "frequency": 11},
    },
    "ont": {
        "consensus_file": "ONT_combined_orthogonal.bed",
        "short_label":    "ONT",
        "long_label":     "Oxford Nanopore",
        "experiment_ids": ["ONT01", "ONT02", "ONT03"],
        "std_cols": {"score": 12, "coverage": 13, "frequency": 14},
    },
}

EXPERIMENTS_SUBDIR = "experiments"
FASTA_FAI        = "hs_rRNAs_NR_046235.fa.fai"
FASTA_FILE       = "hs_rRNAs_NR_046235.fa"
MOD_NAME_MAP_CSV = "data/natural_modifications.csv"

# RNA hg38 tracks — sourced from processed_data/ (output of process_files.py)
HG38_CHROM_SIZES    = os.path.join(SUPP_DIR, "hg38.chrom.sizes")
RNA_HG38_MODALITIES = {
    "Illumina": os.path.join("processed_data", "illumina", "illumina_genome.bed"),
    "MS":       os.path.join("processed_data", "ms",       "ms_genome.bed"),
    "ONT":      os.path.join("processed_data", "ont",      "ont_genome.bed"),
}
RNA_HG38_CONSENSUS  = os.path.join("processed_data", "consensus_genome.bed")
# ─────────────────────────────────────────────────────────────────────────────

COLOR_CODE = {
    "m6A":     "#003f5c",
    "Am":      "#1464a0",
    "Ino":     "#0e8fd4",
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


def load_mod_name_map(csv_path):
    """Read natural_modifications.csv and return {short_name: new_abbrev}."""
    mapping = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            short = row.get("short name", "").strip()
            new   = row.get("new abbrev", "").strip()
            if short and new:
                mapping[short] = new
    return mapping


def mod_to_rgb(mod_type, frequency):
    """
    Convert a mod type + frequency to an RGB string for UCSC.
    Hue = fixed per mod type; saturation scaled by frequency.
    """
    hex_color = COLOR_CODE.get(mod_type.strip(), DEFAULT_COLOR).lstrip("#")
    r, g, b = int(hex_color[0:2], 16) / 255, int(hex_color[2:4], 16) / 255, int(hex_color[4:6], 16) / 255
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s_scaled = s * max(0.0, min(100.0, float(frequency))) / 100.0
    r2, g2, b2 = colorsys.hsv_to_rgb(h, s_scaled, v)
    return f"{int(r2*255)},{int(g2*255)},{int(b2*255)}"


def make_consensus_bed(fixed_bed_paths, output_path):
    """
    Merge processed BED files into a single consensus track.

    Kept as a utility for building consensus from multiple experiment-level BEDs.
    Not called in the main pipeline (consensus files are pre-computed).
    Color: yellow (1/N) → green (N/N) based on observation count.
    """
    total    = len(fixed_bed_paths)
    pos_data = defaultdict(list)

    for bed_path in fixed_bed_paths:
        with open(bed_path) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 9:
                    continue
                chrom    = cols[0]
                start    = cols[1]
                mod_type = cols[3].strip()
                try:
                    freq = float(cols[10]) if len(cols) >= 11 else 0.0
                except ValueError:
                    freq = 0.0
                pos_data[(chrom, start, mod_type)].append((freq, cols))

    rows = []
    for (chrom, start, mod_type), entries in pos_data.items():
        count = len(entries)
        _, best_cols = max(entries, key=lambda x: x[0])
        out_cols     = best_cols[:]
        t            = (count - 1) / (total - 1) if total > 1 else 1.0
        out_cols[8]  = f"{int(255*(1-t))},{int(220+(180-220)*t)},0"
        rows.append(out_cols)

    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols[:11]) + "\n")

    tally = defaultdict(int)
    for entries in pos_data.values():
        tally[len(entries)] += 1
    print(f"    {len(rows)} consensus sites written -> {output_path}")
    for k in sorted(tally):
        print(f"    {tally[k]} sites observed in {k}/{total} samples")


def _get_std_val(cols, idx):
    """Return the value at idx rounded to 2 d.p., or '0' if absent/invalid."""
    if idx is None or idx >= len(cols):
        return "0"
    try:
        return f"{float(cols[idx]):.2f}"
    except (ValueError, TypeError):
        return "0"


def process_bed(input_path, output_path, mod_name_map=None, std_cols=None):
    """
    Read a bedRmod file and write a standardised 14-column BED for bedToBigBed.

    Input column order (bedRmod): chrom chromStart chromEnd name strand score
        thickStart thickEnd itemRgb coverage frequency [extra...]

    Transformations applied:
      - col 0 : remap chrom name via CHROM_MAP
      - col 3 : translate modification name via mod_name_map
      - col 4 : overwritten with rescaled score (0-1000); BED standard position
      - col 5 : set to "+" (rRNA is always forward-strand)
      - col 6/7: thickStart / thickEnd cast to int (may arrive as floats)
      - col 8 : recolored using mod type (hue) and frequency (saturation)
      - col 9 : coverage cast to int (mass spec values can be large floats)
      - cols 11-13: score_std, coverage_std, frequency_std from std_cols indices
      - output is exactly 14 columns, sorted by chrom then position
    """
    rows       = []
    raw_scores = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            cols[0] = CHROM_MAP.get(cols[0], cols[0])

            mod_type = cols[3].strip()
            if mod_name_map:
                mod_type = mod_name_map.get(mod_type, mod_type)
                cols[3] = mod_type

            frequency = float(cols[10]) if len(cols) >= 11 else 100.0

            cols[8] = mod_to_rgb(mod_type, frequency)

            # Score is at col 5 in bedRmod format (col 4 = strand)
            try:
                raw_scores.append(float(cols[5]))
            except (ValueError, IndexError):
                raw_scores.append(0.0)

            # Normalise columns that may arrive as floats
            cols[6] = str(int(float(cols[6])))
            cols[7] = str(int(float(cols[7])))
            if len(cols) > 9:
                try:
                    cols[9] = str(int(float(cols[9])))
                except ValueError:
                    cols[9] = "0"
            if len(cols) >= 11:
                try:
                    cols[10] = f"{float(cols[10]):.2f}"
                except ValueError:
                    pass

            # Extract std values before truncating cols
            sc = std_cols or {}
            s_std = _get_std_val(cols, sc.get("score"))
            c_std = _get_std_val(cols, sc.get("coverage"))
            f_std = _get_std_val(cols, sc.get("frequency"))

            while len(cols) < 11:
                cols.append("0")

            # order: bed9 | score_std | coverage | coverage_std | frequency | frequency_std
            rows.append(cols[:9] + [s_std, cols[9], c_std, cols[10], f_std])

    # Rescale score to 0-1000 → col 4 (BED standard score position)
    if rows:
        min_s       = min(raw_scores)
        max_s       = max(raw_scores)
        score_range = max_s - min_s if max_s != min_s else 1.0
        for i, cols in enumerate(rows):
            cols[4] = str(int((raw_scores[i] - min_s) / score_range * 1000))
            cols[5] = "+"

    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols[:14]) + "\n")

    orig_min = min(raw_scores) if raw_scores else 0
    orig_max = max(raw_scores) if raw_scores else 0
    print(f"    {len(rows)} sites written -> {output_path}")
    print(f"    Original score range: {orig_min:.4f} – {orig_max:.4f} → rescaled to 0–1000")


def process_mrna_bed(input_path, output_path, mod_name_map=None, std_cols=None):
    """
    Like process_bed but for standard BED format (score at col 4, strand at col 5)
    and without CHROM_MAP remapping — used for hg38 mRNA data.
    Outputs 14 columns: bed9+2 (coverage, frequency) + score_std, coverage_std, frequency_std.
    """
    rows       = []
    raw_scores = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            mod_type = cols[3].strip()
            if mod_name_map:
                mod_type = mod_name_map.get(mod_type, mod_type)
                cols[3] = mod_type

            frequency = float(cols[10]) if len(cols) >= 11 else 100.0
            cols[8] = mod_to_rgb(mod_type, frequency)

            # Score is at col 4 in standard BED (col 5 = strand)
            try:
                raw_scores.append(float(cols[4]))
            except (ValueError, IndexError):
                raw_scores.append(0.0)

            cols[6] = str(int(float(cols[6])))
            cols[7] = str(int(float(cols[7])))
            if len(cols) > 9:
                try:
                    cols[9] = str(int(float(cols[9])))
                except ValueError:
                    cols[9] = "0"
            if len(cols) >= 11:
                try:
                    cols[10] = f"{float(cols[10]):.2f}"
                except ValueError:
                    pass

            sc = std_cols or {}
            s_std = _get_std_val(cols, sc.get("score"))
            c_std = _get_std_val(cols, sc.get("coverage"))
            f_std = _get_std_val(cols, sc.get("frequency"))

            while len(cols) < 11:
                cols.append("0")

            # order: bed9 | score_std | coverage | coverage_std | frequency | frequency_std
            rows.append(cols[:9] + [s_std, cols[9], c_std, cols[10], f_std])

    if rows:
        min_s       = min(raw_scores)
        max_s       = max(raw_scores)
        score_range = max_s - min_s if max_s != min_s else 1.0
        for i, cols in enumerate(rows):
            cols[4] = str(int((raw_scores[i] - min_s) / score_range * 1000))
            # col 5 (strand) kept as-is — mRNA sites can be +/-

    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols[:14]) + "\n")

    orig_min = min(raw_scores) if raw_scores else 0
    orig_max = max(raw_scores) if raw_scores else 0
    print(f"    {len(rows)} sites written -> {output_path}")
    print(f"    Original score range: {orig_min:.4f} – {orig_max:.4f} → rescaled to 0–1000")


def write_autosql(out_dir):
    """Write the AutoSql (.as) file for the bed9+5 schema (all technologies)."""
    path = os.path.join(out_dir, "rRNA_mods.as")
    with open(path, "w") as f:
        f.write("""\
table rRNA_mods
"Human rRNA modification sites"
(
string  chrom;         "Chromosome"
uint    chromStart;    "Start position (0-based)"
uint    chromEnd;      "End position"
string  name;          "Modification type"
uint    score;         "Detection score (0-1000)"
char[1] strand;        "Strand"
uint    thickStart;    "Same as chromStart"
uint    thickEnd;      "Same as chromEnd"
uint    itemRgb;       "Item color"
float   score_std;     "Score std"
uint    coverage;      "Coverage (reads)"
float   coverage_std;  "Coverage std"
float   frequency;     "Frequency (%)"
float   frequency_std; "Frequency std (%)"
)
""")
    print(f"    Written {path}")
    return path


def write_rna_hg38_autosql(out_dir):
    """Write the AutoSql (.as) file for hg38 mRNA tracks (bed9+5 schema)."""
    path = os.path.join(out_dir, "rna_hg38_mods.as")
    with open(path, "w") as f:
        f.write("""\
table rna_hg38_mods
"hg38 RNA modification sites"
(
string  chrom;              "Chromosome"
uint    chromStart;         "Start position (0-based)"
uint    chromEnd;           "End position"
string  name;               "Modification type"
uint    score;              "Detection score (0-1000)"
char[1] strand;             "Strand"
uint    thickStart;         "Same as chromStart"
uint    thickEnd;           "Same as chromEnd"
uint    itemRgb;            "Item color"
uint    coverage;           "Coverage (reads)"
float   frequency;          "Frequency (%)"
string  single_letter_code; "Single letter code"
lstring mod_id;             "Modification ID"
string  modality;           "Modality"
)
""")
    print(f"    Written {path}")
    return path


def process_rna_hg38_bed(input_path, output_path, mod_name_map=None, modality=None):
    """
    Process a hg38 mRNA BED file (output of process_files.py) for bigBed conversion.

    Input columns (standard BED format, 13-14 cols):
        chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb
        coverage frequency single_letter_code mod_id [modality]

    If modality is not None, that string is appended to every row (individual files).
    If modality is None, col 13 of each row is used (consensus file).

    Outputs 14 columns: bed9 + coverage, frequency, single_letter_code, mod_id, modality.
    Returns the number of rows written.
    """
    rows = []
    raw_scores = []

    with open(input_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 11:
                continue

            original_name = cols[3].strip()
            # use translated name only for color lookup; keep original as item label
            mod_type = mod_name_map.get(original_name, original_name) if mod_name_map else original_name

            try:
                frequency = float(cols[10])
            except (ValueError, IndexError):
                frequency = 100.0

            cols[8] = mod_to_rgb(mod_type, frequency)

            try:
                raw_scores.append(float(cols[4]))
            except (ValueError, IndexError):
                raw_scores.append(0.0)

            cols[6] = str(int(float(cols[6])))
            cols[7] = str(int(float(cols[7])))
            try:
                cols[9] = str(int(float(cols[9])))
            except (ValueError, IndexError):
                cols[9] = "0"
            try:
                cols[10] = f"{float(cols[10]):.2f}"
            except (ValueError, IndexError):
                pass

            slc      = cols[11].strip() if len(cols) > 11 else ""
            mod_id   = cols[12].strip() if len(cols) > 12 else ""
            row_mod  = modality if modality is not None else (cols[13].strip() if len(cols) > 13 else "")

            rows.append(cols[:9] + [cols[9], cols[10], slc, mod_id, row_mod])

    if rows:
        min_s = min(raw_scores)
        max_s = max(raw_scores)
        score_range = max_s - min_s if max_s != min_s else 1.0
        for i, row in enumerate(rows):
            row[4] = str(int((raw_scores[i] - min_s) / score_range * 1000))

    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row[:14]) + "\n")

    print(f"    {len(rows)} sites written -> {output_path}")
    return len(rows)


def _find_tool(name):
    supp = os.path.join(SUPP_DIR, name)
    return supp if os.path.exists(supp) else name


def scramble_data(dest_dir, sites_per_chrom=80):
    """
    Generate structurally valid but fully randomised BED + FASTA files in
    dest_dir using chromosome names/lengths from data/real/FASTA_FAI.
    One scrambled file is produced per unique consensus_file across all
    TECHNOLOGY_GROUPS (experiments are copies of the consensus, so no extra
    files are needed).
    Short chrom names (18S, 28S, …) are used to match real input format;
    process_bed will remap them via CHROM_MAP.
    """
    os.makedirs(dest_dir, exist_ok=True)
    fai_src = os.path.join(REAL_DIR, FASTA_FAI)
    if not os.path.exists(fai_src):
        print(f"  ERROR: {fai_src} not found — cannot scramble without the real FAI.")
        sys.exit(1)

    # Parse FAI; strip hs_rRNA_ prefix to get short chrom names matching input BED
    chroms = []
    with open(fai_src) as fh:
        for line in fh:
            parts      = line.split("\t")
            short_name = parts[0].replace("hs_rRNA_", "")
            chroms.append((short_name, int(parts[1])))

    mod_types = list(COLOR_CODE.keys())

    # Write a fake FASTA using short names; copy FAI as-is (hs_rRNA_* names)
    with open(os.path.join(dest_dir, FASTA_FILE), "w") as fa:
        for chrom, length in chroms:
            fa.write(f">{chrom}\n")
            seq = "".join(random.choices("ACGT", k=length))
            for i in range(0, length, 60):
                fa.write(seq[i:i+60] + "\n")
    import shutil
    shutil.copy(fai_src, os.path.join(dest_dir, FASTA_FAI))

    seen_files = set()
    for tech_cfg in TECHNOLOGY_GROUPS.values():
        fname = tech_cfg["consensus_file"]
        if fname in seen_files:
            continue
        seen_files.add(fname)
        rows = []
        for chrom, length in chroms:
            positions = sorted(random.sample(range(length - 1), min(sites_per_chrom, length - 1)))
            for pos in positions:
                mod   = random.choice(mod_types)
                score = round(random.uniform(0, 10), 4)
                freq  = round(random.uniform(5, 100), 2)
                cov   = random.randint(10, 500)
                # bedRmod column order: col 4 = strand, col 5 = score
                rows.append([chrom, str(pos), str(pos + 1), mod, "+", str(score),
                              str(pos), str(pos + 1), "0,0,0", str(cov), str(freq)])
        out_path = os.path.join(dest_dir, fname)
        with open(out_path, "w") as fh:
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print(f"    Scrambled {len(rows)} sites -> {out_path}")

    print(f"  Demo data written to {dest_dir}/")


def run_bedtobigbed(bed_path, fai_path, bigbed_path, as_path):
    """Convert a fixed BED to bigBed using bed9+2 schema."""
    cmd    = _find_tool("bedToBigBed")
    result = subprocess.run(
        [cmd, f"-as={as_path}", "-type=bed9+5", bed_path, fai_path, bigbed_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr.strip()}")
        sys.exit(1)
    print(f"    bigBed written -> {bigbed_path}")


def write_hub_txt(out_dir):
    path = os.path.join(out_dir, "hub.txt")
    with open(path, "w") as f:
        f.write("""\
hub        hs_rRNA_modifications
shortLabel Human rRNA Modifications
longLabel  Human rRNA RNA modifications (18S, 28S, 5.8S, 5S)
genomesFile genomes.txt
email      baihe_sun@dfci.harvard.edu
descriptionUrl hubDescription.html
""")
    print(f"    Written {path}")


def write_genomes_txt(out_dir):
    path = os.path.join(out_dir, "genomes.txt")
    with open(path, "w") as f:
        f.write("""\
genome hs_rRNA
trackDb rrna/trackDb.txt
description Human rRNA sequences
twoBitPath rrna/hs_rRNAs_NR_046235.2bit
organism Homo sapiens
defaultPos hs_rRNA_18S:1-1869
orderKey 1
scientificName Homo sapiens
htmlPath hubDescription.html

genome hg38
trackDb hg38/trackDb_hg38.txt
""")
    print(f"    Written {path}")


def get_mod_types(bed_paths):
    """Return a sorted comma-delimited string of unique modification types (col 3)."""
    mod_types = set()
    for path in bed_paths:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) >= 4:
                    mod_types.add(cols[3].strip())
    return ",".join(sorted(mod_types))


def _tech_trackdb_block(tech_key, tech_cfg, consensus_bb, exp_bigbed_names,
                        mod_type_values, priority, html_prefix=""):
    """Return the trackDb stanza text for one technology group.

    html_prefix: prepended to all html fields so paths resolve from the hub root,
                 e.g. "rrna/" when trackDb lives in ucsc_hub/rrna/.
    """
    short  = tech_cfg["short_label"]
    long   = tech_cfg["long_label"]

    subgroup_defs = " ".join(
        f"{eid.lower()}={eid}" for eid in exp_bigbed_names
    )

    subtrack_blocks = []
    for i, (exp_id, bb_url) in enumerate(exp_bigbed_names.items(), start=1):
        subtrack_blocks.append(f"""\
    track {tech_key}_exp_{exp_id}
    parent {tech_key}_experiments_composite
    bigDataUrl {bb_url}
    subGroups sample={exp_id.lower()}
    shortLabel {exp_id}
    longLabel {short} {exp_id} — detected rRNA modifications
    type bigBed 9 +
    itemRgb on
    scoreMin 0
    scoreMax 1000
    filter.frequency 0
    filterByRange.frequency on
    filterLimits.frequency 0:100
    filterLabel.frequency Modification frequency (%)
    filterValues.name {mod_type_values}
    filterLabel.name Modification type
    html {html_prefix}{tech_key}_experiments_composite
    priority {i}

""")

    return f"""\
track {tech_key}_super
superTrack on show
shortLabel {short}
longLabel {short} — Consensus and Supporting Experiments
priority {priority}

track {tech_key}_consensus
parent {tech_key}_super
bigDataUrl {consensus_bb}
shortLabel {short} consensus
longLabel {long}
type bigBed 9 +
itemRgb on
visibility pack
scoreMin 0
scoreMax 1000
filter.frequency 0
filterByRange.frequency on
filterLimits.frequency 0:100
filterLabel.frequency Modification frequency (%)
filterValues.name {mod_type_values}
filterLabel.name Modification type
html {html_prefix}{tech_key}_consensus
priority 1

track {tech_key}_experiments_composite
compositeTrack on
parent {tech_key}_super
shortLabel {short} experiments
longLabel Individual {short} Experiment Tracks
type bigBed 9 +
itemRgb on
visibility hide
subGroup1 sample Sample {subgroup_defs}
dimensions dimX=sample
filterComposite on
html {html_prefix}{tech_key}_experiments_composite
priority 2

{"".join(subtrack_blocks)}"""


def write_trackdb_txt(out_dir, tech_tracks):
    """
    Write trackDb.txt with one supertrack group per technology.

    tech_tracks: list of dicts with keys:
        tech_key, tech_cfg, consensus_bb, exp_bigbed_names, mod_type_values
    """
    path   = os.path.join(out_dir, "trackDb.txt")
    blocks = []
    for priority, entry in enumerate(tech_tracks, start=1):
        blocks.append(_tech_trackdb_block(
            entry["tech_key"],
            entry["tech_cfg"],
            entry["consensus_bb"],
            entry["exp_bigbed_names"],
            entry["mod_type_values"],
            priority,
            html_prefix="rrna/",
        ))
    with open(path, "w") as f:
        f.write("\n".join(blocks) + "\n")
    print(f"    Written {path}")


def write_trackdb_hg38_txt(out_dir, consensus_bb, modality_bbs, all_mod_types):
    """
    Write trackDb_hg38.txt with:
      - SRS000090 VCF track
      - mRNA supertrack containing:
          - consensus track (all modalities, modality shown in popup)
          - composite track with one subtrack per modality with data

    consensus_bb   : bigBed filename for the consensus track, or None
    modality_bbs   : dict {modality_name: bb_filename} for subtracks with data
    all_mod_types  : comma-delimited string of all modification types
    """
    path = os.path.join(out_dir, "trackDb_hg38.txt")

    # Build subgroup/subtrack blocks only for modalities that have data
    modality_keys = {k.lower(): k for k in modality_bbs}
    subgroup_defs = " ".join(f"{k}={v}" for k, v in sorted(modality_keys.items()))
    modality_filter_vals = ",".join(sorted(modality_bbs.keys()))

    subtrack_blocks = []
    for i, (mod_name, bb_name) in enumerate(sorted(modality_bbs.items()), start=1):
        mod_key = mod_name.lower()
        subtrack_blocks.append(f"""\
    track rna_hg38_{mod_key}
    parent rna_hg38_composite
    bigDataUrl {bb_name}
    subGroups modality={mod_key}
    shortLabel RNA {mod_name}
    longLabel RNA Modifications — {mod_name} (hg38)
    type bigBed 9 +
    itemRgb on
    scoreMin 0
    scoreMax 1000
    filter.score 0
    filterByRange.score on
    filterLimits.score 0:1000
    filterLabel.score Score
    filter.frequency 0
    filterByRange.frequency on
    filterLimits.frequency 0:100
    filterLabel.frequency Frequency (%)
    filterValues.name {all_mod_types}
    filterLabel.name Modification type
    html hg38/rna_hg38_composite
    priority {i}

""")

    consensus_block = ""
    if consensus_bb:
        consensus_block = f"""\
track rna_hg38_consensus
parent rna_hg38_super
bigDataUrl {consensus_bb}
shortLabel RNA Consensus
longLabel RNA Modification Consensus — All Modalities (hg38)
type bigBed 9 +
itemRgb on
visibility pack
scoreMin 0
scoreMax 1000
filter.score 0
filterByRange.score on
filterLimits.score 0:1000
filterLabel.score Score
filter.frequency 0
filterByRange.frequency on
filterLimits.frequency 0:100
filterLabel.frequency Frequency (%)
filterValues.name {all_mod_types}
filterLabel.name Modification type
filterValues.modality {modality_filter_vals}
filterLabel.modality Modality
html hg38/rna_hg38_consensus
priority 1

"""

    composite_block = ""
    if modality_bbs:
        composite_block = f"""\
track rna_hg38_composite
compositeTrack on
parent rna_hg38_super
shortLabel RNA by Modality
longLabel RNA Modifications by Modality (hg38)
type bigBed 9 +
itemRgb on
visibility hide
subGroup1 modality Modality {subgroup_defs}
dimensions dimX=modality
filterComposite on
filter.score 0
filterByRange.score on
filterLimits.score 0:1000
filterLabel.score Score
filter.frequency 0
filterByRange.frequency on
filterLimits.frequency 0:100
filterLabel.frequency Frequency (%)
filterValues.name {all_mod_types}
filterLabel.name Modification type
html rna_hg38_composite
priority 2

{"".join(subtrack_blocks)}"""

    with open(path, "w") as f:
        f.write(f"""\
track SRS000090_variants
bigDataUrl https://storage.googleapis.com/srs000090-wgs/SRS000090.two_caller_merged.vcf.gz
shortLabel SRS000090 variants
longLabel SRS000090 WGS merged variant calls (two-caller)
type vcfTabix
visibility squish
maxWindowToDraw 200000

track rna_hg38_super
superTrack on show
shortLabel RNA Modifications
longLabel RNA Modifications — Multi-modality (hg38)
priority 1

{consensus_block}{composite_block}""")
    print(f"    Written {path}")


def write_rna_hg38_htmls(out_dir):
    """Write HTML description files for the hg38 mRNA consensus and composite tracks."""
    with open(os.path.join(out_dir, "rna_hg38_consensus.html"), "w") as f:
        f.write("""\
<h2>RNA Modification Consensus Track (hg38)</h2>
<p><strong>Assembly:</strong> GRCh38 / hg38</p>
<p>Union of all modality calls. The <em>Modality</em> field indicates which
sequencing technology detected each site.</p>

<h3>Fields</h3>
<table border="1" cellpadding="4">
<tr><th>Field</th><th>Description</th></tr>
<tr><td>name</td><td>Short modification name (item label, e.g. m6A)</td></tr>
<tr><td>score</td><td>Detection confidence rescaled 0–1000</td></tr>
<tr><td>coverage</td><td>Read depth at site</td></tr>
<tr><td>frequency</td><td>Modification stoichiometry (%)</td></tr>
<tr><td>single_letter_code</td><td>Single-letter code</td></tr>
<tr><td>modality</td><td>Source technology (Illumina / MS / ONT)</td></tr>
</table>

<h3>Colors</h3>
<p>Hue encodes modification type; saturation scales with frequency.</p>

<h3>Contact</h3>
<p>MISSING</p>
""")
    print("    Written rna_hg38_consensus.html")

    with open(os.path.join(out_dir, "rna_hg38_composite.html"), "w") as f:
        f.write("""\
<h2>RNA Modifications by Modality (hg38)</h2>
<p><strong>Assembly:</strong> GRCh38 / hg38</p>
<p>Individual modality subtracks (Illumina, MS, ONT). Use the subtrack
checkboxes or the <em>Modality</em> filter to show/hide each technology.</p>

<h3>Fields</h3>
<table border="1" cellpadding="4">
<tr><th>Field</th><th>Description</th></tr>
<tr><td>name</td><td>Short modification name (item label, e.g. m6A)</td></tr>
<tr><td>score</td><td>Detection confidence rescaled 0–1000</td></tr>
<tr><td>coverage</td><td>Read depth at site</td></tr>
<tr><td>frequency</td><td>Modification stoichiometry (%)</td></tr>
<tr><td>single_letter_code</td><td>Single-letter code</td></tr>
<tr><td>modality</td><td>Source technology</td></tr>
</table>

<h3>Colors</h3>
<p>Hue encodes modification type; saturation scales with frequency.</p>

<h3>Contact</h3>
<p>MISSING</p>
""")
    print("    Written rna_hg38_composite.html")


def write_track_html(out_dir):
    """Write consensus and composite HTML description files for each technology."""
    for tech_key, tech_cfg in TECHNOLOGY_GROUPS.items():
        short = tech_cfg["short_label"]
        long  = tech_cfg["long_label"]

        with open(os.path.join(out_dir, f"{tech_key}_consensus.html"), "w") as f:
            f.write(f"""\
<h2>{short} Consensus Track</h2>
<p><strong>Long label:</strong> {long}</p>

<h3>Methods</h3>
<table border="1" cellpadding="4">
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Technology</td><td>{short}</td></tr>
<tr><td>QC applied</td><td>MISSING</td></tr>
</table>

<h3>Contact</h3>
<p>MISSING</p>
""")
        print(f"    Written {tech_key}_consensus.html")

        with open(os.path.join(out_dir, f"{tech_key}_experiments_composite.html"), "w") as f:
            f.write(f"""\
<h2>Individual {short} Experiment Tracks</h2>
<p>Individual experiments contributing to the {short} consensus track.</p>

<h3>Methods</h3>
<table border="1" cellpadding="4">
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Technology</td><td>{short}</td></tr>
<tr><td>QC applied</td><td>MISSING</td></tr>
</table>

<h3>Filters</h3>
<p>Use the filter controls to toggle samples or restrict by modification
frequency (%) or modification type.</p>

<h3>Contact</h3>
<p>MISSING</p>
""")
        print(f"    Written {tech_key}_experiments_composite.html")


def write_description_html(out_dir):
    path = os.path.join(out_dir, "hubDescription.html")
    with open(path, "w") as f:
        f.write("""\
<h1>Human rRNA Modifications</h1>
<p>RNA modification tracks for human rRNA sequences (18S, 28S, 5.8S, 5S)
derived from three orthogonal sequencing technologies.</p>

<h2>Track groups</h2>
<table border="1" cellpadding="4">
<tr><th>Group</th><th>Technology</th><th>Description</th></tr>
<tr><td>Illumina</td><td>Short-read sequencing</td><td>Orthogonal consensus modification calls</td></tr>
<tr><td>Mass Spec</td><td>Mass spectrometry</td><td>Consensus modification calls</td></tr>
<tr><td>ONT</td><td>Oxford Nanopore</td><td>Orthogonal consensus modification calls</td></tr>
</table>

<h2>Track colors</h2>
<p>Colors encode modification type (hue) and stoichiometry (saturation).
Vivid = highly modified; pale = partially modified.</p>

<h2>Score</h2>
<p>Detection confidence rescaled to 0–1000 per file independently.
Darker = higher confidence. Scores are not comparable across tracks.</p>

<h2>Modification type color key</h2>
<table border="1" cellpadding="4">
<tr><th>Modification</th><th>Type</th><th>Color family</th></tr>
<tr><td>Am, Cm, Gm, Um, Ym</td><td>2'-O-methylation</td><td>Blue</td></tr>
<tr><td>m6A, m66A, m1A</td><td>N-methyladenosine</td><td>Red / Orange</td></tr>
<tr><td>Y, m1acp3Y</td><td>Pseudouridine</td><td>Teal</td></tr>
<tr><td>m5C, ac4C</td><td>Cytosine modifications</td><td>Green</td></tr>
<tr><td>m7G</td><td>N7-methylguanosine</td><td>Purple</td></tr>
<tr><td>m3U</td><td>3-methyluridine</td><td>Amber</td></tr>
</table>
""")
    print(f"    Written {path}")


def main():
    if "--scramble" in sys.argv:
        print("\n── Generating scrambled demo data ───────────────────────────────────────")
        scramble_data(DEMO_DIR)
        input_dir = DEMO_DIR
    else:
        input_dir = REAL_DIR if DATA_MODE == "real" else DEMO_DIR

    print(f"\n  Using input data from: {input_dir}  (DATA_MODE={DATA_MODE})")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(RRNA_OUTPUT_DIR, exist_ok=True)
    os.makedirs(HG38_OUTPUT_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)

    global COLOR_CODE
    mod_name_map = load_mod_name_map(MOD_NAME_MAP_CSV)
    COLOR_CODE   = {mod_name_map.get(k, k): v for k, v in COLOR_CODE.items()}
    print(f"\n  Loaded {len(mod_name_map)} name mappings from {MOD_NAME_MAP_CSV}")

    fai_path = os.path.join(REF_DIR, FASTA_FAI)
    suffix   = ".demo" if input_dir == DEMO_DIR else ""

    print("\n── Step 1: Write AutoSql file ───────────────────────────────────────────")
    as_path = write_autosql(RRNA_OUTPUT_DIR)

    exp_dir = os.path.join(RRNA_OUTPUT_DIR, EXPERIMENTS_SUBDIR)
    os.makedirs(exp_dir, exist_ok=True)

    tech_tracks = []

    for tech_key, tech_cfg in TECHNOLOGY_GROUPS.items():
        short = tech_cfg["short_label"]
        print(f"\n══ {short} {'═'*(50-len(short))}")

        # ── Consensus ──────────────────────────────────────────────────────────
        src_file      = os.path.join(input_dir, tech_cfg["consensus_file"])
        con_fixed     = os.path.join(SUPP_DIR, f"{tech_key}_consensus{suffix}.fixed.bed")
        con_bb_name   = f"{tech_key}_consensus{suffix}.bigBed"
        con_bb_path   = os.path.join(RRNA_OUTPUT_DIR, con_bb_name)

        std_cols = tech_cfg.get("std_cols")
        print(f"\n  consensus:")
        process_bed(src_file, con_fixed, mod_name_map, std_cols)
        print(f"\n  → bigBed:")
        run_bedtobigbed(con_fixed, fai_path, con_bb_path, as_path)

        mod_type_values = get_mod_types([con_fixed])
        print(f"    Modification types: {mod_type_values}")

        # ── Experiments (simulated copies of consensus) ────────────────────────
        exp_bigbed_names = {}
        for exp_id in tech_cfg["experiment_ids"]:
            exp_fixed   = os.path.join(SUPP_DIR, f"{tech_key}_{exp_id}{suffix}.fixed.bed")
            exp_bb_name = f"{tech_key}_{exp_id}{suffix}.bigBed"
            exp_bb_path = os.path.join(exp_dir, exp_bb_name)

            print(f"\n  experiment {exp_id}:")
            process_bed(src_file, exp_fixed, mod_name_map, std_cols)
            run_bedtobigbed(exp_fixed, fai_path, exp_bb_path, as_path)
            exp_bigbed_names[exp_id] = f"{EXPERIMENTS_SUBDIR}/{exp_bb_name}"

        tech_tracks.append({
            "tech_key":         tech_key,
            "tech_cfg":         tech_cfg,
            "consensus_bb":     con_bb_name,
            "exp_bigbed_names": exp_bigbed_names,
            "mod_type_values":  mod_type_values,
        })

    print("\n── Step 2: Write hub config files ──────────────────────────────────────")
    write_hub_txt(OUTPUT_DIR)
    write_genomes_txt(OUTPUT_DIR)
    write_trackdb_txt(RRNA_OUTPUT_DIR, tech_tracks)
    write_track_html(RRNA_OUTPUT_DIR)
    write_description_html(OUTPUT_DIR)

    print("\n── Step 3: RNA hg38 tracks ──────────────────────────────────────────────")
    if not os.path.exists(HG38_CHROM_SIZES):
        print(f"    Skipping: {HG38_CHROM_SIZES} not found — download with:")
        print(f"    curl -fsSL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -o {HG38_CHROM_SIZES}")
    else:
        rna_as_path = write_rna_hg38_autosql(HG38_OUTPUT_DIR)

        modality_bbs = {}
        all_mod_fixed = []

        for mod_name, src_path in RNA_HG38_MODALITIES.items():
            if not os.path.exists(src_path):
                print(f"    Skipping {mod_name}: {src_path} not found")
                continue
            mod_key  = mod_name.lower()
            fixed    = os.path.join(SUPP_DIR, f"rna_hg38_{mod_key}.fixed.bed")
            bb_name  = f"rna_hg38_{mod_key}.bigBed"
            bb_path  = os.path.join(HG38_OUTPUT_DIR, bb_name)

            print(f"\n  RNA hg38 {mod_name}:")
            count = process_rna_hg38_bed(src_path, fixed, mod_name_map, modality=mod_name)
            if count == 0:
                print(f"    No data — skipping bigBed for {mod_name}")
                continue
            print(f"\n  → bigBed:")
            run_bedtobigbed(fixed, HG38_CHROM_SIZES, bb_path, rna_as_path)
            modality_bbs[mod_name] = bb_name
            all_mod_fixed.append(fixed)

        consensus_bb_name = None
        if os.path.exists(RNA_HG38_CONSENSUS):
            con_fixed         = os.path.join(SUPP_DIR, "rna_hg38_consensus.fixed.bed")
            consensus_bb_name = "rna_hg38_consensus.bigBed"
            con_bb_path       = os.path.join(HG38_OUTPUT_DIR, consensus_bb_name)

            print(f"\n  RNA hg38 consensus:")
            count = process_rna_hg38_bed(RNA_HG38_CONSENSUS, con_fixed, mod_name_map, modality=None)
            if count > 0:
                print(f"\n  → bigBed:")
                run_bedtobigbed(con_fixed, HG38_CHROM_SIZES, con_bb_path, rna_as_path)
                all_mod_fixed.append(con_fixed)
            else:
                print("    No data in consensus — skipping")
                consensus_bb_name = None
        else:
            print(f"    Skipping consensus: {RNA_HG38_CONSENSUS} not found")

        if consensus_bb_name or modality_bbs:
            all_mod_types = get_mod_types(all_mod_fixed) if all_mod_fixed else ""
            print(f"    Modification types: {all_mod_types}")
            write_trackdb_hg38_txt(HG38_OUTPUT_DIR, consensus_bb_name, modality_bbs, all_mod_types)
            write_rna_hg38_htmls(HG38_OUTPUT_DIR)

    fasta_src = os.path.join(input_dir, FASTA_FILE)
    print(f"""
── Done ─────────────────────────────────────────────────────────────────────

Next steps:

  1. Convert FASTA to .2bit (if not done already):
         {SUPP_DIR}/faToTwoBit {fasta_src} {RRNA_OUTPUT_DIR}/hs_rRNAs_NR_046235.2bit

  2. Push to GitHub:
         git add .
         git commit -m "Update hub files"
         git push -u origin main

  3. In UCSC: My Data > Track Hubs > My Hubs
         Paste: https://{GITHUB_USER}.github.io/{GITHUB_REPO}/{OUTPUT_DIR}/hub.txt
""")


if __name__ == "__main__":
    main()

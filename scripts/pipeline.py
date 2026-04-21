"""
pipeline.py

End-to-end pipeline to generate a UCSC assembly hub from rRNA bedRmod files.

Steps:
    1. Recolor BED files — fixed hex color per mod type, saturation encodes frequency
    2. Fix BED files — keep all 11 columns, rescale scores to 0-1000, sort by chrom/position
    3. Write AutoSql (.as) file describing all 11 columns
    4. Convert to bigBed format using UCSC's bedToBigBed with -type=bed9+2 -as=rRNA_mods.as
    5. Generate hub config files — hub.txt, genomes.txt, trackDb.txt

DATA_MODE:
    "real"  — use data/real/ (gitignored, contains actual experimental data)
    "demo"  — use data/demo/ (committed, contains scrambled placeholder data)
              Run with --scramble first to generate demo data from the real FAI.

Usage:
    python scripts/pipeline.py                  # run with current DATA_MODE
    python scripts/pipeline.py --scramble       # generate demo data then run in demo mode

Requirements:
    - bedToBigBed and faToTwoBit in supplementary/ or on PATH
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
               https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit

After running:
    1. Convert FASTA to .2bit:
           supplementary/faToTwoBit data/real/hs_rRNAs_NR_046235.fa ucsc_hub/hs_rRNAs_NR_046235.2bit
       (use data/demo/hs_rRNAs_NR_046235.fa when in demo mode)
    2. Push to GitHub:
           git add .
           git commit -m "Update hub files"
           git push -u origin main
    3. In UCSC: My Data > Track Hubs > My Hubs
           Paste: https://baihesun.github.io/rRNA_dataviz/ucsc_hub/hub.txt
"""

import os
import random
import colorsys
import subprocess
import sys
from collections import defaultdict

# ── Configuration ─────────────────────────────────────────────────────────────
# Set DATA_MODE to "real" to use actual data (gitignored),
# or "demo" to use scrambled placeholder data (committed to GitHub).
DATA_MODE   = "real"

REAL_DIR    = "data/real"
DEMO_DIR    = "data/demo"
SUPP_DIR    = "supplementary"   # intermediate files + binary tools (gitignored)
OUTPUT_DIR  = "ucsc_hub"
GITHUB_USER = "baihesun"
GITHUB_REPO = "rRNA_dataviz"

INPUT_FILES = {
    "ref":      "H.sapiens_ref_bedRmod_All.bed",
    "sample":   "rRNA_mature_bedRmod_Detection_EM_mean2_Log10_1.8_modOnly_Sample_MRI01.bed",
    "filtered": "rRNA_mature_Filtered_MOD_10_MULT_1000_bedRmod_0.99_Allmods_no_m5C.bed",
}
FASTA_FAI  = "hs_rRNAs_NR_046235.fa.fai"
FASTA_FILE = "hs_rRNAs_NR_046235.fa"

# Keys from INPUT_FILES to merge into the consensus track
CONSENSUS_KEYS = ["sample", "filtered"]
# ─────────────────────────────────────────────────────────────────────────────

# Fixed color per modification type (hex)
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
DEFAULT_COLOR = "#808080"  # gray fallback for unknown mod types


def mod_to_rgb(mod_type, frequency):
    """
    Convert a mod type + frequency to an RGB string for UCSC.
    Base color = COLOR_CODE hex (hue + full saturation).
    Saturation is scaled by frequency (0% → gray, 100% → full color).
    """
    hex_color = COLOR_CODE.get(mod_type.strip(), DEFAULT_COLOR).lstrip("#")
    r, g, b = int(hex_color[0:2], 16) / 255, int(hex_color[2:4], 16) / 255, int(hex_color[4:6], 16) / 255
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s_scaled = s * max(0.0, min(100.0, float(frequency))) / 100.0
    r2, g2, b2 = colorsys.hsv_to_rgb(h, s_scaled, v)
    return f"{int(r2*255)},{int(g2*255)},{int(b2*255)}"


def consensus_rgb(count, total):
    """
    Color for the consensus track based on observation count:
      1/N → yellow (255, 220, 0)
      N/N → green  (0, 180, 0)
    Intermediate counts interpolate between the two.
    """
    t = (count - 1) / (total - 1) if total > 1 else 1.0
    r = int(255 * (1 - t))          # 255 → 0
    g = int(220 + (180 - 220) * t)  # 220 → 180
    b = 0
    return f"{r},{g},{b}"


def make_consensus_bed(fixed_bed_paths, output_path):
    """
    Merge processed BED files into a single consensus track.

    For each unique (chrom, start, mod_type) position:
      - Count in how many of the N input files it appears.
      - Color using a sequential red scale (light = rare, dark = universal).
      - Retain the row with the highest frequency; average scores.
    """
    total = len(fixed_bed_paths)
    # (chrom, start, mod_type) -> list of (frequency, cols)
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
        # Use the row with the highest observed frequency
        _, best_cols = max(entries, key=lambda x: x[0])
        out_cols = best_cols[:]
        out_cols[8] = consensus_rgb(count, total)
        rows.append(out_cols)

    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols[:11]) + "\n")

    # Summary
    tally = defaultdict(int)
    for entries in pos_data.values():
        tally[len(entries)] += 1
    print(f"    {len(rows)} consensus sites written -> {output_path}")
    for k in sorted(tally):
        print(f"    {tally[k]} sites observed in {k}/{total} samples")


def process_bed(input_path, output_path):
    """
    Read a bedRmod file and write a cleaned version:
      - Recolor column 9 (itemRgb): hue = mod type, saturation = frequency
      - Rescale score (col 5) to 0-1000 integer range (required by bedToBigBed)
      - Keep all 11 columns (coverage and frequency retained as extra fields)
      - Sort by chrom then start position
    """
    rows = []
    raw_scores = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            mod_type  = cols[3].strip()
            frequency = float(cols[10]) if len(cols) >= 11 else 100.0

            # Recolor col 9
            cols[8] = mod_to_rgb(mod_type, frequency)

            # Record original score for rescaling
            try:
                raw_scores.append(float(cols[4]))
            except ValueError:
                raw_scores.append(0.0)

            # Pad to 11 columns if needed
            while len(cols) < 11:
                cols.append("0")

            rows.append(cols)

    # Rescale scores to 0-1000 (bedToBigBed requires integer scores in col 5)
    if rows:
        min_s = min(raw_scores)
        max_s = max(raw_scores)
        score_range = max_s - min_s if max_s != min_s else 1.0
        for i, cols in enumerate(rows):
            cols[4] = str(int((raw_scores[i] - min_s) / score_range * 1000))

    # Sort by chrom then start position
    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols[:11]) + "\n")

    orig_min = min(raw_scores) if raw_scores else 0
    orig_max = max(raw_scores) if raw_scores else 0
    print(f"    {len(rows)} sites written -> {output_path}")
    print(f"    Original score range: {orig_min:.4f} – {orig_max:.4f} → rescaled to 0–1000")


def write_autosql(out_dir):
    """
    Write the AutoSql (.as) file describing all 11 columns.
    This is required by bedToBigBed when using extra fields beyond col 9.
    """
    path = os.path.join(out_dir, "rRNA_mods.as")
    with open(path, "w") as f:
        f.write("""\
table rRNA_mods
"Human rRNA modification sites"
(
string  chrom;       "Chromosome (rRNA molecule name)"
uint    chromStart;  "Start position (0-based)"
uint    chromEnd;    "End position"
string  name;        "Modification type (e.g. Am, Y, m6A)"
uint    score;       "Detection score rescaled to 0-1000 (higher = more confident)"
char[1] strand;      "Strand (always + for rRNA)"
uint    thickStart;  "Same as chromStart"
uint    thickEnd;    "Same as chromEnd"
uint    itemRgb;     "Color encoding mod type (hue) and frequency (saturation)"
uint    coverage;    "Number of sequencing reads covering this position"
float   frequency;   "Percentage of reads showing the modification (0-100)"
)
""")
    print(f"    Written {path}")
    return path


def _find_tool(name):
    """Return path to a UCSC binary, checking supplementary/ before PATH."""
    supp = os.path.join(SUPP_DIR, name)
    return supp if os.path.exists(supp) else name


def scramble_data(dest_dir, sites_per_chrom=80):
    """
    Generate structurally valid but fully randomised BED + FASTA files in
    dest_dir using chromosome names/lengths from data/real/FASTA_FAI.
    No real positions, scores, or frequencies are preserved.
    """
    os.makedirs(dest_dir, exist_ok=True)
    fai_src = os.path.join(REAL_DIR, FASTA_FAI)
    if not os.path.exists(fai_src):
        print(f"  ERROR: {fai_src} not found — cannot scramble without the real FAI.")
        sys.exit(1)

    # Parse chromosome names and lengths from the real FAI
    chroms = []
    with open(fai_src) as fh:
        for line in fh:
            parts = line.split("\t")
            chroms.append((parts[0], int(parts[1])))

    mod_types = list(COLOR_CODE.keys())

    # Write a fake FASTA + copy FAI (names/lengths only, no real sequence)
    fa_out  = os.path.join(dest_dir, FASTA_FILE)
    fai_out = os.path.join(dest_dir, FASTA_FAI)
    with open(fa_out, "w") as fa:
        for chrom, length in chroms:
            fa.write(f">{chrom}\n")
            seq = "".join(random.choices("ACGT", k=length))
            for i in range(0, length, 60):
                fa.write(seq[i:i+60] + "\n")
    import shutil
    shutil.copy(fai_src, fai_out)

    # Write one scrambled BED per INPUT_FILE
    for _, fname in INPUT_FILES.items():
        out_path = os.path.join(dest_dir, fname)
        rows = []
        for chrom, length in chroms:
            positions = sorted(random.sample(range(length - 1), min(sites_per_chrom, length - 1)))
            for pos in positions:
                mod   = random.choice(mod_types)
                score = round(random.uniform(0, 10), 4)
                freq  = round(random.uniform(5, 100), 2)
                cov   = random.randint(10, 500)
                rows.append([chrom, str(pos), str(pos + 1), mod, str(score),
                              "+", str(pos), str(pos + 1), "0,0,0",
                              str(cov), str(freq)])
        with open(out_path, "w") as fh:
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print(f"    Scrambled {len(rows)} sites -> {out_path}")

    print(f"  Demo data written to {dest_dir}/")


def run_bedtobigbed(bed_path, fai_path, bigbed_path, as_path):
    """Convert a fixed BED file to bigBed using UCSC's bedToBigBed.
    Uses -type=bed9+2 and -as= to retain coverage and frequency as named extra fields."""
    cmd = _find_tool("bedToBigBed")

    result = subprocess.run(
        [cmd, f"-as={as_path}", "-type=bed9+2", bed_path, fai_path, bigbed_path],
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
trackDb trackDb.txt
description Human rRNA sequences
twoBitPath hs_rRNAs_NR_046235.2bit
organism Homo sapiens
defaultPos hs_rRNA_18S:1-1869
orderKey 1
scientificName Homo sapiens
htmlPath hubDescription.html

genome hg38
trackDb trackDb_hg38.txt
""")
    print(f"    Written {path}")


def get_mod_types(bed_paths):
    """
    Scan one or more fixed BED files and return a sorted, pipe-delimited string
    of unique modification type values (column 4) for use in filterValues.name.
    """
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


def write_trackdb_txt(out_dir, bigbed_names, mod_type_values):
    path = os.path.join(out_dir, "trackDb.txt")
    ref_bb, sample_bb, filtered_bb, consensus_bb = bigbed_names
    with open(path, "w") as f:
        f.write(f"""\
track rRNA_reference
bigDataUrl {ref_bb}
shortLabel Reference mods
longLabel Known rRNA modifications
type bigBed 9 +
itemRgb on
visibility pack
scoreMin 0
scoreMax 1000
html rRNA_reference
priority 1

track rRNA_sample_MRI01
bigDataUrl {sample_bb}
shortLabel MRI01 detected
longLabel MRI01 — EM-detected modifications
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
html rRNA_sample_MRI01
priority 2

track rRNA_filtered
bigDataUrl {filtered_bb}
shortLabel Filtered mods
longLabel Filtered modifications
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
html rRNA_filtered
priority 3

track rRNA_consensus
bigDataUrl {consensus_bb}
shortLabel Consensus mods
longLabel Consensus modifications (sample + filtered); color = observation frequency (yellow = 1/2, green = 2/2)
type bigBed 9 +
itemRgb on
visibility pack
scoreMin 0
scoreMax 1000
html rRNA_consensus
priority 4
""")
    print(f"    Written {path}")


TRACK_HTML_TEMPLATE = """\
<h2>{short_label}</h2>
<p><strong>Long label:</strong> {long_label}</p>

<h3>Methods</h3>
<table border="1" cellpadding="4">
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Sequencing technology</td><td>MISSING</td></tr>
<tr><td>QC applied</td><td>MISSING</td></tr>
</table>

<h3>Contact</h3>
<p>MISSING</p>
"""

TRACK_METADATA = {
    "rRNA_reference": {
        "short_label": "Reference mods",
        "long_label":  "Known rRNA modifications",
    },
    "rRNA_sample_MRI01": {
        "short_label": "MRI01 detected",
        "long_label":  "MRI01 — EM-detected modifications",
    },
    "rRNA_filtered": {
        "short_label": "Filtered mods",
        "long_label":  "Filtered modifications",
    },
    "rRNA_consensus": {
        "short_label": "Consensus mods",
        "long_label":  "Consensus modifications (sample + filtered)",
    },
}


def write_track_html(out_dir):
    """Write one HTML description file per track using TRACK_METADATA."""
    for track_name, meta in TRACK_METADATA.items():
        path = os.path.join(out_dir, f"{track_name}.html")
        with open(path, "w") as f:
            f.write(TRACK_HTML_TEMPLATE.format(**meta))
        print(f"    Written {path}")


def write_description_html(out_dir):
    path = os.path.join(out_dir, "hubDescription.html")
    with open(path, "w") as f:
        f.write("""\
<h1>Human rRNA Modifications</h1>
<p>RNA modification tracks for human rRNA sequences (18S, 28S, 5.8S, 5S).</p>
<p>Track colors encode both modification type (hue) and stoichiometry (saturation).
Vivid colors indicate highly modified sites; pale colors indicate partially modified sites.</p>
<h2>Score</h2>
<p>The score (shading intensity) reflects detection confidence, rescaled to 0-1000 per file.
Darker = higher confidence. Note: scores are rescaled independently per file and are not
directly comparable across tracks.</p>
<h2>Consensus track</h2>
<p>The <strong>Consensus</strong> track merges the sample and filtered tracks.
Colors use a sequential red scale based on how often the modification is observed:</p>
<table border="1" cellpadding="4">
<tr><th>Observation</th><th>Color</th></tr>
<tr><td>1 / 2 samples</td><td style="background:#ffdc00">Yellow</td></tr>
<tr><td>2 / 2 samples</td><td style="background:#00b400;color:white">Green</td></tr>
</table>
<h2>Extra fields (visible on click)</h2>
<p>Clicking any site in the browser shows coverage (raw read depth) and frequency
(% of reads showing the modification).</p>
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
    os.makedirs(SUPP_DIR, exist_ok=True)

    fai_path = os.path.join(input_dir, FASTA_FAI)

    print("\n── Step 1 & 2: Recolor and fix BED files ───────────────────────────────")
    fixed_beds   = {}
    bigbed_names = {}

    for key, fname in INPUT_FILES.items():
        input_path  = os.path.join(input_dir, fname)
        base        = os.path.splitext(fname)[0]
        fixed_path  = os.path.join(SUPP_DIR, base + ".fixed.bed")
        bigbed_name = base + ".bigBed"

        print(f"\n  {key}:")
        process_bed(input_path, fixed_path)
        fixed_beds[key]   = fixed_path
        bigbed_names[key] = bigbed_name

    print("\n── Step 2b: Build consensus track ──────────────────────────────────────")
    consensus_fixed   = os.path.join(SUPP_DIR, "consensus.fixed.bed")
    consensus_bb_name = "consensus.bigBed"
    print(f"\n  consensus ({' + '.join(CONSENSUS_KEYS)}):")
    make_consensus_bed(
        [fixed_beds[k] for k in CONSENSUS_KEYS],
        consensus_fixed,
    )
    fixed_beds["consensus"]   = consensus_fixed
    bigbed_names["consensus"] = consensus_bb_name

    print("\n── Step 3: Write AutoSql file ───────────────────────────────────────────")
    as_path = write_autosql(OUTPUT_DIR)

    print("\n── Step 4: Convert to bigBed ────────────────────────────────────────────")
    for key, fixed_path in fixed_beds.items():
        bigbed_path = os.path.join(OUTPUT_DIR, bigbed_names[key])
        print(f"\n  {key}:")
        run_bedtobigbed(fixed_path, fai_path, bigbed_path, as_path)

    print("\n── Step 5: Write hub config files ──────────────────────────────────────")
    mod_type_values = get_mod_types([fixed_beds["sample"], fixed_beds["filtered"]])
    print(f"    Modification types found: {mod_type_values}")
    write_hub_txt(OUTPUT_DIR)
    write_genomes_txt(OUTPUT_DIR)
    write_trackdb_txt(OUTPUT_DIR, [
        bigbed_names["ref"],
        bigbed_names["sample"],
        bigbed_names["filtered"],
        bigbed_names["consensus"],
    ], mod_type_values)
    write_track_html(OUTPUT_DIR)
    write_description_html(OUTPUT_DIR)

    fasta_src = os.path.join(input_dir, FASTA_FILE)
    print(f"""
── Done ─────────────────────────────────────────────────────────────────────

Next steps:

  1. Convert FASTA to .2bit (if not done already):
         {SUPP_DIR}/faToTwoBit {fasta_src} {OUTPUT_DIR}/hs_rRNAs_NR_046235.2bit

  2. Push to GitHub:
         git add .
         git commit -m "Update hub files"
         git push -u origin main

  3. In UCSC: My Data > Track Hubs > My Hubs
         Paste: https://{GITHUB_USER}.github.io/{GITHUB_REPO}/{OUTPUT_DIR}/hub.txt
""")


if __name__ == "__main__":
    main()

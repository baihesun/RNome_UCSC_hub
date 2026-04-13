"""
pipeline.py

End-to-end pipeline to generate a UCSC assembly hub from rRNA bedRmod files.

Steps:
    1. Recolor BED files — hue encodes modification type, saturation encodes frequency
    2. Fix BED files — keep all 11 columns, rescale scores to 0-1000, sort by chrom/position
    3. Write AutoSql (.as) file describing all 11 columns
    4. Convert to bigBed format using UCSC's bedToBigBed with -type=bed9+2 -as=rRNA_mods.as
    5. Generate hub config files — hub.txt, genomes.txt, trackDb.txt

Usage:
    python pipeline.py

Requirements:
    - bedToBigBed in same directory or on PATH
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
    - faToTwoBit in same directory or on PATH
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit

After running:
    1. Convert FASTA to .2bit:
           ./faToTwoBit Test_visualization_data_rRNA/hs_rRNAs_NR_046235.fa ucsc_hub/hs_rRNAs_NR_046235.2bit
    2. Push the entire repo to GitHub (run from repo root):
           cd /Users/baihesun/Downloads/rRNA_dataviz_repo
           git add .
           git commit -m "Update hub files w dna"
           git push -u origin main
    3. In UCSC: My Data > Track Hubs > My Hubs
           Paste: https://baihesun.github.io/rRNA_dataviz/ucsc_hub/hub.txt
"""

import os
import colorsys
import subprocess
import sys
from collections import defaultdict

# ── Edit these paths ──────────────────────────────────────────────────────────
INPUT_DIR   = "Test_visualization_data_rRNA"
OUTPUT_DIR  = "ucsc_hub"
GITHUB_USER = "baihesun"
GITHUB_REPO = "rRNA_dataviz"

INPUT_FILES = {
    "ref":      "H.sapiens_ref_bedRmod_All.bed",
    "sample":   "rRNA_mature_bedRmod_Detection_EM_mean2_Log10_1.8_modOnly_Sample_MRI01.bed",
    "filtered": "rRNA_mature_Filtered_MOD_10_MULT_1000_bedRmod_0.99_Allmods _no_m5C.bed",
}
FASTA_FAI = "hs_rRNAs_NR_046235.fa.fai"

# Keys from INPUT_FILES to merge into the consensus track
CONSENSUS_KEYS = ["sample", "filtered"]
# ─────────────────────────────────────────────────────────────────────────────

# Base hues (0-1) per modification type
# Hue encodes chemical class; saturation encodes frequency
MOD_HUES = {
    # 2'-O-methylation — blue family
    "Am":  0.61,
    "Cm":  0.58,
    "Gm":  0.55,
    "Um":  0.52,
    "Ym":  0.64,
    # N6/N1-methyladenosine — red/orange family
    "m6A":  0.0,
    "m66A": 0.97,
    "m1A":  0.08,
    # Pseudouridine — teal
    "Y":       0.47,
    "m1acp3Y": 0.44,
    # Cytosine modifications — green
    "m5C":  0.33,
    "ac4C": 0.28,
    # Guanosine modifications — purple
    "m7G":  0.78,
    # Uridine modifications — amber
    "m3U":  0.13,
}
DEFAULT_HUE = 0.0  # gray fallback for unknown mod types


def frequency_to_rgb(mod_type, frequency):
    """
    Convert a mod type + frequency value to an RGB string for UCSC.
    Hue        = modification type (chemical class)
    Saturation = frequency (0% = gray, 100% = fully vivid)
    Value      = fixed at 0.85
    """
    hue = MOD_HUES.get(mod_type.strip(), DEFAULT_HUE)
    freq = max(0.0, min(100.0, float(frequency)))
    r, g, b = colorsys.hsv_to_rgb(hue, freq / 100.0, 0.85)
    return f"{int(r*255)},{int(g*255)},{int(b*255)}"


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
            cols[8] = frequency_to_rgb(mod_type, frequency)

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


def run_bedtobigbed(bed_path, fai_path, bigbed_path, as_path):
    """Convert a fixed BED file to bigBed using UCSC's bedToBigBed.
    Uses -type=bed9+2 and -as= to retain coverage and frequency as named extra fields."""
    cmd = "bedToBigBed"
    if os.path.exists("./bedToBigBed"):
        cmd = "./bedToBigBed"

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
""")
    print(f"    Written {path}")


def write_trackdb_txt(out_dir, bigbed_names):
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
priority 4
""")
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
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fai_path = os.path.join(INPUT_DIR, FASTA_FAI)

    print("\n── Step 1 & 2: Recolor and fix BED files ───────────────────────────────")
    fixed_beds   = {}
    bigbed_names = {}

    for key, fname in INPUT_FILES.items():
        input_path  = os.path.join(INPUT_DIR, fname)
        base        = os.path.splitext(fname)[0]
        fixed_path  = os.path.join(OUTPUT_DIR, base + ".fixed.bed")
        bigbed_name = base + ".bigBed"

        print(f"\n  {key}:")
        process_bed(input_path, fixed_path)
        fixed_beds[key]   = fixed_path
        bigbed_names[key] = bigbed_name

    print("\n── Step 2b: Build consensus track ──────────────────────────────────────")
    consensus_fixed  = os.path.join(OUTPUT_DIR, "consensus.fixed.bed")
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
    write_hub_txt(OUTPUT_DIR)
    write_genomes_txt(OUTPUT_DIR)
    write_trackdb_txt(OUTPUT_DIR, [
        bigbed_names["ref"],
        bigbed_names["sample"],
        bigbed_names["filtered"],
        bigbed_names["consensus"],
    ])
    write_description_html(OUTPUT_DIR)

    print(f"""
── Done ─────────────────────────────────────────────────────────────────────

Next steps:

  1. Convert FASTA to .2bit (if not done already):
         ./faToTwoBit {os.path.join(INPUT_DIR, "hs_rRNAs_NR_046235.fa")} {OUTPUT_DIR}/hs_rRNAs_NR_046235.2bit

  2. Push the entire repo to GitHub (run from repo root, not from {OUTPUT_DIR}/):
         cd /Users/{GITHUB_USER}/Downloads/rRNA_dataviz_repo
         git add .
         git commit -m "Update hub files"
         git push -u origin main

  3. In UCSC: My Data > Track Hubs > My Hubs
         Paste: https://{GITHUB_USER}.github.io/{GITHUB_REPO}/{OUTPUT_DIR}/hub.txt
""")


if __name__ == "__main__":
    main()

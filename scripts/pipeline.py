"""
pipeline.py

End-to-end pipeline to generate a UCSC assembly hub from rRNA bedRmod files.

Steps:
    1. Recolor BED files — hue encodes modification type, saturation encodes frequency
    2. Fix BED files — trim to 9 columns, sort by chrom/position
    3. Convert to bigBed format using UCSC's bedToBigBed
    4. Generate hub config files — hub.txt, genomes.txt, trackDb.txt

Usage:
    python pipeline.py

Requirements:
    - bedToBigBed must be in the same directory as the parent repo, or on your PATH
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed
    - faToTwoBit must also be available (same location)
      Download: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/faToTwoBit

After running:
    1. Convert FASTA to .2bit:
           ./faToTwoBit <FASTA_FILE> ucsc_hub/hs_rRNAs_NR_046235.2bit
    2. Push ucsc_hub/ to a public GitHub repo with Pages enabled
    3. In UCSC: My Data > Track Hubs > My Hubs
           Paste: https://<username>.github.io/<repo>/hub.txt
"""

import os
import colorsys
import subprocess
import sys

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
FASTA_FAI  = "hs_rRNAs_NR_046235.fa.fai"
# ─────────────────────────────────────────────────────────────────────────────

# Base hues (0-1) per modification type
# Hue encodes chemical class; saturation will encode frequency
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
DEFAULT_HUE = 0.0  # gray fallback — will be desaturated anyway


def frequency_to_rgb(mod_type, frequency):
    """
    Convert a mod type + frequency value to an RGB string for UCSC.

    Hue    = modification type (chemical class)
    Saturation = frequency (0% mod = fully desaturated/gray, 100% = fully vivid)
    Value  = fixed at 0.85 to keep colors visible but not garish
    """
    hue = MOD_HUES.get(mod_type.strip(), DEFAULT_HUE)
    # Clamp frequency to 0-100 range
    freq = max(0.0, min(100.0, float(frequency)))
    saturation = freq / 100.0
    value = 0.85
    r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
    return f"{int(r*255)},{int(g*255)},{int(b*255)}"


def rescale_score(scores):
    """Rescale a list of float scores to 0-1000 integer range for UCSC."""
    min_s = min(scores)
    max_s = max(scores)
    score_range = max_s - min_s if max_s != min_s else 1.0
    return [str(int((s - min_s) / score_range * 1000)) for s in scores]


def process_bed(input_path, output_path):
    """
    Read a bedRmod file and write a cleaned version:
      - Recolor column 9 (itemRgb) using hue=mod type, saturation=frequency
      - Rescale column 5 (score) to 0-1000
      - Trim to 9 columns (drop coverage and frequency cols used internally)
      - Sort by chrom (col 1) then start position (col 2)
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

            mod_type = cols[3].strip()

            # Use frequency from col 11 if available, else default to 100
            frequency = float(cols[10]) if len(cols) >= 11 else 100.0

            # Recolor col 9 based on mod type + frequency
            cols[8] = frequency_to_rgb(mod_type, frequency)

            try:
                raw_scores.append(float(cols[4]))
            except ValueError:
                raw_scores.append(0.0)

            rows.append(cols)

    # Rescale scores
    if rows:
        rescaled = rescale_score(raw_scores)
        for i, cols in enumerate(rows):
            cols[4] = rescaled[i]

    # Trim to 9 columns and sort
    rows = [cols[:9] for cols in rows]
    rows.sort(key=lambda c: (c[0], int(c[1])))

    with open(output_path, "w") as out:
        for cols in rows:
            out.write("\t".join(cols) + "\n")

    print(f"    {len(rows)} sites written -> {output_path}")


def run_bedtobigbed(bed_path, fai_path, bigbed_path):
    """Convert a fixed BED file to bigBed using UCSC's bedToBigBed."""
    # Look for bedToBigBed in current directory or PATH
    cmd = "bedToBigBed"
    if os.path.exists("./bedToBigBed"):
        cmd = "./bedToBigBed"

    result = subprocess.run(
        [cmd, bed_path, fai_path, bigbed_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr.strip()}")
        sys.exit(1)
    print(f"    bigBed written -> {bigbed_path}")


def write_hub_txt(out_dir):
    path = os.path.join(out_dir, "hub.txt")
    with open(path, "w") as f:
        f.write(f"""\
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
    ref_bb, sample_bb, filtered_bb = bigbed_names
    with open(path, "w") as f:
        f.write(f"""\
track rRNA_reference
bigDataUrl {ref_bb}
shortLabel Reference mods
longLabel Known H. sapiens rRNA modifications (all types, curated)
type bigBed 9
itemRgb on
visibility pack
priority 1

track rRNA_sample_MRI01
bigDataUrl {sample_bb}
shortLabel MRI01 detected
longLabel MRI01 — EM-detected modifications 
type bigBed 9
itemRgb on
visibility pack
priority 2

track rRNA_filtered
bigDataUrl {filtered_bb}
shortLabel Filtered mods
longLabel Filtered modifications 
type bigBed 9
itemRgb on
visibility pack
priority 3
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
<h2>Modification type color key</h2>
<table>
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
    fixed_beds = {}
    bigbed_names = {}

    for key, fname in INPUT_FILES.items():
        input_path = os.path.join(INPUT_DIR, fname)
        base = os.path.splitext(fname)[0]
        fixed_path  = os.path.join(OUTPUT_DIR, base + ".fixed.bed")
        bigbed_name = base + ".bigBed"
        bigbed_path = os.path.join(OUTPUT_DIR, bigbed_name)

        print(f"\n  {key}:")
        process_bed(input_path, fixed_path)
        fixed_beds[key]   = fixed_path
        bigbed_names[key] = bigbed_name

    print("\n── Step 3: Convert to bigBed ────────────────────────────────────────────")
    for key, fixed_path in fixed_beds.items():
        bigbed_path = os.path.join(OUTPUT_DIR, bigbed_names[key])
        print(f"\n  {key}:")
        run_bedtobigbed(fixed_path, fai_path, bigbed_path)

    print("\n── Step 4: Write hub config files ──────────────────────────────────────")
    write_hub_txt(OUTPUT_DIR)
    write_genomes_txt(OUTPUT_DIR)
    write_trackdb_txt(OUTPUT_DIR, [
        bigbed_names["ref"],
        bigbed_names["sample"],
        bigbed_names["filtered"],
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

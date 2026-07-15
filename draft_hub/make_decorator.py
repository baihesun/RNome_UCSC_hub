"""
make_decorator.py — build a decorator bigBed for the HRP draft sequence hub.

Reads the polyA consensus tiered bed and the HRP PSL file, intersects them,
and writes a decorator bigBed that marks each modification position on every
overlapping PSL alignment item, colored by modification type.

Run from the repo root:
    python draft_hub/make_decorator.py

Output: draft_hub/hg38/HRP_decorator.bigBed
"""

import bisect
import os
import subprocess
import sys

# ── Paths ──────────────────────────────────────────────────────────────────────

PSL_FILE      = "draft_sequence/HRP.psl"
MOD_BED       = "final_data/polyA_RNA_hg38/polyA_RNA_hg38_consensus_tiered.bed"
CHROM_SIZES   = "utils/hg38.chrom.sizes"
BEDTOBIGBED   = "utils/bedToBigBed"
OUT_DIR       = "draft_hub/hg38"
OUT_BED       = os.path.join(OUT_DIR, "HRP_decorator.bed")
OUT_AS        = os.path.join(OUT_DIR, "decorator.as")
OUT_BIGBED    = os.path.join(OUT_DIR, "HRP_decorator.bigBed")

# ── Modification colors ────────────────────────────────────────────────────────

COLOR_CODE = {
    "m6A": "#003f5c", "Am": "#1464a0", "Ino": "#0e8fd4", "I": "#0e8fd4",
    "m1A": "#00b4d8", "m66A": "#00d4b4", "m1acp3Y": "#00c875",
    "Y":   "#90be6d", "Ym":  "#c9b84c", "Um": "#ffa600", "m3U": "#ff7c43",
    "Gm":  "#ff6361", "m7G": "#dd7371", "Cm": "#ef476f", "m5C": "#d45087",
    "ac4C":"#7b2d8b",
}
DEFAULT_HEX = "#808080"

# ── AutoSQL schema ─────────────────────────────────────────────────────────────

DECORATOR_AS = """\
table decorator
"RNA modification decorators for HRP sequence alignment"
(
string   chrom;              "Chromosome"
uint     chromStart;         "Start position of modification"
uint     chromEnd;           "End position of modification"
string   name;               "Modification name (e.g. m6A, m5C)"
uint     score;              "Score: 1000=tier1, 500=tier2"
char[1]  strand;             "Strand (+ or -)"
uint     thickStart;         "Same as chromStart"
uint     thickEnd;           "Same as chromEnd"
uint     color;              "Primary RGB color packed as uint"
int      blockCount;         "Number of blocks (always 1)"
int[blockCount] blockSizes;  "Block size (modification width)"
int[blockCount] chromStarts; "Block start relative to chromStart (always 0)"
string   decoratedItem;      "Parent PSL item: chr:tStart-tEnd:qName"
string   style;              "Decoration style: block"
string   fillColor;          "Fill color as hex string"
string   glyph;              "Glyph type (ignored for block style)"
)
"""

# ── Helpers ────────────────────────────────────────────────────────────────────

def hex_to_uint(hexcode):
    """Convert #RRGGBB to packed RGB uint."""
    h = hexcode.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return (r << 16) | (g << 8) | b


def load_valid_chroms(sizes_path):
    chroms = set()
    with open(sizes_path) as fh:
        for line in fh:
            chroms.add(line.split()[0])
    return chroms


def load_psl_items(psl_path, valid_chroms):
    """
    Parse PSL file and return per-chrom sorted list of
    (tStart, tEnd, qName) tuples.  Skips the 5-line PSL header.
    """
    items = {}  # chrom → sorted list of (tStart, tEnd, qName)
    with open(psl_path) as fh:
        for _ in range(5):
            next(fh, None)  # skip psLayout header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 17:
                continue
            t_name  = cols[13]
            t_start = int(cols[15])
            t_end   = int(cols[16])
            q_name  = cols[9]
            if t_name not in valid_chroms:
                continue
            items.setdefault(t_name, []).append((t_start, t_end, q_name))

    for chrom in items:
        items[chrom].sort()
    print(f"  Loaded {sum(len(v) for v in items.values()):,} PSL items "
          f"across {len(items)} chromosomes")
    return items


def find_overlapping_psl(psl_items, chrom, start, end):
    """Return list of (tStart, tEnd, qName) that overlap [start, end)."""
    chrom_items = psl_items.get(chrom, [])
    if not chrom_items:
        return []
    starts = [x[0] for x in chrom_items]
    # All items whose tStart < end and tEnd > start
    idx = bisect.bisect_left(starts, end)  # first item starting at or after end
    results = []
    for i in range(idx - 1, -1, -1):
        t_start, t_end, q_name = chrom_items[i]
        if t_end <= start:
            break
        results.append((t_start, t_end, q_name))
    return results


def load_mod_sites(bed_path):
    """Read consensus tiered bed → list of (chrom, start, end, name, score, strand)."""
    sites = []
    with open(bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name, score, strand = cols[3], cols[4], cols[5]
            sites.append((chrom, start, end, name, score, strand))
    return sites


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading chromosome sizes...")
    valid_chroms = load_valid_chroms(CHROM_SIZES)

    print("Loading PSL items...")
    psl_items = load_psl_items(PSL_FILE, valid_chroms)

    print("Loading modification sites...")
    mod_sites = load_mod_sites(MOD_BED)
    print(f"  {len(mod_sites):,} modification sites")

    print("Intersecting modifications × PSL items...")
    rows = []
    no_overlap = 0
    for chrom, start, end, name, score, strand in mod_sites:
        hex_color = COLOR_CODE.get(name, DEFAULT_HEX)
        color_uint = hex_to_uint(hex_color)
        block_size = end - start  # always 1 for these sites

        overlapping = find_overlapping_psl(psl_items, chrom, start, end)
        if not overlapping:
            no_overlap += 1
            continue

        for t_start, t_end, q_name in overlapping:
            decorated_item = f"{chrom}:{t_start}-{t_end}:{q_name}"
            rows.append([
                chrom,
                str(start),
                str(end),
                name,
                score,
                strand,
                str(start),          # thickStart
                str(end),            # thickEnd
                str(color_uint),     # color as uint
                "1",                 # blockCount
                f"{block_size},",    # blockSizes (trailing comma = UCSC array)
                "0,",                # chromStarts
                decorated_item,
                "block",
                hex_color,           # fillColor
                ".",                 # glyph (ignored for block)
            ])

    print(f"  {len(rows):,} decorator entries "
          f"({no_overlap} mod sites had no overlapping PSL item)")

    print(f"Sorting and writing {OUT_BED} ...")
    rows.sort(key=lambda r: (r[0], int(r[1])))
    with open(OUT_BED, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"Writing AutoSQL schema → {OUT_AS}")
    with open(OUT_AS, "w") as fh:
        fh.write(DECORATOR_AS)

    print(f"Running bedToBigBed → {OUT_BIGBED}")
    cmd = [
        BEDTOBIGBED,
        "-type=bed12+4",
        "-tab",
        f"-as={OUT_AS}",
        OUT_BED,
        CHROM_SIZES,
        OUT_BIGBED,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] {result.stderr.strip()}")
        sys.exit(1)

    size_mb = os.path.getsize(OUT_BIGBED) / 1e6
    print(f"  Done → {OUT_BIGBED} ({size_mb:.1f} MB)")
    print(f"\nAdd to draft_hub/hg38/trackDb.txt:")
    print(f"  decorator.default.bigDataUrl HRP_decorator.bigBed")


if __name__ == "__main__":
    main()

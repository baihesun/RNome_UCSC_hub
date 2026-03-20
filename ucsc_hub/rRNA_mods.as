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

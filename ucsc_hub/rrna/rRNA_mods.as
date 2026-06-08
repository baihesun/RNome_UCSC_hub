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

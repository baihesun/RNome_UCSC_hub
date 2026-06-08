table rna_mods_v2
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

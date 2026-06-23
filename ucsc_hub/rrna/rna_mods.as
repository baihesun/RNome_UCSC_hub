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

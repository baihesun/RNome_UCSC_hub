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

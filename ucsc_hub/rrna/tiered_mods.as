table tiered_mods
"Tiered cross-platform consensus RNA modification sites"
(
string  chrom;      "Chromosome or RNA name"
uint    chromStart; "Start position"
uint    chromEnd;   "End position"
string  name;       "Modification short name"
uint    score;      "Display score (1000=tier1, 500=tier2)"
char[1] strand;     "Strand (+ or -)"
uint    thickStart; "Same as chromStart"
uint    thickEnd;   "Same as chromEnd"
uint    itemRgb;    "Display color (R,G,B) — hue by modification type, vivid=tier1/pale=tier2"
string  tier;       "Confidence tier: tier1 = strongest cross-platform evidence, tier2 = supported by >=2 datasets"
)

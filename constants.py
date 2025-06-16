RESTRICTION_SITES = {
    "BsmBI": ["CGTCTC", "GAGACG"],
    "BsaI": ["GGTCTC", "GAGACC"],
    "PmeI": ["GTTTAAAC"],
    "NotI": ["GCGGCCGC"]
}

OVERHANGS = {
                    "Promoter": ("GCATCGTCTCATCGGTCTCAAACG", "GGAAGATGAGACCTGAGACGGCAT"),
                    "N-terminal tag": ("GCATCGTCTCATCGGTCTCAAAGATG", "TCCGGTATGTGAGACCTGAGACGGCAT"),
                    "CDS": ("GCATCGTCTCATCGGTCTCATATG", "GGATCCTGAGACCTGAGACGGCAT"),
                    "C-terminal tag": ("GCATCGTCTCATCGGTCTCAATCCTCAGGT", "TCCGGTGGCTGAGACCTGAGACGGCAT"),
                    "Terminator": ("GCATCGTCTCATCGGTCTCATGGCTAATGA", "GGGCTGTGAGACCTGAGACGGCAT")
                }
NO_TAG_OVERHANGS = {
                    "Promoter": ("GCATCGTCTCATCGGTCTCAAACG", "TCCGGTATGTGAGACCTGAGACGGCAT"),
                    "Terminator": ("GCATCGTCTCATCGGTCTCAATCCTAATGA", "GGGCTGTGAGACCTGAGACGGCAT")
                }

AA_TO_CODON = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['CTT', 'TTG', 'TTA', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*': ['TAA', 'TAG', 'TGA'],
}

CODON_TO_AA = {codon: aa for aa, codons in AA_TO_CODON.items() for codon in codons}

HELPER_RULES = [
    ("pL6", ["Pos1", "Pos2", "Pos3", "Pos4", "Pos5"]),
    ("pL56", ["Pos1", "Pos2", "Pos3", "Pos4"]),
    ("pL456", ["Pos1", "Pos2", "Pos3"]),
    ("pL3456", ["Pos1", "Pos2"]),
]

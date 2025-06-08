import re
restriction_sites = {
    "BsmBI": ["CGTCTC", "GAGACG"],
    "BsaI": ["GGTCTC", "GAGACC"],
    "PmeI": ["GTTTAAAC"]
}

amino_acid_to_codon = {
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
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
}

codon_to_aa = {codon: aa for aa, codons in amino_acid_to_codon.items() for codon in codons}

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def remove_stop_codon(dna):
    if dna[-3:] in ['TAA', 'TAG', 'TGA']:
        return dna[:-3]
    return dna

def optimize_dna_sequence(dna_seq):
    max_attempts = 100
    for _ in range(max_attempts):
        modified = False
        for enzyme, patterns in restriction_sites.items():
            for site in patterns:
                for strand in [dna_seq, reverse_complement(dna_seq)]:
                    match = re.search(site, strand)
                    if match:
                        start = match.start()
                        pos = (start // 3) * 3
                        codon = dna_seq[pos:pos+3]
                        if codon in codon_to_aa:
                            aa = codon_to_aa[codon]
                            alternatives = [c for c in amino_acid_to_codon[aa] if c != codon]
                            if alternatives:
                                dna_seq = dna_seq[:pos] + alternatives[0] + dna_seq[pos+3:]
                                modified = True
                        break
        if not modified:
            break
    return dna_seq

def count_re_sites(dna_seq):
    dna_seq = dna_seq.upper()

    counts = {}
    for enzyme, sites in restriction_sites.items():
        count = sum(len(re.findall(site, dna_seq)) for site in sites)
        counts[enzyme] = count
    return counts

def highlight_dna_and_protein(dna_seq, width=81):
    from html import escape

    dna_seq = dna_seq.upper().replace("\n", "").replace(" ", "")
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }

    sites = restriction_sites

    highlight_pos = set()
    for patterns in sites.values():
        for pattern in patterns:
            start = 0
            while True:
                start = dna_seq.find(pattern, start)
                if start == -1:
                    break
                highlight_pos.update(range(start, start + len(pattern)))
                start += 1

    html_output = ""
    for i in range(0, len(dna_seq), width):
        line = dna_seq[i:i + width]
        
        # Highlight DNA
        highlighted_dna = ""
        for j, nt in enumerate(line):
            idx = i + j
            if idx in highlight_pos:
                highlighted_dna += f"<span style='background-color:yellow'>{escape(nt)}</span>"
            else:
                highlighted_dna += escape(nt)

        # Translate DNA and space-align
        highlighted_protein = ""
        for j in range(0, len(line) - 2, 3):
            codon = line[j:j + 3]
            idxs = [i + j, i + j + 1, i + j + 2]
            aa = codon_table.get(codon, "X")
            spacer = "  "  # two spaces between AAs
            if any(x in highlight_pos for x in idxs):
                highlighted_protein += f"<span style='background-color:orange'>{aa}</span>{spacer}"
            else:
                highlighted_protein += f"{aa}{spacer}"

        html_output += (
            f"<div style='font-family:monospace; white-space:pre;'>"
            f"{highlighted_dna}\n{highlighted_protein.strip()}</div>\n"
        )

    return html_output

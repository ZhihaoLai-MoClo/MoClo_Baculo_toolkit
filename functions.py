import re
import requests
from constants import RESTRICTION_SITES,CODON_TO_AA,AA_TO_CODON

def BsmBI_dg(sequence):
    disequence= sequence + sequence
    fwd_site = "CGTCTC"
    rev_site = "GAGACG"
    fwd_index = disequence.find(fwd_site)
    cutfwd = disequence[fwd_index +7:]
    rev_index = cutfwd.find(rev_site)
    return cutfwd[:rev_index -5]

def BsaI_dg(sequence):
    disequence= sequence + sequence
    fwd_site = "GGTCTC"
    rev_site = "GAGACC"
    fwd_index = disequence.find(fwd_site)
    cutfwd = disequence[fwd_index +7:]
    rev_index = cutfwd.find(rev_site)
    return cutfwd[:rev_index -5]

def BsmBI_dg_rev(sequence):
    disequence= sequence + sequence
    fwd_site = "CGTCTC"
    rev_site = "GAGACG"
    rev_index = disequence.find(rev_site)
    cutrev = disequence[rev_index -7:]
    rev_index = cutrev.find(fwd_site)
    return cutrev[:rev_index +5]

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def rm_stop_codon(dna):
    while len(dna) >= 3: 
        if dna[-3:] in ['TAA', 'TAG', 'TGA']:
            dna = dna[:-3]  
        else:
            break  
    return dna

def rm_start_codon(dna):
    if dna[:3] in ['ATG']:
        return dna[3:]
    return dna

def translate_and_find_stops(dna_seq):
    dna_seq = re.sub(r'\s+', '', dna_seq.upper())  # Clean up input
    protein_seq = ""
    stop_positions = []

    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = CODON_TO_AA.get(codon, 'X')
        protein_seq += aa
        if aa == '*':
            codon_index = i // 3
            stop_positions.append(codon_index)

    return protein_seq, stop_positions


def optimize_dna_sequence(dna_seq):
    max_attempts = 100
    for _ in range(max_attempts):
        modified = False
        for enzyme, patterns in RESTRICTION_SITES.items():
            for site in patterns:
                match = re.search(site, dna_seq)
                if match:
                    start = match.start()
                    pos = (start // 3) * 3
                    codon = dna_seq[pos:pos+3]
                    if codon in CODON_TO_AA:
                        aa = CODON_TO_AA[codon]
                        alternatives = [c for c in AA_TO_CODON[aa] if c != codon]
                        if alternatives:
                            dna_seq = dna_seq[:pos] + alternatives[0] + dna_seq[pos+3:]
                            modified = True
                    break
        if not modified:
            break
    for _ in range(max_attempts):
        modified = False
        for enzyme, patterns in RESTRICTION_SITES.items():
            for site in patterns:
                match = re.search(site, dna_seq)
                if match:
                    start = match.start()
                    pos = (start // 3) * 3
                    codon = dna_seq[pos+3:pos+6]
                    if codon in CODON_TO_AA:
                        aa = CODON_TO_AA[codon]
                        alternatives = [c for c in AA_TO_CODON[aa] if c != codon]
                        if alternatives:
                            dna_seq = dna_seq[:pos+3] + alternatives[0] + dna_seq[pos+6:]
                            modified = True
                    break
        if not modified:
            break
    return dna_seq



def count_re_sites(dna_seq):
    dna_seq = dna_seq.upper()

    counts = {}
    for enzyme, sites in RESTRICTION_SITES.items():
        count = sum(len(re.findall(site, dna_seq)) for site in sites)
        counts[enzyme] = count
    return counts

def highlight_dna_and_protein(dna_seq, width=81):
    from html import escape

    dna_seq = dna_seq.upper().replace("\n", "").replace(" ", "")
    sites = RESTRICTION_SITES

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
            aa = CODON_TO_AA.get(codon, "X")
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

def get_ccds_link(uniprot_id):
    base_id = uniprot_id.split('-')[0]
    url = f"https://rest.uniprot.org/uniprotkb/{base_id}.json"
    res = requests.get(url)

    if res.status_code != 200:
        return None, f"Failed to retrieve UniProt entry for {base_id}."

    data = res.json()
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref["database"] == "CCDS":
            ccds_id = xref["id"]
            link = f"https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA={ccds_id}"
            return link, None

    return None, "No CCDS entry found in UniProt cross-references."

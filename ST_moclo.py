import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import StringIO
import requests
from xml.etree import ElementTree
import re
from rm_site import restriction_sites, amino_acid_to_codon, codon_to_aa, remove_stop_codon, optimize_dna_sequence, count_re_sites, highlight_dna_and_protein

st.set_page_config(page_title="MoClo-YTK Assembly Tool", layout="wide")
st.title("MoClo Cloning Platform")

HELPER_RULES = [
    ("pL6", ["Pos1", "Pos2", "Pos3", "Pos4", "Pos5"]),
    ("pL56", ["Pos1", "Pos2", "Pos3", "Pos4"]),
    ("pL456", ["Pos1", "Pos2", "Pos3"]),
    ("pL3456", ["Pos1", "Pos2"]),
]

def BsmBI_dg(sequence):
    disequence= sequence + sequence
    fwd_site = "CGTCTC"
    rev_site = "GAGACG"
    fwd_index = disequence.find(fwd_site)
    cutfwd = disequence[fwd_index +6:]
    rev_index = cutfwd.find(rev_site)
    return cutfwd[:rev_index -5]

def BsaI_dg(sequence):
    disequence= sequence + sequence
    fwd_site = "GGTCTC"
    rev_site = "GAGACC"
    fwd_index = disequence.find(fwd_site)
    cutfwd = disequence[fwd_index +6:]
    rev_index = cutfwd.find(rev_site)
    return cutfwd[:rev_index -5]

page = st.sidebar.selectbox("Select Mode", [
    "Level 1 Assembly", 
    "Level 2 Assembly", 
    "Generate new CDS",
    "Restriction Site Remover",
    "Input New Sequence",
    "CSV editing and FSTA download",
    ])


uploaded_file = st.file_uploader("📄 Upload CSV (name, type, sequence, no_tag_sequence, parts_used)", type="csv")

if uploaded_file:
    df = pd.read_csv(uploaded_file)

    if "parts_used" not in df.columns:
        df["parts_used"] = ""

    if page == "Restriction Site Remover":
        st.header("Restriction Site Remover for CDS")
        st.markdown("Enter a coding sequence below. The tool will remove BsaI, BsmBI, and PmeI sites.")
        dna_input = st.text_area("Input DNA Sequence (ATGC only)", height=200).strip().upper()
        dna_input = re.sub(r'\s+', '', dna_input.upper())

        if dna_input:
            site_counts = count_re_sites(dna_input)
            st.subheader("Restriction Site Count in Original Sequence")
            for enzyme, count in site_counts.items():
                st.write(f"🔹 {enzyme}: {count} site(s)")
        if dna_input[:3] not in ['ATG']:
            st.text(f"No ATG found")
        if dna_input[-3:] in ['TAA', 'TAG', 'TGA']:
            st.text(f"STOP codon not removed!")
        
        if st.button("Show restriction sites"):
            result = highlight_dna_and_protein(dna_input)
            st.markdown(result, unsafe_allow_html=True)
        
        if st.button("Visualise optimised DNA"):
            if not dna_input or not re.fullmatch(r'[ATGC]+', dna_input):
                st.error("Please enter a valid DNA sequence containing only A, T, G, C.")
            else:
                st.text(f"Original sequence length: {len(dna_input)} bp")
                dna_cleaned = remove_stop_codon(dna_input)
                dna_optimized = optimize_dna_sequence(dna_cleaned)
                st.session_state.dna_optimized = dna_optimized
                st.success("Optimized DNA Sequence (Restriction sites removed)")
                st.code(dna_optimized, language="text")
                st.text(f"Optimized sequence length: {len(dna_optimized)} bp")
                site_counts = count_re_sites(dna_optimized)
                st.subheader("Restriction Site Count again")
                for enzyme, count in site_counts.items():
                    st.write(f"🔹 {enzyme}: {count} site(s)")
                optimised_result = highlight_dna_and_protein(dna_optimized)
                st.markdown(optimised_result, unsafe_allow_html=True)
        new_name = st.text_input("Part Name")
        selected_type = st.selectbox("Part Type", ["N-terminal tag", "CDS", "C-terminal tag"])
        if st.button("Add optimised sequence to CSV file"):
            dna_cleaned = remove_stop_codon(dna_input)
            dna_optimized = optimize_dna_sequence(dna_cleaned)
            if dna_optimized.startswith("ATG"):
                dna_optimized = dna_optimized[-3:]
                overhangs = {
                    "N-terminal tag": ("TCGGTCTCAAAGATG", "TCCGGTATGTGAGACC"),
                    "CDS": ("TCGGTCTCATATG", "GGATCCTGAGACC"),
                    "C-terminal tag": ("TCGGTCTCAATCCTCAGGT", "TCCGGTGGCTGAGACC"),
                }
                left_oh, right_oh = overhangs.get(selected_type, ("", ""))
                final_seq = left_oh + dna_optimized + right_oh
                new_row = pd.DataFrame([{
                    "name": new_name,
                    "type": selected_type,
                    "sequence": final_seq,
                    "no_tag_sequence": "",
                    "parts_used": ""
                }])
            df = pd.concat([df, new_row], ignore_index=True)
            csv = df.to_csv(index=False)
            st.success("New part added with overhangs")
            st.download_button("Download Updated CSV", csv, "updated_parts.csv", "text/csv")
   

    if page == "Level 1 Assembly":
        st.header("Level 1 Assembly")
        num_constructs = st.slider("Number of constructs", 1, 10, 2)
        types = ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"]
        backbones = df[df['type'] == 'Backbone']['name'].tolist()

        constructs = []
        for i in range(num_constructs):
            st.subheader(f"Construct {i+1}")
            construct = {}
            construct['name'] = st.text_input(f"Construct name {i+1}", key=f"name_{i}")
            part_names = []
            seq_parts = []
            start = 1
            parts_used = []

            for t in types:
                options = df[df['type'] == t]['name'].tolist()
                options = ["None"] + options if t in ["N-terminal tag", "C-terminal tag"] else options
                sel = st.selectbox(f"{t} [{i+1}]", options, key=f"{t}_{i}")
                if sel != "None":
                    row = df[df['name'] == sel].iloc[0]
                    seq = BsaI_dg(row['sequence'])
                elif t == "N-terminal tag":
                    row = df[df['type'] == 'Promoter'].iloc[0]
                    seq = BsaI_dg(row['no_tag_sequence'])
                    sel = f"{row['name']}_no_tag"
                elif t == "C-terminal tag":
                    row = df[df['type'] == 'Terminator'].iloc[0]
                    seq = BsaI_dg(row['no_tag_sequence'])
                    sel = f"{row['name']}_no_tag"
                else:
                    continue
                end = start + len(seq) - 1
                parts_used.append(f"{sel}({start}-{end})")
                seq_parts.append(seq)
                start = end + 1
                part_names.append(sel)

            backbone = st.selectbox(f"Backbone for construct {i+1}", backbones, key=f"bb_{i}")
            row = df[df['name'] == backbone].iloc[0]
            seq = BsaI_dg(row['sequence'])
            end = start + len(seq) - 1
            parts_used.append(f"{backbone}({start}-{end})")
            seq_parts.append(seq)

            
            construct['sequence'] = BsmBI_dg(''.join(seq_parts))
            construct['parts_used'] = '|'.join(parts_used)
            construct['type'] = backbone
            construct['no_tag_sequence'] = ""
            constructs.append(construct)

        if st.button("Assemble Constructs"):
            new_df = pd.DataFrame(constructs)
            st.dataframe(new_df)
            all_df = pd.concat([df, new_df[df.columns]], ignore_index=True)
            csv = all_df.to_csv(index=False)
            st.download_button("Download Updated CSV", csv, "assembled_constructs.csv", "text/csv")

    elif page == "Level 2 Assembly":
        st.header("Level 2 Assembly")
        pos_cols = [f"Pos{i}" for i in range(1, 7)]
        lv2bb_options = df[df["type"].str.contains("Lv2BB")]["name"].tolist()
        filtered_df = df[df["type"].str.contains("Pos|Lv2BB") | df["name"].isin([h for h, _ in HELPER_RULES])]
        num_lv2 = st.slider("Number of constructs", 1, 5, 2)
        constructs = []

        for i in range(num_lv2):
            st.subheader(f"Lv2 Construct {i+1}")
            c = {"Construct Name": st.text_input(f"Name [{i+1}]", key=f"lv2_name_{i}")}
            for pos in pos_cols:
                options = ["None"] + filtered_df[filtered_df["type"].str.contains(pos)]["name"].tolist()
                c[pos] = st.selectbox(f"{pos} [{i+1}]", options, key=f"{pos}_{i}")
            c["Lv2BB"] = st.selectbox(f"Lv2BB [{i+1}]", lv2bb_options, key=f"Lv2BB_{i}")
            constructs.append(c)

        if st.button("Assemble Lv2 Constructs"):
            output = []
            for c in constructs:
                seq_parts = []
                parts_with_coords = []
                curr_pos = 1
                for pos in pos_cols:
                    if c[pos] != "None":
                        row = df[df["name"] == c[pos]].iloc[0]
                        seq = row["sequence"]
                        seq_parts.append(seq)
                        end = curr_pos + len(seq) - 1
                        parts_with_coords.append(f"{c[pos]}({curr_pos}-{end})")
                        curr_pos = end + 1

                for helper, required in HELPER_RULES:
                    if all(c[p] != "None" for p in required):
                        row = df[df["name"] == helper]
                        if not row.empty:
                            row = row.iloc[0]
                            seq = row["sequence"]
                            seq_parts.append(seq)
                            end = curr_pos + len(seq) - 1
                            parts_with_coords.append(f"{helper}({curr_pos}-{end})")
                            curr_pos = end + 1
                        break

                lv2bb = df[df["name"] == c["Lv2BB"]].iloc[0]
                seq = BsmBI_dg(lv2bb["sequence"])
                seq_parts.append(seq)
                end = curr_pos + len(seq) - 1
                parts_with_coords.append(f"{lv2bb['name']}({curr_pos}-{end})")
                
                full_seq = ''.join(seq_parts)
                part_list = "|".join(parts_with_coords)
                output.append({"name": c["Construct Name"] or f"Lv2_Construct_{constructs.index(c)+1}", "type": c["Lv2BB"], "sequence": full_seq, "no_tag_sequence": "", "parts_used": part_list})

            out_df = pd.DataFrame(output)
            full_df = pd.concat([df, out_df[df.columns]], ignore_index=True)
            st.dataframe(out_df)
            csv_out = StringIO()
            full_df.to_csv(csv_out, index=False)
            st.download_button("Download Updated CSV", csv_out.getvalue(), "assembled_level2.csv", "text/csv")

            export_txt = StringIO()
            for idx, row in out_df.iterrows():
                export_txt.write(f">{row['name']}\n{row['sequence']}\n")
            st.download_button("Download FASTA-style File", export_txt.getvalue(), "assembled_level2_sequences.txt", "text/plain")

            for idx, row in out_df.iterrows():
                st.subheader(f"Map: {row['name']}")
                fig, ax = plt.subplots(figsize=(8, 0.5))
                start = 0
                for part in row['parts_used'].split('|'):
                    name, coords = part.split('(')
                    end = int(coords.strip(')').split('-')[1])
                    width = end - start
                    patch = mpatches.Rectangle((start, 0), width, 0.5, label=name, edgecolor='black')
                    ax.add_patch(patch)
                    ax.text(start + width/2, 0.2, name, ha='center', va='center', fontsize=8)
                    start = end
                ax.set_xlim(0, start)
                ax.set_ylim(0, 1)
                ax.axis('off')
                st.pyplot(fig)

    elif page == "CSV editing and FSTA download":
        st.header("Edit Sequences & Batch Download FASTA")
        st.subheader("Edit Sequences")
        edited_df = st.data_editor(df, num_rows="dynamic", use_container_width=True)

        # Multi-select for download
        st.subheader("Select Sequences to Download as FASTA")
        selected_names = st.multiselect("Choose entries to export", edited_df["name"].tolist())

        if selected_names:
            fasta_output = StringIO()
            for _, row in edited_df[edited_df["name"].isin(selected_names)].iterrows():
                fasta_output.write(f">{row['name']}\n{row['sequence']}\n")

            st.download_button("Download FASTA", fasta_output.getvalue(), "selected_sequences.fasta", "text/plain")
        csv_output = edited_df.to_csv(index=False)
        st.download_button(" Download Edited CSV", csv_output, "edited_sequences.csv", "text/csv")
    
    if page == "Input New Sequence":
        st.header("Add New Sequence")
        new_name = st.text_input("Part Name")
        new_seq = st.text_area("DNA Sequence (no stop codon)")
        selected_type = st.selectbox("Part Type", ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"])

        if st.button("Add Overhangs"):
            overhangs = {
                "Promoter": ("TCGGTCTCAAACG", "GGAAGATGAGACC"),
                "N-terminal tag": ("TCGGTCTCAAAGATG", "TCCGGTATGTGAGACC"),
                "CDS": ("TCGGTCTCATATG", "GGATCCTGAGACC"),
                "C-terminal tag": ("TCGGTCTCAATCCTCAGGT", "TCCGGTGGCTGAGACC"),
                "Terminator": ("TCGGTCTCATGGCTAATGA", "GGGCTGTGAGACC")
            }
            no_tag_overhangs = {
                "Promoter": ("TCGGTCTCAAACG", "TCCGGTATGTGAGACC"),
                "Terminator": ("TCGGTCTCAATCCTAATGA", "GGGCTGTGAGACC")
            }

            left_oh, right_oh = overhangs.get(selected_type, ("", ""))
            if new_seq.endswith("TAA") or new_seq.endswith("TAG") or new_seq.endswith("TGA"):
                new_seq = new_seq[:-3]  # remove stop codon
            if new_seq.startswith("ATG"):
                new_seq = new_seq[3:]
            final_seq = left_oh + new_seq + right_oh

            if selected_type in ["Promoter", "Terminator"]:
                final_notag_seq = no_tag_overhangs[selected_type][0] + new_seq + no_tag_overhangs[selected_type][1]
            else:
                final_notag_seq = ""
            
            new_row = pd.DataFrame([{
                    "name": new_name,
                    "type": selected_type,
                    "sequence": final_seq,
                    "no_tag_sequence": final_notag_seq,
                    "parts_used": ""
                }])
            df = pd.concat([df, new_row], ignore_index=True)
            csv = df.to_csv(index=False)
            st.success("New part added with overhangs")
            st.download_button("Download Updated CSV", csv, "updated_parts.csv", "text/csv")

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

if page == "Generate new CDS":
    st.title("🔗 UniProt Isoform to CCDS Link")
    uniprot_id = st.text_input("Enter UniProt Isoform ID (e.g., P68871-1)")

    if st.button("Find DNA and protein sequence"):
        with st.spinner("Fetching CCDS link..."):
            ccds_link, error = get_ccds_link(uniprot_id)
            if error:
                st.error(error)
            else:
                st.success("CCDS link found!")
                st.markdown(f"[Open CCDS Entry for {uniprot_id}]({ccds_link})")


    

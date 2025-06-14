import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import StringIO
import re
import requests
import random
from datetime import datetime

#from xml.etree import ElementTree
import re
from functions import restriction_sites, amino_acid_to_codon, codon_to_aa, overhangs, no_tag_overhangs, BsaI_dg, BsmBI_dg, BsmBI_dg_rev, get_ccds_link, translate_and_find_stops, rm_stop_codon,rm_start_codon, optimize_dna_sequence, count_re_sites, highlight_dna_and_protein

st.set_page_config(page_title="MoClo-YTK Assembly Tool", layout="wide")
st.title("MoClo Cloning Platform")
st.markdown("""
        <style>
        .aligned-paragraph {
            width: 100%;
            font-family: 'Lucida Handwriting', monospace;
            font-size: 16px;
            line-height: 1.6;
            text-align: justify;
            text-indent: 0; /* No indent on first line */
        }
        .aligned-paragraph::first-line {
            text-align: left;
        }
        .aligned-paragraph span {
            display: block;
            text-align: right;
        }
        </style>

        <p class="aligned-paragraph">
            MoClo is hard to learn, but easy to build plasmids at the bench.
            <span>---- Prof. Chen Davidovich</span>
        </p>
        """, unsafe_allow_html=True)

    

HELPER_RULES = [
    ("pL6", ["Pos1", "Pos2", "Pos3", "Pos4", "Pos5"]),
    ("pL56", ["Pos1", "Pos2", "Pos3", "Pos4"]),
    ("pL456", ["Pos1", "Pos2", "Pos3"]),
    ("pL3456", ["Pos1", "Pos2"]),
]


page = st.sidebar.selectbox("Select function", [
    "Level 1 Assembly", 
    "Level 2 Assembly", 
    "Generate level 0 parts",
    "Tables editor",
    
    ])

uploaded_file = st.file_uploader("Upload CSV (name, type, sequence, no_tag_sequence, parts_used)", type="csv")
try:
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        if "parts_used" not in df.columns:
            df["parts_used"] = ""
    else:
        df = pd.read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ_XIjOkjMVCARB0vxfz4fSQoBF5IdtRBtbKCK9IJ4vTYEgAk-RshAwvJ3uU_q4zM9cxvm-GKSDYx_Q/pub?output=csv")
        st.header (f"Example file loaded, please upload your data to use it")
except Exception as e:
    st.error(f"Failed to load example file: {e}")
    st.stop()


if page == "Generate level 0 parts":
    
        
    Subpage = st.sidebar.selectbox("Select sequence type", [
        "Coding (ORFs, Tags)",  "Non-coding (Promoters, Terminators)",
    ])
    if Subpage == "Coding (ORFs, Tags)":
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.subheader ("UniProt isoform DNA squence finder")
        st.sidebar.markdown("UniProt isoform DNA squence finder will help you find the DNA sequence of a protein using its uniprot ID. Make sure you enter the ID including the isoform number.")
        st.header("UniProt isoform DNA squence finder")
        uniprot_id = st.text_input("Enter UniProt ISOFORM ID (e.g., Q92833-1 not Q92833)")

        if st.button("Find DNA and protein sequence"):
            with st.spinner("Fetching CCDS link..."):
                ccds_link, error = get_ccds_link(uniprot_id)
                if error:
                    st.error(error)
                else:
                    st.success("CCDS link found! Please double check!")
                    st.markdown(f"[Open CCDS entry for {uniprot_id}]({ccds_link})")
            
        st.header("Restriction sites remover")
        st.sidebar.subheader ("Restriction sites remover")
        st.sidebar.markdown("Restriction sites remover tool will remove stop codon at the end, BsaI, BsmBI, NotI and PmeI retriction sites by silent mutations." \
        "At the end of this process, you have an option to add the sequence directly to the CSV file, but only the BsmBI and BsaI sites will be added. " \
        "Backbone sequences are not added to reduce file size."
        )
        
        dna_input = st.text_area("Input DNA Sequence (ATGC only)", height=200).strip().upper()
        dna_input = re.sub(r'\s+', '', dna_input.upper())
        dna_input = rm_stop_codon(dna_input)

        if dna_input[:3] not in ['ATG']:
            st.text(f"No ATG found!")

        if dna_input:
            site_counts = count_re_sites(dna_input)
            st.subheader("Restriction Site Count in Original Sequence")
            for enzyme, count in site_counts.items():
                st.write(f"🔹 {enzyme}: {count} site(s)")
            
            protein, stops = translate_and_find_stops(dna_input)
            if len(stops) > 1:
                st.error(f"Internal stop codon found at position {stops}")
                st.text(protein)

        if st.button("Show restriction sites"):
                result = highlight_dna_and_protein(dna_input)
                st.markdown(result, unsafe_allow_html=True)
        if st.button("Visualise optimised DNA"):
            if not dna_input or not re.fullmatch(r'[ATGC]+', dna_input):
                st.error("Please enter a valid DNA sequence containing only A, T, G, C.")
            else:
                st.text(f"Original sequence length: {len(dna_input)} bp")
                dna_cleaned = rm_stop_codon(dna_input)
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
            dna_cleaned = rm_stop_codon(dna_input)
            dna_optimized = optimize_dna_sequence(dna_cleaned)
            dna_optimized = rm_start_codon(dna_optimized)

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
            st.download_button("Download updated CSV of Level 0 parts", csv, "updated_parts.csv", "text/csv")
        
    elif Subpage == "Non-coding (Promoters, Terminators)":
        st.header("Add new sequence")
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown("This tool can be used to add parts directly. It will check the restriction site but wont try to remove it. " \
        "For promoter and terminator, it will also add no tag sequence. " \
        "At the end of this process, you have an option to add the sequence directly to the CSV file, but only the BsmBI and BsaI sites will be added. " \
        "Backbone sequences are not added to reduce file size.")
        new_name = st.text_input("Part Name")
        new_seq = st.text_area("DNA sequence")
        if new_seq:
            site_counts = count_re_sites(new_seq)
            st.subheader("Restriction sites count")
            for enzyme, count in site_counts.items():
                st.write(f"🔹 {enzyme}: {count} site(s)")
                
        selected_type = st.selectbox("Part Type", ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"])

        if st.button("Add Overhangs"):


            left_oh, right_oh = overhangs.get(selected_type, ("", ""))
            if selected_type in ["N-terminal tag", "CDS", "C-terminal tag"]:
                    st.text(f"Only stop codons will be removed")
                    new_seq = rm_start_codon(new_seq) #added back in OH
                    new_seq = rm_stop_codon(new_seq)
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
            st.download_button("Download updated CSV of level 0 parts", csv, "updated_parts.csv", "text/csv")
   
if page == "Level 1 Assembly": 
        st.header("Level 1 Assembly")
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool will help you contruct many level 1 plasmids all at once. " \
        "It will atomatically select the no_tag sequence for promoter and terminator when no N- or C-terminal tags are added."
        "At the end of this process the level 1 plasmid will be digested by BsmBI and stored in the CSV file to reduce file size.")
        num_constructs = st.number_input("Number of constructs", min_value=1, max_value=10)
        types = ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"]
        constructs = []

        for i in range(num_constructs):
            st.subheader(f"Construct {i+1}")
            construct = {}
            construct['name'] = st.text_input(f"Construct name {i+1}", key=f"name_{i}")
            part_names = []
            seq_parts = []
            start = 1
            parts_used = []

            cols = st.columns([1,2, 2, 2,1,1])
            
            sel_pro = cols[0].selectbox(f"Promoter [{i+1}]", df[df['type'] == "Promoter"]['name'].tolist(), key=f"Promoter_{i}")
            sel_N_tag = cols[1].selectbox(f"N-terminal tag [{i+1}]", ["None"] + df[df['type'] == "N-terminal tag"]['name'].tolist(), key=f"Ntag_{i}")
            sel_CDS = cols[2].selectbox(f"CDS [{i+1}]", df[df['type'] == "CDS"]['name'].tolist(), key=f"CDS_{i}")
            sel_C_tag = cols[3].selectbox(f"C-terminal tag [{i+1}]", ["None"] + df[df['type'] == "C-terminal tag"]['name'].tolist(), key=f"Ctag_{i}")
            sel_ter = cols[4].selectbox(f"Terminator [{i+1}]", df[df['type'] == "Terminator"]['name'].tolist(), key=f"Terminator_{i}")
            sel_BB = cols[5].selectbox(f"Backbone [{i+1}]", df[df['type'] == 'Backbone']['name'].tolist(), key=f"bb_{i}")

            for t in types:
                if t == "Promoter":
                    row = df[df['name'] == sel_pro].iloc[0]
                    if sel_N_tag != "None":
                        seq = BsaI_dg(row['sequence'])
                        name = row['name']
                    else:
                        seq = BsaI_dg(row['no_tag_sequence'])
                        name = row['name'] + "_no_tag"

                elif t == "N-terminal tag":
                    if sel_N_tag == "None":
                        continue
                    row = df[df['name'] == sel_N_tag].iloc[0]
                    seq = BsaI_dg(row['sequence'])
                    name = sel_N_tag

                elif t == "CDS":
                    row = df[df['name'] == sel_CDS].iloc[0]
                    seq = BsaI_dg(row['sequence'])
                    name = sel_CDS

                elif t == "C-terminal tag":
                    if sel_C_tag == "None":
                        continue
                    row = df[df['name'] == sel_C_tag].iloc[0]
                    seq = BsaI_dg(row['sequence'])
                    name = sel_C_tag

                elif t == "Terminator":
                    row = df[df['name'] == sel_ter].iloc[0]
                    if sel_C_tag != "None":
                        seq = BsaI_dg(row['sequence'])
                        name = row['name']
                    else:
                        seq = BsaI_dg(row['no_tag_sequence'])
                        name = row['name'] + "_no_tag"

                else:
                    continue

                end = start + len(seq) - 1
                parts_used.append(f"{name}({start}-{end})")
                seq_parts.append(seq)
                part_names.append(name)
                start = end + 1

            backbone = sel_BB
            row = df[df['name'] == backbone].iloc[0]
            seq = BsaI_dg(row['sequence'])
            end = start + len(seq) - 1
            parts_used.append(f"{backbone}({start}-{end})")
            seq_parts.append(seq)

            construct['sequence'] = ''.join(seq_parts)
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
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool will help you contruct level 2 plasmids and show linear maps of them. When finding level 1 plasmids you can use the simpified map to check if they are the one you want.")
        pos_cols = [f"Pos{i}" for i in range(1, 7)]
        lv2bb_options = df[df["type"].str.contains("Lv2BB")]["name"].tolist()
        filtered_df = df[df["type"].str.contains("Pos|Lv2BB") | df["name"].isin([h for h, _ in HELPER_RULES])]
        num_lv2 = st.number_input("Number of constructs", min_value=1, max_value=10)
        constructs = []

        for i in range(num_lv2):
            st.subheader(f"Lv2 Construct {i+1}")
            c = {"Construct Name": st.text_input(f"Name [{i+1}]", key=f"lv2_name_{i}")}
            for pos in pos_cols:
                options = ["None"] + filtered_df[filtered_df["type"].str.contains(pos)]["name"].tolist()
                c[pos] = st.selectbox(f"{pos} [{i+1}]", options, key=f"{pos}_{i}")
                if c[pos] != "None":
                    row = filtered_df[filtered_df["name"] == c[pos]]
                    parts_used_info = row.iloc[0]["parts_used"]
                    re_parts_used_info = re.sub(r'\([^)]*\)', '', parts_used_info)
                    st.text(f"{pos} construct linear map: {re_parts_used_info}")
            c["Lv2BB"] = st.selectbox(f"Lv2BB [{i+1}]", lv2bb_options, key=f"Lv2BB_{i}")
            constructs.append(c)

        if st.button("Assemble Lv2 constructs"):
            output = []
            for c in constructs:
                seq_parts = []
                parts_with_coords = []
                curr_pos = 1
                for pos in pos_cols:
                    if c[pos] != "None":
                        row = df[df["name"] == c[pos]].iloc[0]
                        seq = BsmBI_dg(row["sequence"])
                        seq_parts.append(seq)
                        end = curr_pos + len(seq) - 1
                        parts_with_coords.append(f"{c[pos]}({curr_pos}-{end})")
                        curr_pos = end + 1

                for helper, required in HELPER_RULES:
                    if all(c[p] != "None" for p in required):
                        row = df[df["name"] == helper]
                        if not row.empty:
                            row = row.iloc[0]
                            seq = BsmBI_dg(row["sequence"])
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
            st.dataframe(out_df.drop(columns=["no_tag_sequence"], errors="ignore"))
            csv_out = StringIO()
            full_df.to_csv(csv_out, index=False)
            st.download_button("Download updated CSV", csv_out.getvalue(), "assembled_level2.csv", "text/csv")

            export_txt = StringIO()
            for idx, row in out_df.iterrows():
                export_txt.write(f">{row['name']}\n{row['sequence']}\n")
            st.download_button("Download FASTA file of Level 2 plasmids", export_txt.getvalue(), "assembled_level2_sequences.fasta", "text/plain")

            
            color_map = {}

            def get_random_color():
                return "#" + ''.join(random.choices('0123456789ABCDEF', k=6))

            for idx, row in out_df.iterrows(): 
                st.subheader(f"Map: {row['name']}")
                fig, ax = plt.subplots(figsize=(8, 0.5))
                start = 0
                for part in row['parts_used'].split('|'):
                    name, coords = part.split('(')
                    end = int(coords.strip(')').split('-')[1])
                    width = end - start
                    
                    
                    if name not in color_map:
                        color_map[name] = get_random_color()
                    color = color_map[name]
                    patch = mpatches.Rectangle((start, 0), width, 0.5, color=color, edgecolor='black')
                    ax.add_patch(patch)
                    ax.text(start + width / 2, 0.6, name, ha='center', va='bottom', fontsize=8, rotation=45)

                    start = end

                ax.set_xlim(0, start)
                ax.set_ylim(0, 1.5)
                ax.axis('off')
                st.pyplot(fig)


elif page == "Tables editor":
        st.title("Edit Sequences & Batch Download FASTA")
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool allow you to edit the CSV file directly and export selcted parts as FSTA files.")
        st.sidebar.title (":warning: Warning! :warning:")
        st.sidebar.markdown ("There is no fail safe mechanism build in for this tool. " \
        "Please be caucious when you edit sequence. Level 0 and level 1 plasmids are stored as digested form to reduce file size. " \
        "The backbone sequences will be added upon FSTA download." \
        "The plasmids in the demo files are complete plasmid sequecne from our paper." \
        "The demo plasmid sequences are avalible from addgene. Do no try to download the sequence here! ")
        st.subheader("Edit Sequences")
        edited_df = st.data_editor(df, num_rows="dynamic", use_container_width=True)

        
        st.subheader("Select Sequences to Download as FASTA")
        selected_names = st.multiselect("Choose entries to export", edited_df["name"].tolist())

        if selected_names:
            fasta_output = StringIO()
            for _, row in edited_df[edited_df["name"].isin(selected_names)].iterrows():
                fasta_output.write(f">{row['name']}\n{row['sequence']}\n")

            st.download_button("Download selected sequences as FASTA", fasta_output.getvalue(), "selected_sequences.fasta", "text/plain")
        csv_output = edited_df.to_csv(index=False)
        st.download_button("Download Edited CSV", csv_output, "edited_sequences.csv", "text/csv")


bottom_sidebar = st.sidebar.empty()
bottom_sidebar.markdown("""---  
**Useful Links**  
- [Check the code on Github](https://github.com/ZhihaoLai-MoClo/MoClo_Baculo_toolkit/)  
- [Cite our MoClo Baculo paper](https://www.sciencedirect.com/science/article/pii/S0022283625000099)  
- [Get MoClo Baculo toolkit from Addgene](https://www.addgene.org/kits/davidovich-moclo-baculo/)  
- [Get protocol](https://www.davidovich-lab.com/moclobaculo)  
- [Davidovich Lab](https://www.davidovich-lab.com/)  
""")

api_url = f"https://api.github.com/repos/ZhihaoLai-MoClo/MoClo_Baculo_toolkit/commits"
try:
    response = requests.get(api_url)
    response.raise_for_status()
    commits = response.json()
    latest_commit_date = commits[0]["commit"]["committer"]["date"]
    # Format nicely
    date_obj = datetime.strptime(latest_commit_date, "%Y-%m-%dT%H:%M:%SZ")
    st.sidebar.markdown("Made by Zhihao Lai")
    st.sidebar.markdown(f"**Last updated on GitHub:** {date_obj.strftime('%Y-%m-%d %H:%M')} UTC")
except Exception as e:
    st.sidebar.warning("Could not fetch update date on GitHub.")
    st.sidebar.text(f"Error: {e}")

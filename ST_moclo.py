import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import StringIO
import re
import requests
import random
from datetime import datetime
import re
from functions import BsaI_dg, BsmBI_dg, get_ccds_link, translate_and_find_stops, rm_stop_codon,rm_start_codon, optimize_dna_sequence, count_re_sites, highlight_dna_and_protein
from constants import OVERHANGS, NO_TAG_OVERHANGS, HELPER_RULES
from PIL import Image



st.set_page_config(page_title="MoClo Assembly Platform", page_icon="🧬", layout="wide")
st.title("MoClo Assembly Platform")
st.markdown("""
    <style>
    .quote-block {
        width: 100%;
        font-family: 'Times New Roman';
        font-size: 24px;
        line-height: 1.8;
        margin-bottom: 1em;
    }
    .quote-attribution {
        text-align: right;
        font-size: 20px;
        font-style: italic;
    }
    </style>

    <div class="quote-block">
        MoClo is hard to learn, but easy to build plasmids at the bench.
        <div class="quote-attribution">— Prof. Chen Davidovich</div>
    </div>
""", unsafe_allow_html=True)
    
page = st.sidebar.selectbox("Select function", [
    "Main page",
    "Generate level 0 parts",
    "Level 1 Assembly", 
    "Level 2 Assembly", 
    "Tables editor",
    ])

if page == "Main page":

    st.subheader ("What is MoClo Baculo toolkit")
    st.text("The modular cloning (MoClo) Baculo toolkit is designed for constructing multigene expression vectors for the baculovirus system and, through compatibility with the MoClo Yeast toolkit, also for yeast. Vector construction by MoClo Baculo utilises Golden Gate assembly, which does not require PCR, primers or the sequencing of intermediate products. For users of the baculovirus system, MoClo Baculo is compatible with the Bac-to-Bac and biGBac systems.")
    st.subheader("Where can I find the protocol")
    st.markdown("A step-by-step protocol can be found on [Davidovich Lab website](https://www.davidovich-lab.com/moclobaculo).")
    st.subheader("How does Golden Gate Assembly work")
    st.markdown("There is detailed introduction on [Snapgene](https://www.snapgene.com/guides/golden-gate-assembly)")
    st.subheader("Tutorial video")
    st.text("If you are only interested in learning how this webtool works, while happy to skip all the general introduction to MoClo Baculo, then skip to minute 24 in the video.")
    st.video("https://youtu.be/g7LITscOQyM")

    
    api_url = f"https://api.github.com/repos/ZhihaoLai-MoClo/MoClo_Baculo_toolkit/commits"
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        commits = response.json()
        latest_commit_date = commits[0]["commit"]["committer"]["date"]
        date_obj = datetime.strptime(latest_commit_date, "%Y-%m-%dT%H:%M:%SZ")
        
        st.sidebar.markdown(f"**Tool last updated on GitHub:** {date_obj.strftime('%Y-%m-%d %H:%M')} UTC")
    except Exception as e:
        st.sidebar.warning("Could not fetch update date on GitHub.")
        st.sidebar.text(f"Error: {e}")
    st.subheader("Start by uploading your CSV file or try our example file.")

uploaded_file = st.file_uploader("Upload CSV with columns (name, type, sequence, no_tag_sequence, parts_used)", type="csv")
try:
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        if "parts_used" not in df.columns:
            df["parts_used"] = ""
        required_columns = {'name', 'sequence', 'type'}
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            st.error(f"Missing required columns: {', '.join(missing_columns)}")
    else:
        df = pd.read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ_XIjOkjMVCARB0vxfz4fSQoBF5IdtRBtbKCK9IJ4vTYEgAk-RshAwvJ3uU_q4zM9cxvm-GKSDYx_Q/pub?output=csv")
        st.subheader (f"Example file loaded")
except Exception as e:
    st.error(f"Failed to load example file: {e}")
    st.stop()

if "cached_df" not in st.session_state:
    st.session_state.cached_df = df.copy()
if st.button("Reset cache to the uploaded file"): 
    st.session_state.cached_df = df.copy()



elif page == "Generate level 0 parts":
    lv0bb_row = df[df["name"] == "pM.000"].iloc[0]
    lv0bb_seq = BsmBI_dg(lv0bb_row['sequence'])
        
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
        "At the end of this process, you have an option to add the lv0 plasmid directly to the CSV file."
        )
        
        dna_input = st.text_area("Input DNA Sequence (ATGC only)", height=200).strip().upper()
        dna_input = re.sub(r'\s+', '', dna_input.upper())
        if not dna_input or not re.fullmatch(r'[ATGC]+', dna_input):
                st.error("Please enter a valid DNA sequence containing only A, T, G, C.")
        dna_input = rm_stop_codon(dna_input)

        if dna_input[:3] not in ['ATG']:
            st.error(f"No ATG found!")

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

        if new_name in df["name"].values:
            st.warning("Name already exists. Please choose a unique name.")
        elif new_name in set(st.session_state.cached_df["name"].values):
            st.warning("Name already exists as cached plasmids")

        selected_type = st.selectbox("Part Type", ["N-terminal tag", "CDS", "C-terminal tag"])
        if st.button("Add optimised sequence to CSV file"):
            dna_cleaned = rm_stop_codon(dna_input)
            dna_optimized = optimize_dna_sequence(dna_cleaned)
            dna_optimized = rm_start_codon(dna_optimized)

            left_oh, right_oh = OVERHANGS.get(selected_type, ("", ""))
            final_seq = left_oh + dna_optimized + right_oh
            lv0_seq = []
            lv0_seq.append(BsmBI_dg(final_seq))
            lv0_seq.append(lv0bb_seq)

            new_row = pd.DataFrame([{
                        "name": new_name,
                        "type": selected_type,
                        "sequence": ''.join(lv0_seq),
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
        st.sidebar.markdown("This tool can be used to add parts directly. It will check the restriction sites but will not remove them. " \
        "For promoter and terminator, it will also add no tag sequence. " \
        "At the end of this process, you have an option to add the sequence directly to the CSV file.")
        new_name = st.text_input("Part Name")

        if new_name in df["name"].values:
            st.warning("Name already exists. Please choose a unique name.")
        elif new_name in set(st.session_state.cached_df["name"].values):
            st.warning("Name already exists as cached plasmids")

        new_seq = st.text_area("DNA sequence")
        if new_seq:
            site_counts = count_re_sites(new_seq)
            st.subheader("Restriction sites count")
            for enzyme, count in site_counts.items():
                st.write(f"🔹 {enzyme}: {count} site(s)")
                
        selected_type = st.selectbox("Part Type", ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"])

        if st.button("Add Overhangs"):
            left_oh, right_oh = OVERHANGS.get(selected_type, ("", ""))
            if selected_type in ["N-terminal tag", "CDS", "C-terminal tag"]:
                    st.text(f"Only stop codons will be removed")
                    new_seq = rm_start_codon(new_seq) #added back in OH
                    new_seq = rm_stop_codon(new_seq)
            final_seq = left_oh + new_seq + right_oh
            
            lv0_seq = []
            lv0_notag_seq = []
            lv0_seq.append(BsmBI_dg(final_seq))
            lv0_seq.append(lv0bb_seq)

            if selected_type in ["Promoter", "Terminator"]:
                final_notag_seq = NO_TAG_OVERHANGS[selected_type][0] + new_seq + NO_TAG_OVERHANGS[selected_type][1]
                lv0_notag_seq.append(BsmBI_dg(final_notag_seq))
                lv0_notag_seq.append(lv0bb_seq)
            else:
                lv0_notag_seq = ""
        
            new_row = pd.DataFrame([{
                        "name": new_name,
                        "type": selected_type,
                        "sequence": ''.join(lv0_seq),
                        "no_tag_sequence": ''.join(lv0_notag_seq),
                        "parts_used": ""
                    }])
            df = pd.concat([df, new_row], ignore_index=True)
            csv = df.to_csv(index=False)
            st.success("New part added with overhangs")
            st.download_button("Download updated CSV of level 0 parts", csv, "updated_parts.csv", "text/csv")
   
elif page == "Level 1 Assembly": 
        st.header("Level 1 Assembly")
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool will help you contruct many level 1 plasmids all at once. " \
        "It will automatically select the right promoter and terminator when no N- or C-terminal tags are added."
        )
        num_constructs = st.number_input("Number of constructs", min_value=1, max_value=10)
        types = ["Promoter", "N-terminal tag", "CDS", "C-terminal tag", "Terminator"]
        constructs = []

        for i in range(num_constructs):
            st.subheader(f"Construct {i+1}")
            construct = {}
            new_name=st.text_input(f"Construct name {i+1}", key=f"name_{i}")
            if new_name in df["name"].values:
                st.warning("Name already exists. Please choose a unique name.")
            elif new_name in set(st.session_state.cached_df["name"].values):
                st.warning("Name already exists as cached plasmids")
            construct['name'] = new_name
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
            sel_BB = cols[5].selectbox(f"Backbone [{i+1}]", df[df['type'] == "Lv1BB"]['name'].tolist(), key=f"bb_{i}")

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

            row = df[df['name'] == sel_BB].iloc[0]
            bb_seq = BsaI_dg(row['sequence'])
            end = start + len(bb_seq) - 1
            parts_used.append(f"{sel_BB}({start}-{end})")
            seq_parts.append(bb_seq)

            construct['sequence'] = ''.join(seq_parts)
            construct['parts_used'] = '|'.join(parts_used)
            construct['type'] = sel_BB
            construct['no_tag_sequence'] = ""
            constructs.append(construct)
        
        cache_lv1 = st.checkbox("Cache your level 1 plasmids for level 2 assembly.")

        if st.button("Assemble Constructs"):
            new_df = pd.DataFrame(constructs)
            st.dataframe(new_df)
            all_df = pd.concat([df, new_df[df.columns]], ignore_index=True)
            if cache_lv1:
                st.session_state.cached_df = pd.concat([st.session_state.cached_df, new_df], ignore_index=True)
                st.success("Constructs added to cache.")
            csv = all_df.to_csv(index=False)
            st.download_button("Download Updated CSV", csv, "assembled_level1_constructs.csv", "text/csv")
        
        

elif page == "Level 2 Assembly":
        st.header("Level 2 Assembly")
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool will help you contruct level 2 plasmids and show linear maps of them. When finding level 1 plasmids you can use the simpified map to check if they are the one you want.")
        read_cache = st.checkbox("Read plasmids from cache for level 2 assembly")
        if read_cache:
            df = st.session_state.cached_df
        else:
            df=df
        
        pos_cols = [f"Pos{i}" for i in range(1, 7)]
        lv2bb_options = df[df["type"].str.contains("Lv2BB")]["name"].tolist()
        filtered_df = df[df["type"].str.contains("Pos|Lv2BB") | df["name"].isin([h for h, _ in HELPER_RULES])]
        num_lv2 = st.number_input("Number of constructs", min_value=1, max_value=10)
        constructs = []

        for i in range(num_lv2):
            st.subheader(f"Lv2 Construct {i+1}")
            new_name=st.text_input(f"Name [{i+1}]", key=f"lv2_name_{i}")
            if new_name in df["name"].values:
                st.warning("Name already exists. Please choose a unique name.")
            elif new_name in set(st.session_state.cached_df["name"].values):
                st.warning("Name already exists as cached plasmids")
            c = {"Construct Name": new_name}
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
            st.download_button("Download FASTA file of Level 2 plasmids", export_txt.getvalue(), "assembled_level2_constructs.fasta", "text/plain")

            
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
        read_cache = st.checkbox("show and edit cached plasmids")
        if read_cache:
            df = st.session_state.cached_df
        else:
            df=df
        st.sidebar.title (":bulb: Tips! :bulb:")
        st.sidebar.markdown ("This tool allow you to edit the CSV file directly and export selcted parts as FSTA files.")
        st.sidebar.title (":warning: Warning! :warning:")
        st.sidebar.markdown ("There is no fail safe mechanism build in for this tool. " \
        "Please be caucious when you edit sequence.")
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
st.sidebar.markdown("Made by Zhihao Lai using Streamlit")

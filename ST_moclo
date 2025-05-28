import streamlit as st
import pandas as pd
from io import StringIO

st.set_page_config(page_title="MoClo Assembly Tool", layout="wide")
st.title("MoClo Cloning Platform")

# Helper rules (most specific first)
HELPER_RULES = [
    (["Pos1", "Pos2", "Pos3", "Pos4", "Pos5"], "pL6"),
    (["Pos1", "Pos2", "Pos3", "Pos4"], "pL56"),
    (["Pos1", "Pos2", "Pos3"], "pL456"),
    (["Pos1", "Pos2"], "pL3456"),
]

REQUIRED_PARTS = ["Promoter", "CDS", "Terminator"]
OPTIONAL_PARTS = ["N-terminal tag", "C-terminal tag"]

page = st.sidebar.selectbox("Select Mode", ["Level 1 Golden Gate Assembly", "Level 2 Golden Gate Assembly"])

uploaded_file = st.file_uploader("Upload CSV with columns: name, type, sequence, no_tag_sequence", type="csv")

if uploaded_file:
    df = pd.read_csv(uploaded_file)
    if not {"name", "type", "sequence", "no_tag_sequence"}.issubset(df.columns):
        st.error("CSV must contain: name, type, sequence, no_tag_sequence")
        st.stop()

    if page == "Level 1 Golden Gate Assembly":
        st.header("Level 1 Golden Gate Assembly")
        part_options = {
            part: df[df["type"].str.lower() == part.lower()]["name"].tolist()
            for part in REQUIRED_PARTS + OPTIONAL_PARTS
        }
        backbone_options = df[df["type"].str.lower() == "backbone"]["name"].tolist()

        num_constructs = st.slider("Number of constructs", 1, 10, 3)
        construct_inputs = []

        for i in range(num_constructs):
            st.markdown(f"### Construct {i+1}")
            row = {}
            row["Construct Name"] = st.text_input(f"Name [{i+1}]", key=f"name_{i}")
            for part in REQUIRED_PARTS:
                row[part] = st.selectbox(f"{part} [{i+1}]", part_options[part], key=f"{part}_{i}")
            for part in OPTIONAL_PARTS:
                row[part] = st.selectbox(f"{part} [{i+1}]", ["None"] + part_options[part], key=f"{part}_{i}")
            row["Backbone"] = st.selectbox(f"Backbone [{i+1}]", backbone_options, key=f"Backbone_{i}")
            construct_inputs.append(row)

        if st.button("Assemble Constructs"):
            output = []
            for i, row in enumerate(construct_inputs):
                seq_parts = []
                used_parts = []

                promoter_row = df[df["name"] == row["Promoter"]].iloc[0]
                seq_parts.append(promoter_row["no_tag_sequence"] if row["N-terminal tag"] == "None" else promoter_row["sequence"])
                used_parts.append(row["Promoter"])

                if row["N-terminal tag"] != "None":
                    ntag_row = df[df["name"] == row["N-terminal tag"]].iloc[0]
                    seq_parts.append(ntag_row["sequence"])
                    used_parts.append(row["N-terminal tag"])

                cds_row = df[df["name"] == row["CDS"]].iloc[0]
                seq_parts.append(cds_row["sequence"])
                used_parts.append(row["CDS"])

                if row["C-terminal tag"] != "None":
                    ctag_row = df[df["name"] == row["C-terminal tag"]].iloc[0]
                    seq_parts.append(ctag_row["sequence"])
                    used_parts.append(row["C-terminal tag"])
                    terminator_row = df[df["name"] == row["Terminator"]].iloc[0]
                    seq_parts.append(terminator_row["sequence"])
                else:
                    terminator_row = df[df["name"] == row["Terminator"]].iloc[0]
                    seq_parts.append(terminator_row["no_tag_sequence"])
                used_parts.append(row["Terminator"])

                backbone_row = df[df["name"] == row["Backbone"]].iloc[0]
                seq_parts.append(backbone_row["sequence"])
                used_parts.append(row["Backbone"])

                full_seq = "".join(seq_parts)
                construct_name = row["Construct Name"] or f"Construct_{i+1}"
                part_list = "|".join(used_parts)

                output.append({
                    "name": construct_name,
                    "type": row["Backbone"],
                    "sequence": full_seq,
                    "no_tag_sequence": "",
                    "parts_used": part_list
                })

            out_df = pd.DataFrame(output)
            full_df = pd.concat([df, out_df[df.columns]], ignore_index=True)
            st.success("Constructs assembled")
            st.dataframe(out_df)
            csv_out = StringIO()
            full_df.to_csv(csv_out, index=False)
            st.download_button("Download Updated CSV", csv_out.getvalue(), "assembled_level1.csv", "text/csv")

    elif page == "Level 2 Golden Gate Assembly":
        st.header("Level 2 Golden Gate Assembly (BsmBI)")
        pos_cols = [f"Pos{i}" for i in range(1, 7)]
        lv2bb_options = df[df["type"].str.contains("Lv2BB")]["name"].tolist()

        filtered_df = df[df["type"].str.contains("Pos|Lv2BB") | df["name"].isin([h for _, h in HELPER_RULES])]

        num_lv2 = st.slider("Number of constructs", 1, 5, 2)
        constructs = []

        for i in range(num_lv2):
            st.markdown(f"### Lv2 Construct {i+1}")
            construct = {"Construct Name": st.text_input(f"Name [{i+1}]", key=f"lv2_name_{i}")}
            for pos in pos_cols:
                options = ["None"] + filtered_df[filtered_df["type"].str.contains(pos)]["name"].tolist()
                construct[pos] = st.selectbox(f"{pos} [{i+1}]", options, key=f"{pos}_{i}")
            construct["Lv2BB"] = st.selectbox(f"Lv2BB [{i+1}]", lv2bb_options, key=f"Lv2BB_{i}")
            constructs.append(construct)

        if st.button("⚡ Assemble Lv2 Constructs"):
            output = []
            for c in constructs:
                seq_parts = []
                used_parts = []

                # Add Pos1–Pos6
                for pos in pos_cols:
                    val = c[pos]
                    if val != "None":
                        row = df[df["name"] == val].iloc[0]
                        seq_parts.append(row["sequence"])
                        used_parts.append(val)

                # Add most specific matching helper plasmid
                for required, helper in HELPER_RULES:
                    if all(c.get(p, "None") != "None" for p in required):
                        row = df[df["name"] == helper]
                        if not row.empty:
                            seq_parts.append(row.iloc[0]["sequence"])
                            used_parts.append(helper)
                        break

                # Add Lv2BB
                lv2bb_row = df[df["name"] == c["Lv2BB"]].iloc[0]
                seq_parts.append(lv2bb_row["sequence"])
                used_parts.append(c["Lv2BB"])

                full_seq = "".join(seq_parts)
                name = c["Construct Name"] or f"Lv2_Construct_{constructs.index(c)+1}"
                part_list = "|".join(used_parts)

                output.append({
                    "name": name,
                    "type": c["Lv2BB"],
                    "sequence": full_seq,
                    "no_tag_sequence": "",
                    "parts_used": part_list
                })

            out_df = pd.DataFrame(output)
            full_df = pd.concat([df, out_df[df.columns]], ignore_index=True)
            st.success("Level 2 Constructs Assembled")
            st.dataframe(out_df)
            csv_out = StringIO()
            full_df.to_csv(csv_out, index=False)
            st.download_button("Download Updated CSV", csv_out.getvalue(), "assembled_level2.csv", "text/csv")

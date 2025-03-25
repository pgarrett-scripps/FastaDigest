import streamlit as st
import pandas as pd
import peptacular as pt
import plotly.express as px
import re

st.set_page_config(page_title="Fasta-Digest", page_icon="🦠", layout="wide")

def main():
    with st.expander("Common Proteases"):
        st.write(pt.PROTEASES)
    with st.sidebar:
        st.title("FASTA Digest 🦠")

        st.caption("This tool digests proteins in a FASTA file using user-defined protease configurations.")

        uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta"])

        df = pd.DataFrame([
            {"regex": "(?=C)", "order": 1, "mc": 0, "complete": True, "semi": False},
            {"regex": "(?=D)", "order": 1, "mc": 0, "complete": True, "semi": False},
            {"regex": "(?<=[KR])", "order": 2, "mc": 0, "complete": True, "semi": False},
        ])

        peptide_regex_filter = st.text_input("Peptide Regex Filter", value="[M]",
                                             help="Filter digested peptides based on regex pattern")

        st.caption("Sequence Digestion Parameters")
        digestion_parameters = st.data_editor(
            df,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "regex": st.column_config.TextColumn(required=True, help="Protease Regex"),
                "order": st.column_config.NumberColumn(help="Order of digestion", min_value=1, step=1, default=1),
                "mc": st.column_config.NumberColumn(help="Number of missed cleavages", min_value=0, step=1, default=0),
                "semi": st.column_config.CheckboxColumn(help="Semi-enzymatic digestion", default=False),
                "complete": st.column_config.CheckboxColumn(help="Complete digestion", default=True)
            }
        )

        c1, c2 = st.columns(2)
        min_len = c1.number_input("Minimum Peptide Length", min_value=1, value=7, step=1)
        max_len = c2.number_input("Maximum Peptide Length", min_value=1, value=30, step=1)
        percentage_mode = st.checkbox("Show histograms as percentage of total peptides")

        if not st.button("Run", type="primary", use_container_width=True):
            st.stop()

    if not uploaded_file:
        st.warning("Please upload a FASTA file.")
        st.stop()

    with st.spinner("Processing..."):
        proteins = pt.parse_fasta(uploaded_file)

        if not proteins:
            st.error("No proteins found in the uploaded FASTA file.")
            st.stop()

        digestion_parameters = digestion_parameters.sort_values(by="order")

        enzyme_configs = {}
        for _, row in digestion_parameters.iterrows():
            order = row['order']

            if order not in enzyme_configs:
                enzyme_configs[order] = pt.EnzymeConfig(
                    regex=[row['regex']],
                    missed_cleavages=row['mc'],
                    semi_enzymatic=row['semi'],
                    complete_digestion=row['complete']
                )
            else:
                # If the order already exists, we need to update the existing EnzymeConfig
                # Append values to the existing EnzymeConfig object
                enzyme_configs[order].regex.append(row['regex'])
                enzyme_configs[order].missed_cleavages += row['mc']

                # Ensure complete digestion is True if any row for the same order is True
                enzyme_configs[order].semi_enzymatic = (
                        enzyme_configs[order].semi_enzymatic or row['semi']
                )

                # Ensure complete digestion is True if any row for the same order is True
                enzyme_configs[order].complete_digestion = (
                        enzyme_configs[order].complete_digestion or row['complete']
                )


        # amke list
        enzyme_configs = [config for order, config in enzyme_configs.items()]

        # show table of enzyme configs

        st.subheader("Enzyme Configurations")
        st.dataframe(pd.DataFrame(enzyme_configs), use_container_width=True)

        all_digested_peptides, all_proteins = [], []
        for header, sequence in proteins:
            sequence_annot = pt.parse(sequence)
            spans = list(pt.sequential_digest(sequence=sequence, enzyme_configs=enzyme_configs, min_len=min_len, max_len=max_len, return_type='span'))
            peptides = [sequence_annot.slice(span[0], span[1]).serialize() for span in spans]

            # num of cleavage sites to include
            cleavage_len = 2

            nterm_cleavages = [(max(span[0]-cleavage_len, 0), span[0], span[2]) for span in spans]
            peptide_nterms = [sequence_annot.slice(span[0], min(len(sequence_annot), span[0]+cleavage_len)).serialize() for span in spans]
            nterm_seqs = [sequence_annot.slice(nterm[0], nterm[1]).serialize() for nterm in nterm_cleavages]
            n_term_cleavage_sts = [f'{nseq}.{nterm}' for nterm, nseq in zip(peptide_nterms, nterm_seqs)]

            cterm_cleavages = [(span[1], min(span[1]+cleavage_len, len(sequence_annot)), span[2]) for span in spans]
            peptide_cterms = [sequence_annot.slice(max(span[1]-cleavage_len, 0), span[1]).serialize() for span in spans]
            cterm_seqs = [sequence_annot.slice(cterm[0], cterm[1]).serialize() for cterm in cterm_cleavages]
            c_term_cleavage_sts = [f'{cterm}.{cseq}' for cterm, cseq in zip(peptide_cterms, cterm_seqs)]

            filtered_peptides = []
            for peptide, n, c in zip(peptides, n_term_cleavage_sts, c_term_cleavage_sts):

                pass_regex_filter = not peptide_regex_filter or bool(re.search(peptide_regex_filter, peptide))
                if pass_regex_filter:
                    filtered_peptides.append(peptide)

                try:
                    neutral_mass = pt.mass(peptide)
                except pt.AmbiguousAminoAcidError:
                    neutral_mass = None

                all_digested_peptides.append(
                    {
                     "peptide": peptide,
                     "n_cleavage_site": n,
                     "c_cleavage_site": c,
                     'is_filtered': pass_regex_filter,
                     "neutral_mass": neutral_mass,
                     "peptide_len": pt.sequence_length(peptide),
                     "KRH_count": peptide.count("K") + peptide.count("R") + peptide.count("H"),
                     "fasta_header": header,
                     })

            protein_cov = pt.percent_coverage(sequence, peptides)
            filtered_cov = pt.percent_coverage(sequence, filtered_peptides)

            all_proteins.append({"fasta_header": header,
                                 'peptides': len(peptides),
                                 'filtered_peptides': len(filtered_peptides),
                                 'coverage': round(protein_cov * 100, 2),
                                 'filtered_coverage': round(filtered_cov * 100, 2),
                                 })

    df_peptides = pd.DataFrame(all_digested_peptides)
    df_peptides['charge_range'] = [f"{max(c, 1)}-{c+2}" for c in df_peptides['KRH_count']]
    df_proteins = pd.DataFrame(all_proteins)

    t1, t2, t3 = st.tabs(["Peptides", "Proteins", "Plots"])

    with t1:
        st.subheader("Peptide Statistics")
        c1, c2 = st.columns(2)
        c1.metric("Total Peptides", len(df_peptides))
        c2.metric("Filtered Peptides", df_peptides['is_filtered'].sum())

        st.dataframe(df_peptides, hide_index=True)
        st.download_button(
            label="Download Peptide Data as CSV",
            data=df_peptides.to_csv(index=False).encode('utf-8'),
            file_name='peptides.csv',
            mime='text/csv',
            on_click='ignore',
            use_container_width=True,
            type="primary"
        )

    with t2:
        st.subheader("Protein Statistics")
        c1, c2, c3 = st.columns(3)
        c1.metric("Total Proteins", len(df_proteins))
        c2.metric("Proteins with Peptides", df_proteins['peptides'].sum())
        c3.metric("Proteins with Filtered Peptides", df_proteins['filtered_peptides'].sum())
        st.dataframe(df_proteins, hide_index=True)
        st.download_button(
            label="Download Protein Data as CSV",
            data=df_proteins.to_csv(index=False).encode('utf-8'),
            file_name='proteins.csv',
            mime='text/csv',
            on_click='ignore',
            use_container_width=True,
            type="primary"
        )

    with t3:
        st.subheader("Peptide Distributions")
        histnorm = "percent" if percentage_mode else None

        # Plot Mass distribution
        fig_mass = px.histogram(df_peptides, x="neutral_mass", title="Peptide Mass Distribution", nbins=50, histnorm=histnorm)
        st.plotly_chart(fig_mass, use_container_width=True)

        # Plot Length distribution
        fig_length = px.histogram(df_peptides, x="peptide_len", title="Peptide Length Distribution", nbins=50,
                                  histnorm=histnorm)
        st.plotly_chart(fig_length, use_container_width=True)

        # Plot KRH count distribution
        fig_krh = px.histogram(df_peptides, x="KRH_count", title="KRH Count Distribution", nbins=20, histnorm=histnorm)
        st.plotly_chart(fig_krh, use_container_width=True)


if __name__ == "__main__":
    main()

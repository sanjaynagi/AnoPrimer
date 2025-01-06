import allel
import gget
import malariagen_data
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def retrieve_data_resource(species):
    assert species in [
        "gambiae_sl",
        "funestus",
    ], f"species {species} not recognised, please use 'gambiae_sl' or 'funestus'"
    if species == "gambiae_sl":
        data_resource = malariagen_data.Ag3(
            url="gs://vo_agam_release_master_us_central1", pre=True
        )
    elif species == "funestus":
        data_resource = malariagen_data.Af1(
            url="gs://vo_afun_release_master_us_central1", pre=True
        )
    return data_resource


def check_my_oligo(
    sequence, sample_sets="3.0", sample_query=None, width=700, height=400
):
    """
    Align a sequence to AgamP3, retrieve ag3 frequencies in this region and plot.
    Only works with An.gambiae_sl for now.
    """

    print("Aligning sequence to AgamP3 genome with BLAT")
    blat_df = gget.blat(sequence=sequence, seqtype="DNA", assembly="anoGam3")
    if blat_df is None:
        print(f"No hit for {sequence}")
        return

    contig, start, end = blat_df.loc[0, ["chromosome", "start", "end"]]
    contig = contig.replace("chr", "")
    region_span = f"{contig}:{start}-{end}"
    print("plotting frequencies in ag3 data")

    fig = plot_sequence_frequencies(
        data_resource=malariagen_data.Ag3(
            url="gs://vo_agam_release_master_us_central1", pre=True
        ),
        region=region_span,
        sample_sets=sample_sets,
        sample_query=sample_query,
        width=width,
        height=height,
    )
    return fig


def plot_sequence_frequencies(
    data_resource, region, sample_sets=None, sample_query=None, width=700, height=400
):
    """Retrieve frequencies"""

    snps = data_resource.snp_calls(
        region=region, sample_sets=sample_sets, sample_query=sample_query
    )
    ref_alt_arr = snps["variant_allele"].compute().values.astype(str)
    freq_arr = (
        allel.GenotypeArray(snps["call_genotype"]).count_alleles().to_frequencies()
    )
    pos = snps["variant_position"].compute().values
    df = pd.DataFrame(
        {
            "position": pos,
            "base": ref_alt_arr[:, 0],
            "alt_frequency": freq_arr[:, 1:].sum(axis=1),
        }
    )  # Make dataframe for plotting
    df["base_pos"] = df["base"] + "_" + df["position"].astype(str)
    # Get the frequency of each base and store as data frame
    freq_df = _get_base_freqs(_addZeroCols(freq_arr), ref_alt_arr).filter(like="freq")

    data = pd.concat([df, freq_df], axis=1)

    fig = _plotly_frequencies(
        data=data,
        region=region,
        sample_sets=sample_sets,
        sample_query=sample_query,
        width=width,
        height=height,
    )
    return fig


def _plotly_frequencies(
    data, region, sample_sets, sample_query=None, width=700, height=400, save=False
):
    import plotly.graph_objects as go

    hover_template = "<br>".join(
        [
            "Base / Position: %{customdata[4]}",
            "Total Alternate freq: %{y}",
            "A_freq: %{customdata[0]}",
            "C_freq: %{customdata[1]}",
            "G_freq: %{customdata[2]}",
            "T_freq: %{customdata[3]}",
        ]
    )
    # Color scatterpoints blue if segregating SNP
    color = [-1 if v == 0 else 1 if v > 0 else 0 for v in data["alt_frequency"]]
    colorscale = [[0, "lightgray"], [0.5, "lightgray"], [1, "dodgerblue"]]

    fig = go.Figure(
        go.Scatter(
            x=data["position"],
            y=data["alt_frequency"],
            customdata=data[["A_freq", "C_freq", "G_freq", "T_freq", "base_pos"]],
            hovertemplate=hover_template,
            mode="markers",
            marker=dict(
                size=14,
                color=color,
                colorscale=colorscale,
                line=dict(width=2, color="black"),
            ),
            marker_symbol="circle",
        )
    )
    # Set xticks to be the REF allele
    fig.update_xaxes(
        tickmode="array",
        tickangle=0,
        tickvals=data["position"].to_list(),
        ticktext=data["base"].to_list(),
    )
    fig.update_yaxes(
        tickmode="array",
        tickvals=np.arange(0, 1, 0.2),
        range=[0, 1],
        title="Alternate allele frequency",
    )
    # Set plot title
    if sample_query is not None:
        title_text = f"{region} | {sample_sets} | {sample_query} | allele frequencies"
    else:
        title_text = f"{region} | {sample_sets} | allele frequencies"

    fig.update_layout(
        height=height,
        width=width,
        title_text=title_text,
        title_x=0.5,
        template="simple_white",
        showlegend=False,
    )
    fig.show()
    return fig


#### utility functions ####


def round_floats_in_df(df, decimal_places=1):
    import numpy as np

    def round_if_float(val):
        if isinstance(val, (int, np.integer)):
            return val
        try:
            float_val = float(val)
            if float_val.is_integer():
                return int(float_val)
            return round(float_val, decimal_places)
        except (ValueError, TypeError):
            return val

    return df.applymap(round_if_float)


def extract_trailing_digits(string):
    import re

    match = re.search(r"\d+$", string)
    if match:
        return match.group(0)
    else:
        return None


def _retrieve_span(primer_df, gdna_pos, oligo, assay_type, pair):
    primer_loc = primer_df.loc[f"primer_{oligo}", str(pair)][0]
    primer_loc = primer_loc + 1 if oligo == "reverse" else primer_loc
    primer_size = primer_df.loc[f"primer_{oligo}", str(pair)][1]

    if any(item in assay_type for item in ["gDNA", "probe"]):
        pos_arr = gdna_pos
    elif assay_type == "cDNA primers":
        pos_arr = np.array([])
        exon_spans = np.array(_consecutive(gdna_pos)) + 1
        for span in exon_spans:
            pos_arr = np.append(pos_arr, np.arange(span[0], span[1] + 1)).astype(int)

    if oligo in ["forward", "probe"]:
        pos = pos_arr[primer_loc : primer_loc + primer_size]
    elif oligo == "reverse":
        pos = np.flip(pos_arr[primer_loc - primer_size : primer_loc])

    return pos


def _get_primer_arrays(
    species, contig, gdna_pos, sample_sets, assay_type, sample_query=None
):
    """
    Load genotype data from Ag3/Af1 resource and return allele frequencies
    at entire input sequence region
    """
    data_resource = retrieve_data_resource(species=species)

    if any(item in assay_type for item in ["gDNA", "probe"]):
        span_str = f"{contig}:{gdna_pos.min()}-{gdna_pos.max()}"
        ds_snps = data_resource.snp_calls(
            region=span_str, sample_sets=sample_sets, sample_query=sample_query
        )  # get genotypes
        ref_alt_arr = ds_snps["variant_allele"].compute().values
        geno = ds_snps["call_genotype"]
        freq_arr = allel.GenotypeArray(geno).count_alleles().to_frequencies()
        pos_arr = gdna_pos
    elif assay_type == "cDNA primers":
        freq_arr = []
        ref_alt_arr = []
        pos_arr = np.array([])
        exon_spans = np.array(_consecutive(gdna_pos)) + 1
        for span in exon_spans:
            span_str = f"{contig}:{span[0]}-{span[1]}"
            ds_snps = data_resource.snp_calls(
                region=span_str, sample_sets=sample_sets, sample_query=sample_query
            )  # get genotypes
            ref_alts = ds_snps["variant_allele"]
            geno = ds_snps["call_genotype"]
            freqs = (
                allel.GenotypeArray(geno).count_alleles().to_frequencies()
            )  # calculate allele frequencies
            freqs = _addZeroCols(freqs)
            freq_arr.append(freqs)
            ref_alt_arr.append(ref_alts)
            pos_arr = np.append(pos_arr, np.arange(span[0], span[1] + 1)).astype(int)
        freq_arr = np.concatenate(freq_arr)
        ref_alt_arr = np.concatenate(ref_alt_arr)

    return (freq_arr, ref_alt_arr.astype("U13"), pos_arr)


def _get_primer_alt_frequencies(
    species, primer_df, gdna_pos, pair, sample_sets, assay_type, contig, sample_query
):
    """
    Find the genomic locations of pairs of primers, and runs span_to_freq
    to get allele frequencies at those locations
    """

    oligos, _ = _return_oligo_list(assay_type)
    base_freqs, ref_alt_arr, pos_arr = _get_primer_arrays(
        species=species,
        contig=contig,
        gdna_pos=gdna_pos,
        sample_sets=sample_sets,
        assay_type=assay_type,
        sample_query=sample_query,
    )

    freq_arr = base_freqs[:, 1:].sum(axis=1)

    di = {}
    for oligo in oligos:
        primer_loc = primer_df.loc[f"primer_{oligo}", str(pair)][0]
        primer_loc = primer_loc + 1 if oligo == "reverse" else primer_loc
        primer_size = primer_df.loc[f"primer_{oligo}", str(pair)][1]
        if oligo in ["forward", "probe"]:
            freq = freq_arr[primer_loc : primer_loc + primer_size]
            base_freqs_arr = base_freqs[primer_loc : primer_loc + primer_size, :]
            ref = ref_alt_arr[primer_loc : primer_loc + primer_size, 0]
            ref_alt = ref_alt_arr[primer_loc : primer_loc + primer_size, :]
            pos = pos_arr[primer_loc : primer_loc + primer_size]
        elif oligo == "reverse":
            freq = np.flip(freq_arr[primer_loc - primer_size : primer_loc])
            base_freqs_arr = base_freqs[primer_loc - primer_size : primer_loc, :]
            base_freqs_arr = np.flip(base_freqs_arr, axis=0)
            ref = ref_alt_arr[primer_loc - primer_size : primer_loc, 0]
            ref = np.array(list(_rev_complement("".join(ref))), dtype=str)
            ref_alt = ref_alt_arr[primer_loc - primer_size : primer_loc, :]
            ref_alt = _complement(np.flip(ref_alt, axis=0))
            pos = np.flip(pos_arr[primer_loc - primer_size : primer_loc])

        df = pd.DataFrame(
            {"position": pos, "base": ref, "alt_frequency": freq}
        )  # Make dataframe for plotting
        df["base_pos"] = df["base"] + "_" + df["position"].astype(str)
        assert df.shape[0] == primer_size, "Wrong size primers"

        freq_df = _get_base_freqs(_addZeroCols(base_freqs_arr), ref_alt).filter(
            like="freq"
        )
        df = pd.concat([df, freq_df], axis=1)
        di[oligo] = df
    return di


def _plotly_primers(
    primer_df,
    res_dict,
    name,
    assay_type,
    sample_sets,
    target,
    out_dir=None,
):
    oligos, _ = _return_oligo_list(assay_type)
    if len(oligos) == 2:
        plt_title = ["Forward primer", "Reverse primer"]
    elif len(oligos) == 3:
        plt_title = ["Forward primer", "Reverse primer", "Probe"]
    elif len(oligos) == 1:
        plt_title = ["Probe"]

    title_list = []
    for pair in primer_df:
        for oligo in plt_title:
            title_list.append(f"{oligo} {pair}")

    hover_template = "<br>".join(
        [
            "Base / Position: %{customdata[4]}",
            "Total Alternate freq: %{y}",
            "A_freq: %{customdata[0]}",
            "C_freq: %{customdata[1]}",
            "G_freq: %{customdata[2]}",
            "T_freq: %{customdata[3]}",
        ]
    )

    fig = make_subplots(
        rows=len(primer_df.columns),
        cols=len(oligos),
        subplot_titles=title_list,
        horizontal_spacing=0.03,
        vertical_spacing=0.08,
    )
    fig.update_annotations(font_size=13)
    for idx, oligo in enumerate(oligos):
        idx = idx + 1
        for pair in primer_df:
            row_i = int(pair)

            color = [
                -1 if v == 0 else 1 if v > 0 else 0
                for v in res_dict[pair][oligo]["alt_frequency"]
            ]
            colorscale = [[0, "lightgray"], [0.5, "lightgray"], [1, "dodgerblue"]]

            # Write text to plot for Tm, GC, span, and 3/5'
            tm = np.round(primer_df.loc[f"primer_{oligo}_tm", pair], 2)
            gc = np.round(primer_df.loc[f"primer_{oligo}_gc_percent", pair], 2)
            span = f"{int(res_dict[pair][oligo]['position'].min())}-{int(res_dict[pair][oligo]['position'].max())}"

            for col in res_dict[pair][oligo].columns:
                if col.endswith("freq"):
                    res_dict[pair][oligo][col] = np.round(res_dict[pair][oligo][col], 2)

            fig.add_trace(
                go.Scatter(
                    x=res_dict[pair][oligo]["base_pos"],
                    y=res_dict[pair][oligo]["alt_frequency"],
                    customdata=res_dict[pair][oligo][
                        ["A_freq", "C_freq", "G_freq", "T_freq", "base_pos"]
                    ],
                    hovertemplate=hover_template,
                    mode="markers",
                    marker=dict(
                        size=14,
                        color=color,
                        colorscale=colorscale,
                        line=dict(width=2, color="black"),
                    ),
                    marker_symbol="circle",
                ),
                row=row_i,
                col=idx,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[pair][oligo]["base_pos"][0],
                y=0.8,
                text="5'",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[pair][oligo]["base_pos"].to_numpy()[-1],
                y=0.8,
                text="3'",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[pair][oligo]["base_pos"].to_numpy()[4],
                y=0.92,
                text=span,
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[pair][oligo]["base_pos"].to_numpy()[-7],
                y=0.92,
                text=f"GC={gc}",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[pair][oligo]["base_pos"].to_numpy()[-3],
                y=0.92,
                text=f"TM={tm}",
                showarrow=False,
            )

            fig.update_xaxes(
                row=row_i,
                col=idx,
                tickmode="array",
                tickvals=res_dict[pair][oligo]["base_pos"],
                ticktext=res_dict[pair][oligo]["base"],
                tickangle=0,
                mirror=True,
            )
            if idx > 1:
                fig.update_yaxes(
                    row=row_i,
                    col=idx,
                    range=[0, 1],
                    tickvals=np.arange(0, 1, 0.2),
                    showticklabels=False,
                    mirror=True,
                )
            else:
                fig.update_yaxes(
                    row=row_i,
                    col=idx,
                    tickvals=np.arange(0, 1, 0.2),
                    range=[0, 1],
                    mirror=True,
                )

            if ((row_i % 2) == 0) & (idx == 1):
                fig.update_yaxes(row=row_i, col=idx, title="Alternate allele frequency")

    if any(item in assay_type for item in ["gDNA"]):
        title_text = f"{name} primer pairs | {sample_sets} | target {target} bp"
    elif assay_type == "probe":
        title_text = f"{name} probe | {sample_sets} | target {target} bp"
    elif assay_type == "cDNA primers":
        title_text = f"{name} primer pairs | {sample_sets} | target {target}"

    # fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)
    fig.update_layout(
        height=220 * len(primer_df.columns),
        width=500 * len(oligos),
        title_text=title_text,
        title_x=0.5,
        template="simple_white",
        showlegend=False,
    )
    if out_dir:
        fig.write_html(f"{name}_{assay_type}_snps.html")
        fig.write_image(f"{name}_{assay_type}_snps.pdf")
    fig.show()


def _get_gDNA_locs(gff, contig, start, end):
    locgff = gff.query(
        f"contig == '{contig}' & type == 'exon' & start < {end} & end > {start}"
    )
    min_ = locgff.start.min() - 100
    max_ = locgff.end.max() + 100
    genegff = gff.query(
        f"contig == '{contig}' & type == 'gene' & start < {end} & end > {start}"
    )
    return (locgff, min_, max_, genegff)


def _get_qPCR_locs(gff, contig, transcript):
    # Load geneset (gff)
    locgff = gff.query(f"Parent == '{transcript}' & type == 'exon'")
    min_ = locgff.start.min() - 200
    max_ = locgff.end.max() + 200
    genegff = gff.query(
        f"contig == '{contig}' & type == 'gene' & start > {min_} & end < {max_}"
    )
    return (locgff, min_, max_, genegff)


def _return_oligo_list(assay_type):
    if assay_type == "probe":
        oligos = ["probe"]
        row_start = 9
    elif any(item == assay_type for item in ["gDNA primers", "cDNA primers"]):
        oligos = ["forward", "reverse"]
        row_start = 11
    elif assay_type == "gDNA primers + probe":
        oligos = ["forward", "reverse", "probe"]
        row_start = 12
    return (oligos, row_start)


def _convert_results_dict_naming(primer_dict):
    k = {}
    for key in primer_dict.keys():
        if "LEFT" in key:
            nkey = key.replace("LEFT", "forward")
        elif "RIGHT" in key:
            nkey = key.replace("RIGHT", "reverse")
        elif "INTERNAL" in key:
            nkey = key.replace("INTERNAL", "probe")
        else:
            nkey = key
        k[nkey.lower()] = primer_dict[key]
    return k


def _complement(x):
    if x == "A":
        return "T"
    elif x == "C":
        return "G"
    elif x == "G":
        return "C"
    elif x == "T":
        return "A"


_complement = np.vectorize(_complement)


def _rev_complement(seq):
    BASES = "NRWSMBDACGTHVKSWY"
    return "".join([BASES[-j] for j in [BASES.find(i) for i in seq][::-1]])


def _consecutive(data, stepsize=1):
    arr = np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)
    arr = [[a.min(), a.max()] for a in arr]
    return arr


def _addZeroCols(freqs):
    freqlength = freqs.shape[1]
    needed = 4 - freqlength
    if needed > 0:
        for i in range(needed):
            freqs = np.column_stack([freqs, np.repeat(0, freqs.shape[0])])
    return freqs


def _get_base_freqs(freqs, ref_alt_array):
    assert freqs.shape == ref_alt_array.shape, "Shape of arrays is different"
    freq_df = pd.DataFrame(ref_alt_array)
    for i_base in range(4):
        for i in range(freqs.shape[0]):
            base = ref_alt_array[i, i_base]
            freq_df.loc[i, f"{base}_freq"] = freqs[i, i_base]
    return freq_df


def load_params_json(param_type):
    import requests

    url = f"https://raw.githubusercontent.com/sanjaynagi/AnoPrimer/refs/heads/main/AnoPrimer/params/{param_type}.json"
    # Send a GET request to the URL
    response = requests.get(url)
    # Raise an exception for bad status codes
    response.raise_for_status()
    # Parse the JSON content
    data = response.json()
    return data

import allel
import gget
import malariagen_data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import primer3
import seaborn as sns
from matplotlib import patches
from plotly.subplots import make_subplots

ag3 = malariagen_data.Ag3(url="gs://vo_agam_release/", pre=True)


def prepare_gDNA_sequence(
    target_loc,
    amplicon_size_range,
    genome_seq,
    assay_name,
    assay_type,
    probe_exclude_region_size=20,
):
    """
    Extracts sequence of interest from genome sequence

    PARAMETERS
    ----------
    target_loc : int
        The target location of the SNP in the genome sequence
    amplicon_size_range : list
        The minimum and maximum size of the amplicon to design primers for
    genome_seq : dask.array.core.Array
        The genome sequence from ag3.genome_sequence()
    assay_name : str
        The name of the assay
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    """
    target = int(target_loc)
    # Set up range for the input sequence, we'll take the max range
    # of the amplicon size and add that either side of the target SNP
    amp_size_range = int(np.max(amplicon_size_range))
    start = target - amp_size_range
    end = target + amp_size_range
    # join array into be one character string, and store the positions
    # of these sequences for later
    target_sequence = "".join(genome_seq[start : end - 1].compute().astype(str))
    gdna_pos = np.arange(start, end).astype(int) + 1
    print(f"The target sequence is {len(target_sequence)} bases long")

    # We need the target snp indices within the region of interest
    target_loc_primer3 = int(np.where(gdna_pos == target)[0])
    target_loc_primer3 = [target_loc_primer3, 10]
    print(f"the target snp is {target_loc_primer3[0]} bp into our target sequence")

    seq_parameters = {
        "SEQUENCE_ID": assay_name,
        "SEQUENCE_TEMPLATE": target_sequence,
        "SEQUENCE_TARGET": target_loc_primer3,
        "GENOMIC_TARGET": target,
    }

    if "probe" in assay_type:
        seq_parameters["SEQUENCE_INTERNAL_EXCLUDED_REGION"] = [
            [1, target_loc_primer3[0] - probe_exclude_region_size],
            [
                target_loc_primer3[0] + probe_exclude_region_size,
                len(target_sequence)
                - (target_loc_primer3[0] + probe_exclude_region_size),
            ],
        ]

    return (target_sequence, gdna_pos, seq_parameters)


def prepare_cDNA_sequence(transcript, genome_seq, assay_name, cDNA_exon_junction):
    """
    Extract exonic sequence for our transcript and record exon-exon junctions

    PARAMETERS
    ----------
    transcript : str
        The AGAP identifier of the transcript to design primers for
    genome_seq : dask.array.core.Array
        The genome sequence from ag3.genome_sequence()
    assay_name : str
        The name of the assay
    cDNA_exon_junction : bool
        If True, cDNA primers will be designed to span an exon-exon junction. We strongly recommend for qPCR purposes.
    """

    # subset gff to your gene
    gff = ag3.geneset()
    gff = gff.query(f"type == 'exon' & Parent == '{transcript}'")
    # Get fasta sequence for each of our exons, and remember gDNA position
    seq = dict()
    gdna_pos = dict()
    for i, exon in enumerate(zip(gff.start, gff.end)):
        seq[i] = "".join(np.array(genome_seq)[exon[0] - 1 : exon[1]].astype(str))
        gdna_pos[i] = np.arange(exon[0] - 1, exon[1])
    # concatenate exon FASTAs into one big sequence
    gdna_pos = np.concatenate(list(gdna_pos.values()))
    target_mRNA_seq = "".join(seq.values())

    # Get list of exon junction positions
    exon_junctions = np.array(np.cumsum(gff.end - gff.start))[:-1]
    exon_sizes = np.array(gff.end - gff.start)[:-1]
    exon_junctions_pos = [ex + gff.iloc[i, 3] for i, ex in enumerate(exon_sizes)]
    print(f"Exon junctions for {transcript}:", exon_junctions, exon_junctions_pos, "\n")
    seq_parameters = {
        "SEQUENCE_ID": assay_name,
        "SEQUENCE_TEMPLATE": target_mRNA_seq,
        "GENOMIC_TARGET": transcript,
    }

    if cDNA_exon_junction:
        seq_parameters["SEQUENCE_OVERLAP_JUNCTION_LIST"] = list(
            map(int, exon_junctions)
        )

    return (target_mRNA_seq, gdna_pos, seq_parameters)


def prepare_sequence(
    target,
    assay_type,
    assay_name,
    genome_seq,
    amplicon_size_range,
    cDNA_exon_junction=True,
):
    """
    Prepare the sequence for primer3, depending on cDNA or gDNA input type

    PARAMETERS
    ----------
    target : str
        The target to design primers for. For gDNA primers, this should be a contig:position string,
        for example '2L:28545767'. For cDNA primers, this should be an AGAP identifier.
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    assay_name : str
        The name of the assay
    genome_seq : dask.array.core.Array
        The genome sequence from ag3.genome_sequence()
    amplicon_size_range : list
        The minimum and maximum size of the amplicon to design primers for
    cDNA_exon_junction : bool
        If True, cDNA primers will be designed to span an exon-exon junction. We strongly recommend for qPCR purposes.
    """

    if any(item in assay_type for item in ["gDNA", "probe"]):
        # genomic DNA
        target_sequence, gdna_pos, seq_parameters = prepare_gDNA_sequence(
            target_loc=target,
            amplicon_size_range=amplicon_size_range,
            genome_seq=genome_seq,
            assay_name=assay_name,
            assay_type=assay_type,
        )
    elif assay_type == "cDNA primers":
        # quantitative PCR
        target_sequence, gdna_pos, seq_parameters = prepare_cDNA_sequence(
            transcript=target,
            genome_seq=genome_seq,
            assay_name=assay_name,
            cDNA_exon_junction=cDNA_exon_junction,
        )

    return (target_sequence, gdna_pos, seq_parameters)


def primer_params(
    assay_type,
    primer_parameters=None,
    n_primer_pairs=None,
    amplicon_size_range=None,
    generate_defaults=False,
):
    """
    adds necessary parameters depending on assay_type, or can
    generate the default parameters

    PARAMETERS
    ----------
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    primer_parameters : dict
        A dictionary of primer3 parameters to use for primer design. If 'default' is specified, default
        primer3 parameters will be generated.
    n_primer_pairs : int
        The number of primer pairs to design.
    amplicon_size_range : list
        The minimum and maximum size of the amplicon to design primers for
    generate_defaults : bool
        If True, default primer3 parameters will be generated.
    """

    if generate_defaults:
        primer_parameters = {
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_TASK": "generic",
            "PRIMER_MIN_SIZE": 17,
            "PRIMER_MAX_SIZE": 24,
            "PRIMER_OPT_TM": 60.0,
            "PRIMER_MIN_TM": 57.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 30.0,
            "PRIMER_MAX_GC": 70.0,
            # this parameter is the minimum distance between successive pairs.
            # If 1, it means successive primer pairs could be identical bar one base shift
            "PRIMER_MIN_THREE_PRIME_DISTANCE": 3,
            # Probe size preferences if selected, otherwise ignored
            "PRIMER_INTERNAL_OPT_SIZE": 16,
            "PRIMER_INTERNAL_MIN_SIZE": 10,
            "PRIMER_INTERNAL_MAX_SIZE": 22,
            # Probe Tm considerations are quite relaxed, assumed that LNAs will be used
            # later to affect TM
            "PRIMER_INTERNAL_MIN_TM": 45,
            "PRIMER_INTERNAL_MAX_TM": 65,
            # Extra primer3 parameters can go here
            # In the same format as above
        }

    primer_parameters["PRIMER_NUM_RETURN"] = n_primer_pairs
    primer_parameters["PRIMER_PRODUCT_SIZE_RANGE"] = amplicon_size_range

    if assay_type == "gDNA primers + probe":
        primer_parameters["PRIMER_PICK_INTERNAL_OLIGO"] = 1
        primer_parameters["PRIMER_PICK_RIGHT_PRIMER"] = 1
        primer_parameters["PRIMER_PICK_LEFT_PRIMER"] = 1
    elif assay_type == "probe":
        primer_parameters["PRIMER_PICK_INTERNAL_OLIGO"] = 1
        primer_parameters["PRIMER_PICK_RIGHT_PRIMER"] = 0
        primer_parameters["PRIMER_PICK_LEFT_PRIMER"] = 0
    return primer_parameters


def primer3_run_statistics(primer_dict, assay_type):
    """
    Prints out primer3 run statistics from the primer3 results dictionary.

    PARAMETERS
    ----------
    primer_dict : dict
        The primer3 results dictionary returned by primer3.designPrimers()
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    """
    _, row_start = _return_oligo_list(assay_type)
    primer_dict = _convert_results_dict_naming(primer_dict)
    # Convert the dict into a pandas dataframe
    primer_df = pd.DataFrame.from_dict(primer_dict.items())
    # Rename the columns
    primer_df = primer_df.rename(columns={0: "parameter", 1: "value"})
    explanations_df = primer_df.iloc[
        :row_start, :
    ]  # Take the first 7 rows which are general
    # Loop through each row and print information
    for (
        idx,
        row,
    ) in explanations_df.iterrows():
        print(row["parameter"], " : ", row["value"], "\n")


def primer3_to_pandas(primer_dict, assay_type):
    """
    Convert primer3 results to pandas dataframe

    PARAMETERS
    ----------
    primer_dict : dict
        The primer3 results dictionary returned by primer3.designPrimers()
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'

    RETURNS
    -------
    primer_df : pandas.DataFrame
        A pandas DataFrame containing the primer sequences and their information.
    """
    oligos, row_start = _return_oligo_list(assay_type)
    # Convert the dict into a pandas dataframe
    primer_dict = _convert_results_dict_naming(primer_dict)
    primer_df = pd.DataFrame.from_dict(primer_dict.items())
    # Rename the columns
    primer_df = primer_df.rename(columns={0: "parameter", 1: "value"})
    # Create a column which is primer pair #, and a column for primer
    # parameter which does not contain primer pair #
    primer_df = primer_df.iloc[row_start:, :].copy()
    primer_df["primer_pair"] = primer_df["parameter"].str.extract("([0-9][0-9]|[0-9])")
    primer_df["parameter"] = primer_df["parameter"].str.replace(
        "(_[0-9][0-9]|_[0-9])", "", regex=True
    )

    # Put the different primer pairs in different columns
    primer_df = primer_df.pivot(
        index="parameter", columns="primer_pair", values="value"
    )

    # Get a list of the rows we need
    primer_span = [f"primer_{oligo}" for oligo in oligos]
    required_info = ["sequence", "TM", "GC_PERCENT"]
    required_info = [
        p + "_" + y for y in required_info for p in primer_span
    ] + primer_span
    required_info = (
        required_info + ["primer_PAIR_PRODUCT_SIZE"]
        if assay_type != "probe"
        else required_info
    )
    required_info = [string.lower() for string in required_info]

    # Subset data frame
    primer_df = primer_df.loc[required_info, np.arange(primer_df.shape[1]).astype(str)]
    return primer_df


def plot_primer_ag3_frequencies(
    primer_df,
    gdna_pos,
    contig,
    sample_sets,
    assay_type,
    seq_parameters,
    out_dir=None,
    sample_query=None,
):
    """
    Loop through n primer pairs, retrieving frequency data and plot allele frequencies

    PARAMETERS
    ----------
    primer_df : pandas.DataFrame
        A pandas DataFrame containing the primer sequences and their information, returned by primer3_to_pandas()
    gdna_pos : numpy.ndarray
        The genomic positions of the target sequence
    contig : str
        The contig of the target sequence, e.g. '2L'
    sample_sets : str or list of str, optional
        Can be a sample set identifier (e.g., "AG1000G-AO") or a list of
        sample set identifiers (e.g., ["AG1000G-BF-A", "AG1000G-BF-B"]) or a
        release identifier (e.g., "3.0") or a list of release identifiers.
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    seq_parameters : dict
        A dictionary of parameters for primer3, returned by prepare_sequence()
    out_dir : str, optional
        The directory to write output files to. If not specified, outputs will not be written to file.
    sample_query: str, optional
        A pandas query string which will be evaluated against the sample
        metadata e.g., "taxon == 'coluzzii' and country == 'Burkina Faso'" to subset
        the genotypes to a specific set of samples.
    """

    if sample_query is not None:
        print(f"Subsetting allele frequencies to {sample_query}")

    name = seq_parameters["SEQUENCE_ID"]
    target = seq_parameters["GENOMIC_TARGET"]

    res_dict = {}
    # Loop through each primer pair and get the frequencies of alternate alleles, storing in dict
    for i in range(len(primer_df.columns)):
        res_dict[i] = _get_primer_alt_frequencies(
            primer_df, gdna_pos, i, sample_sets, assay_type, contig, sample_query
        )

    # Plot data with plotly
    _plotly_primers(
        primer_df=primer_df,
        res_dict=res_dict,
        name=name,
        assay_type=assay_type,
        sample_sets=sample_sets,
        target=target,
        out_dir=out_dir,
    )

    return res_dict


def plot_primer_locs(
    primer_res_dict,
    primer_df,
    contig,
    seq_parameters,
    assay_type,
    legend_loc="best",
    out_dir=None,
):
    """
    Plot the position of the primer sets in relation to any nearby exons

    PARAMETERS
    ----------
    primer_res_dict : dict
        A dictionary containing the primer3 results, returned by primer3.designPrimers()
    primer_df : pandas.DataFrame
        A pandas DataFrame containing the primer sequences and their information, returned by primer3_to_pandas()
    contig : str
        The contig of the target sequence, e.g. '2L'
    seq_parameters : dict
        A dictionary of parameters for primer3, returned by prepare_sequence()
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    legend_loc : str, optional
        The location of the legend in the plot. Can be 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center' or 'center'.
    out_dir : str, optional
        The directory to write output files to. If not specified, outputs will not be written to file.
    """
    oligos, _ = _return_oligo_list(assay_type)
    assay_name = seq_parameters["SEQUENCE_ID"]
    # Load geneset (gff)
    gff = ag3.geneset()
    if any(item in assay_type for item in ["gDNA", "probe"]):
        start = seq_parameters["GENOMIC_TARGET"] - 500
        end = seq_parameters["GENOMIC_TARGET"] + 500
        locgff, min_, max_, genegff = _get_gDNA_locs(gff, contig, start, end)
        min_ = np.min([min_, start])
        max_ = np.max([max_, end])
    elif assay_type == "cDNA primers":
        locgff, min_, max_, genegff = _get_qPCR_locs(
            gff, contig, seq_parameters["GENOMIC_TARGET"]
        )

    if locgff.empty:
        print("No exons in close proximity for loc plot")
        return

    fig, ax = plt.subplots(1, 1, figsize=[16, 4])
    # configure axes
    if min_ in ["inf", "NaN"]:
        min_ = start
    if max_ in ["inf", "NaN"]:
        max_ = end

    ax.set_xlim(min_, max_)
    ax.set_ylim(-0.5, 1.5)
    ax.ticklabel_format(useOffset=False)
    ax.axhline(0.5, color="k", linewidth=0.7, linestyle="--")
    sns.despine(ax=ax, left=True, bottom=False)
    # ax.set_yticks(ticks=[0.2,1.2], size=20)#labels=['- ', '+']
    ax.tick_params(
        top=False, left=False, right=False, labelleft=False, labelbottom=True
    )
    ax.tick_params(axis="x", which="major", labelsize=13)
    ax.set_ylabel("Exons")
    ax.set_xlabel(f"Chromosome {contig} position", fontdict={"fontsize": 14})
    # Add rectangles for exons one at a time
    for _, exon in locgff.iterrows():
        ex_start, ex_end = exon[["start", "end"]]
        e_name = exon["Name"][-2:]
        strand = exon["strand"]
        if strand == "+":
            rect = patches.Rectangle(
                (ex_start, 0.55),
                ex_end - ex_start,
                0.3,
                linewidth=3,
                edgecolor="none",
                facecolor="grey",
                alpha=0.9,
            )
            ax.text((ex_start + ex_end) / 2, 0.65, e_name)
        else:
            rect = patches.Rectangle(
                (ex_start, 0.45),
                ex_end - ex_start,
                -0.3,
                linewidth=3,
                edgecolor="none",
                facecolor="grey",
                alpha=0.9,
            )
            ax.text((ex_start + ex_end) / 2, 0.3, e_name)
        ax.add_patch(rect)

    tot_genes = genegff.shape[0]
    for i, gene in genegff.reset_index(drop=True).iterrows():
        start, end = gene[["start", "end"]]
        diff = np.diff([min_, max_])
        interval = diff / tot_genes + 1
        name_point = min_ + (interval * i + 1)
        strand = gene["strand"]
        if strand == "+":
            rect = patches.Rectangle(
                (start, 0.55),
                end - start,
                0.3,
                linewidth=3,
                edgecolor="black",
                facecolor="none",
            )
            ax.text(
                name_point, 0.95, s=gene["ID"], fontdict={"fontsize": 12}, weight="bold"
            )
        else:
            rect = patches.Rectangle(
                (start, 0.45),
                end - start,
                -0.3,
                linewidth=3,
                edgecolor="black",
                facecolor="none",
            )
            ax.text(
                name_point, -0.3, s=gene["ID"], fontdict={"fontsize": 12}, weight="bold"
            )
        ax.add_patch(rect)

    pal = sns.color_palette("Set2", len(primer_df.columns))
    handles, labels = ax.get_legend_handles_labels()
    for pair in primer_df:
        pair = int(pair)
        for oligo in oligos:
            lower, upper = (
                primer_res_dict[pair][oligo]["position"].min(),
                primer_res_dict[pair][oligo]["position"].max(),
            )

            if oligo == "forward":
                plt.arrow(
                    lower,
                    0.8 + (2 / (10 - (pair))),
                    upper - lower,
                    0,
                    width=0.03,
                    length_includes_head=True,
                    color=pal[pair],
                )
            elif oligo == "reverse":
                plt.arrow(
                    upper,
                    0.8 + (2 / (10 - (pair))),
                    lower - upper,
                    0,
                    width=0.03,
                    length_includes_head=True,
                    color=pal[pair],
                )
            elif oligo == "probe":
                ax.axhline(y=0.8 + (2 / (10 - (pair))), xmin=lower, xmax=upper)
                line = plt.Line2D(
                    (lower, upper),
                    (0.8 + (2 / (10 - (pair))), 0.8 + (2 / (10 - (pair)))),
                    lw=2.5,
                    color=pal[pair],
                )
                ax.add_line(line)
            # manually define a new patch
        patch = patches.Patch(color=pal[pair], label=f"pair {pair}")
        # handles is a list, so append manual patch
        handles.append(patch)
    # plot the legend
    plt.legend(handles=handles, loc=legend_loc)
    if out_dir:
        fig.savefig(f"{out_dir}/{assay_name}_primer_locs.png", dpi=300)


def gget_blat_genome(primer_df, assay_type, assembly="anoGam3"):
    """
    Aligns primers to the AgamP3 genome with BLAT.

    PARAMETERS
    ----------
    primer_df : pandas.DataFrame
        A pandas DataFrame containing the primer sequences and their information, returned by primer3_to_pandas()
    assay_type : str
        The type of assay, either 'gDNA primers' or 'gDNA primers + probe', 'cDNA primers'
        or 'probe'
    assembly : str, optional
        The genome assembly to use with Blat. 'anoGam3' is Anopheles gambiae. Please see the gget documentation for more information.
    """
    oligos, _ = _return_oligo_list(assay_type=assay_type)

    pair_dict = {}
    for primer_pair in primer_df:
        oligo_list = []
        for oligo in oligos:
            seq = primer_df[primer_pair].loc[f"primer_{oligo}_sequence"]
            blat_df = gget.blat(sequence=seq, seqtype="DNA", assembly=assembly)
            if blat_df is None:
                print(f"No hit for {oligo} - pair {primer_pair}")
                continue
            blat_df.loc[:, "primer"] = f"{oligo}_{primer_pair}"
            oligo_list.append(blat_df.set_index("primer"))

        if oligo_list:
            pair_dict[primer_pair] = pd.concat(oligo_list)
        elif not oligo_list:
            continue

    if pair_dict:
        return pd.concat(pair_dict)
    else:
        print("No HITs found for these primer pairs")


def check_and_split_target(target, assay_type):
    # split contig from target
    if target.startswith("AGAP"):
        assert (
            assay_type == "cDNA primers"
        ), "an AGAP identifier is specified, but the assay type is not cDNA primers. Please provide a contig:position identifier for gDNA primers."
        gff = ag3.geneset()
        assert (
            target in gff["ID"].to_list()
        ), f"requested target {target} not in ag3 transcript set"
        contig = gff.query(f"ID == '{target}'")["contig"].unique()[0]
        return (contig, target)
    else:
        assert isinstance(
            target, str
        ), "For genomic DNA the target should be a string, such as '2L:28545767'"
        contig, target = target.split(":")
        assert contig in [
            "2L",
            "2R",
            "3L",
            "3R",
            "X",
        ], "target contig not recognised, should be 2L, 2R, 3L, 3R or X"
        return (contig, int(target))


def designPrimers(
    assay_type,
    assay_name,
    min_amplicon_size,
    max_amplicon_size,
    n_primer_pairs,
    target,
    primer_parameters,
    sample_sets,
    sample_query=None,
    out_dir=None,
    cDNA_exon_junction=True,
):
    """
    Run whole AgamPrimer workflow to design primers/probes with in one function

    Parameters
    ----------
    assay_type: {'gDNA primers', 'cDNA primers', 'gDNA primers + probe', 'probe}, str
        The type of oligos we wish to design. If 'gDNA primers' or 'cDNA primers' are specified,
        then only primers will be designed. If 'gDNA primers + probe' is specified, then primers
        and a probe will be designed. If 'probe' is specified, then only an internal probe will
        be designed.
    assay_name : str
        A name to give the assay, used for naming output files.
    min_amplicon_size : int
        The minimum size of the amplicon we wish to design primers for.
    max_amplicon_size : int
        The maximum size of the amplicon we wish to design primers for. cDNA primers for gene
        expression assays should be designed with a max amplicon size of ~120.
    n_primer_pairs : int
        The number of primer pairs to design.
    target : str
        The target to design primers for. For gDNA primers, this should be a contig:position string,
        for example '2L:28545767'. For cDNA primers, this should be an AGAP identifier.
    primer_parameters : dict or 'default'
        A dictionary of primer3 parameters to use for primer design. If 'default' is specified, default
        primer3 parameters will be generated.
    sample_sets : str or list of str, optional
        Can be a sample set identifier (e.g., "AG1000G-AO") or a list of
        sample set identifiers (e.g., ["AG1000G-BF-A", "AG1000G-BF-B"]) or a
        release identifier (e.g., "3.0") or a list of release identifiers.
    sample_query: str, optional
        A pandas query string which will be evaluated against the sample
        metadata e.g., "taxon == 'coluzzii' and country == 'Burkina Faso'".
    out_dir : str, optional
        The directory to write output files to. If not specified, outputs will not be written to file.
    cDNA_exon_junction : bool, optional
        If True, cDNA primers will be designed to span an exon-exon junction. We strongly recommend
        that this is set to True. In the case of gDNA primers, this parameter is ignored.

        Returns
        -------
        primer_df : pandas.DataFrame
            A pandas DataFrame containing the primer sequences and their information.
        blat_df : pandas.DataFrame
            A pandas DataFrame containing the BLAT alignment information for the primers.
    """

    # check target is valid for assay type and find contig
    oligos, _ = _return_oligo_list(assay_type)
    contig, target = check_and_split_target(target=target, assay_type=assay_type)
    amplicon_size_range = [[min_amplicon_size, max_amplicon_size]]

    # adds some necessary parameters depending on assay type
    if primer_parameters == "default":
        primer_parameters = primer_params(
            primer_parameters=None,
            assay_type=assay_type,
            n_primer_pairs=n_primer_pairs,
            amplicon_size_range=amplicon_size_range,
            generate_defaults=True,
        )
    else:
        primer_parameters = primer_params(
            primer_parameters=primer_parameters,
            assay_type=assay_type,
            n_primer_pairs=n_primer_pairs,
            amplicon_size_range=amplicon_size_range,
            generate_defaults=False,
        )
    # load genome sequence
    genome_seq = ag3.genome_sequence(region=contig)
    print(f"Our genome sequence for {contig} is {genome_seq.shape[0]} bp long")

    target_sequence, gdna_pos, seq_parameters = prepare_sequence(
        target=target,
        assay_type=assay_type,
        assay_name=assay_name,
        genome_seq=genome_seq,
        amplicon_size_range=amplicon_size_range,
        cDNA_exon_junction=cDNA_exon_junction,
    )

    # run primer3
    primer_dict = primer3.designPrimers(
        seq_args=seq_parameters, global_args=primer_parameters
    )

    if assay_type != "probe":
        # check if primer3 has returned any primers
        if int(extract_trailing_digits(primer_dict["PRIMER_PAIR_EXPLAIN"])) == 0:
            print(
                f"No primers found for {assay_name}. For cDNA primers, this is more likely to occur if the target contains only one exon-exon junction. see troubleshooting below for more information. We suggest relaxing the primer parameters \n"
            )
            print(primer3_run_statistics(primer_dict, assay_type))
            return (None, None)

    # AgamPrimer.primer3_run_statistics(primer_dict, assay_type)
    primer_df = primer3_to_pandas(primer_dict=primer_dict, assay_type=assay_type)

    # plot frequencies of alleles in primer pairs
    results_dict = plot_primer_ag3_frequencies(
        primer_df=primer_df,
        gdna_pos=gdna_pos,
        contig=contig,
        sample_sets=sample_sets,
        sample_query=sample_query,
        assay_type=assay_type,
        seq_parameters=seq_parameters,
        out_dir=out_dir,
    )

    # plot primer locations on genome
    plot_primer_locs(
        primer_res_dict=results_dict,
        primer_df=primer_df,
        assay_type=assay_type,
        contig=contig,
        seq_parameters=seq_parameters,
        legend_loc="lower left",
        out_dir=out_dir,
    )

    # add genomic span to primer_df
    spans = []
    for oligo in oligos:
        for i in range(len(primer_df.columns)):
            start = int(results_dict[i][oligo]["position"].min())
            end = int(results_dict[i][oligo]["position"].max())
            span = f"{contig}:{start}-{end}"
            spans.append(span)

    spans = np.array(spans).reshape(len(oligos), len(primer_df.columns))
    spans_df = pd.DataFrame(spans)
    spans_df.index = [f"{o}_genomic_span" for o in oligos]
    spans_df.columns = primer_df.columns
    primer_df = pd.concat([primer_df, spans_df], axis=0).drop(
        index=[f"primer_{o}" for o in oligos]
    )

    # check primers for specificity against the genome and write to file
    blat_df = gget_blat_genome(primer_df, assay_type, assembly="anoGam3")

    # write primer3 and blat output to file
    if out_dir:
        primer_df.to_csv(f"{out_dir}/{assay_name}.{assay_type}.tsv", sep="\t")
        primer_df.to_excel(f"{out_dir}/{assay_name}.{assay_type}.xlsx")
        blat_df.to_csv(f"{out_dir}/{assay_name}.{assay_type}.blat.tsv", sep="\t")

    return (primer_df, blat_df)


def extract_trailing_digits(string):
    import re

    match = re.search(r"\d+$", string)
    if match:
        return match.group(0)
    else:
        return None


def _get_primer_arrays(contig, gdna_pos, sample_sets, assay_type, sample_query=None):
    if any(item in assay_type for item in ["gDNA", "probe"]):
        span_str = f"{contig}:{gdna_pos.min()}-{gdna_pos.max()}"
        snps = ag3.snp_calls(
            region=span_str, sample_sets=sample_sets, sample_query=sample_query
        )  # get genotypes
        ref_alt_arr = snps["variant_allele"].compute().values
        geno = snps["call_genotype"]
        freq_arr = allel.GenotypeArray(geno).count_alleles().to_frequencies()
        pos_arr = gdna_pos
    elif assay_type == "cDNA primers":
        freq_arr = []
        ref_alt_arr = []
        pos_arr = np.array([])
        exon_spans = np.array(_consecutive(gdna_pos)) + 1
        for span in exon_spans:
            span_str = f"{contig}:{span[0]}-{span[1]}"
            snps = ag3.snp_calls(
                region=span_str, sample_sets=sample_sets, sample_query=sample_query
            )  # get genotypes
            ref_alts = snps["variant_allele"]
            geno = snps["call_genotype"]
            freqs = (
                allel.GenotypeArray(geno).count_alleles().to_frequencies()
            )  # calculate allele frequencies
            freqs = _addZeroCols(freqs)
            freq_arr.append(freqs)
            ref_alt_arr.append(ref_alts)
            pos_arr = np.append(pos_arr, np.arange(span[0], span[1] + 1).astype(int))
        freq_arr = np.concatenate(freq_arr)
        ref_alt_arr = np.concatenate(ref_alt_arr)

    return (freq_arr, ref_alt_arr.astype("U13"), pos_arr)


def _get_primer_alt_frequencies(
    primer_df, gdna_pos, pair, sample_sets, assay_type, contig, sample_query
):
    """
    Find the genomic locations of pairs of primers, and runs span_to_freq
    to get allele frequencies at those locations
    """

    oligos, _ = _return_oligo_list(assay_type)
    base_freqs, ref_alt_arr, pos_arr = _get_primer_arrays(
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
        for i in primer_df:
            i = int(i)
            row_i = i + 1

            color = [
                -1 if v == 0 else 1 if v > 0 else 0
                for v in res_dict[i][oligo]["alt_frequency"]
            ]
            colorscale = [[0, "lightgray"], [0.5, "lightgray"], [1, "dodgerblue"]]

            tm = np.round(primer_df.loc[f"primer_{oligo}_tm", str(i)], 2)
            gc = np.round(primer_df.loc[f"primer_{oligo}_gc_percent", str(i)], 2)
            span = f"{int(res_dict[i][oligo]['position'].min())}-{int(res_dict[i][oligo]['position'].max())}"
            # Write text to plot for Tm, GC, span, and 3/5'

            fig.add_trace(
                go.Scatter(
                    x=res_dict[i][oligo]["base_pos"],
                    y=res_dict[i][oligo]["alt_frequency"],
                    customdata=res_dict[i][oligo][
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
                x=res_dict[i][oligo]["base_pos"][0],
                y=0.8,
                text="5'",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[i][oligo]["base_pos"].to_numpy()[-1],
                y=0.8,
                text="3'",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[i][oligo]["base_pos"].to_numpy()[4],
                y=0.92,
                text=span,
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[i][oligo]["base_pos"].to_numpy()[-7],
                y=0.92,
                text=f"GC={gc}",
                showarrow=False,
            )
            fig.add_annotation(
                row=row_i,
                col=idx,
                x=res_dict[i][oligo]["base_pos"].to_numpy()[-3],
                y=0.92,
                text=f"TM={tm}",
                showarrow=False,
            )

            fig.update_xaxes(
                row=row_i,
                col=idx,
                tickmode="array",
                tickvals=res_dict[i][oligo]["base_pos"],
                ticktext=res_dict[i][oligo]["base"],
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
        height=200 * len(primer_df.columns),
        width=500 * len(oligos),
        title_text=title_text,
        title_x=0.5,
        template="simple_white",
        showlegend=False,
    )
    if out_dir:
        fig.write_html(f"{name}_{assay_type}.html")
        fig.write_image(f"{name}_{assay_type}.pdf")
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

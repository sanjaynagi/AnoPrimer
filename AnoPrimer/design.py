import numpy as np
import pandas as pd
import primer3

from .evaluate import AnoPrimerResults
from .utils import (
    _convert_results_dict_naming,
    _return_oligo_list,
    extract_trailing_digits,
    retrieve_data_resource,
)


def design_primers(
    species,
    assay_type,
    assay_name,
    target,
    min_amplicon_size,
    max_amplicon_size,
    n_primer_pairs,
    primer_parameters,
    cDNA_exon_junction=True,
):
    """
    Run whole AnoPrimer workflow to design primers/probes

    Parameters
    ----------
    species: str
        The species to design primers for. Either 'gambiae_sl' or 'funestus'
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
    cDNA_exon_junction : bool, optional
        If True, cDNA primers will be designed to span an exon-exon junction. We strongly recommend
        that this is set to True. In the case of gDNA primers, this parameter is ignored.

        Returns
        -------
        AnoPrimerResults
            An AnoPrimerResults object containing the results of the primer design.
    """

    data_resource = retrieve_data_resource(species=species)

    # check target is valid for assay type and find contig
    contig, target = check_and_split_target(
        species=species, target=target, assay_type=assay_type
    )
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
    genome_seq = data_resource.genome_sequence(region=contig)
    print(f"Our genome sequence for {contig} is {genome_seq.shape[0]} bp long")

    seq_parameters = prepare_sequence(
        species=species,
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

    # AnoPrimer.primer3_run_statistics(primer_dict, assay_type)
    primer_df = primer3_to_pandas(primer_dict=primer_dict, assay_type=assay_type)

    return AnoPrimerResults(
        species=species,
        data_resource=data_resource,
        contig=contig,
        assay_type=assay_type,
        assay_name=assay_name,
        target=target,
        primer_df=primer_df,
        seq_parameters=seq_parameters,
        primer_parameters=primer_parameters,
    )


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
        "GENOMIC_DNA_POSITIONS": gdna_pos,
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

    return seq_parameters


def prepare_cDNA_sequence(
    species, transcript, genome_seq, assay_name, cDNA_exon_junction
):
    """
    Extract exonic sequence for our transcript and record exon-exon junctions

    PARAMETERS
    ----------
    species: str
        The species to design primers for. Either 'gambiae_sl' or 'funestus'
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
    data_resource = retrieve_data_resource(species)

    gff = data_resource.geneset()
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
        "GENOMIC_DNA_POSITIONS": gdna_pos,
    }

    if cDNA_exon_junction:
        seq_parameters["SEQUENCE_OVERLAP_JUNCTION_LIST"] = list(
            map(int, exon_junctions)
        )

    return seq_parameters


def prepare_sequence(
    species,
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
    species: str
        The species to design primers for. Either 'gambiae_sl' or 'funestus'
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
        seq_parameters = prepare_gDNA_sequence(
            target_loc=target,
            amplicon_size_range=amplicon_size_range,
            genome_seq=genome_seq,
            assay_name=assay_name,
            assay_type=assay_type,
        )
    elif assay_type == "cDNA primers":
        # quantitative PCR
        seq_parameters = prepare_cDNA_sequence(
            species=species,
            transcript=target,
            genome_seq=genome_seq,
            assay_name=assay_name,
            cDNA_exon_junction=cDNA_exon_junction,
        )

    return seq_parameters


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


def check_and_split_target(species, target, assay_type):
    data_resource = retrieve_data_resource(species=species)

    # split contig from target
    if target.startswith("AGAP") or target.startswith("LOC"):
        assert (
            assay_type == "cDNA primers"
        ), "an AGAP/AFUN identifier is specified, but the assay type is not cDNA primers. Please provide a contig:position identifier for gDNA primers."
        gff = data_resource.geneset()
        assert (
            target in gff["ID"].to_list()
        ), f"requested target {target} not in ag3/af1 transcript set"
        contig = gff.query(f"ID == '{target}'")["contig"].unique()[0]
        return (contig, target)
    else:
        assert isinstance(
            target, str
        ), "For genomic DNA the target should be a string, such as '2L:28545767'"
        contig, target = target.split(":")
        if species == "gambiae_sl":
            assert contig in [
                "2L",
                "2R",
                "3L",
                "3R",
                "X",
                "2RL",
                "3RL",
            ], "target contig not recognised, should be 2L, 2R, 3L, 3R, 2RL, 3RL, X"
        elif species == "funestus":
            assert contig in [
                "2RL",
                "3RL",
                "X",
            ], "target contig not recognised, should be 2RL, 3RL or X"
        return contig, int(target)

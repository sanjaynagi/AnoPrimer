import gget
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .utils import (
    _get_gDNA_locs,
    _get_primer_alt_frequencies,
    _get_qPCR_locs,
    _plotly_primers,
    _retrieve_span,
    _return_oligo_list,
    retrieve_data_resource,
    round_floats_in_df,
)


class AnoPrimerResults:
    """
    A class to represent the results of primer design and provide methods for analysis and visualization.

    This class encapsulates the results of primer design for a specific assay and provides
    methods to analyze and visualize these results, including SNP frequencies in ag3/af1, primer locations,
    and genome alignment with gget (BLAT).

    Attributes:
        species (str): The species for which primers were designed.
        data_resource: The data resource object used for retrieving genomic data.
        contig (str): The contig or chromosome on which the primers are located.
        target (str): The target sequence or location for primer design.
        assay_type (str): The type of assay (e.g., 'gDNA primers', 'cDNA primers').
        assay_name (str): The name of the assay.
        df (pd.DataFrame): DataFrame containing the designed primers and their properties.
        seq_parameters (dict): Parameters used for sequence preparation.
        primer_parameters (dict): Parameters used for primer design with primer3.
        target_sequence (str): The target DNA sequence.
        gdna_pos (np.array): Array of genomic DNA positions of the target input sequence.
    """

    def __init__(
        self,
        species,
        data_resource,
        assay_type,
        assay_name,
        contig,
        target,
        primer_df,
        seq_parameters,
        primer_parameters,
    ):
        """
        Initialize the AnoPrimerResults object with the results of primer design.

        Args:
            species (str): The species for which primers were designed.
            data_resource: The data resource object used for retrieving genomic data.
            assay_type (str): The type of assay (e.g., 'gDNA primers', 'cDNA primers').
            assay_name (str): The name of the assay.
            contig (str): The contig or chromosome on which the primers are located.
            target (str): The target sequence or location for primer design.
            primer_df (pd.DataFrame): DataFrame containing the designed primers and their properties.
            seq_parameters (dict): Parameters used for sequence preparation.
            primer_parameters (dict): Parameters used for primer design.
        """
        self.species = species
        self.data_resource = data_resource
        self.contig = contig
        self.target = target
        self.assay_type = assay_type
        self.assay_name = assay_name

        self.seq_parameters = seq_parameters
        self.primer_parameters = primer_parameters

        # Extract additional attributes from seq_parameters
        self.target_sequence = seq_parameters.get("SEQUENCE_TEMPLATE")
        self.gdna_pos = np.array(seq_parameters.get("GENOMIC_DNA_POSITIONS"))

        self.df = primer_df
        self.df = round_floats_in_df(self.add_spans_to_df(), decimal_places=2)

    def evaluate_primers(
        self,
        sample_sets,
        sample_query=None,
        out_dir=None,
        legend_loc="best",
        assembly="anoGam3",
    ):
        """
        Evaluate the designed primers by plotting SNP frequencies and primer locations.

        This method generates plots of SNP frequencies for each primer pair and the locations
        of the primers in relation to nearby exons on the genome. It also performs a BLAT alignment
        of the primers to the genome using the gget library.

        Args:
            sample_sets (str or list): Sample set identifier(s) to use for frequency calculation.
            sample_query (str): A pandas query string to filter the samples.
            out_dir (str): Directory to save the output files.
            legend_loc (str, optional): Location of the legend in the plot. Default is "best".
            assembly (str, optional): The genome assembly to use with BLAT. Default is "anoGam3".
        """
        # Plot SNP frequencies for each primer pair
        self.plot_primer_snp_frequencies(
            sample_sets=sample_sets,
            sample_query=sample_query,
            out_dir=out_dir,
        )

        # Plot primer locations in relation to nearby exons
        self.plot_primer_locs(
            legend_loc=legend_loc,
            out_dir=out_dir,
        )

        # Perform BLAT alignment of primers to the genome
        if self.species == "gambiae_sl":
            blat_df = self.gget_blat_genome(assembly=assembly)
            if out_dir is not None and blat_df is not None:
                blat_df.to_csv(f"{out_dir}/{self.assay_name}_blat_results.csv")

    def add_spans_to_df(self):
        df = self.df
        oligos, _ = _return_oligo_list(assay_type=self.assay_type)

        oligo_spans = {}
        for oligo in oligos:
            spans = []
            for pair in df:
                pos = _retrieve_span(
                    df,
                    gdna_pos=self.gdna_pos,
                    oligo=oligo,
                    assay_type=self.assay_type,
                    pair=pair,
                )
                span = f"{self.contig}:{pos.min()}-{pos.max()}"
                spans.append(span)

            oligo_spans[oligo] = pd.Series(
                spans, name=f"primer_{oligo}_span", index=self.df.columns
            )
            df = pd.concat([df, oligo_spans[oligo].to_frame().T])

        return df

    def summarise_metadata(self, sample_sets=None, sample_query=None):
        """
        Retrieve a summary of metadata for samples in the ag3/af1 resource.

        This method creates a pivot table summarizing the sample counts by sample set,
        year, country, and taxon.

        Args:
            sample_sets (str or list, optional): Sample set identifier(s) to filter the data.
            sample_query (str, optional): A pandas query string to filter the samples.

        Returns:
            pd.DataFrame: A pivot table summarizing the sample metadata.
        """
        # Retrieve sample metadata based on the provided filters
        df_samples = self.data_resource.sample_metadata(
            sample_sets=sample_sets, sample_query=sample_query
        )

        # Create a pivot table to summarize the data
        pivot_country_year_taxon = df_samples.pivot_table(
            index=["sample_set", "year", "country"],
            columns=["taxon"],
            values="sample_id",
            aggfunc="count",
            fill_value=0,
        )

        return pivot_country_year_taxon

    def plot_primer_snp_frequencies(
        self,
        sample_sets,
        sample_query=None,
        out_dir=None,
    ):
        """
        Plot SNP frequencies for each primer pair.

        This method retrieves allele frequency data for each primer pair and generates
        a plot using the Plotly library.

        Args:
            sample_sets (str or list): Sample set identifier(s) to use for frequency calculation.
            sample_query (str, optional): A pandas query string to filter the samples.
            out_dir (str, optional): Directory to save the output files. If None, files are not saved.

        Returns:
            dict: A dictionary containing the frequency data for each primer pair.
        """
        name = self.assay_name
        target = self.seq_parameters["GENOMIC_TARGET"]

        if sample_query is not None:
            print(f"Subsetting allele frequencies to {sample_query}")

        # Loop through each primer pair and get the frequencies of alternate alleles
        res_dict = {}
        for pair in self.df:
            res_dict[pair] = _get_primer_alt_frequencies(
                species=self.species,
                primer_df=self.df,
                gdna_pos=self.gdna_pos,
                pair=pair,
                assay_type=self.assay_type,
                contig=self.contig,
                sample_sets=sample_sets,
                sample_query=sample_query,
            )

        # Generate the plot using Plotly
        _plotly_primers(
            primer_df=self.df,
            res_dict=res_dict,
            name=name,
            assay_type=self.assay_type,
            sample_sets=sample_sets,
            target=target,
            out_dir=out_dir,
        )

        return res_dict

    def plot_primer_locs(
        self,
        legend_loc="best",
        out_dir=None,
    ):
        """
        Plot the positions of primer sets in relation to nearby exons.

        This method generates a matplotlib figure showing the locations of primers
        and nearby exons on the genome.

        Args:
            legend_loc (str, optional): Location of the legend in the plot. Default is "best".
            out_dir (str, optional): Directory to save the output files. If None, files are not saved.
        """
        # Retrieve necessary data and parameters
        data_resource = retrieve_data_resource(species=self.species)
        oligos, _ = _return_oligo_list(self.assay_type)
        assay_name = self.assay_name
        exon_id_col = "Name" if self.species == "gambiae_sl" else "ID"

        # Load geneset (gff) and get relevant locations
        gff = data_resource.geneset()
        if any(item in self.assay_type for item in ["gDNA", "probe"]):
            start = self.seq_parameters["GENOMIC_TARGET"] - 500
            end = self.seq_parameters["GENOMIC_TARGET"] + 500
            locgff, min_, max_, genegff = _get_gDNA_locs(gff, self.contig, start, end)
            min_ = np.min([min_, start])
            max_ = np.max([max_, end])
        elif self.assay_type == "cDNA primers":
            locgff, min_, max_, genegff = _get_qPCR_locs(
                gff, self.contig, self.seq_parameters["GENOMIC_TARGET"]
            )
            start, end = min_, max_

        if locgff.empty:
            print("No exons in close proximity for loc plot")
            return

        # remove duplicate exons from different transcripts
        locgff = locgff.drop_duplicates("Name")

        # Set up the plot
        fig, ax = plt.subplots(1, 1, figsize=[16, 4])
        self._configure_plot_axes(ax, min_, max_, start, end)
        # Plot exons, genes, primer spans
        self._plot_exons(ax, locgff, exon_id_col)
        self._plot_genes(ax, genegff, min_, max_)
        handles = self._plot_primers(ax, oligos, locgff)

        # Add legend and save if out_dir is provided
        plt.legend(handles=handles, loc=legend_loc)
        if out_dir:
            fig.savefig(
                f"{out_dir}/{assay_name}_primer_locs.png", dpi=300, bbox_inches="tight"
            )

    def _configure_plot_axes(self, ax, min_, max_, start, end):
        """Helper method to configure plot axes."""
        if min_ in ["inf", "NaN"]:
            min_ = start
        if max_ in ["inf", "NaN"]:
            max_ = end
        ax.set_xlim(min_, max_)
        ax.set_ylim(-0.5, 1.5)
        ax.ticklabel_format(useOffset=False)
        ax.axhline(0.5, color="k", linewidth=0.7, linestyle="--")
        sns.despine(ax=ax, left=True, bottom=False)
        ax.tick_params(
            top=False, left=False, right=False, labelleft=False, labelbottom=True
        )
        ax.tick_params(axis="x", which="major", labelsize=13)
        ax.set_ylabel("Exons")
        ax.set_xlabel(f"Chromosome {self.contig} position", fontdict={"fontsize": 14})

    def _plot_exons(self, ax, locgff, exon_id_col):
        """Helper method to plot exons."""
        for _, exon in locgff.iterrows():
            ex_start, ex_end = exon[["start", "end"]]
            e_name = exon[exon_id_col][-2:]
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

    def _plot_genes(self, ax, genegff, min_, max_):
        """Helper method to plot genes."""
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
                    name_point,
                    0.95,
                    s=gene["ID"],
                    fontdict={"fontsize": 12},
                    weight="bold",
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
                    name_point,
                    -0.3,
                    s=gene["ID"],
                    fontdict={"fontsize": 12},
                    weight="bold",
                )
            ax.add_patch(rect)

    def _plot_primers(self, ax, oligos, locgff):
        """Helper method to plot primers."""

        def _generate_primer_pair_positions(num_pairs, start=1, end=1.45):
            if num_pairs == 1:
                return [start]
            step = (end - start) / (num_pairs - 1)
            return [start + i * step for i in range(num_pairs)]

        def _is_exon_exon_junction(start, end, locgff):
            """Check if the primer spans an exon-exon junction."""
            exons = locgff[(locgff["start"] <= end) & (locgff["end"] >= start)]
            return len(exons) > 1

        pal = sns.color_palette("Set2", len(self.df.columns))
        handles, labels = ax.get_legend_handles_labels()
        pair_ypos = _generate_primer_pair_positions(len(self.df.columns))

        for pair in self.df:
            pair = int(pair)
            pair_idx = pair - 1  # python based indexing
            for oligo in oligos:
                oligo_pos = _retrieve_span(
                    primer_df=self.df,
                    gdna_pos=self.gdna_pos,
                    pair=pair,
                    oligo=oligo,
                    assay_type=self.assay_type,
                )
                lower, upper = oligo_pos.min(), oligo_pos.max()

                # dashed lines for exon junction spanning primers
                exonexonjunction = _is_exon_exon_junction(lower, upper, locgff)
                linestyle = "dotted" if exonexonjunction else "solid"

                if oligo == "forward":
                    plt.arrow(
                        lower,
                        pair_ypos[pair_idx],
                        upper - lower,
                        0,
                        width=0.03,
                        length_includes_head=True,
                        color=pal[pair_idx],
                        linestyle=linestyle,
                    )
                elif oligo == "reverse":
                    plt.arrow(
                        upper,
                        pair_ypos[pair_idx],
                        lower - upper,
                        0,
                        width=0.03,
                        length_includes_head=True,
                        color=pal[pair_idx],
                        linestyle=linestyle,
                    )
                elif oligo == "probe":
                    ax.axhline(y=pair_ypos[pair_idx], xmin=lower, xmax=upper)
                    line = plt.Line2D(
                        xdata=(lower, upper),
                        ydata=(pair_ypos[pair_idx], pair_ypos[pair_idx]),
                        linewidth=2,
                        linestyle=(0, (1, 1)),
                        color=pal[pair_idx],
                    )
                    ax.add_line(line)

            patch = patches.Patch(color=pal[pair_idx], label=f"pair {pair}")
            handles.append(patch)

        return handles

    def gget_blat_genome(self, assembly="anoGam3"):
        """
        Align primers to the genome using BLAT.

        This method uses the gget library to perform BLAT alignment of the designed primers
        against the specified genome assembly.

        Args:
            assembly (str, optional): The genome assembly to use with BLAT. Default is "anoGam3".

        Returns:
            pd.DataFrame or None: A DataFrame containing BLAT alignment results, or None if no hits found.
        """
        oligos, _ = _return_oligo_list(assay_type=self.assay_type)

        pair_dict = {}
        for primer_pair in self.df:
            oligo_list = []
            for oligo in oligos:
                seq = self.df[primer_pair].loc[f"primer_{oligo}_sequence"]
                blat_df = gget.blat(sequence=seq, seqtype="DNA", assembly=assembly)
                if blat_df is None:
                    print(f"No hit for {oligo} - pair {primer_pair}")
                    continue
                blat_df.loc[:, "primer"] = f"{oligo}_{primer_pair}"
                oligo_list.append(blat_df.set_index("primer"))

            if oligo_list:
                pair_dict[primer_pair] = pd.concat(oligo_list)

        if pair_dict:
            return pd.concat(pair_dict)
        else:
            print("No HITs found for these primer pairs")
            return None

    def to_csv(self, file_path, **kwargs):
        self.df.to_csv(file_path, **kwargs)

    def to_excel(self, file_path, **kwargs):
        self.df.to_excel(file_path, **kwargs)

    def __str__(self):
        return f"AnoPrimerResults for {self.assay_name} ({self.assay_type})"

    def __repr__(self):
        return self.__str__()

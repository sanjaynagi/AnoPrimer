<div align="center">

[<img src="https://github.com/sanjaynagi/AgamPrimer/blob/main/graphics/AgamPrimer_logo.png?raw=True" width="400"/>](https://github.com/sanjaynagi/AgamPrimer/blob/main/graphics/AgamPrimer_logo.png?raw=True)   

[![DOI](https://zenodo.org/badge/503315581.svg)](https://zenodo.org/badge/latestdoi/503315581)
[![Execute notebook](https://github.com/sanjaynagi/AgamPrimer/workflows/Execute%20notebook/badge.svg)](https://github.com/sanjaynagi/AgamPrimer/actions?query=workflow:"Execute+notebook")
[![GitHub release](https://img.shields.io/github/release/sanjaynagi/AgamPrimer?include_prereleases=&sort=semver&color=blue)](https://github.com/sanjaynagi/AgamPrimer/releases/)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)

</div>


Primer design in *Anopheles gambiae s.l* taking into account SNP variation in primer binding sites, using primer3-py and malariagen_data. Supports genomic DNA primers, cDNA primers and hybridisation probes:    


<div align="center">


[![Custom badge](https://img.shields.io/endpoint?color=white&logo=Google%20Colab&url=https%3A%2F%2Fraw.githubusercontent.com%2Fsanjaynagi%2FAgamPrimer%2Fmain%2Fgraphics%2Fbadge.json)](https://colab.research.google.com/github/sanjaynagi/AgamPrimer/blob/main/notebooks/AgamPrimer-long.ipynb)   


 [![Custom badge](https://img.shields.io/endpoint?color=red&logo=Google%20Colab&url=https%3A%2F%2Fraw.githubusercontent.com%2Fsanjaynagi%2FAgamPrimer%2Fmain%2Fgraphics%2Fbadge-short.json)](https://colab.research.google.com/github/sanjaynagi/AgamPrimer/blob/main/notebooks/AgamPrimer-short.ipynb)   

</div>



#### Release notes

- 0.5.12 - for cDNA primers, designing over exon-exon junctions is now optional due to `cDNA_exon_junction` argument. Docstring added to `designPrimers()`. 
- 0.5.11 - feedback on troubleshooting now provided in `designPrimers()` function if primer design fails
- 0.5.10 - region string provided instead of contig 
- 0.5.9 - 'cDNA primers' replaces 'qPCR primers' throughout.
- 0.5.8 - Two notebooks, short and long, including all in one function. can now check the primers for specificity with gget blat implementation. Most versions >0.4.0  were for a workshop and github actions CI.
- 0.4.0 - plot_primer_ag3_frequencies() now uses plotly to make an interactive plot, where one can hover over primer bases, which returns the exact frequency for each base.
- 0.3.4 - Fix bug in get_gDNA_sequence() and in plot_primer_locs()
- 0.3.3 - Introduced feature to enable a sample query (e.g subset to a specific species within a cohort)
- 0.3.2 - minor fix to remove printing of variable in `plot_primer_ag3_frequencies()`
- 0.3.1 - minor fix to bug introduced in 0.3.0 to `plot_primer()` function
- 0.3.0 - Support for probe design added, functions restructured to accommodate this

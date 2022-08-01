[<img src="https://github.com/sanjaynagi/AgamPrimer/blob/main/graphics/AgamPrimer_logo.png?raw=True" width="400"/>](https://github.com/sanjaynagi/AgamPrimer/blob/main/graphics/AgamPrimer_logo.png?raw=True)    
[![DOI](https://zenodo.org/badge/503315581.svg)](https://zenodo.org/badge/latestdoi/503315581)


Primer Design in *Anopheles gambiae s.l* taking into account SNP variation in primer binding sites, using primer3-py and malariagen_data. Supports:

- genomic DNA primers
- quantitative PCR primers (cDNA, designed to span exon-exon junctions)
- probe design

------------>------------>------------>------------>   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sanjaynagi/AgamPrimer/blob/main/Primer-Design-in-Anopheles-gambiae.ipynb)    <------------<------------<------------<------------

#### Release notes

- 0.3.2 - minor fix to remove printing of variable in `plot_primer_ag3_frequencies()`
- 0.3.1 - minor fix to bug introduced in 0.3.0 to `plot_primer()` function 
- 0.3.0 - Support for probe design added, functions restructured to accommodate this

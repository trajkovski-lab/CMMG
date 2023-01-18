# Comprehensive Mouse Microbiome Genome (CMMG) catalog

## Citation

[Publication](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009947)

[Pre-print](https://doi.org/10.1101/2021.03.18.435958)

## Get Started


### Quantify a mouse or human megagenome using a host-specific genome collection

We created Kreken2 + braken databases for CMMG and the Unified human gastrointestinal genome collection (UHGG) with a manually curated taxonomy. They are availalb at https://ezmeta.unige.ch/CMMG/Kraken2db/

We reccomend to use our snakemake pipeline for efficient use of Kraken2 + Braken https://github.com/SilasK/Kraken

### Analysis

Once you have the Kraken2 quantification you can use this [jupyter notebook](https://colab.research.google.com/github/trajkovski-lab/CMMG/blob/main/notebooks/Analyze-cold-adapted-microbiota.ipynb) that shows how CMMG can be used for the functional and taxonomic analysis of mouse metagenome data.


![Volcanoplot](Figures/Volcanoplot.svg) ![barplot](Figures/barplot_Butyrate.svg)


The same code can be used for the analysis of **human metagenome** data, by using our [functionally annotions](https://ezmeta.unige.ch/CMMG) of UHHG.

#### Code and source data for:

_Douglas GM, Shapiro BJ. 2024. Pseudogenes act as a neutral reference for detecting selection in prokaryotic pangenomes. Nat Ecol Evol 1–11. doi:[10.1038/s41559-023-02268-6](https://doi.org/10.1038/s41559-023-02268-6)_

The published paper requires is behind a paywall, but you can read an earlier version of the manuscript as a [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.05.17.541134v2). In addition, this [blog post](https://communities.springernature.com/posts/using-degenerating-genes-to-understand-the-evolution-of-rare-intact-genes-across-bacteria) is a more accessible description of our work.

Please feel free to [open an issue](https://github.com/gavinmdouglas/pangenome_pseudogene_null/issues) if you have any questions.


#### Repository structure:

* `display_source` - Source data for each display item, included for convenience.

* `scripts`

	* `analyses` - R scripts to generate reported models and to run key statistical tests

	* `broad_analysis_processing` - R scripts for analyzing and processing files for broad pangenome analysis

	* `display` - R scripts for generating display items

	* `indepth_analysis_processing` - R scripts for analyzing and processing files for indepth pangenome analysis

	* `preprocessing` - R and Python scripts for preprocessing raw and intermediate files

	* `sanity_checks` - Quick scripts for running checks of specific results reported in the manuscript (i.e., regenerating results with independent code).

	* `text_results` - Small scripts for computing values reported in text (intended to be run on datafiles distributed on FigShare repository)


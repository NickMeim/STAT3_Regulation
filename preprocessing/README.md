## Algorithms for the pre-processing of the raw data
1. distance_scores.R: Function calculating the GSEA score distance[^1][^2] between samples, as adopted by this [repository](https://github.com/BioSysLab/deepSIBA) which is part of the deepSIBA publication[^3].
2. preprocessingDistances.R: Loading and pre-processing data to calculate in the end the GSEA distances between a subset of the data for quality evaluation purposes later.
3.  enrichment_calculations.R: It is a global function that can perform gene set enrichment analysis (GSEA) of various prior knowledge sets such as: KEGG pathways, GO terms, Transcription factors (based on some regulatory network) and MSIG database. For more info see the script.

The produced data are used to create the figures stored in the figures folder and the results in the results folder.

## References
[^1]: F. Iorio, R. Bosotti, E. Scacheri, V. Belcastro, P. Mithbaokar, R. Ferriero, L. Murino, R. Tagliaferri, N. Brunetti-Pierri and A. Isacchi, Proc. Natl. Acad. Sci. U. S. A., 2010, 107, 14621–14626
[^2]: F. Li, Y. Cao, L. Han, X. Cui, D. Xie, S. Wang and X. Bo, OMICS: J. Integr. Biol., 2013, 17, 116–118
[^3]: Fotis, C., Meimetis, N., Sardis, A., & Alexopoulos, L. G. (2021). DeepSIBA: chemical structure-based inference of biological alterations using deep learning. Molecular Omics, 17(1), 108-120.

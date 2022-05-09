## Post processing of therapeutic discovery results
1.	pathToSTAT3.py : Script to validate that a path from source gene candidates to the STAT3 target exists in the human intracellular network.
2.	valGEO.R : Validating therapeutic candidates in GEO
3.	distance_scores.R: Function calculating the GSEA score distance[^1][^2] between samples, as adopted by this [repository](https://github.com/BioSysLab/deepSIBA) which is part of the deepSIBA publication[^3].
4.	enrichment_calculations.R: It is a global function that can perform gene set enrichment analysis (GSEA) of various prior knowledge sets such as: KEGG pathways, GO terms, Transcription factors (based on some regulatory network) and MSIG database. For more info see the script.
5.	STAT3_profiling.R: It contains code for builiding a profile for STAT3 knockdown based on the shRNA signatures in the data used (see data folder). It performs pathway analysis to infer a representative KEGG pathway profile, TF activity inference, signaling network inference using the CARNIVAL[^4] R package and finally it identifies the candidate genes to be tested, as well as those that have a similar KEGG pathway signature with STAT3 knockdown in the end.
6.	CompoundScreening.R : It contains the code for screening drugs that induce a similar effect as STAT3 knockdown.

## References
[^1]: F. Iorio, R. Bosotti, E. Scacheri, V. Belcastro, P. Mithbaokar, R. Ferriero, L. Murino, R. Tagliaferri, N. Brunetti-Pierri and A. Isacchi, Proc. Natl. Acad. Sci. U. S. A., 2010, 107, 14621–14626
[^2]: F. Li, Y. Cao, L. Han, X. Cui, D. Xie, S. Wang and X. Bo, OMICS: J. Integr. Biol., 2013, 17, 116–118
[^3]: Fotis, C., Meimetis, N., Sardis, A., & Alexopoulos, L. G. (2021). DeepSIBA: chemical structure-based inference of biological alterations using deep learning. Molecular Omics, 17(1), 108-120.
[^4]: Liu, Anika, et al. "From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL." NPJ systems biology and applications 5.1 (2019): 1-10.


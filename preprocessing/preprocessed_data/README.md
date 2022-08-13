## Preprocessed data
1. dupl_mean_dist.rds : RDS file to be opened in R, containing the GSEA distances (see overview page with references) between a subset of data that have duplicate signatures.
2. FilteredOmnipath.tsv: Text file (tsv format) containing the interactions of the OmniPath resource[^1] used to see if we can access STAT3 from another gene via a signalling cascade.
Duplicate signatures are samples that are derived by testing the same drug on the same cell-type,with same dose and for the same treatment duration.
3. PKN-Model.tsv: Trimmed prior knowledge network derived from OmniPath[^1].
4. MaintainInter.tsv: Interactions forced to maintain during trimming.
5. targetd_tfs.tsv: Transcription factors that are targetd by some shRNA knock-down.

## References
[^1]: Türei, D., Valdeolivas, A., Gul, L., Módos, D., Ceccarelli, F., Palacio, N., Ivanova, O., Klein, M., Gábor, A., &amp; Ölbei, M. (2016). OmniPath: intra- &amp; intercellular signaling knowledge. Retrieved May 3, 2022.


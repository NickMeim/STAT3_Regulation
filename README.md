## Overview
This study aims to computationally identify therapeutic desings for cancer e.g. :
1. Reversing NK-cells dysregulation
2. Regulating STAT3 activity

The current repository contains code for:
1. Initial evaluation of the quality and preprocessing of the data.
2. Comparison of the transcriptomic signatures via the GSEA score distance function[^1][^2] as adopted by this [repository](https://github.com/BioSysLab/deepSIBA)[^3]
3. Machine learning algorithms to manipulate transcriptomic profiles and make disease-relevent predictions (still under development).

## Data
The transcriptomic signatures (level 5 profiles) of the L1000 CMap resource[^4] are used for this study, together with data from the KEGG database[^5] (accessed via the Bioconductor resource[^6]).

**Details on how to access these data can be found in the data folder**, but generally the main resources can be accessed [here](https://clue.io/data/CMap2020?fbclid=IwAR1Uc379nDYELH8lYU9MPI9TiAT3054_55g72Ymbgm7FAW7WZnPD3YBCXeI#LINCS2020)

## Folder structure
1. **data** : Folder that should contain the raw data of the study.
2. **figures** : Folder containing the scripts to produce the figures of the study (also containing the figures).
3. **learning (still under development)** : Folder containing machine and deep learning algorithms.
4. **preprocessing** : Folder containing scripts to pre-process the raw data and evaluate their quality.
	*preprocessed_results : Here the pre-processed data to be used in subsequent analysis are stored.
5. **results** : Here the results of a subsequent analysis should be stored.

## Instalation
The study utilizes multiple resources from the Python and R programming languages.

**R dependencies**:

```bash
# clone the source code on your directory
$ git clone https://github.com/NickMeim/final_project_440.git
```


## References
[^1]: F. Iorio, R. Bosotti, E. Scacheri, V. Belcastro, P. Mithbaokar, R. Ferriero, L. Murino, R. Tagliaferri, N. Brunetti-Pierri and A. Isacchi, Proc. Natl. Acad. Sci. U. S. A., 2010, 107, 14621–14626
[^2]: F. Li, Y. Cao, L. Han, X. Cui, D. Xie, S. Wang and X. Bo, OMICS: J. Integr. Biol., 2013, 17, 116–118
[^3]: Fotis, C., Meimetis, N., Sardis, A., & Alexopoulos, L. G. (2021). DeepSIBA: chemical structure-based inference of biological alterations using deep learning. Molecular Omics, 17(1), 108-120.
[^4]: Subramanian, Aravind, et al. "A next generation connectivity map: L1000 platform and the first 1,000,000 profiles." Cell 171.6 (2017): 1437-1452.
[^5]: Kanehisa, Minoru. "The KEGG database." (2002): 91-101.
[^6]: Gentleman, Robert C., et al. "Bioconductor: open software development for computational biology and bioinformatics." Genome biology 5.10 (2004): 1-16.


## Overview
This study aims to computationally identify therapeutic desings for cancer e.g. :
1. Reversing NK-cells dysregulation
2. Regulating STAT3 activity

The current repository contains code for:
1. Initial evaluation of the quality and preprocessing of the data.
2. Comparison of the transcriptomic signatures via the GSEA score distance function[^1][^2] as adopted by this [repository](https://github.com/BioSysLab/deepSIBA) which is part of the deepSIBA publication[^3].
3. Machine learning algorithms to manipulate transcriptomic profiles and make disease-relevent predictions (still under development).

## Data
The transcriptomic signatures (level 5 profiles) of the L1000 CMap resource[^4] are used for this study, together with data from the KEGG database[^5] (accessed via the Bioconductor resource[^6]).

The transcriptomic profiles were generated by measuring 978 important (landmark) genes in cancer with a luminex bead-based assay and computationally infering the rest[^4]. 

**Details on how to access these data can be found in the data folder**, but generally the main resources can be accessed [here](https://clue.io/data/CMap2020?fbclid=IwAR1Uc379nDYELH8lYU9MPI9TiAT3054_55g72Ymbgm7FAW7WZnPD3YBCXeI#LINCS2020)

## Folder structure
1. **data** : Folder that should contain the raw data of the study.
2. **figures** : Folder containing the scripts to produce the figures of the study (also containing the figures).
3. **learning (still under development)** : Folder containing machine and deep learning algorithms.
4. **preprocessing** : Folder containing scripts to pre-process the raw data and evaluate their quality.
	* **preprocessed_data** : Here the pre-processed data to be used in subsequent analysis are stored.
5. **results** : Here the results of a subsequent analysis should be stored.

## Installation
The study utilizes multiple resources from the Python and R programming languages.

**R dependencies**: 

If you want, you can run the setup_rlibs.R script in Rstudio or any R command line environment. 

This will install **the latest versions of the libraries used (not the specific ones that were used in the study)**.

**Please keep in mind that this will update some packages which may cause conflicts in your setup.**

You can check the list bellow and manually install your preferences without running setup_rlibs.R (which you can open however for guidance).

**Important Note:**
* **This installation was performed in a WINDOWS enviroment.** 
* **For a Linux installation there might be needed some manual installation of external dependencies (especially) for tidyverse. Please check libraries' documentation online**
* ** For a MAC installation we engourage checking online**

In a quick overview the following R libraries and versions (**although any version of the following libraries is appropriate**) were/are used:
1. [R](https://cran.r-project.org/bin/windows/base/) version 4.1.2
2. [tidyverse](https://www.tidyverse.org/packages/) 1.3.1
3. [BiocManager](https://www.bioconductor.org/install/) 1.30.16
4. [cmapR](https://bioconductor.org/packages/release/bioc/html/cmapR.html) 1.4.0
5. [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) 3.13.0
6. [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html) 2.36.0
7. [doFuture](https://cran.r-project.org/web/packages/doFuture/index.html) 0.12.0
8. [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) 1.8.2
9. [ggplot2](https://ggplot2.tidyverse.org/) 3.3.5
10. [ggpubr](https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0) 0.4.0
11. [GeneExpressionSignature](https://www.bioconductor.org/packages/release/bioc/html/GeneExpressionSignature.html) 1.38.0

**Python dependencies**: 
First install conda (anaconda) environment in your computer and then you can use the commands **in a bash-terminal** after the list of libraries.

**Important Note:**
* **Pytorch GPU installation CHANGES according to your NVIDIA GPU and cuda version. Check the [pytorch installation guide](https://pytorch.org/get-started/locally/) here for more information.**
* **This installation was performed in a WINDOWS enviroment. For other environments check libraries' documentation** 

In a quick overview the following Python libraries and versions (**although different versions are POSSIBLY also appropriate**) were/are used:
1. [python](https://www.python.org/downloads/) 3.8.8
2. [seaborn](https://seaborn.pydata.org/installing.html) 0.11.2 (version does not matter for this library)
3. [numpy](https://numpy.org/install/) 1.20.3 (version does not matter for this library)
4. [pandas](https://pandas.pydata.org/docs/getting_started/install.html) 1.3.5 (version does not matter for this library)
5. [matplotlib](https://anaconda.org/conda-forge/matplotlib) 3.5.1 (version does not matter for this library)
6. [scipy](https://anaconda.org/anaconda/scipy) 1.7.3
7. CUDA Version: 11.1 - 11.5 **(Very important for pytorch GPU installation)**
8. [pytorch](https://pytorch.org/get-started/locally/) **(Important to download a compatible version for your system. Click on pytorch to view more information on its original website)**
9. GPU: NVIDIA GeForce RTX 3060

```bash
# After installing anaconda create a conda environment:
conda create -n myenv python=3.8.8
conda install numpy
conda install pandas
conda install -c conda-forge matplotlib
conda install seaborn
conda install -c anaconda scipy
# This is hardware-specific
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
```


## References
[^1]: F. Iorio, R. Bosotti, E. Scacheri, V. Belcastro, P. Mithbaokar, R. Ferriero, L. Murino, R. Tagliaferri, N. Brunetti-Pierri and A. Isacchi, Proc. Natl. Acad. Sci. U. S. A., 2010, 107, 14621–14626
[^2]: F. Li, Y. Cao, L. Han, X. Cui, D. Xie, S. Wang and X. Bo, OMICS: J. Integr. Biol., 2013, 17, 116–118
[^3]: Fotis, C., Meimetis, N., Sardis, A., & Alexopoulos, L. G. (2021). DeepSIBA: chemical structure-based inference of biological alterations using deep learning. Molecular Omics, 17(1), 108-120.
[^4]: Subramanian, Aravind, et al. "A next generation connectivity map: L1000 platform and the first 1,000,000 profiles." Cell 171.6 (2017): 1437-1452.
[^5]: Kanehisa, Minoru. "The KEGG database." (2002): 91-101.
[^6]: Gentleman, Robert C., et al. "Bioconductor: open software development for computational biology and bioinformatics." Genome biology 5.10 (2004): 1-16.


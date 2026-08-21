# Van Mooy Lab Global Ocean Lipidome Dataset
This repository contains untargeted lipid from surface ocean samples across the globe. It was originally created to host data from [Holm *et al.* 2022](https://www.science.org/doi/10.1126/science.abn7455), but has been updated with new data since. As of the latest update, it includes annotated chemical data from 1228 lipid species across 930 samples from the surface ocean.

Since version 1.0.2 data is now offered at three levels of processing to aid the user in their own investigation. Latest data can be found in [/Data/csv/Joined/v110/](https://github.com/hholm/OceanLipidome/tree/main/Data/csv/Joined/v110/). Including sample metadata this comprises four files `.csv` files per version. See  [Holm_2022_quantify](https://github.com/hholm/OceanLipidome/blob/main/Scripts/Holm_2022_quantify.Rmd)  for code.

**global_Suspended_metadata_final**: Sample metadata.
**global_suspended_peakarea**: No corrections or post-processing applied. Units are in peakarea.  
**global_suspended_pg_oncolumn**: Quantified using standard curves and blank subtracted with field blanks. Units are in picograms.  
**global_suspended_pg_perLiter**: Quantified using standard curves, blank subtracted with field blanks, recovery corrected with DNP-PE and volume of seawater filtered. Units of picograms per liter.

Convensions are as described below and generally match the format of [LOBSTAHS](https://github.com/vanmooylipidomics/LOBSTAHS) output. If you use data from this dataset in your work please cite the original publication introducing it:
 >  H. C. Holm, H. F. Fredricks, S. M. Bent, D. P. Lowenstein, J. E. Ossolinski, K. W. Becker, W. M. Johnson, K. Schrage, B. A. S. V. Mooy, Global ocean lipidomes show a universal relationship between temperature and lipid unsaturation. _Science_ **376**, 1487–1491 (2022). 
 >  [10.1126/science.abn7455](doi.org/10.1126/science.abn7455)

If you have any questions or inquiries about this, please do not hesitate to contact me at hholm@whoi.edu. 

Version 1.1.0 - 8/21/2026 - [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7035947.svg)](https://doi.org/10.5281/zenodo.7035947)

## Screened Data
#### Sample Naming Convention
Sample filenames follow the following naming scheme: "S","B", or "QC" indicates whether each file is a field sample, blank, or a pooled quality control sample. These are numbered and labeled with the cruise they were sampled on. The final sequence beginning with "QE00-" denotes the unique file number of the original file.</p>

### Lipid Annotation Data
All lipid annotations including quantified mass for each sample are available in **Holm_global_suspended_final.csv** or in the similarly labeled tab in **Holm_Combined_Data.xlsx**. Additionally, lipids annotation information such as chemical formula, retention time, adduct mass, and annotations are included in the preliminary columns of this sheet. Lipid annotations are rows, samples are column. All lipid units are in picograms on-column. 

#### Lipid Annotation Metadata
- Lipid metadata is contained in columns A-I, MS/MS validations are in columns J-M, adduct validations in columns N-S, followed by samples T-AOG.
- Raw MS file with relevant MS/MS scans for compound are denoted in Column J with information about observed fragments in Column K. Whether a DHA or EPA moiety was observed was denoted in column L and M.
- Whether proper adducts were seen between positive and negative ion mode is denoted in column N and O. Adduct summaries are given in P-S with LOBSTAHS codes. 

### Sample Metadata
- Compiled environmental data collected with each lipid sample can be found in **Holm_Global_Suspended_metadata.csv** or in the similarly labeled tab in **Holm_Combined_Data.xlsx**. Units are noted in the second row of the spreadsheet. Unavailable measurements are marked with "NA"; measurements below the relevant detection limit are marked with "bdl".

#### Additional Data
- TAG quantification data and modeled output from %EPA maps are available in **Holm_TAG_Curve_Compare.csv** and **EPA_projections.csv**

## Code
- Code for lipid annotation screening, quantification and dataset merger, and figure generation are available respectively (**Holm_2022_prepOrbidata_v3.R**, **Holm_2022_quantify.RMD**, **Holm_2022_figure_generation.RMD**)

## Raw Files Access
Raw MS files which are too large to be hosted on GitHub (241.05 GB) are hosted on MetaboLights in the following repository: https://www.ebi.ac.uk/metabolights/MTBLS2838/files

## Associated Publications
Data from this collection has been featured in the following publications. If you used data from this dataset in published works please message me to be added. 

>  1. H. C. Holm, H. F. Fredricks, S. M. Bent, D. P. Lowenstein, J. E. Ossolinski, K. W. Becker, W. M. Johnson, K. Schrage, B. A. S. V. Mooy,
> Global ocean lipidomes show a universal relationship between
> temperature and lipid unsaturation. _Science_ **376**, 1487–1491
> (2022). [10.1126/science.abn7455](doi.org/10.1126/science.abn7455)
> 
>  2. W. Liu, H. C. Holm, J. S. Lipp, H. F. Fredricks, B. A. S. Van Mooy, K.-U. Hinrichs, Unraveling plankton adaptation in global oceans
> through the untargeted analysis of lipidomes. _Science Advances_
> **11**, eads4605 (2025). [10.1126/sciadv.ads4605](doi.org/10.1126/sciadv.ads4605)
>  3. Daniel Lowenstein, Henry Holm, Helen Fredricks et al. Lipid stoichiometry and biomarkers reflect microbial acclimation and
> nutrient stress across the Atlantic Ocean, 19 November 2025, PREPRINT
> (Version 1) available at Research Square
> [10.21203/rs.3.rs-8022751/v1](https://doi.org/10.21203/rs.3.rs-8022751/v1)
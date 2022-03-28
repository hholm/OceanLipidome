# Van Mooy Lab Global Ocean Lipidome Dataset
This repository contains the applicable data and code for the results of the paper "Global ocean lipidomes show a universal relationship between temperature and lipid unsaturation" by Holm et al. 2022. If you have any questions or inquiries about this, please do not hesitate to contact me at hholm@whoi.edu. 

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

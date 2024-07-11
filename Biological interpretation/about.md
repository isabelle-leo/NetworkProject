## Inputs
### RSS data 
- rss_data_full.txt (contains RSS value/type per peptide)

### Annotation files (info_data)
- CellAtlasAnnotationsVsSubCellBarCode.txt (Sub-cellular location)
- DisorderAnnotation.txt (Disorder fractions)
- IsoformAnnotation.txt (Isoform information by gene symbol/Uniprot ID)
- interpro_description_all_meltome 1.txt (Domain description; start/end; by gene symbol)

## Scripts
### RSS analysis
- Fast_comparsion.R (fast and convenient way to generate result tables for stat tests)
  - First, run the RSS analysis script, and assign rss_data_conventional_info (or just run the whole script).
  - Then, go to **line 61**, assign **group_interest** (e.g. disorder_quantile) and **alter_locations** (the information of interest [ioi] e.g. 'Q4').
  - Finally, press *shift-cmd-enter* to run the whole script and get the following results:
   - Kruskal-Wallis rank sum test
   - p-value table of Wilcox test
   - p-value table of Kolmogorov-Smirnov test
   - Information of alter and null group
   - (OPTIONAL) kurtosis value
- disorder_RSS.R
- subcelluar_RSS.R

### Functions
- Plotting_functions.R



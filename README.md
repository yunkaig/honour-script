# honour-script
This repository contains scripts built in-house for my Honours project.

## Tree of files:\
# honour-script/\
# ├── DawnRank\
# │   ├── DawnRankAnalysis.main.R           - main script for DawnRank analysis\
# │   └── dawnrank.visualisation.R          - script for DawnRank visualisation\
# ├── DriverNet\
# │   ├── DriverNetAnalysis.main.R          - main script for DriverNet analysis\
# │   └── drivernet.visualisation.R         - script for DriverNet visualisation\
# ├── GSEA\
# │   └── gsea.analysis.R                   - Downstream analysis of GSEA results\
# ├── general                           - scripts with functions used for more than one type of analysis\
# │   ├── PlottingFunctions.R               - plotting heatmaps and barplots for driver gene prediction\
# │   ├── inf.graph.fix.R                   - checking and fixing influence graph (gene interaction network)\
# │   ├── influence.graph.R                 - building the initial influence graph (gene-gene interaction, which gene ID matching)\
# │   ├── influence.graph.ppi.R             - building protein-protein interaction network (data from Dr Samuel Lee)\
# │   ├── influence.graph.ppi.and.ggi.R     - combining ppi and ggi to get the gene interaction network\
# │   ├── mut.and.exp.matrix.R              - patient data pre-processing to generate mutation and expression matrices\
# │   ├── parsePathways.R                   - parse the gene set file (.gmt) downloaded from MSigDB and helper functions\
# │   └── validation.mut.R                  - sample ID matching and analysis using the larger validation cohort (n = 190)\
# └── singscore\
#     ├── singscore.R                       - main script for singscore analysis\
#     └── singscore.visualisation.R         - visualisation of singscore results\


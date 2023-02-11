# Replication Files for Single Proxy Control (Tchetgen Tchetgen, Park, Richardson, 2023) 

This Github repository contains replication files for [Tchetgen Tchetgen, Park, Richardson (2023)](https://fill.later "SPC").


## Data

The dataset is used in Universal Difference-in-Differences for Causal Inference in Epidemiology [(Tchetgen Tchetgen, Park, Richardson, 2023)](https://arxiv.org/abs/2302.00840 "ORECEPI")
The source of the dataset are given below:
* Pre- and Post-treatment Outcomes, Treatment, and log population: zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY
* log population density and proportion of female: https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados


## Code

* SPC_Analysis.R replicates the main analysis in Section 4 and the sensitivity analysis in the Appendix.
* SPC_Sensitivity_Cluster.R replicates the grid-search to find the sensitivity parameters (recommended to use a cluster computing system).
* SPC_Ft.R contains functions used in SPC_Analysis.R and SPC_Sensitivity_Cluster.R

## References
Tchetgen Tchetgen, Park, Richardson (2023 **Universal Difference-in-Differences for Causal Inference in Epidemiology**, _arXiv:2302.00840_ [[link](https://arxiv.org/abs/2302.00840 "ORECEPI")]
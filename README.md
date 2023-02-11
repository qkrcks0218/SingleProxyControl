# Replication Files for Single Proxy Control (Tchetgen Tchetgen, Park, Richardson, 2023) 

This Github repository contains replication files for [Tchetgen Tchetgen, Park, Richardson (2023)](https://fill.later "SPC").


## Data

The dataset contains birth rate information from 603 municipalities in two states of Brazil, Pernambuco and Rio Grande do Sul. 
Municipality-level birth rates were measured in 2014 and 2016, before and after the 2015 Zika virus outbreak.
More details on the source of the dataset are given below:
* Pre- and Post-treatment Outcomes, Treatment, and log population: zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY. See Taddeo, Amorim, Aquino (2022) (https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "Zika_Brazil_2022") for details.
* log population density and proportion of female: https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados
* The dataset is also used in Universal Difference-in-Differences for Causal Inference in Epidemiology [(Tchetgen Tchetgen, Park, Richardson, 2023)](https://arxiv.org/abs/2302.00840 "ORECEPI")



## Code

* SPC_Analysis.R replicates the main analysis in Section 4 and the sensitivity analysis in the Appendix.
* SPC_Sensitivity_Cluster.R replicates the grid-search to find the sensitivity parameters (recommended to use a cluster computing system).
* SPC_Ft.R contains functions used in SPC_Analysis.R and SPC_Sensitivity_Cluster.R

## References

Taddeo, Amorim, Aquino (2022) **Causal Measures Using Generalized Difference-in-difference Approach with Nonlinear Models**, _Statistics and Its Interface_, 15(4):399-413 [[link] (https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "Zika_Brazil_2022")]


https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns
Tchetgen Tchetgen, Park, Richardson (2023) **Universal Difference-in-Differences for Causal Inference in Epidemiology**, _arXiv:2302.00840_ [[link](https://arxiv.org/abs/2302.00840 "ORECEPI")]


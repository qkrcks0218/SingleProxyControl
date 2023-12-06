# Replication Files for Single Proxy Control (Park, Richardson, Tchetgen Tchetgen, 2023) 

This Github repository contains replication files for [Single Proxy Control (Park, Richardson, Tchetgen Tchetgen, 2023)](https://arxiv.org/abs/2302.06054 "SPC").


## Data

The dataset contains birth rate information from 603 municipalities in two states of Brazil, Pernambuco and Rio Grande do Sul. 
Municipality-level birth rates were measured in 2013, 2014 and 2016, before and after the 2015 Zika virus outbreak.
More details on the source of the dataset are given below:
* Pre- and Post-treatment Outcomes, Treatment, and log population: zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY. 
See [Taddeo, Amorim, Aquino (2022)](https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "ZB") for details.
* log population density and proportion of female: https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados
* The dataset is also used in Universal Difference-in-Differences for Causal Inference in Epidemiology [(Tchetgen Tchetgen, Park, Richardson, 2023)](https://arxiv.org/abs/2302.00840 "ORECEPI")



## Code

* 0.SPC_Ft.R and 0.Functions_MM.R contain functions used in the rest R-files.
* 1.SPC_PARA_CV_Parallel.R replicates estimates obtained from parametric estimators in the Appendix (recommended to use a cluster computing system). 
* 2.SPC_PARA_Summary.R summarizes estimates obtained from 1.SPC_PARA_CV_Parallel.R.
* 3.SPC_PARA_Sensitivity_Parallel.R replicates the sensitivity analysis in the Appendix (recommended to use a cluster computing system). 
* 4.SPC_PARA_Sensitivity_Summary.R summarizes estimates obtained from 3.SPC_PARA_Sensitivity_Parallel.R.
* 5.SPC_NP_Parallel.R replicates estimates obtained from a semiparametric estimator in the main paper (recommended to use a cluster computing system). 
* 6.SPC_Summary.R summarizes parametric and semiparametric estimates and replicates the paper in the main paper.

## References

Taddeo, Amorim, Aquino (2022) **Causal Measures Using Generalized Difference-in-difference Approach with Nonlinear Models**, _Statistics and Its Interface_, 15(4):399-413 [[link](https://www.intlpress.com/site/pub/pages/journals/items/sii/content/vols/0015/0004/a001/index.php?mode=ns "ZB")]

Tchetgen Tchetgen, Park, Richardson (2023) **Universal Difference-in-Differences for Causal Inference in Epidemiology**, _arXiv:2302.00840_ [[link](https://arxiv.org/abs/2302.00840 "ORECEPI")]

Park, Richardson, Tchetgen Tchetgen (2023) **Single Proxy Control**, _arXiv:2302.06054_ [[link](https://arxiv.org/abs/2302.06054 "SPC")]


# mahogany-r
In an interesting paper, [Chimeli & Soares (AER, 2017)](https://www.aeaweb.org/articles?id=10.1257/app.20160055) use difference-in-differences to investigate whether regulations on mahogany extraction increased homicides rates and violence levels in the Brazilian Amazon region. They conclude in favor of a causal link between illegality of mahogany trade and homicides and find the results to be consistent across several robustness checks. The authors shared the database used in their research, as well as the Stata code for their main results (available [here](https://www.openicpsr.org/openicpsr/project/113679/version/V1/view)).

In this repository, I try (and succeed) to reproduce the authors' results in R. I believe this may be useful to those interested in learning to implement diff-in-diff in R. In addition, here you can find some experiments in which I try to implement propensity score estimators using the same data. I also provide my comments on the results (in Portuguese), as submitted to my Professors at the end of a master's course on Policy Evaluation.

Here is a description of all the available files:
* `AEJApp-2016-0055-Stata_Dataset.dta`: database used in Chimeli & Soares (2017), as provided by the authors
* `replicating-results.R`: R script reproducing the paper's main findings
* `ps-matching.R`: R script implementing some propensity score estimators using the mahogany data; the script also checks common support and balance, assesses matching quality and implements diff-in-diff regressions using matched samples
* `df_mahogany.rds`: diff-in-diff-ready database created from the database provided by the authors; this file is generated by `replicating-results.R`
* `df_probit_unmatched.dta`: database with propensity scores for implementing kernel matching in Stata; this file is generated by `ps-matching.R`
* `kernel-matching-stata.do`: Stata script to perform kernel matching using `psmatch2`
* `df_probit_matched.dta`: database with kernel matching results; this file is generated by `kernel-matching-stata.do`
* `ps-helpers.R`: R script with custom functions used in `ps-matching.R`
* `benchmark-models.rData`: `felm` objects containing estimated models corresponding to columns 1-4 of table 2 from Chimeli & Soares (2017); these models are estimated using `replicating-results.R` and are stored so they can be compared to matched-DD models in `ps-matching.R`
* `paper-aps-andre-luis-menegatti.pdf`: my comments on the results, in Portuguese

**Note:** This was my first time doing econometrics in R. Thus, the code could obviously DRYer and more elegant. On the other hand, I think that its current form will be more useful to begginers than a less repetitive version.

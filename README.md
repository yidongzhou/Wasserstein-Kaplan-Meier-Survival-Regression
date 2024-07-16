# Wasserstein-Kaplan-Meier Survival Regression

This repository contains codes necessary to replicate **Zhou and Müller (2024+)**: “Wasserstein-Kaplan-Meier Survival Regression.” The `wkm` function in the `code` folder, short for Wasserstein-Kaplan-Meier Survival Regression, is designed for analyzing survival data from heterogeneous populations. This function requires a data frame (`df`) containing columns for survival time (`time`), censoring indicator (`censor`), and other covariates. The key purpose of this function is to estimate survival functions for different subgroups defined by the covariates. Control options (`optns`) can be specified to adjust the behavior of the estimation process, including setting lower and upper bounds of the support of the measure. The output is a `wkm` object containing quantile functions (`qf`) corresponding to the provided data, the domain grid of these quantile functions (`qfsupp`), the input data frame (`df`), and the specified control options (`optns`).

## Folder Structure

The folder structure of this repo is as follows:

| Folder      | Detail                                                                                                    |
|:------------|:----------------------------------------------------------------------------------------------------------|
| code        | R scripts for the proposed approach                                                                       |
| data        | Data files                                                                                                |
| mdn         | Python implementation for Han et al. (2022 PMLR), adapted from https://github.com/XintianHan/Survival-MDN |
| sim         | R scripts for simulations                                                                                 |

## code

| Data file  | Detail                                       |
|:-----------|:---------------------------------------------|
| wkm.R      | Wasserstein-Kaplan-Meier Survival Regression |
| wkmSim.R   | Simplified version of `wkm` for simulations  |

## data

| Data file              | Detail                                                                |
|:-----------------------|:----------------------------------------------------------------------|
| se1.RData, se2.RData   | Simulation results for subsection 5.2                                 |
| sem.RData, sel.RData   | Simulation results for subsection S.2.1 of the Supplementary Material |
| seb.RData              | Simulation results for subsection S.2.2 of the Supplementary Material |
| sew.RData              | Simulation results for subsection S.2.3 of the Supplementary Material |
| ser1.RData, ser2.RData | Simulation results for subsection S.2.4 of the Supplementary Material |
| seh1.RData, seh2.RData | Simulation results for subsection S.2.5 of the Supplementary Material |

## mdn

| Data file | Detail                                                                                               |
|:----------|:-----------------------------------------------------------------------------------------------------|
| sim.py    | Python scripts to replicate Han et al. (2022 PMLR) in Subsection S.2.3 of the Supplementary Material |

## sim

R scripts to replicate simulation results in subsection 5.2 of the main text and Section S.2 of the Supplementary Material.

| Data file | Detail                                                 |
|:----------|:-------------------------------------------------------|
| sim.R     | Subsection 5.2                                         |
| sim2.R    | Subsection S.2.1 of the Supplementary Material         |
| simb.R    | Subsection S.2.2 of the Supplementary Material         |
| simmdn.R  | Subsection S.2.3 of the Supplementary Material         |
| simrsf.R  | Subsection S.2.4 of the Supplementary Material         |
| simh.R    | Subsection S.2.5 of the Supplementary Material         |
| vis.R     | Scripts to obtain the tables and boxplots in the paper |

## Report Errors

To report errors, please contact <ydzhou@ucdavis.edu>. Comments and suggestions are welcome.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-imbensxu" class="csl-entry">

Zhou, Y. and Müller, H.G., 2023. Wasserstein-Kaplan-Meier Survival Regression.

</div>

</div>



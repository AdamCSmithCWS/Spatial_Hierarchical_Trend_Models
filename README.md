# Spatial_Hierarchical_Trend_Models

Code and data to support [Smith et al. 2023](https://doi.org/10.1093/ornithapp/duad056) describing spatially explicit models for long-term monitoring data of birds in North America.

Adam C Smith, Allison Binley, Lindsay Daly, Brandon PM Edwards, Danielle Ethier, Barbara Frei, David Iles, Timothy D Meehan, Nicole L Michel, Paul A Smith

A pre-print is here: <https://doi.org/10.32942/X2088D>

## Abstract

Population trend estimates form the core of avian conservation assessments in North America and indicate important changes in the state of the natural world. The models used to estimate these trends would be more efficient and informative for conservation if they explicitly considered the spatial locations of the monitoring data. We created spatially explicit versions of some standard status and trend models applied to long-term monitoring data for birds across North America. We compared the spatial models to simpler non-spatial versions of the same models, fitting them to simulated data and real data from three broad-scale monitoring programs: the North American Breeding Bird Survey (BBS), the Christmas Bird Count (CBC), and a collection of programs we refer to as Migrating Shorebird Surveys (MSS). All the models generally reproduced the simulated trends and population trajectories when there were many data, and the spatial models performed better when there were fewer data and in locations where the local trends differed from the range-wide means. When fit to real data, the spatial models revealed interesting spatial patterns in trend, such as increases along the Appalachian Mountains for the Eastern Whip-poor-will, that were much less apparent in results from the non-spatial versions. The spatial models also had higher out-of-sample predictive accuracy than the non-spatial models for a selection of species using BBS data. The spatially explicit sharing of information allows fitting the models with much smaller strata, allowing for finer-grained patterns in trends. Spatially informed trends will facilitate more locally relevant conservation, highlight areas of conservation successes and challenges, and help generate and test hypotheses about the spatially dependent drivers of population change.

## Data

All of the data required to reproduce the analyses and results in this paper are included in the *Data* folder. This includes:

-   saved data objects from the R-package `bbsBayes2` that include the data from the North American Breeding Bird Survey for all species modeled here. Note: additional BBS data from the same 2022-data-release used in this paper can be downloaded using the package `bbsBayes2` , using the following arguments: `bbsBayes2::fetch_bbs_data(release = 2022)`.

-   the Christmas Bird Count data used in the analyses of American Dipper in this paper.

-   the Migrating Shorebird Survey data used in the analyses of Red Knot trends in this paper.

-   the simulated data used to demonstrate the models.

-   the complete list of approximated population trends for CBC and BBS data used in the prior predictive simulations.

## Scripts

There are numbered R-scripts that include all of the code required to reproduce the results in this study. The numbering reflects an approximate sequence because some later scripts require results generated in the earlier scripts.

-   The scripts numbered 1 through 3, create, model, and summarize the simulated data using the models. Figures 1 through 3 are created using these scripts.

-   The scripts numbered 4a through 4c prepare and fit the models to the example data from the BBS, CBC, and MSS surveys.

-   Script 5 confirms the convergence of the models fit to the real data.

-   Script 6 uses results from previous scripts to generate many of the figures in the main paper.

-   Script 7 runs the full cross-validation analyses used in the paper and creates Figure 7.

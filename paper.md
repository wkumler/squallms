---
title: 'squallms: Squashing qualms with speedy quality assurance via lasso labeling for untargeted mass spectrometry data'
tags:
    - R
    - mass spectrometry
    - quality control
    - shiny application
authors:
    - name: William Kumler
      orcid: 0000-0002-5022-8009
      affiliation: 1
      corresponding: true
    - name: Anitra E. Ingalls
      orcid: 0000-0003-1953-7329
      affiliation: 1
affiliations:
    - name: University of Washington School of Oceanography, Seattle, WA 98195, USA
      index: 1
date: '18 July 2024'
bibliography: paper.bib
---

# Summary

Distinguishing between high and low quality mass features remains a major bottleneck in untargeted mass spectrometry. Here, we present the R package squallms that streamlines the process of chromatographic feature annotation by extracting useful metrics from the raw data, labeling groups of similar features interactively, and building a logistic classification model. The package is available via Bioconductor (https://bioconductor.org/packages/devel/bioc/html/squallms.html) and interfaces with other mass spectrometry packages there such as XCMS [@Smith2006]. We expect this functionality to significantly reduce the effort required to explore and curate increasingly common large-scale chromatography-based mass spectrometry datasets.

# Statement of need

Chromatographic mass spectrometry is a powerful way to characterize the molecular composition of chemical and biological samples. The development of data-driven algorithms for untargeted mass spectrometry has massively increased the amount of information obtainable from these datasets, but these algorithms are still plagued by false positive (noise) features [@Gika2019; @Pirttila2022; @Gloaguen2022]. Though novel metrics are actively being developed to quantify the quality of features extracted by untargeted software [@Kantz2019; @Kumler2023], they are difficult to translate across datasets due to differences in instrument setup and sensitivity. Additionally, removing features via sequential thresholding fails to make use of the multivariate advantage and requires the use of thresholds set arbitrarily by the user [@Olivieri2008].

A potential alternative is the use of a multivariate logistic model for the annotation of chromatographic feature quality. However, creating labeled datasets for training such models is often dull and time consuming because it requires the manual review of many individual chromatograms by an expert. Instead, we propose that many similar features can be annotated simultaneously if organized in such a way that an expert can denote an entire group of features as good or bad with a single interaction. These annotations can then be used to train a logistic regression model so that an entire dataset can be quality annotated with confidence quickly and easily [@Kantz2019].

One way to organize features in such a way that similar ones are near to each other is by interpolating a multi-file extracted ion chromatogram to a shared set of retention times and treating each combination of file and retention time as a dimension in a principal component analysis (\autoref{fig:EIC}).

![Coercion of a typical extracted ion chromatogram (upper plot) to a file-by-time matrix. A single high-quality chromatographic feature is shown in all six samples after normalization both in "side view" with intensity on the y axis (upper plot) and "top down" with intensity corresponding to fill color (lower plot).\label{fig:EIC}](joss_fig1.png)

Features that are visually similar tend to cluster tightly in the first few principal components while noise features scatter randomly in this space (\autoref{fig:PCA}).

![Subset of chromatographic features shown both individually in top down view and in principal component space. The upper plot shows eighteen chromatographic features and their associated ID beginning with the letters FT. Features 134, 137, 138, 141, and 144 all appear to be high quality and look very similar with a bright green central stripe and dark purple edges. The two final small multiples (PC1 and PC2) show the first and second principal components for comparison. The lower plot shows how the eighteen features fall in principal component space, with the five high-quality features noted above all clustering to the left side of the plot, separately from the other features.\label{fig:PCA}](joss_fig2.png)

However, the location of this high-quality cluster in principal component space is not easily predictable because feature quality is not necessarily the largest source of variance in all datasets. This necessitated the development of an interactive feature selection method that we have implemented in the squallms package using R's shiny and plotly libraries [@Chang2024; @Sievert2020]. The main interface is shown below in \autoref{fig:ss} for all features, not just the subset used in the demo above.

![Interactive group labeling dashboard of squallms.\label{fig:ss}](joss_fig3.png)

The package also contains a manual training tool, also developed using the shiny package, for rapid single feature annotation by binding quality assessments to specific keys for ease of annotation. Finally, the package includes a basic logistic regression model designed to accept the feature labels and feature quality metrics which produces a best-fit estimate of each feature's quality as not all features are likely to be labeled using either the aggregate or the manual methods.

# Acknowledgements

This work was supported by the University of Washington eScience Institute through their Data Science Incubator program. We would like to thank Bryna Hazelton and Valentina Staneva for thoughtful discussion on algorithm design and package development. Funding was provided by the Simons Foundation (SCOPE Award ID 329108 to AEI, SF Award ID 385428 to AEI).

# References



# Project Description

The content of this repository provides an excerpt from a larger scientific work in the field of molecular structural biochemistry focusing on a protein-tRNA complex responsible for the formation of a specific post-transcriptional transfer RNA modification. Here, the focus is on the statistical modeling of this interaction.

<p>
    <img src="/pictures/project2.png" alt="intro">
    <em>Figure 1.1 Project at a glance: The stoichiometry of a ThiI-tRNA complex is investigated in order to elucidate whether dimeric ThiI binds one or two tRNA substrates simultaneously. This is accomplished by fitting a set of five binding models all derived from appropriate underlying 1:1 or 1:2 binding equilibria to a fluorescence spectroscopic titration profile.</em>
</p>


### Research Background
ThiI is a transfer RNA modifying enzyme responsible for the post-transcriptional sulfur incorporation into uridine at position eight of cytosolic tRNAs. Such a ThiI-tRNA complex is investigated in this work aimed at unraveling its true underlying stoichiometry.

### Methods

A tryptophan fluorescence spectroscopic titration experiment was carried out. An analytic work flow was established including ordinary non linear least squares regression, parameter identifiability analysis, model evaluation procedures such as KL divergence analysis and cross validation as well as generic and robust ANOVA approaches. On the basis of this work flow it is possible to objectively infer the stoichiometry from ThiI-tRNA binding profiles by statistical comparison of a set of variably parameterized binding models all derived from the appropriate underlying 1:1 (ThiI:tRNA) or 1:2 binding equilibria (fig. 1.1).

### Major Results

There is a cumulative weight of evidence of 96.8% that the statistically best approximating binding model is among the set describing binding stoichiometries of 1:2. This finding is supported by contrast analysis in the context of a generic ANOVA approach revealing that on average the one-site binding model yields statistically higher cross validation scores than the two-site binding models.

### Conclusions

Overall results from the complete analysis workflow provide evidence that the underlying stoichiometry of complexes between ThiI and its full-length tRNA substrates is indeed a 1:2 stoichiometry. However, it has to be stressed that the results of this work are just first insights based on a preliminary experimental set up and data set. Ultimately, more biological replications are going to be needed based on an improved experimental set up concerning the applied concentration range, and the number of data points. This improvement is going to allow for a more robust direct comparison of binding models but also for a meta-analytic statistical analysis of CV-score distributions. 

## Analysis Workflow
The complete analysis workflow of the project has already been published more detailed while giving more introductory background as well as greater discussion sections on [RPubs](https://rpubs.com/DeTwes/modeling_thii_trna_interaction). 

However, in an attempt to additionally give an overview in a concise way, each of the six analysis steps (see fig. 1.2) were allocated to individual Rmarkdown-files which are deposited in this repository. Also, to each Rmd-file belongs a static website published on Rpubs. In what follows are  short descriptions of the Rmd-file contents.

<p>
    <img src="/pictures/analysis_overview.png" alt="intro">
    <em>Figure 1.2 Overview of the analysis workflow: Binding models were fit by ordinary non-linear least squares regression and validated by maximum likelihood profiling as part of a parameter identifiability analysis. Binding models were then evaluated based on KL divergence analysis taking into account model complexity and by cross-validation taking into account the predictive power. The distribution of cross-validation scores were further investigated by means of ANOVA and robust ANOVA.</em>
</p>

### 01-nls_regression.Rmd ([Rpubs1.1](https://rpubs.com/DeTwes/NLS-Regression))

This Rmd-file covers the measurement set up of the tryptophan fluoresence spectroscopic titration experiment as well as  data collection and data structures. In addition, the set of binding models is described and the non linear least squares regression approach to fit the ThiI-tRNA interaction profile. 
### 02-ml_profiling.Rmd ([Rpubs1.2](https://rpubs.com/DeTwes/Identifiability-Analysis))
This Rmd-file represents the parameter identifiability analysis. This likelihood approach was used in a qualitative way in order to validate whether the currently applied measurement set up and quality, respectively, as well as the respective binding model parameterizations are appropriate for the sets of parameter to be uniquely identifiable with finite precision. 

### 03-modelSel_AIC.Rmd ([Rpubs1.3](https://rpubs.com/DeTwes/KL-Divergence))

In this Rmd-file the first binding model evaluation step is presented that is based on the Akaike information Criterion (AIC). The AIC efficiently trades off goodness-of-fit as quantified by the RMSE of the residuals and model complexity as quantified by the number of parameter. Based on the AIC the relative strength of support for each of the binding models in the set could be assessed.

### 04-modelSel_CV.Rmd ([Rpubs1.4](https://rpubs.com/DeTwes/completeCV))

In this Rmd-file the second binding model evaluation step is presented that is based on complete 4-fold cross validation (CV). CV offers the possibility to simulate future data sets. In this way, dependencies of binding model performances on variations in the data set as would inevitably be caused by future titration experiments can be uncovered.

### 05-CV_ANOVA.Rmd ([Rpubs1.5](https://rpubs.com/DeTwes/CV-ANOVA))

Model evaluation based on a complete 4-fold CV generated CV-score distributions for each binding model. This allowed for addition of a second dimension to this analysis by further investigating the central tendencies of those distributions by means of analysis of variance (ANOVA) which is presented in this Rmd-file.

### 06-CV_robANOVA.Rmd ([Rpubs1.6](https://rpubs.com/DeTwes/robANOVA))

In this Rmd-file a robust ANOVA approach is covered intended to validate the outcomes of the generic ANOVA approach. Specifically, the ANOVA model was made robust against the assumption of normally distributed residuals by using trimmed CV-score distributions. This file also contains a summary of the whole analysis workflow.
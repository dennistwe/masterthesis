# Project Description

The content of this repository provides an excerpt from a larger scientific work in the field of molecular structural biochemistry focusing on a protein-tRNA complex responsible for the formation of a specific post-transcriptional transfer RNA modification.

### Research Background
ThiI is a transfer RNA modifying enzyme responsible for the post-transcriptional sulfur incorporation into uridine at position eight of cytosolic tRNAs. Such a ThiI-tRNA complex is investigated in this work aimed at unraveling its true underlying stoichiometry.

### Methods

A tryptophan fluorescence spectroscopic titration experiment was carried out. An analytic work flow was established including ordinary non linear least squares regression, parameter identifiability analysis, model evaluation procedures such as KL divergence analysis and cross validation as well as generic and robust ANOVA approaches. On the basis of this work flow it is possible to objectively infer the stoichiometry from ThiI-tRNA binding profiles by statistical comparison of a set of variably parameterized binding models all derived from the appropriate underlying 1:1 (ThiI:tRNA) or 1:2 binding equilibria.

### Major Results

### Conclusions
<p>
    <img src="/pictures/project2.png" alt="intro">
    <em>**Figure 1.1|Project at a glance:** The stoichiometry of a ThiI-tRNA complex is investigated in order to elucidate whether dimeric ThiI binds one or two tRNA substrates simultaneously. This is accomplished by fitting a set of five binding models all derived from appropriate underlying 1:1 or 1:2 binding equilibria to a fluorescence spectroscopic titration profile.</em>
</p>
Normal text would continue here...Normal text would continue here...Normal text would continue here...Normal text would continue here...

## Analysis Workflow
Each individual Rmd file covers one particular analysis step as depicted in figure 1.2. Also, to each Rmd file belongs a static website deposited on Rpubs.

<p>
    <img src="/pictures/analysis_overview.png" alt="intro">
    <em>Shown is a glimpse into the whole analysis workflow. </em>
</p>

### 01-nls_regression.Rmd ([Rpubs1](https://rpubs.com/DeTwes))
In dieser file werden zum einen die Messmethodik und der erzeugte Datensatz vorgestellt. Zum anderen wird der Satz von 5 Bindungsmodellen und deren Zusammenhänge erklärt. Anschließend wird eine nicht-lineare regressionsanalyse durchgeführt, welche initial die beschreibende Stärke der Modelle beleuchtet.

### 02-ml_profiling.Rmd
Diese File beinhaltet die Validierung der Bindungsmodelle mittels Parameter Identifiability analysis. Diese maximum likelihood methodik wurde in einem qualitativen ansatz verwendet um zu ermitteln ob der Informationsgehalt in dem Bindungsprofil sowie die Parameterisierung der Modelle es erlauben, dass modellspezifische Params uniquely and with finite precision identifizierbar sind. 

### 03-modelSel_AIC.Rmd
In dieser File wird die Bindungsmodell-Evaluierung mittels der sog. Kullback-Leibler Divergence Analyse durchgeführt. Dabei wird die informationsth. KL-Distanz durch das Akaike Informationskriterium (AIC) approximiert.

### 04-modelSel_CV.Rmd

## Repository Organization

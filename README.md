# Project Description
<p>
    <img src="/pictures/project2.png" alt="intro">
    <em>image_caption</em>
</p>
Normal text would continue here...Normal text would continue here...Normal text would continue here...Normal text would continue here...

# Analysis Workflow
Each individual Rmd file covers one particular analysis step as depicted in figure 1.2. Also, to each Rmd file belongs a static website deposited on Rpubs.

<p>
    <img src="/pictures/analysis_overview.png" alt="intro">
    <em>Shown is a glimpse into the whole analysis workflow. </em>
</p>

## 01-nls_regression.Rmd ([Rpubs1](https://rpubs.com/DeTwes))
In dieser file werden zum einen die Messmethodik und der erzeugte Datensatz vorgestellt. Zum anderen wird der Satz von 5 Bindungsmodellen und deren Zusammenhänge erklärt. Anschließend wird eine nicht-lineare regressionsanalyse durchgeführt, welche initial die beschreibende Stärke der Modelle beleuchtet.

## 02-ml_profiling.Rmd
Diese File beinhaltet die Validierung der Bindungsmodelle mittels Parameter Identifiability analysis. Diese maximum likelihood methodik wurde in einem qualitativen ansatz verwendet um zu ermitteln ob der Informationsgehalt in dem Bindungsprofil sowie die Parameterisierung der Modelle es erlauben, dass modellspezifische Params uniquely and with finite precision identifizierbar sind. 

## 03-modelSel_AIC.Rmd
In dieser File wird die Bindungsmodell-Evaluierung mittels der sog. Kullback-Leibler Divergence Analyse durchgeführt. Dabei wird die informationsth. KL-Distanz durch das Akaike Informationskriterium (AIC) approximiert.

## 04-modelSel_CV.Rmd

# Repository Organization

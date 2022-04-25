# gMCStool

<div align="center">
    <img src="www/logo_gMCStool_2.png" alt="Welcome" style="width: 75%"/>
    <p> #### *"A tool for discovering essential genes in cancer metabolism"* </p>
</div>

<br>

## Welcome!

We present gMCStool, a user-friedly web-tool to predict essential genes in the latest reconstruction of the human metabolism, [Human1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7331181/), freely available in [GitHub](https://github.com/SysBioChalmers/Human-GEM). 
We have chosen the game 'jenga' as the icon of our tool, because the jenga tower falls when you take out the most vulnerable part of the structure. In our case, we are searching for the most vulnerable parts of the metabolism of cancer cells in order to disrupt cellular proliferation. 
In order to fin metabolic vulnerabilities in cancer cells, we employ the concept of genetic Minimal Cut Sets ([gMCSs](https://academic.oup.com/bioinformatics/article/35/3/535/5056753)), a metabolic network-based approach to synthetic lethality, and RNA-seq data. 
<br>

gMCStool has been developed to achieve the following analysis:
<ol>
    <li>Find essential genes for cancer metabolism based on RNA-seq data.</li>
    <li>Predict putative companion biomarkers for predicted essential genes.</li>
    <li>Identify essential tasks or essential metabolites associated with predicted essential genes.</li>
</ol>

## Run it!
gMCStool has been developed using R and Shiny. You can clone the repository in your own computer and run locally simply by typing:

```{r eval=F, echo=T}
shiny::runApp('./app.R')
```

The code is prepared to adjust the amount of available RAM and number of cores, based on the information of the PC. Furthermore, there are two variables to activate real-time tables in the tabs 4 and 5. A full tutorial is available in the help tab or in [TouTube](https://www.youtube.com/watch?v=k9OM2lhDHEU).


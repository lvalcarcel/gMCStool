# gMCStool

<div align="center">
    <img src="www/logo_gMCStool_2.png" alt="Welcome" style="width: 75%"/>
    <p> <b> <i> A tool for discovering essential genes in cancer metabolism </i> </b> </p>
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


gMCStool has been developed using R and Shiny. It can be launched either on the public web application of gMCStool, on a local machine using R and RShiny or within Docker. 

### Web application

link: [gMCStool](http://biotecnun.unav.es:8080/)

### Local

Clone the repository

```git
git clone https://gitlab.com/lvalcarcel/gmcstool.git
```

Run RShiny

```r
shiny::runApp('app.R')
```

The code is prepared to adjust the amount of available RAM and number of cores, based on the information of the PC. Furthermore, there are two variables to activate real-time tables in the tabs 4 and 5. 

#### Run gMCStool in a Docker App

Run the following command:

```git
sudo docker run -p 3838:3838 lvalcarcel/gmcstool R -e 'shiny::runApp("/root/gMCStool")'
```

Now, you can point your browser to [http://localhost:3838](http://localhost:3838) and use gMCStool!



## Quick start

There is a complete manual of [PDF](https://github.com/lvalcarcel/gMCStool/blob/385356bd9d6aee0cc6cc23c96431155729ba36c4/instructions/gMCStool_help_tab.pdf), in the help tab of the tool or in or in [YouTube](https://www.youtube.com/watch?v=k9OM2lhDHEU).

<br>


## Cite Us


*Luis V. Valcárcel, Edurne San José-Enériz, Raquel Ordoñez, Iñigo Apaolaza, Ana Valcárcel, Leire Garate, Jesús San Miguel, Antonio Pineda-Lucena, Xabier Agirre, Felipe Prósper, Francisco J. Planes.. "gMCStool App"*


Companion code to the paper: \
**_Tumour lineage shapes BRCA-mediated phenotypes_**

Code in this repository can be used to recreate the essential parts of the main figures in the manuscript. All code is in R.

### Instructions
Install required R packages:
```r
install.packages(c('knitr', 'plyr', 'tidyverse', 'patchwork', 'binom', 'ggrepel', 'grid', 'readxl', 'ggbeeswarm', 'here', 'survival', 'survminer'))
```

Get supplementary data from manuscript and download germline data from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001858.v1.p1).

Place these in the `data` folder, with the names `Taylor_SItables.xlsx`, `germline_mutations.maf`, and `germline_cnvs.txt`, or manually enter the file paths in the `prerequisites.R` script in the `r` directory.

An HTML file containing the figures can be generated from the command-line for each main figure as such:
```shell
R -e "rmarkdown::render('r/figure-1.Rmd', output_file = 'figure-1.html')"
```

### Citation
- URL: https://www.nature.com/articles/s41586-019-1382-1
- DOI: 10.1038/s41586-019-1382-1

### Contact
E-mail any questions to [jonssonp@mskcc.org](mailto:jonssonp@mskcc.org?subject=[GitHub]%20BRCA%20paper).

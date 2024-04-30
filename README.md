# A genome-wide DNA methylation study integrated with gene expression data reveals novel insights in the epigenetic landscape of Systemic Sclerosis

To unravel DNA methylation abnormalities associated with Systemic Sclerosis (SSc) we performed an epigenome-wide association study (EWAS) using whole blood samples from 179 SSc patients and 241 unaffected individuals. This analysis yielded 525 differentially methylated positions (DMPs), enriched in immune-related pathways, highlighting integrins as a potential therapeutic option for SSc. Analysis of transcription factors revealed the relevance of the myeloid CEBP family in the epigenetic and transcriptomic drift of SSc. Integration with gene expression data from the same individuals revealed 842 significant correlations between DMPs and gene expression levels, highlighting neutrophils as a key cell type in SSc pathogenesis. Overall, this study provides a comprehensive overview of epigenetic changes in SSc and their consequences on gene expression, underscoring the pivotal role of myeloid cells in the pathogenesis of the disease and uncovering new potential biomarkers and therapeutic options.

## EWAS pipeline
Below is the necessary code to conduct the different analyses performed in this study, divided into different sections, following the order of the manuscript. 
### Principal component analysis

```

```
### Type 3 ANOVA analysis for covariates

```

```
### Infering cell composition from EpiDish R package
Since it is a study of whole blood, changes in the cell composition could lead to false positives, so it is important to consider potential differences in cell counts. Thus, estimated blood cells proportions were inferred from the methylation beta values using [epiDish](https://bioconductor.org/packages/release/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html) R package. This analysis was carried out through the robust partial correlations (RPC) method with 100,000 iterations using [Reinius et al (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3405143/) as the reference panel.
```
library(EpiDish)
data(centDHSbloodDMC.m)

B <- lumi::m2beta(Mvalues_merge)

cell_counts <- epidish(B, ref.m = centDHSbloodDMC.m, method = "RPC", maxit = 100000)
cells <- cell_counts$estF
cells <- cells[Data_SE$Sample_Name,]
table(row.names(cells)==Data_SE$Sample_Name)


DATA <- cbind(Data_SE, cells)
```
### Remove CpGs associated with sex, age, and cell composition

```

```
### Differential analysis with DECO

```

```
### DMPs location analysis

```

```
### Epigenetic clock analysis

```

```
### Transcription factor activity analysis

```

```
### Differential expression analysis and integration with DNA methylation data

```

```
### Gene ontology enrichment of eQTMs

```

```
### Serum protein cytokines analysis

```

```

**NOTE:** [HOMER](http://homer.ucsd.edu/homer/) and [GREAT](http://great.stanford.edu/public/html/) analysis are not detailed here as they were performed through their corresponding APIs.  

To cite this paper:

> Martinez-Lopez, J., Estupiñan-Moreno, E., Ortiz-Fernandez, L. et al., A genome-wide DNA methylation study integrated with gene expression data reveals novel insights in the epigenetic landscape of Systemic Sclerosis, _JOURNAL TBD_, 2024, DOI:

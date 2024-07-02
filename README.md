# A genome-wide DNA methylation study integrated with gene expression data reveals novel insights in the epigenetic landscape of Systemic Sclerosis

To unravel DNA methylation abnormalities associated with Systemic Sclerosis (SSc) we performed an epigenome-wide association study (EWAS) using whole blood samples from 179 SSc patients and 241 unaffected individuals. This analysis yielded 525 differentially methylated positions (DMPs), enriched in immune-related pathways, highlighting integrins as a potential therapeutic option for SSc. Analysis of transcription factors revealed the relevance of the myeloid CEBP family in the epigenetic and transcriptomic drift of SSc. Integration with gene expression data from the same individuals revealed 842 significant correlations between DMPs and gene expression levels, highlighting neutrophils as a key cell type in SSc pathogenesis. Overall, this study provides a comprehensive overview of epigenetic changes in SSc and their consequences on gene expression, underscoring the pivotal role of myeloid cells in the pathogenesis of the disease and uncovering new potential biomarkers and therapeutic options.

## EWAS pipeline
Below is the necessary code to conduct the different analyses performed in this study, divided into different sections, following the order of the manuscript. 
### Principal component analysis (PCA)
To identify undesired biasing effects, a principal component analysis (PCA) was carried out. Samples deviating 4 standard deviations from the cluster centroid were removed from further analysis. 
```
M <-as.matrix(M_SE)

B <- lumi::m2beta(M)

B[B==0] <- 0.0001
B[B> 0.98999999] <- 0.989
B <- na.omit(B)

M_PCA <- t(lumi::beta2m(B))

PCA <- prcomp(x=M_PCA, retx = T, center = T, scale. = T, rank.=10)

summary(PCA)

#CALCULATE OUTLIERS
PCs <- data.frame(PCA["x"])
#Cuadramos ambos DataFrames
PCs <- PCs[DATA_SE$Sample_Name,]

outliers<-apply(PCs, 2, function(x) which(abs(x - mean(x)) > (4 * sd(x)) ))

table(row.names(PCs) == DATA_SE$Sample_Name)
PCA <- cbind(DATA_SE, PCs)
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
### Type 3 ANOVA analysis for covariates and remove CpGs associated with sex, age, and cell composition
In order to avoid any potential bias, the effect of sex, age and cell composition was calculated through a type III Anova test with R. Additionally, CpGs associated with sex, age and cell composition were estimated through regression analysis with [limma](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/intro.html). The resulting significant CpGs from this analysis with an estimated false discovery rate (FDR<0.05) were removed from further analysis.
```
#ANOVA
library(car)
DNA_TOT$Sex <- as.factor(DNA_TOT$Sex)
DNA_TOT$Disease2 <- as.numeric(factor(DNA_TOT$Disease, levels=c('SSc', 'CTRL'), labels=c(1,0)))
COV <- Anova(lm(Disease2~Sex+Age+B+NK+CD4T+CD8T+Mono+Neutro, data=DNA_TOT), type='III') 

#limma
Gender <- factor(DATA_SE$Sex)
Age <- as.numeric(DATA_SE$Age)
Cells <- data.frame(DATA_SE[,c(13,14,15,16,17,18)])

M_SE <- M_SE[,DATA_SE$Sample_Name] 
table(DATA_SE$Sample_Name==colnames(M_SE))

mod <-model.matrix(~0+Gender+Age+Cells$Neutro) #We used Neutrophils as covariate as cell types proportions correlate between them and it is the main cell type in whole blood
colnames(mod) <- make.names(colnames(mod))
contrast_matrix <- makeContrasts(Comp="GenderF-GenderM", levels=mod)

ft <-  lmFit(M_SE, mod)
ft2 <- contrasts.fit(ft, contrast_matrix)
ft3 <-eBayes(ft2)
summary(decideTests(ft3))

SEGE <- topTable(ft3, n=nrow(M_SE), adjust.method = "BH")
```
### Differential analysis with DECO
[DECO](https://github.com/fjcamlab/deco), a bioinformatic tool specifically designed for analyzing heterogeneous cohorts, was employed to identify differentially methylated positions (DMPs). Through a randomized subsampling approach and an established false discovery rate (FDR) threshold of less than 0.05, this tool performs an iterative differential analysis with 100 iterations. Differential analysis conducted in each iteration is executed using limma, and the resulting p-values from all iterations are combined following a Chi-square distribution with degrees of freedom equal to twice the number of iterations. This tool also labels the resulting significant CpGs depending on the overlapping of the methylation signal across SSc patients and HD cohorts. Hence, DECO identifies if these differentially methylated CpGs are impaired in both groups or they affect only a subset of patients or the whole of one group.
```
bpparam <- MulticoreParam(8, RNGseed = 7739465, progressbar = TRUE) #Para hacerlo en memoria compartida

table(Data_SE$Sample_Name==colnames(M_DECO))

Data_DECO <- Data_SE$Disease
names(Data_DECO) <- Data_SE$Sample_Name

deco1 <- decoRDA(data=M_DECO, classes=Data_DECO, iterations = 100, r=30, q.val = 0.01, bpparam = bpparam)

table(is.infinite(c(deco1$data, deco1$results)))
na.omit(deco1)

deco2 <-  decoNSCA(sub = deco1, v=50, method=NULL, bpparam = bpparam, k.control=NULL, k.case = NULL, samp.perc = 0.05, rep.thr = 1) 

hMatrixSSc <- NSCAcluster(deco2)$Case$NSCA$h
hMatrixCTRL <- NSCAcluster(deco2)$Control$NSCA$h

sig_CPGs_PRE <- deco2@featureTable[,c(1,2,3,8,9,11,12,23)]


sig_CPGs <- sig_CPGs_PRE[!sig_CPGs_PRE$ID%in%row.names(SEGE[SEGE$adj.P.Val<0.05,]),] #We remove CpGs associated with covariates
DMPs <- sig_CPGs[sig_CPGs$Profile=="Minority",]
```
### DMPs annotation
DMPs, and CpGs included in the differential analysis, were annotated using the illumina annotations for the hg19 genomic assembly.
```
require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") #Get hg19 annotations

annot_CPGs <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot_CPGs <- annot_CPGs[match(row.names(sig_CPGs),annot_CPGs$Name),]
table(row.names(sig_CPGs)==row.names(annot_CPGs))

annotated_CPGs <- data.frame(cbind(sig_CPGs,annot_CPGs))
```
### Epigenetic clock analysis
EAA was inferred through the online [DNA Methylation Age Calculator](https://dnamage.clockfoundation.org/), utilizing DNA methylation B values as input. EAA was determined as the residuals from the regression between the actual age of the patient or CTRL and the epigenetic age calculated from Hannum or Horvath reference models. The quantiles from the IEAA of both clocks were used to categorize individuals into accelerated and normal aging groups. Differences in the frequencies of individuals with accelerated epigenetic aging in each SSc and CTRL group were evaluated through a Chi-square test, accompanied by odds ratio (OR) calculation.
```
#SEE IF EA IS WELL CALCULATED
HANNUM <- cor.test(DNA_TOT$Age, DNA_TOT$DNAmAgeHannum, method = "spearman")
HORVATH <- cor.test(DNA_TOT$Age, DNA_TOT$DNAmAge, method = "spearman")

HANNUM
HORVATH

summary(DNA_TOT$IEAA) #We observe the quartiles by focusing on the intrinsic epigenetic accelerated age of Horvath epigenetic clock
#Those with a difference greater than 2.5506 years are epigenetically aging faster

summary(DNA_TOT$IEAA.Hannum) #We observe the quartiles by focusing on the intrinsic epigenetic accelerated age of Hannum epigenetic clock
#Those with a difference greater than 1.82928 years are epigenetically aging faster


#This loop classifies patients in different groups depending wether they are normally aging or rapidly aging in terms of Epigenetic Age
i <- 1
for(i in 1:nrow(DNA_TOT)){
  if (DNA_TOT[i, 'AgeAccelerationResidualHannum']>2.5506){
    DNA_TOT[i, 'CategoryHannum'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'AgeAccelerationResidualHannum']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryHannum'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'AgeAccelerationResidual']>2.5506){
    DNA_TOT[i, 'CategoryResiduals'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'AgeAccelerationResidual']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryResiduals'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'DNAmTLAdjAge']>2.5506){
    DNA_TOT[i, 'CategoryTLs'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'DNAmTLAdjAge']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryTLs'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'IEAA']>2.5506){
    DNA_TOT[i, 'CategoryIEAA'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'IEAA']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryIEAA'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'IEAA.Hannum']>2.5506){
    DNA_TOT[i, 'CategoryIEAAHannum'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'IEAA.Hannum']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryIEAAHannum'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'EEAA']>2.5506){
    DNA_TOT[i, 'CategoryEEAA'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'EEAA']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryEEAA'] <- 'Decelerated'
  #}
  if (DNA_TOT[i, 'AgeAccelPheno']>2.5506){
    DNA_TOT[i, 'CategoryEAAPheno'] <- 'Accelerated'
  }
  #if (DNA_TOT[i, 'AgeAccelPheno']<(-1.84981)){
  #  DNA_TOT[i, 'CategoryEAAPheno'] <- 'Decelerated'
  #}
} #You have to repeat the table for Hannum threshold


Hannum <- chisq.test(table(DNA_TOT$CategoryHannum, DNA_TOT$Disease))
Hannum
Residuals <- chisq.test(table(DNA_TOT$CategoryResiduals, DNA_TOT$Disease))
Residuals
IEAA <- chisq.test(table(DNA_TOT$CategoryIEAA, DNA_TOT$Disease))
IEAA
IEAAHannum <- chisq.test(table(DNA_TOT$CategoryIEAAHannum, DNA_TOT$Disease))
IEAAHannum
EEAA <- chisq.test(table(DNA_TOT$CategoryEEAA, DNA_TOT$Disease))
EEAA

#Odds ratio calculation
Hannum$OR <- (Hannum$observed[1,2]*(Hannum$observed[2,1]))/(Hannum$observed[1,1]*(Hannum$observed[2,2]))
Residuals$OR <- (Residuals$observed[1,2]*(Residuals$observed[2,1]))/(Residuals$observed[1,1]*(Residuals$observed[2,2]))
IEAA$OR <- (IEAA$observed[1,2]*(IEAA$observed[2,1]))/(IEAA$observed[1,1]*(IEAA$observed[2,2]))
IEAAHannum$OR <- (IEAAHannum$observed[1,2]*(IEAAHannum$observed[2,1]))/(IEAAHannum$observed[1,1]*(IEAAHannum$observed[2,2]))
EEAA$OR <- (EEAA$observed[1,2]*(EEAA$observed[2,1]))/(EEAA$observed[1,1]*(EEAA$observed[2,2]))

#OR confidence intervals calculations
#CI 95%(OR)=OR×exp(±z×sqrt(1/a+1/b+1/c+1/d)
#Being z the z value for a qnorm distribution for 0.975 (0.025 per each tail)

z_value <- qnorm(0.975) # For a 95% confidence interval
#CI95 - LOWER
Hannum$ORCI[1] <- Hannum$OR*exp(-z_value*sqrt(1/Hannum$observed[1,2]+1/(Hannum$observed[2,1])+1/Hannum$observed[1,1]+1/(Hannum$observed[2,2])))
Residuals$ORCI[1] <- Residuals$OR*exp(-z_value*sqrt(1/Residuals$observed[1,2]+1/(Residuals$observed[2,1])+1/Residuals$observed[1,1]+1/(Residuals$observed[2,2])))
IEAA$ORCI[1] <- IEAA$OR*exp(-z_value*sqrt(1/IEAA$observed[1,2]+1/(IEAA$observed[2,1])+1/IEAA$observed[1,1]+1/(IEAA$observed[2,2])))
IEAAHannum$ORCI[1] <- IEAAHannum$OR*exp(-z_value*sqrt(1/IEAAHannum$observed[1,2]+1/(IEAAHannum$observed[2,1])+1/IEAAHannum$observed[1,1]+1/(IEAAHannum$observed[2,2])))
EEAA$ORCI[1] <- EEAA$OR*exp(-z_value*sqrt(1/EEAA$observed[1,2]+1/(EEAA$observed[2,1])+1/EEAA$observed[1,1]+1/(EEAA$observed[2,2])))

#CI95 - UPPER
Hannum$ORCI[2] <- Hannum$OR*exp(z_value*sqrt(1/Hannum$observed[1,2]+1/(Hannum$observed[2,1])+1/Hannum$observed[1,1]+1/(Hannum$observed[2,2])))
Residuals$ORCI[2] <- Residuals$OR*exp(z_value*sqrt(1/Residuals$observed[1,2]+1/(Residuals$observed[2,1])+1/Residuals$observed[1,1]+1/(Residuals$observed[2,2])))
IEAA$ORCI[2] <- IEAA$OR*exp(z_value*sqrt(1/IEAA$observed[1,2]+1/(IEAA$observed[2,1])+1/IEAA$observed[1,1]+1/(IEAA$observed[2,2])))
IEAAHannum$ORCI[2] <- IEAAHannum$OR*exp(z_value*sqrt(1/IEAAHannum$observed[1,2]+1/(IEAAHannum$observed[2,1])+1/IEAAHannum$observed[1,1]+1/(IEAAHannum$observed[2,2])))
EEAA$ORCI[2] <- EEAA$OR*exp(z_value*sqrt(1/EEAA$observed[1,2]+1/(EEAA$observed[2,1])+1/EEAA$observed[1,1]+1/(EEAA$observed[2,2])))
```
### Differential expression analysis
Individuals deviating>4 standard deviation from cluster centroid in PCA analysis or that presented atypical gene expression counts in comparison to the mean values of the remaining samples were excluded from further analysis. Low expressed genes were removed (CPM<1) and the remaining genes were normalized by the trimmed mean of M-values (TMM) method. Differential expression analyses were carried out with a quasi-likelihood negative binomial generalized log-linear model implemented in  the [edgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/intro.html) package. All genes with a FDR<0.05 were considered statistically significant and defined as differentially expressed genes (DEGs).Sex, age and cell composition were included as covariates in the analysis.
```
#PCA
PCA <- prcomp(t(COUNTS_SE), rank. = 10)
plot(PCA$x)
outliers<-apply(PCA$x, 2, function(x) which(abs(x - mean(x)) > (4 * sd(x)) ))
outliers <- unlist(outliers)
toremove <- c()
for(i in 1:length(outliers)){toremove[i] <- outliers[i]}
toremove

#Remove outliers
COUNTS_SE <- COUNTS_SE[,-toremove]
DESIGN_SE <- DESIGN_SE[-toremove,]


###Diferential expression

Disease <- DESIGN_SE[,1] #Extract patients and controls metadatadata
table(names(Disease) == colnames(COUNTS_SE) )#Check the order
rownames(COUNTS_SE) <- sub("\\..*$", "", rownames(COUNTS_SE)) #Make it the same names

#Annotate genes
data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]
rm(annotEnsembl63)
y <- DGEList(counts=COUNTS_SE, genes=annot[rownames(COUNTS_SE),])
#Filtering and normalization
isexpr <- filterByExpr(y, group=Disease, min.count=1) #Filter out genes with less than 1 CPM (threshold)
table(isexpr) #How many
hasannot <- rowSums(is.na(y$genes))==0 #Remove genes without annotations
y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE] #Recompute the librery size
dim(y)
barplot(y$samples$lib.size*1e-6, names=1:ncol(y$counts), ylab="Library size (millions)")

#TMM Normalization
y <- calcNormFactors(y)
head(y$samples) 

#Take out the covariates to include them in the model
Sex <- DESIGN_SE[,2]
Age  <- DESIGN_SE[,3]
Cell <- DESIGN_SE[,13]

desig <- model.matrix(~Sex+Age+Cell+Disease)
y <- estimateDisp(y, desig, robust=TRUE) #Estimate dispersion
plotBCV(y)
fit <- glmQLFit(y, desig, robust=TRUE)
plotQLDisp(fit)
#Differential expression
qlf <- glmQLFTest(fit)
topTags(qlf,n=25)
summary(decideTests(qlf))
DEGs <- topTags(qlf,n=nrow(y$genes))$table
DEG <- DEGs[DEGs$FDR<0.05,]
```
### Transcription factor activity analysis
To estimate the activity of TFs from differential expression data, we utilized human TFs data from the [collecTRI](https://github.com/saezlab/CollecTRI) package. decoupleR package was used to infer their activity in regards to differential expression results following the TF activity inference in [bulk RNA-seq guidelines](https://saezlab.github.io/decoupleR/articles/tf_bk.html#loading-packages). TFs with an activity score greater than 1 in absolute value and p-value<0.05 were considered significant. 
```
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)

## ----"load data"------------------------------------------------------------------------------------------------------
inputs_dir <- system.file("extdata", package = "decoupleR")
data <- readRDS(file.path(inputs_dir, "bk_data.rds"))

## ----"counts"---------------------------------------------------------------------------------------------------------
# Remove NAs and set row names
counts <- data$counts %>%
  dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x)) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix()
head(counts)

## ----"design"---------------------------------------------------------------------------------------------------------
design <- data$design
design

## ----"deg"------------------------------------------------------------------------------------------------------------
DEGs <- read.delim("~/PRECISESADS/Meth/JML_methylation/eQTMs/DEGs.tsv", header=FALSE)
DEGs$t <- -log10(DEGs$V4)*DEGs$V5 
colnames(DEGs) <- c('Number','ID','P.Value.NotAdjusted','P.Value','logFC','F','t')
DEGs$ID <- make.names(DEGs$ID)
deg <- DEGs[,c(2,5,7,4)] %>% 
  filter(!ID%in%c('AC012652.1', 'AC138035.1', 'AL592284.1', 'X5S_rRNA', 'X7SK', 'Y_RNA'))%>%
  filter(!is.na(t)) %>% 
  column_to_rownames(var = 'ID') %>%
  as.matrix()
head(deg)

## ----"collectri"------------------------------------------------------------------------------------------------------
net <- get_collectri(organism='human', split_complexes=FALSE)
net

## ----"contrast_umean", message=FALSE----------------------------------------------------------------------------------
# Run umean
contrast_acts <- run_ulm(mat=deg[, 't', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

SIGNIF.TFs <- contrast_acts[contrast_acts$p_value<0.05&abs(contrast_acts$score)>1,]
```
### eQTM analysis - Integration with DNA methylation data
Expression quantitative trait methylation (eQTMs) were calculated through a Pearson correlation test between a DMP and a differentially expressed gene (DEG). This integrative approach was applied using the [MatrixEQTL](https://github.com/andreyshabalin/MatrixEQTL) R package. A maximum distance of 1 Mb between DMP and DEG was defined. Those eQTMs with a FDR<0.05 were considered significant. The pipeline was obtained as from Dr. Laura C. Terron-Camero [github repository](https://github.com/lterroncamero/GCA_CD14_Analysis/blob/main/Transcriptome/Integration.R).
```
## Definition of variables ####

dir<-"~/PRECISESADS/Meth/JML_methylation/eQTMs/"

# RNA-Seq data
Exp<- cpm(y, log = TRUE)
Exp

#Set the genes names
orden <- which(rownames(Exp)%in%rownames(y$genes)) #Check the order
rownames(Exp) <-  y$genes$Symbol

boxplot(Exp, main = "CPM Values by Sample", xlab = "Sample", ylab = "log(CPM) - TMM")
Exp
#write.table(Exp, "TMM.normalized.expressiondata.csv", col.names = T, row.names = T, sep = "\t", quote = F)
Exp <- "TMM.normalized.expressiondata.csv"


#Follow this guide to transform gtf into bed files 
# https://sciberg.com/resources/bioinformatics-scripts/converting-gtf-files-into-bed-files

Exp_loc<-fread("hg19.annotation.bed")
Exp_loc<-data.frame(Exp_loc[,c(7,1,2,3)])
colnames(Exp_loc) <- c("geneid","chr","start","end")
#write.table(Exp_loc, "Exp_loc.csv", col.names = T, row.names=F, sep = "\t", quote = F)
Exp_loc <- "Exp_loc.csv"

# Methylation data
#write.table(M_QTM, "M.Values.eQTMs.csv", col.names = T, row.names = T, sep = "\t", quote = F)
Meth_id<-"M.Values.eQTMs.csv"
require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annot_CPGs <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot_CPGs <- data.frame(annot_CPGs@listData)
annot_CPGs_study <- annot_CPGs[which(annot_CPGs$Name%in%rownames(M_SE)),]
write.table(annot_CPGs_study[,c("Name", "chr", "pos")], "Meth_loc.csv", col.names = T, row.names = T, sep = "\t", quote = F)
Meth_loc<-"Meth_loc.csv"

#write.table(t(as.matrix(DESIGN_SE[,c("disSSc", "sexMale", "age", "neutros")])), "covariates.csv", col.names = T, row.names = T, sep = "\t",quote = F)
covariates<-"covariates.csv"

pvOutputThreshold_cis<-0.01
pvOutputThreshold_trans<-0.0001

# Distance
dist_cis=1000000

## Location of the package with the data files.
base.dir = dir

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
Meth_file_name = paste(base.dir, Meth_id, sep="");
meth_location_file_name = paste(base.dir, Meth_loc, sep="");

# Gene expression file name
expression_file_name = paste(base.dir, Exp , sep="");
gene_location_file_name = paste(base.dir, Exp_loc , sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, covariates , sep="");

# Output file name
output_file_name_cis = tempfile(tmpdir = dir);
output_file_name_pericoeldelospalotes = tempfile(tmpdir= dir);

# Only associations significant at this level will be saved

##Load Meth data ####

meth = SlicedData$new();
meth$fileDelimiter = "\t";      # the TAB character
meth$fileOmitCharacters = "NA"; # denote missing values;
meth$fileSkipRows = 1;          # one row of column labels
meth$fileSkipColumns = 1;       # one column of row labels
meth$fileSliceSize = 2000;      # read file in slices of 2,000 rows
meth$LoadFile(Meth_file_name);

## Load gene expression data ####

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates ####

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

methpos<-read.table(meth_location_file_name, header = T,stringsAsFactors = F, sep = "\t")
methpos<-methpos[,1:3]
colnames(methpos)<-c("snpid","chr","pos")
genespos <-read.table(gene_location_file_name, header = T, stringsAsFactors = F, sep = "\t")
colnames(genespos)<-c("geneid","chr","start","end")

#Genos and phenos sort
nb.cols <- length(meth$columnNames)

same.order<-(sum(colnames(gene) == colnames(meth)) == nb.cols) & (sum(colnames(cvrt) == 
                                                                        colnames(meth)) == nb.cols)
if(!same.order){
  message("re-order the columns of the phenotype object")
  new.col.orderG <- rep(NA, nb.cols)
  new.col.orderS <- rep(NA, nb.cols)
  new.col.orderC <- rep(NA, nb.cols)
  coln.go <- meth$columnNames
  coln.gn <- meth$columnNames
  coln.p <- gene$columnNames
  coln.c <- cvrt$columnNames
  
  for(i in 1:nb.cols){
    #colocar gene en base a snp ordenado
    new.col.orderG[i] <- which(coln.p == coln.gn[i])
    new.col.orderS[i] <- which(coln.go == coln.gn[i])
    new.col.orderC[i] <- which(coln.c == coln.gn[i])
    
  }
  #colocar en base a snp ordenado
  gene$ColumnSubsample(new.col.orderG)
  meth$ColumnSubsample(new.col.orderS)
  cvrt$ColumnSubsample(new.col.orderC)
}

me = Matrix_eQTL_main(
  snps = meth,
  gene = gene,
  cvrt= cvrt,
  output_file_name = output_file_name_pericoeldelospalotes,
  pvOutputThreshold = pvOutputThreshold_trans,
  useModel = modelLINEAR,
  #errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = methpos,
  genepos = genespos,
  cisDist = dist_cis,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

#load("eQTMs_cisYtransDEF.RData")

#Beta
me$cis$eqtls$beta_se=me$cis$eqtls$beta/me$cis$eqtls$statistic
me$trans$eqtls$beta_se=me$trans$eqtls$beta/me$trans$eqtls$statistic

#R2 for cis
dfFull<-as.numeric(me$param$dfFull)
tstat_cis<-me$cis$eqtls$statistic
r_cis=tstat_cis/sqrt(dfFull+(tstat_cis^2))
R2_cis<-r_cis^2
me$cis$eqtls$R2<-R2_cis


## Generate trans
trans1<-merge(me$trans$eqtls,methpos,by.x="snps",by.y="snpid")
names(trans1)[names(trans1) == 'chr'] <- 'meth.chr'
names(trans1)[names(trans1) == 'pos'] <- 'meth.pos'
#
trans<-merge(trans1,genespos,by.x="snps",by.y="geneid")
names(trans)[names(trans) == 'chr'] <- 'gene.chr'
names(trans)[names(trans) == 'left'] <- 'gene.start'
names(trans)[names(trans) == 'right'] <- 'gene.end'
#
write.table(trans1,paste(dir,"/trans_results.tsv",sep=""),sep="\t",row.names = F)

## Generate cis
me$cis$eqtls$FDR<-p.adjust(me$cis$eqtls$pvalue,method="fdr",n=me$cis$ntests)

cis1<-merge(me$cis$eqtls,methpos,by.x="snps",by.y="snpid")
names(cis1)[names(cis1) == 'chr'] <- 'meth.chr'
names(cis1)[names(cis1) == 'pos'] <- 'meth.pos'

cis<-merge(cis1,genespos,by.x="gene",by.y="geneid")
names(cis)[names(cis) == 'chr'] <- 'gene.chr'
names(cis)[names(cis) == 'left'] <- 'gene.start'
names(cis)[names(cis) == 'right'] <- 'gene.end'

write.table(na.omit(cis),paste(dir,"/cis_results.txt",sep=""),sep="\t",row.names = F)

#Relaciones cis que son de este estudio
#cis <- read.table("cis_results.txt", header = T)

cis_results <- read.delim("~/PRECISESADS/Meth/JML_methylation/eQTMs/cis_results.txt")
cis_top <- cis_results[cis_results$FDR<0.05&cis_results$gene%in%DEG$Symbol&cis_results$snps%in%DMPs$ID,]

```
### Gene ontology enrichment of eQTMs
For the eQTMs, GO enrichment analyses were carried out with the [EnrichR](https://github.com/wjawaid/enrichR) R package. GO Biological terms with an adjusted p-value<0.05 and with a minimum count of 6 genes were considered significant. 
```
pista <- enrichr(cis_top$Gene, databases = "GO_Biological_Process_2021")
pista$GO_Biological_Process_2021$ID <- gsub(".*\\((GO:\\d+)\\).*", "\\1", pista$GO_Biological_Process_2021$Term) #


CIS.ENRICHMENT <- data.frame(Name = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Term"],
                   FDR = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Overlap"],
                   PValue = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"P.value"],
                   OR = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Odds.Ratio"],
                   Genes=pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Genes"])
                   
CIS.ENRICHMENT.TOP <- CIS.ENRICHMENT[CIS.ENRICHMENT$NumGenes>6,]
```
### Serum protein cytokines analysis
To analyze the association between DMPs, DEGs, and serum proteins in SSc patients, Spearman correlations were performed between protein levels assessed with different assays and the beta values and normalized gene counts, respectively. The p-values of these correlations were adjusted using the Benjamini-Hochberg method for both DMPs and DEGs analyses, with correlations having a FDR<0.05 were considered significant. 
```
#Import eQTMs
EQTMs <- read.delim("~/JML_methylation/eQTMs/SupplementaryMaterial - Supplementary Table 13 - eQTMs - DEF.tsv")
CpGeQTMlist <- c(unique(EQTMs$CpG)) #Take out the unique CpGs name
GENeQTMlist <- c(paste(unique(EQTMs$Gene), "_T", sep = '')) #Paste a _T to distinguish between some proteins and genes

colnames(CYTOKINES)[3:99]
NOCALL <- c(3:99)
NOCALL <- NOCALL[-c(14,18,25,26,28,33,38,41:45,49,54,69)] #Columns with proteins data 

#Correlations with DNA methylation data
CORRELATIONS.CPG <- matrix(NA, rep(0,), nrow=length(CpGeQTMlist)*length(NOCALL), ncol=7) #Create an empty matrix
CORRELATIONS.CPG <- data.frame(CORRELATIONS.CPG)
colnames(CORRELATIONS.CPG) <- c('CpG','Protein', 'n', 't', 'df', 'pvalue', 'cor')
contador <- 1
for(i in 1:length(CpGeQTMlist)){
    for(j in NOCALL){
    DUMB <- cor.test(CYTOKINES[,j], CYTOKINES[,CpGeQTMlist[i]])
    CORRELATIONS.CPG[contador,1] <- eQTMlist[i] 
    CORRELATIONS.CPG[contador,2] <- colnames(CYTOKINES)[j] 
    CORRELATIONS.CPG[contador,3] <- min(c(length(na.omit(CYTOKINES[,j])), length(na.omit(CYTOKINES[,CpGeQTMlist[i]]))))
    CORRELATIONS.CPG[contador,4] <- DUMB$statistic
    CORRELATIONS.CPG[contador,5] <- DUMB[["parameter"]][["df"]]
    CORRELATIONS.CPG[contador,6] <- DUMB$p.value
    CORRELATIONS.CPG[contador,7] <- DUMB[["estimate"]][["cor"]]
        
    contador <- contador+1
  }
  
}

CORRELATIONS.CPG <- na.omit(CORRELATIONS.CPG) 
CORRELATIONS.CPG$FDR <- p.adjust(CORRELATIONS.CPG$pvalue, method = 'fdr')
table(CORRELATIONS.CPG$FDR<0.05)
SIG.CORRELATIONS.CPG <- CORRELATIONS.CPG[which(CORRELATIONS.CPG$FDR<0.05),]

#Do the exact same with gene expression data instead of DNA methylation data
CORRELATIONS.GEN <- matrix(NA, rep(0,), nrow=length(GENeQTMlist)*length(NOCALL), ncol=7)
CORRELATIONS.GEN <- data.frame(CORRELATIONS.GEN)
colnames(CORRELATIONS.GEN) <- c('Gene','Protein','n', 't', 'df', 'pvalue', 'cor') #, '95CIlow','95CIhigh'
contador <- 1
for(i in 1:length(GENeQTMlist)){
  #print(paste(i, ' - Comienza el análisis de ', eQTMlist[i]))
  for(j in NOCALL){
    DUMB <- cor.test(CYTOKINES[,j], CYTOKINES[,GENeQTMlist[i]])
    CORRELATIONS.GEN[contador,1] <- GENeQTMlist[i] 
    CORRELATIONS.GEN[contador,2] <- colnames(CYTOKINES)[j] 
    CORRELATIONS.GEN[contador,3] <- min(c(length(na.omit(CYTOKINES[,j])), length(na.omit(CYTOKINES[,GENeQTMlist[i]]))))
    CORRELATIONS.GEN[contador,4] <- DUMB$statistic
    CORRELATIONS.GEN[contador,5] <- DUMB[["parameter"]][["df"]]
    CORRELATIONS.GEN[contador,6] <- DUMB$p.value
    CORRELATIONS.GEN[contador,7] <- DUMB[["estimate"]][["cor"]]
    #CORRELATIONS[contador,7] <- DUMB[["conf.int"]][1] Quito los Confint porque da fallo a veces
    #CORRELATIONS[contador,8] <- DUMB[["conf.int"]][2] 
    
    contador <- contador+1
  }
  #print(paste(i, ' - Ha terminado el análisis de ', eQTMlist[i]))
}

CORRELATIONS.GEN$Gene <- gsub('_T', '', CORRELATIONS.GEN$Gene)
CORRELATIONS.GEN <- na.omit(CORRELATIONS.GEN) 
CORRELATIONS.GEN$FDR <- p.adjust(CORRELATIONS.GEN$pvalue, method = 'fdr')
table(CORRELATIONS.GEN$FDR<0.05)
SIG.CORRELATIONS.GEN <- CORRELATIONS.GEN[which(CORRELATIONS.GEN$FDR<0.05),]
```

**NOTE:** [HOMER](http://homer.ucsd.edu/homer/) and [GREAT](http://great.stanford.edu/public/html/) analysis are not detailed here as they were performed through their corresponding APIs or web resources.  

To cite this paper:

> Martinez-Lopez, J., Estupiñan-Moreno, E., Ortiz-Fernandez, L. et al., A genome-wide DNA methylation study integrated with gene expression data reveals novel insights in the epigenetic landscape of Systemic Sclerosis, _JOURNAL TBD_, 2024, DOI:

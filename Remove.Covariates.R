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

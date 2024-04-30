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

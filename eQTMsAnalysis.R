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


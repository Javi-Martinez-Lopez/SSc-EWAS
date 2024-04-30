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

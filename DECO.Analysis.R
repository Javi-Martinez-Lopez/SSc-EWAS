bpparam <- MulticoreParam(8, RNGseed = 7739465, progressbar = TRUE) #Para hacerlo en memoria compartida

table(Data_SE$Sample_Name==colnames(M_DECO))

Data_DECO <- Data_SE$Disease
names(Data_DECO) <- Data_SE$Sample_Name

deco1 <- decoRDA(data=M_DECO, classes=Data_DECO, iterations = 100, r=30, q.val = 0.01, bpparam = bpparam)

table(is.infinite(c(deco1$data, deco1$results)))
na.omit(deco1)

deco2 <-  decoNSCA(sub = deco1, v=50, method=NULL, bpparam = bpparam, k.control=NULL, k.case = NULL, samp.perc = 0.05, rep.thr = 1) #Da un fallo Biocparallel que tendré que solucionar en un futuro.

hMatrixSSc <- NSCAcluster(deco2)$Case$NSCA$h
hMatrixCTRL <- NSCAcluster(deco2)$Control$NSCA$h

sig_CPGs_PRE <- deco2@featureTable[,c(1,2,3,8,9,11,12,23)]


sig_CPGs <- sig_CPGs_PRE[!sig_CPGs_PRE$ID%in%row.names(SEGE[SEGE$adj.P.Val<0.05,]),] #We remove CpGs associated with covariates
DMPs <- sig_CPGs[sig_CPGs$Profile=="Minority",]

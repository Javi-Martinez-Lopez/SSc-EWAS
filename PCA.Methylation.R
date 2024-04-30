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

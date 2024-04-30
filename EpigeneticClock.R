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

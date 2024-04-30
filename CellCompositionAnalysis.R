library(EpiDish)
data(centDHSbloodDMC.m)

B <- lumi::m2beta(Mvalues_merge)

cell_counts <- epidish(B, ref.m = centDHSbloodDMC.m, method = "RPC", maxit = 100000)
cells <- cell_counts$estF
cells <- cells[Data_SE$Sample_Name,]
table(row.names(cells)==Data_SE$Sample_Name)


DATA <- cbind(Data_SE, cells)

pista <- enrichr(cis_top$Gene, databases = "GO_Biological_Process_2021")
pista$GO_Biological_Process_2021$ID <- gsub(".*\\((GO:\\d+)\\).*", "\\1", pista$GO_Biological_Process_2021$Term) #


CIS.ENRICHMENT <- data.frame(Name = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Term"],
                   FDR = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Overlap"],
                   PValue = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"P.value"],
                   OR = pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Odds.Ratio"],
                   Genes=pista$GO_Biological_Process_2021[pista$GO_Biological_Process_2021$Adjusted.P.value<0.05,"Genes"])
                   
CIS.ENRICHMENT.TOP <- CIS.ENRICHMENT[CIS.ENRICHMENT$NumGenes>6,]

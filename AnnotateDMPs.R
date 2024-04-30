require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") #Get hg19 annotations

annot_CPGs <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot_CPGs <- annot_CPGs[match(row.names(sig_CPGs),annot_CPGs$Name),]
table(row.names(sig_CPGs)==row.names(annot_CPGs))

annotated_CPGs <- data.frame(cbind(sig_CPGs,annot_CPGs))

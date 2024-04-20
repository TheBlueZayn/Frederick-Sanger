library(dplyr)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
getwd()
setwd("E:/hackio_genomics_workshop/stage3/AMR")
AMR_files <- list.files("E:/hackio_genomics_workshop/stage3/AMR")
AMR_sep.files <- lapply(AMR_files, read.csv)
AMR_all = bind_rows(AMR_sep.files)

AMR_trim = AMR_all %>% select(Name,X..Identity.to.reference.sequence,Gene.symbol)

matrixAMR=acast(AMR_trim,Name~Gene.symbol, value.var ="X..Identity.to.reference.sequence", fun.aggregate = NULL)
AMRdf=as.data.frame(matrixAMR)
colnames(AMRdf)
AMRdf_trim=AMRdf[, c("catA1", "catA2", "cmlA5", "floR", "sul1", "sul2", "dfrA19", "dfrA23","blaOXA-10", "blaSCO-1","blaTEM-1","blaCTX-M-15","blaSHV-12","ere(A)","mph(A)","qnrA1","qnrB2")]

setwd("E:/hackio_genomics_workshop/stage3/plasmidtyping")
#read csv files generated from abricate tool using plasmidfinder database

plasmidtyping_files <- list.files("E:/hackio_genomics_workshop/stage3/plasmidtyping")
plasmidtyping_sep.files <- lapply(plasmidtyping_files, read.csv)

#merge all the files
plasmidtyping_all = bind_rows(plasmidtyping_sep.files)

plasmidtyp_trim = plasmidtyping_all %>% select(X.FILE,X.IDENTITY,GENE)
matrixplasmidtyp=acast(plasmidtyp_trim,X.FILE~GENE, value.var ="X.IDENTITY", fun.aggregate = NULL)
plasmidtyp_df=as.data.frame(matrixplasmidtyp)
colnames(plasmidtyp_df)
plasmidtyp_df_trim=plasmidtyp_df[, c("Col(BS512)_1","IncA/C2_1","IncHI2A_1","IncHI2_1","IncI1_1_Alpha","IncQ1_1","IncY_1")]

##serotyping
getwd()
setwd("E:/hackio_genomics_workshop/stage3/serotyping")
serotyping_files <- list.files("E:/hackio_genomics_workshop/stage3/serotyping")
serotyping_sep.files <- lapply(serotyping_files, read_tsv)
serotyping_all = cbind.data.frame(serotyping_sep.files)
merged_df <- do.call(rbind, serotyping_sep.files)
sero_phenotype=merged_df[,c(1,9)]
sero_phenotype=sero_phenotype %>% separate(`Sample name`, c("sample_name","delete"))
sero_phenotype = sero_phenotype[,c(1,3)]
sample_names=unlist(sero_phenotype$sample_name)
rownames(sero_phenotype) = sample_names
sero_phenotype=sero_phenotype[-1]
rownames(sero_phenotype) = sample_names
Heatmap(sero_phenotype) 

sero_genes=merged_df[,c(1,4,5,6)]
rownames(sero_genes) = sample_names
sero_genes=sero_genes[,c(2,3,4)]
Heatmap(sero_genes)


##heatmap plotting
Heatmap(AMRdf_trim,name = "AMR",col=col_fun, column_title = "AMR_genes") + Heatmap(plasmidtyp_df_trim,name = "replicongenes", col=col_fun, column_title = "Plasmid_Typing") + Heatmap(sero_phenotypename = "insilico_serotype", col=col_fun, column_title = "Serotyping")
column_ha = HeatmapAnnotation(AMRgenes = runif(17))
col_hp = HeatmapAnnotation(replicongenes = runif(7))
Heatmap(AMRdf_trim, name = "AMR", col=col_fun, column_title = "AMR_genes", row.names(sample_names)) + Heatmap(plasmidtyp_df_trim, name = "replicongenes", col=col_fun, column_title = "Plasmid_Typing") + Heatmap(sero_phenotype, name = "insilico_serotype", col=col_fun, column_title = "Serotyping") + Heatmap(sero_genes, name = "serotype_genes", col=col_fun)

dim(AMRdf_trim)
dim(plasmidtyp_df_trim)

col_fun=colorRamp2(c(0, 1, 2), c("red", "grey","yellow"))
seq(-3,3)
col_fun(seq(-3,3))
#finalheatmap
Heatmap(sero_phenotype) + Heatmap(plasmidtyp_df_trim) + Heatmap(AMRdf_trim) + Heatmap(sero_genes)

write.csv(AMRdf_trim, "E:/hackio_genomics_workshop/AMR_data_for_heatmap.csv")
write.csv(plasmidtyp_df_trim, "E:/hackio_genomics_workshop/plasmidtyping_data_for_heatmap.csv")

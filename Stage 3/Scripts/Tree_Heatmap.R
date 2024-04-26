#load the required packages
library(dplyr)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggtree)


###phylogenetic tree generation###

#set working directory
setwd("E:/hackio_genomics_workshop/stage3/tree")
# Load data 
tree <- read.tree("RAxML_result.raxml.bestTree")
#Plot the tree
plot(tree, show.tip.label = F)

##rownames of AMR and replicon genes must correspond to the tip labels of the tree
#get the tip labels
writeClipboard(tree$tip.label)

#make the tip labels in a vector to reorder the AMR and plasmidtyping genes according to them
tip_lable <- c("ERR9516298",
               "ERR9516208",
               "SRR8294664",
               "ERR9516234",
               "ERR9516178",
               "ERR9516213",
               "ERR3635344",
               "ERR9516202",
               "SRR8159842",
               "SRR11903203",
               "ERR9516265",
               "ERR9516228",
               "ERR9516335",
               "ERR9516226",
               "SRR13221299",
               "ERR9516307",
               "ERR9516190",
               "ERR9516242",
               "ERR9516356",
               "ERR9516264",
               "ERR9516292",
               "ERR9516269",
               "ERR9516319",
               "ERR9516309",
               "ERR9516306",
               "SRR11927771",
               "ERR9516192",
               "SRR8159839",
               "SRR8327695",
               "SRR7514404",
               "ERR9516341",
               "ERR9516211",
               "ERR9516191",
               "ERR9516287",
               "ERR9855643",
               "ERR9516186",
               "ERR9516233",
               "ERR9516215",
               "ERR9516204",
               "ERR9516260",
               "ERR9516194",
               "ERR9516279",
               "ERR9516349",
               "ERR9516338",
               "ERR9516266",
               "ERR9516189",
               "SRR11927768",
               "ERR9516354",
               "ERR9516219",
               "SRR11927779")



##AMR matix generation
#set the working directory to where the csv files of the AMR genes each sample is located
getwd()
setwd("E:/hackio_genomics_workshop/stage3/AMR")
#load the csv files
AMR_files <- list.files("E:/hackio_genomics_workshop/stage3/AMR")
#read the csv files, separatly
AMR_sep.files <- lapply(AMR_files, read.csv)
#merge the csv files into a one dataframe
AMR_all = bind_rows(AMR_sep.files)
#select the rows the will be used: sample name, gene symbol, identity to ref seq corresponds to the presence of the gene of not
AMR_trim = AMR_all %>% select(Name,X..Identity.to.reference.sequence,Gene.symbol)
#make a dataframe to let the gene names as colnames and the values are the identity to ref
#where 1 is 100% identity corresponding to presence and 0 is N/A corresponding to absence
AMR_ordered=reshape2::dcast(AMR_trim,Name~Gene.symbol, value.var ="X..Identity.to.reference.sequence", fun.aggregate = NULL)
#make the rownames of AMR_ordered dataframe the sample names
sample_names = AMR_ordered$Name
rownames(AMR_ordered) = sample_names
colnames(AMR_ordered)

#extract the genes that are in phenicol class
phenicol= AMR_ordered[, c("catA1", "catA2", "cmlA5", "floR","Name")]
phenicol = phenicol %>%
  slice(match(tip_lable, Name))
phenicol=phenicol[,c("catA1", "catA2", "cmlA5", "floR")]

#extract the genes that are in sulfonamide class
sulphonamide = AMR_ordered[, c("sul1", "sul2","Name")]
sulphonamide = sulphonamide %>%
  slice(match(tip_lable, Name))
sulphonamide=sulphonamide[,c("sul1", "sul2")]

#extract the genes that are in trimethoprim class
trimethoprim = AMR_ordered[, c("dfrA19", "dfrA23", "Name")]
trimethoprim = trimethoprim %>%
  slice(match(tip_lable, Name))
trimethoprim=trimethoprim[,c("dfrA19", "dfrA23")]

#extract the genes that are in betalactam class
betalactam = AMR_ordered[, c("blaOXA-10", "blaSCO-1", "blaTEM-1", "blaCTX-M-15","blaSHV-12","Name")]
betalactam = betalactam %>%
  slice(match(tip_lable, Name))
betalactam=betalactam[,c("blaOXA-10", "blaSCO-1", "blaTEM-1", "blaCTX-M-15","blaSHV-12")]


#extract the genes that are in macrolide class
macrolide = AMR_ordered[, c("ere(A)","mph(A)","Name")]
macrolide = macrolide %>%
  slice(match(tip_lable, Name))
macrolide=macrolide[,c("ere(A)","mph(A)")]

#extract the genes that are in quinolone class
quinolone = AMR_ordered[,c("qnrA1","qnrB2","Name")]
quinolone = quinolone %>%
  slice(match(tip_lable, Name))
quinolone=quinolone[,c("qnrA1","qnrB2")]


##plasmidtyping##
setwd("E:/hackio_genomics_workshop/stage3/plasmidtyping")
#read csv files generated from abricate tool using plasmidfinder database

plasmidtyping_files <- list.files("E:/hackio_genomics_workshop/stage3/plasmidtyping")
plasmidtyping_sep.files <- lapply(plasmidtyping_files, read.csv)

#merge all the files
plasmidtyping_all = bind_rows(plasmidtyping_sep.files)
#select the columns that contain the name of the sample, replicon genes, and the %identity to ref 
plasmidtyp_trim = plasmidtyping_all %>% select(X.FILE,X.IDENTITY,GENE)
#reorder the plasmidtyping trimmed file to let the sample name rownames and replicon genes the colnames
plasmidtyp_ordered=reshape2::dcast(plasmidtyp_trim,X.FILE~GENE, value.var ="X.IDENTITY", fun.aggregate = NULL)
colnames(plasmidtyp_ordered)

#make the rownames of ordered plasmidtyping dataframe the sample names
sample_names_pt = plasmidtyp_ordered$X.FILE
rownames(plasmidtyp_ordered) = sample_names_pt
#select only the specified replicon genes that will be plotted in the heatmap
plasmidtyp_ordered_selected=plasmidtyp_ordered[, c("X.FILE","Col(BS512)_1","IncA/C2_1","IncHI2A_1","IncHI2_1","IncI1_1_Alpha","IncQ1_1","IncY_1")]



f_plasmidtyp_ordered_selected = plasmidtyp_ordered_selected %>%
  slice(match(tip_lable, X.FILE))
f_plasmidtyp_ordered_selected=f_plasmidtyp_ordered_selected[, c(2,3,4,5,6,7,8)]

##plotting the summary data
setwd("E:/hackio_genomics_workshop/stage3/summary_data")
summary_data = read.csv("E:/hackio_genomics_workshop/stage3/summary_data/summary_data.csv")
rownames(summary_data) = tip_lable

#choose the colours of the heatmap
col_fun=colorRamp2(c(0, 1, 2), c("white", "blue","blue"))
seq(-3,3)
col_fun(seq(-3,3))


#draw the heatmap 
ht1 = Heatmap(phenicol, row_order = tip_lable, name = "phenicol",col=col_fun, column_title = "phenicol", column_title_rot = 90, border_gp = gpar(col = "black"),show_heatmap_legend = FALSE)
ht2 = Heatmap(sulphonamide,row_order = tip_lable, name="sulphonamide", col=col_fun, column_title = "sulphonamide", column_title_rot = 90,border_gp = gpar(col = "black"), show_heatmap_legend = FALSE)
ht3 = Heatmap(trimethoprim, row_order = tip_lable, name = "trimethoprim", col = col_fun, column_title = "trimethoprim", column_title_rot = 90,border_gp = gpar(col = "black"),show_heatmap_legend = FALSE)
ht4 = Heatmap(betalactam, row_order = tip_lable, name = "betalactam", col = col_fun, column_title = "betalactam", column_title_rot = 90,border_gp = gpar(col = "black"),show_heatmap_legend = FALSE)
ht5 = Heatmap(macrolide, row_order = tip_lable, name = "macrolide", col = col_fun, column_title = "macrolide", column_title_rot = 90,border_gp = gpar(col = "black"),show_heatmap_legend = FALSE)
ht6 = Heatmap(quinolone, row_order = tip_lable, name = "quinolone", col = col_fun, column_title = "quinolone", column_title_rot = 90,border_gp = gpar(col = "black"),show_heatmap_legend = FALSE)

country_factor <- factor(summary_data$SampleIsolationCountry)
date = anno_points(summary_data$SampleIsolationDate)
source_factor <- factor(summary_data$SampleSource)

country_date_source= rowAnnotation(Isolation_country=country_factor, Sample_source=source_factor, Date=anno_points(summary_data$SampleIsolationDate))

ht7 = Heatmap(f_plasmidtyp_ordered_selected, row_order = tip_lable, name = "replicon genes", col = col_fun, column_title = "replicon genes", column_title_rot = 90,border_gp = gpar(col = "black"), right_annotation = country_date_source, show_heatmap_legend = FALSE)


#draw all heatmaps together

ht_list=ht1+ht2+ht3+ht4+ht5+ht6+ht7

#add ligands
lgd_list = Legend(labels = c("present", "absent"), title = "AMR and replicon genes", pch = 15, type = "points",
                  legend_gp = gpar(col = c("blue","white")))

draw(ht_list, annotation_legend_list = lgd_list)

#AMRdf_trim=AMRdf[, c("catA1", "catA2", "cmlA5", "floR", "sul1", "sul2", "dfrA19", "dfrA23","blaOXA-10", "blaSCO-1","blaTEM-1","blaCTX-M-15","blaSHV-12","ere(A)","mph(A)","qnrA1","qnrB2")]

















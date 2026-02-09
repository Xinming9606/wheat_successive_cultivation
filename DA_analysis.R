rm(list = ls())
setwd("C:/Users/xuxm/Desktop/Adele_amp_analysis")

library(readxl)
library(tidyverse)
library(nloptr)
library(phyloseq)
library(stringr)
library(ggpubr)
library(magrittr)
library(qwraps2)
library(pander)
library(microbiome)
panderOptions('table.caption.prefix', NULL)
panderOptions('table.continues', NULL)
panderOptions('table.emphasize.rownames', FALSE)

source("ancom_bc.R")

# 1. Feature Table Construction

otu_table <- read.table("asv_table.tsv", # Loading otutable
                      sep="\t", 
                      header=TRUE, 
                      stringsAsFactors=TRUE, 
                      row.names = 1, 
                      comment.char ="")

taxonomy <- read.table("asv_taxa.tsv",  # Loading otutable
                      sep="\t", 
                      header=TRUE, 
                      stringsAsFactors=TRUE, 
                      row.names = 1, 
                      comment.char ="",
                      fill=TRUE)

metadata <- read.table("metadata_combi.tsv", # Loading Metadata
                       sep="\t", 
                       header=TRUE, 
                       stringsAsFactors=TRUE, 
                       row.names = 1, 
                       comment.char ="" )
metadata$SampleID <- rownames(metadata)

#aggregate your taxonomy on different levels
#genus level
ASV_ID <- taxonomy$genus
otu_table <- cbind(ASV_ID,otu_table)
genus_table=otu_table%>%group_by(ASV_ID)%>%
  summarise_all(sum)
genus_table <- genus_table[-78,]
genus_table=as.data.frame(genus_table)
# 2. Differential Abundance (DA) Analysis

##Uninoculated Resilient vs Susceptible
# Subset meta data
meta.data=metadata%>%filter(lib_response%in%c("uninoculated_resilient", "uninoculated_susceptible"))
meta.data$response=as.character(meta.data$response)
pander(table(meta.data$response))
# Subset OTU table
obs.abn=genus_table
rownames(obs.abn)=obs.abn$ASV_ID
obs.abn=obs.abn[, -1]
selected_cols <- rownames(meta.data)
obs.abn <- obs.abn[, selected_cols, drop = FALSE]


# Pre-Processing
feature.table=obs.abn; sample.var="SampleID"; group.var="response"; 
zero.cut=0.90; lib.cut=1000; neg.lb=TRUE
pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, group.var, 
                                      zero.cut, lib.cut, neg.lb)
feature.table=pre.process$feature.table
library.size=pre.process$library.size
group.name=pre.process$group.name
group.ind=pre.process$group.ind
struc.zero=pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05

# Run ANCOM-BC
start_time <- Sys.time()
out=ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
             tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

# Test result
res.ANCOM_BC=data.frame(genus=rownames(out$feature.table), out$res, 
                        struc.zero[rownames(out$feature.table), ], 
                        row.names = NULL, stringsAsFactors=FALSE, check.names = FALSE)
alpha.adj=0.05/nrow(res.ANCOM_BC)
critic.val=qnorm(1-alpha.adj/2)
res.ANCOM_BC=res.ANCOM_BC%>%transmute(genus,
                                      `log fold change (resilient - susceptible)` = -`mean.difference (susceptible - resilient)`,
                                      se = `se (susceptible - resilient)`,
                                      ci.lo.adj=`log fold change (resilient - susceptible)`-critic.val*se, 
                                      ci.up.adj=`log fold change (resilient - susceptible)`+critic.val*se, 
                                      p.val, q.val, diff.abn, 
                                      `structural.zero (resilient)`, 
                                      `structural.zero (susceptible)`)
res.ANCOM_BC=res.ANCOM_BC%>%arrange(q.val)
write_csv(res.ANCOM_BC, "C:/Users/xuxm/Desktop/Adele_amp_analysis/genus_uninoculated_resilient_vs_susceptible.csv")
summary(res.ANCOM_BC$`log fold change (resilient - susceptible)` > 0)

#Library resilient vs susceptible 
meta.data=metadata%>%filter(lib_response%in%c("library_resilient", "library_susceptible"))
meta.data$response=as.character(meta.data$response)
pander(table(meta.data$response))
# Subset OTU table
obs.abn=genus_table
rownames(obs.abn)=obs.abn$ASV_ID
obs.abn=obs.abn[, -1]
selected_cols <- rownames(meta.data)
obs.abn <- obs.abn[, selected_cols, drop = FALSE]

# Pre-Processing
feature.table=obs.abn; sample.var="SampleID"; group.var="response"; 
zero.cut=0.90; lib.cut=1000; neg.lb=TRUE
pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, group.var, 
                                      zero.cut, lib.cut, neg.lb)
feature.table=pre.process$feature.table
library.size=pre.process$library.size
group.name=pre.process$group.name
group.ind=pre.process$group.ind
struc.zero=pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05

# Run ANCOM-BC
start_time <- Sys.time()
out=ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
             tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

# Test result
res.ANCOM_BC=data.frame(genus=rownames(out$feature.table), out$res, 
                        struc.zero[rownames(out$feature.table), ], 
                        row.names = NULL, stringsAsFactors=FALSE, check.names = FALSE)
alpha.adj=0.05/nrow(res.ANCOM_BC)
critic.val=qnorm(1-alpha.adj/2)
res.ANCOM_BC=res.ANCOM_BC%>%transmute(genus,
                                      `log fold change (resilient - susceptible)` = -`mean.difference (susceptible - resilient)`,
                                      se = `se (susceptible - resilient)`,
                                      ci.lo.adj=`log fold change (resilient - susceptible)`-critic.val*se, 
                                      ci.up.adj=`log fold change (resilient - susceptible)`+critic.val*se, 
                                      p.val, q.val, diff.abn, 
                                      `structural.zero (resilient)`, 
                                      `structural.zero (susceptible)`)
res.ANCOM_BC=res.ANCOM_BC%>%arrange(q.val)
write_csv(res.ANCOM_BC, "C:/Users/xuxm/Desktop/Adele_amp_analysis/genus_library_resilient_vs_susceptible.csv")
summary(res.ANCOM_BC$`log fold change (resilient - susceptible)` > 0)
#Mode   FALSE    TRUE 


##nostress library vs nostress uninoculated
# Subset meta data
meta.data=metadata%>%filter(lib_response%in%c("uninoculated_nostress", "library_nostress"))
meta.data$lib.unin=as.character(meta.data$lib.unin)
pander(table(meta.data$lib.unin))
# Subset OTU table
obs.abn=genus_table
rownames(obs.abn)=obs.abn$ASV_ID
obs.abn=obs.abn[, -1]
selected_cols <- rownames(meta.data)
obs.abn <- obs.abn[, selected_cols, drop = FALSE]


# Pre-Processing
feature.table=obs.abn; sample.var="SampleID"; group.var="lib.unin"; 
zero.cut=0.90; lib.cut=1000; neg.lb=TRUE
pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, group.var, 
                                      zero.cut, lib.cut, neg.lb)
feature.table=pre.process$feature.table
library.size=pre.process$library.size
group.name=pre.process$group.name
group.ind=pre.process$group.ind
struc.zero=pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05

# Run ANCOM-BC
start_time <- Sys.time()
out=ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
             tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

# Test result
res.ANCOM_BC=data.frame(genus=rownames(out$feature.table), out$res, 
                        struc.zero[rownames(out$feature.table), ], 
                        row.names = NULL, stringsAsFactors=FALSE, check.names = FALSE)
alpha.adj=0.05/nrow(res.ANCOM_BC)
critic.val=qnorm(1-alpha.adj/2)
res.ANCOM_BC=res.ANCOM_BC%>%transmute(genus,
                                      `log fold change (library - uninoculated)` = -`mean.difference (uninoculated - library)`,
                                      se = `se (uninoculated - library)`,
                                      ci.lo.adj=`log fold change (library - uninoculated)`-critic.val*se, 
                                      ci.up.adj=`log fold change (library - uninoculated)`+critic.val*se, 
                                      p.val, q.val, diff.abn, 
                                      `structural.zero (library)`, 
                                      `structural.zero (uninoculated)`)
res.ANCOM_BC=res.ANCOM_BC%>%arrange(q.val)
write_csv(res.ANCOM_BC, "C:/Users/xuxm/Desktop/Adele_amp_analysis/genus_nostress_library_vs_uninoculated.csv")
summary(res.ANCOM_BC$`log fold change (library - uninoculated)` > 0)


##resilient library vs resilient uninoculated
# Subset meta data
meta.data=metadata%>%filter(lib_response%in%c("uninoculated_resilient", "library_resilient"))
meta.data$lib.unin=as.character(meta.data$lib.unin)
pander(table(meta.data$lib.unin))
# Subset OTU table
obs.abn=genus_table
rownames(obs.abn)=obs.abn$ASV_ID
obs.abn=obs.abn[, -1]
selected_cols <- rownames(meta.data)
obs.abn <- obs.abn[, selected_cols, drop = FALSE]


# Pre-Processing
feature.table=obs.abn; sample.var="SampleID"; group.var="lib.unin"; 
zero.cut=0.90; lib.cut=1000; neg.lb=TRUE
pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, group.var, 
                                      zero.cut, lib.cut, neg.lb)
feature.table=pre.process$feature.table
library.size=pre.process$library.size
group.name=pre.process$group.name
group.ind=pre.process$group.ind
struc.zero=pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05

# Run ANCOM-BC
start_time <- Sys.time()
out=ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
             tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

# Test result
res.ANCOM_BC=data.frame(genus=rownames(out$feature.table), out$res, 
                        struc.zero[rownames(out$feature.table), ], 
                        row.names = NULL, stringsAsFactors=FALSE, check.names = FALSE)
alpha.adj=0.05/nrow(res.ANCOM_BC)
critic.val=qnorm(1-alpha.adj/2)
res.ANCOM_BC=res.ANCOM_BC%>%transmute(genus,
                                      `log fold change (library - uninoculated)` = -`mean.difference (uninoculated - library)`,
                                      se = `se (uninoculated - library)`,
                                      ci.lo.adj=`log fold change (library - uninoculated)`-critic.val*se, 
                                      ci.up.adj=`log fold change (library - uninoculated)`+critic.val*se, 
                                      p.val, q.val, diff.abn, 
                                      `structural.zero (library)`, 
                                      `structural.zero (uninoculated)`)
res.ANCOM_BC=res.ANCOM_BC%>%arrange(q.val)
write_csv(res.ANCOM_BC, "C:/Users/xuxm/Desktop/Adele_amp_analysis/genus_resilient_library_vs_uninoculated.csv")
summary(res.ANCOM_BC$`log fold change (library - uninoculated)` > 0)


##susceptible library vs susceptible uninoculated
# Subset meta data
meta.data=metadata%>%filter(lib_response%in%c("uninoculated_susceptible", "library_susceptible"))
meta.data$lib.unin=as.character(meta.data$lib.unin)
pander(table(meta.data$lib.unin))
# Subset OTU table
obs.abn=genus_table
rownames(obs.abn)=obs.abn$ASV_ID
obs.abn=obs.abn[, -1]
selected_cols <- rownames(meta.data)
obs.abn <- obs.abn[, selected_cols, drop = FALSE]


# Pre-Processing
feature.table=obs.abn; sample.var="SampleID"; group.var="lib.unin"; 
zero.cut=0.90; lib.cut=1000; neg.lb=TRUE
pre.process=feature_table_pre_process(feature.table, meta.data, sample.var, group.var, 
                                      zero.cut, lib.cut, neg.lb)
feature.table=pre.process$feature.table
library.size=pre.process$library.size
group.name=pre.process$group.name
group.ind=pre.process$group.ind
struc.zero=pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name=group.name; grp.ind=group.ind; adj.method="bonferroni"
tol.EM=1e-5; max.iterNum=100; perNum=1000; alpha=0.05

# Run ANCOM-BC
start_time <- Sys.time()
out=ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, 
             tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

# Test result
res.ANCOM_BC=data.frame(genus=rownames(out$feature.table), out$res, 
                        struc.zero[rownames(out$feature.table), ], 
                        row.names = NULL, stringsAsFactors=FALSE, check.names = FALSE)
alpha.adj=0.05/nrow(res.ANCOM_BC)
critic.val=qnorm(1-alpha.adj/2)
res.ANCOM_BC=res.ANCOM_BC%>%transmute(genus,
                                      `log fold change (library - uninoculated)` = -`mean.difference (uninoculated - library)`,
                                      se = `se (uninoculated - library)`,
                                      ci.lo.adj=`log fold change (library - uninoculated)`-critic.val*se, 
                                      ci.up.adj=`log fold change (library - uninoculated)`+critic.val*se, 
                                      p.val, q.val, diff.abn, 
                                      `structural.zero (library)`, 
                                      `structural.zero (uninoculated)`)
res.ANCOM_BC=res.ANCOM_BC%>%arrange(q.val)
write_csv(res.ANCOM_BC, "C:/Users/xuxm/Desktop/Adele_amp_analysis/genus_susceptible_library_vs_uninoculated.csv")
summary(res.ANCOM_BC$`log fold change (library - uninoculated)` > 0)

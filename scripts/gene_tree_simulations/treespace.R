library(treespace)
library(RColorBrewer) 
library(ggplot2)
library(reshape2)

cp_tree <- read.tree("cp_tree_renamed.tre")
sp_tree <- read.tree("MPEST_accession.tre")
gene_trees <- read.tree("gene_trees_renamed.tre")


############################ Make data frame with tips assigned to categories #####################################################

# Create a data frame with tip names as rows and no columns
df <- data.frame(sp_tree$tip.label, row.names = TRUE)
# Add a column with category information to the data frame
df$species <- ifelse(grepl("flava",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. flava",
                     ifelse(grepl("Sminor",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. minor",
                            ifelse(grepl("jonesii|rubra|alabamensis|leucophylla|alata",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. rubra complex",
                                   ifelse(grepl("psittacina",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. psittacina",
                                          ifelse(grepl("purpurea|rosea",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. purpurea",
                                                 ifelse(grepl("oreophila",rownames(df), ignore.case = TRUE,fixed = FALSE), "S. oreophila",
                                                        "Outgroup"))))))
tip_labels <- rownames(df)
d <- cbind(tip_labels, data.frame(df, row.names=NULL))
d <- d %>% select (species, tip_labels)
df <- d

###############################

#trees need to be rooted
cp_tree <- root(cp_tree, outgroup = "Dcalifornica", resolve.root = TRUE)
sp_tree <- root(sp_tree, outgroup = "Dcalifornica")
for (i in 1:length(gene_trees)) {gene_trees[[i]] <- midpoint.root(gene_trees[[i]])}

#calculate tree distance
relatedTreeDist(gene_trees,d)

#make collapsed tree
ref_tree <- makeCollapsedTree(sp_tree, d, warnings = TRUE)

plot(ref_tree)

treeConcordance(ref_tree,sp_tree,df)
treeConcordance(ref_tree,cp_tree,df)
treeConcordance(ref_tree,gene_trees[[1]],df)

refTreeDist(refTree = ref_tree, trees = gene_trees, df)

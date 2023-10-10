library(ape)
library(phytools)
library(tidyverse)
library(phangorn)

setwd("~/Research/plastome paper/data")

# read in species tree
sp_tree <- read.tree("MPEST_accession.tre")

# read in simulated trees
simtrees <- read.tree("simtrees.newick")

######################################### Rename gene trees##################################################
# read all gene trees and put them in a multiphylo object
filenames <- list.files("./genetrees_IQ/", pattern="*", full.names=TRUE)
gene_trees <- lapply(filenames, read.tree)
class(gene_trees)<-"multiPhylo"

#make a list of IDs for the gene trees
gene_tree_ID <- substr(filenames, 16,nchar(filenames)-22)

# remove underscores and typos in tip labels from all gene trees
rename.tips.phylo <- function(tree) {
  tree$tip.label <- str_remove_all(tree$tip.label,'[_]')
  tree$tip.label <- str_replace_all(tree$tip.label,"SS","S")
  tree$tip.label <- str_replace_all(tree$tip.label,"DcalifornicaOR","Dcalifornica")
  tree$tip.label <- str_replace_all(tree$tip.label,"SpurpureapurpureaWI1","SpurpureapurpureaWI")
  tree <- drop.tip(tree, c("DcalifornicaUN2","DcalifornicaUN1"))
  return(tree)
}

gene_trees_renamed <- lapply(gene_trees, rename.tips.phylo)
class(gene_trees_renamed)<-"multiPhylo"

write.tree(gene_trees_renamed, file = "gene_trees_renamed.tre")

################################## Rename  my cp tree ##################################

cp_tree <- read.tree(file = "08202021.treefile")

#assign all tip labels starting with "m" to a vector
cp_tips_to_drop <- cp_tree$tip.label[grepl("^m",cp_tree$tip.label)]

# drop all tips starting with "m"
cp_tree <- drop.tip(cp_tree, cp_tips_to_drop)

# remove my sample names (e.g. j004)
cp_tree$tip.label <- str_sub(cp_tree$tip.label, 6,-8)

cp_tree <- drop.tip(cp_tree, c("DarlingtoniaUN2","DarlingtoniaUN1","PsittacinaUK"))

# add S, make second letter lowercase
cp_tree$tip.label <- paste0("S", tolower(substr(cp_tree$tip.label, 1, nchar(cp_tree$tip.label)-4)), substr(cp_tree$tip.label, nchar(cp_tree$tip.label)-3, nchar(cp_tree$tip.label)))


# replace naming conventions
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "rosea", "purpureavenosaburkii")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "jonesii", "rubrajonesii")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "alabamensis", "rubraalabamensis")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "montana", "venosamontana")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "SpurpureapurpureaWI1", "SpurpureapurpureaWI")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "Sflavarubricorpa", "Sflavarubricorpora")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "SdarlingtoniaOR", "Dcalifornica")
cp_tree$tip.label <- str_replace_all(cp_tree$tip.label, "SheliamphoraVE", "Hminor")

setdiff(cp_tree$tip.label, sp_tree$tip.label)
setdiff(sp_tree$tip.label, cp_tree$tip.label)

write.tree(cp_tree, file = "cp_tree_renamed.tre")

#################################### Check gene tree tips #######################################################
# initialize empty vector
GT_tip_labels <- c()
tip_label_list <- list()

# add tip labels from every gene tree to vector
for (tree in gene_trees_renamed) {
  tip_label_list <- append(tip_label_list,list(tree$tip.label))
  GT_tip_labels <- c(GT_tip_labels, tree$tip.label)
}

# preserved samples across all gene trees
Reduce(intersect, tip_label_list)

# histogram of number of samples
hist(lengths(tip_label_list))

# number of tips in species tree
length(sp_tree$tip.label)

#only unique tip labels
GT_tip_labels_unique <- unique(GT_tip_labels)
sort(GT_tip_labels_unique)

setdiff(GT_tip_labels_unique, sp_tree$tip.label)
setdiff(sp_tree$tip.label, GT_tip_labels_unique)

################################# Prune trees ###############################################

# read in trees
#cp_tree <- read.tree("cp_tree_pruned.tre")
#sp_tree <- read.tree("MPEST_accession.tre")


# taxa in cp_tree not in sp_tree
setdiff(cp_tree$tip.label, sp_tree$tip.label)

# taxa in sp_tree not in cp_tree
setdiff(sp_tree$tip.label, cp_tree$tip.label)

# vector of tips to prune
species <- c("DcalifornicaUN2", "DcalifornicaUN1")

# prune tips
pruned.tree <- drop.tip(cp_tree, species)

# prune everything except tips
pruned.tree<-drop.tip(cp_tree,cp_tree$tip.label[-match(species, cp_tree$tip.label)])

cp_tree <- pruned.tree

# rename tip
cp_tree$tip.label[1] <- "Dcalifornica"

# write tree to file
write.tree(cp_tree, file="cp_tree_pruned2.tre")

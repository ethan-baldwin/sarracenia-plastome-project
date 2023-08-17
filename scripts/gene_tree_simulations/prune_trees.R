library(ape)

setwd("~/Research/plastome paper/data")

# read in trees
cp_tree <- read.tree("cp_tree_pruned.tre")
sp_tree <- read.tree("MPEST_accession.tre")

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

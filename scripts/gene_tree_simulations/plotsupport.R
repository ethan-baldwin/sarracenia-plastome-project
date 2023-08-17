library(ape)
library(phytools)
library(ggplot2)
library(dplyr)
library(ggtree)
setwd("~/Research/plastome paper/data")

# read trees
sp_tree <- read.tree("MPEST_accession.tre")
simtrees <- read.tree("simtrees.newick")
cp_tree <- read.tree("cp_tree_pruned2.tre")

#pp <- prop.part(simtrees)

# calculate clade frequencies
clad_sp <- prop.clades(sp_tree, simtrees, rooted = TRUE)
clad_cp <- prop.clades(cp_tree, simtrees, rooted = TRUE)


# attach frequencies to trees (replacing bootstraps)
sp_tree$node.label <- clad_sp
cp_tree$node.label <- clad_cp

#convert NA to 0
sp_tree$node.label[is.na(sp_tree$node.label)] <- 0
cp_tree$node.label[is.na(cp_tree$node.label)] <- 0

#convert to percentages
sp_tree$node.label <- (sp_tree$node.label)/10
cp_tree$node.label <- (cp_tree$node.label)/10


# plot trees
ggtree(sp_tree, branch.length = "none") + 
  geom_tiplab() +
  geom_nodelab(aes(label=label),
                    size=4,
                    hjust = 1.3,
                    vjust = -0.3) +
  scale_x_continuous(limits = c(0.0, 25))

ggtree(cp_tree, branch.length = "none") + 
  geom_tiplab() +
  geom_nodelab(aes(label=label),
               size=4,
               hjust = 1.2,
               vjust = -0.25) +
  scale_x_continuous(limits = c(0.0, 18))


# plot.phylo(sp_tree, type = "phylogram", use.edge.length = FALSE, font = 1)
# drawSupportOnEdges(boot)
# nodelabels(clad)


# layout(1)
# par(mar = rep(2, 4))
# plot(sp_tree, main = "Bipartition vs. Clade Support Values")
# drawSupportOnEdges(boot)
# nodelabels(clad)
# legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
#        pt.bg = c("green", "lightblue"), pt.cex = 2.5)

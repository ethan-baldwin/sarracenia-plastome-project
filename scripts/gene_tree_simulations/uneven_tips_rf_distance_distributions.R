library(phangorn)
library(phytools)
library(tidyverse)
library(ggpubr)
library(TreeDist)
library(forcats)
library(stringr)
library(ggplot2)

############################### Functions ###################################
# drop tips that aren't in both trees, then calculate generalized information theory based tree distance
rf.dist.tip.changes <- function(tree1,tree2) {
  diff1 <- setdiff(tree1$tip.label,tree2$tip.label)
  diff2 <- setdiff(tree2$tip.label,tree1$tip.label)
  tree1_comp <- drop.tip(tree1,diff1)
  tree2_comp <- drop.tip(tree2,diff2)
  TreeDistance(tree1_comp,tree2_comp)
}

# drop tips that aren't in both trees, then visualize tree distance
visualize.tip.changes <- function(Dist.metric,tree1,tree2) {
  diff1 <- setdiff(tree1$tip.label,tree2$tip.label)
  diff2 <- setdiff(tree2$tip.label,tree1$tip.label)
  tree1_comp <- drop.tip(tree1,diff1)
  tree2_comp <- drop.tip(tree2,diff2)
  VisualizeMatching(Dist.metric,tree1_comp,tree2_comp,Plot = TreeDistPlot)
}

# loop over a set of gene trees and run rf.dist.tip.changes against a reference tree
multi_RF <- function(ref_tree,multitree) {
  # make a variable name from the input. e.g. if sp_tree is the input, vname <- "sp_tree_dist"
  vname <- paste(deparse(substitute(multitree)),deparse(substitute(ref_tree)),"dist", sep="_")
  # assign vname to the global environment
  assign(vname, c(), envir = .GlobalEnv)
  # calculate Robinson-Fould distances for each tree in multitree vs. the ref_tree
  for (i in 1:length(multitree)) {
    x <- rf.dist.tip.changes(ref_tree, multitree[[i]])
    new_vector <- append(get(vname) ,x)
    assign(vname, new_vector, envir = .GlobalEnv)
  }
}

# rename tips in tree that dendropy appended a "_1" to
rename_tips <- function(tree) {
  for (i in 1:length(tree)){
    tree[[i]]$tip.label <- str_remove(tree[[i]]$tip.label,"_1")}
  return(tree)
}

######################################################

setwd("~/Research/plastome paper/data")

# read in trees
cp_tree <- read.tree("cp_tree_renamed.tre") # chloroplast tree
sp_tree <- read.tree("MPEST_accession.tre") # species tree (Stephens et al 2015)
gene_trees <- read.tree("gene_trees_renamed.tre") # gene trees (Stephens et al 2015)
simtrees1 <- read.tree("simtrees.1.newick") # simulated gene trees, branch lengths unscaled
simtrees2 <- read.tree("simtrees.2.newick") # simulated gene trees, branch lengths scaled by 2
simtrees4 <- read.tree("simtrees.4.newick") # simulated gene trees, branch lengths scaled by 4

# remove underscore from simtrees/organellar trees
simtrees1 <- rename_tips(simtrees1)
simtrees2 <- rename_tips(simtrees2)
simtrees4 <- rename_tips(simtrees4)

# compare species tree and chloroplast tree to all 
multi_RF(cp_tree,gene_trees)
multi_RF(sp_tree,gene_trees)
multi_RF(cp_tree,simtrees1)
multi_RF(sp_tree,simtrees1)
multi_RF(cp_tree,simtrees2)
multi_RF(sp_tree,simtrees2)
multi_RF(cp_tree,simtrees4)
multi_RF(sp_tree,simtrees4)

################################ Poster plot ######################################
# get RF dist between sp_tree and cp_tree
sp2cp<-rf.dist.tip.changes(sp_tree,cp_tree)

sp_sim4_distro <- data.frame(Tree_distance=simtrees4_sp_tree_dist, Tree = "Species tree")

p <- ggplot(sp_sim4_distro, aes(x=Tree_distance)) +
  geom_histogram(bins=80,alpha=0.3, color = "black", size = 0.05,fill="blue",position = "identity") +
  geom_vline(aes(xintercept=sp2cp),
              color="red", size=0.5)+
  xlab("Information theoretic Robinson-Foulds distance")+
  ggtitle("Distance from species tree")+
  theme(axis.title = element_text(size = 6),
        title = element_text(size=8),
        # axis.title.y = element_blank(),
        axis.text = element_text(size = 5),
        panel.background = element_rect(fill = "white", color = "grey"))
  # theme_classic()

#p + expand_limits(x=1)
p + expand_limits(x = c(0, 1)) + coord_cartesian(ylim=c(6.5,130))
ggsave("../Figures/Final figures/Figure4.svg", width=80, height = 55, units = "mm", dpi = 300)


################################ Plot distributions #####################################
# get RF dist between sp_tree and cp_tree
sp2cp<-rf.dist.tip.changes(sp_tree,cp_tree)

# make a data frame from the distribution of RF distances
sp_gene_distro <- data.frame(Tree_distance=gene_trees_sp_tree_dist, Tree = "Species tree", ref="Gene trees")
cp_gene_distro <- data.frame(Tree_distance=gene_trees_cp_tree_dist, Tree = "Plastome tree",ref="Gene trees")
sp_sim1_distro <- data.frame(Tree_distance=simtrees1_sp_tree_dist, Tree = "Species tree",ref="Simulated trees (branch_lengths*1)")
cp_sim1_distro <- data.frame(Tree_distance=simtrees1_cp_tree_dist, Tree = "Plastome tree",ref="Simulated trees (branch_lengths*1)")
sp_sim2_distro <- data.frame(Tree_distance=simtrees2_sp_tree_dist, Tree = "Species tree",ref="Simulated trees (branch_lengths*2)")
cp_sim2_distro <- data.frame(Tree_distance=simtrees2_cp_tree_dist, Tree = "Plastome tree",ref="Simulated trees (branch_lengths*2)")
sp_sim4_distro <- data.frame(Tree_distance=simtrees4_sp_tree_dist, Tree = "Species tree",ref="Simulated trees (branch_lengths*4)")
cp_sim4_distro <- data.frame(Tree_distance=simtrees4_cp_tree_dist, Tree = "Plastome tree",ref="Simulated trees (branch_lengths*4)")

distro <- rbind(sp_sim1_distro,cp_sim1_distro,sp_sim2_distro,cp_sim2_distro,sp_sim4_distro,cp_sim4_distro)

# plot all distributions together
ggplot(distro, aes(x=Tree_distance,fill=Tree,color=Tree)) +
  geom_histogram(bins=30,alpha=0.3,position = "identity") +
  facet_wrap(~fct_rev(ref),nrow=3,)+
  geom_vline(aes(xintercept=sp2cp),
             color="blue", linetype="dashed", size=1)+
  xlab("Normalized Robinson-Foulds distance")+
  theme(legend.title=element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.background = element_rect(colour = "black", fill = "white"))

ggsave("3branchlengths.svg")

# plot only gene tree simulations
ggplot(distro_sim, aes(x=Tree_distance,fill=Tree,color=Tree)) +
  geom_histogram(bins=30,alpha=0.3,position = "identity") +
  geom_vline(aes(xintercept=sp2cp),
             color="blue", linetype="dashed", size=1)+
  theme(legend.title=element_blank(),
        # text = element_text(size = 20),
        # axis.text = element_text(size = 10),
        # axis.title = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.background = element_rect(colour = "black", fill = "white"))

ggplot(distro_organellar, aes(x=Tree_distance,fill=Tree,color=Tree)) +
  geom_histogram(bins=30,alpha=0.3,position = "identity") +
  geom_vline(aes(xintercept=sp2cp),
             color="blue", linetype="dashed", size=1)+
  theme(legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.background = element_rect(colour = "black", fill = "white"))

#compare individual gene trees w/ simulated trees
single_gene_tree <- gene_trees[[5]]
multi_RF(single_gene_tree,simtrees_subset)

sg_sim_distro <- data.frame(RF=simtrees_subset_single_gene_tree_dist, Tree = "Single gene vs. Simulated gene trees")
distro_sg_plastome <- rbind(cp_sim_distro,sg_sim_distro)

ggplot(distro_sg_plastome, aes(x=RF,fill=Tree,color=Tree)) +
  geom_histogram(bins=30,fill="white",alpha=0.5,position = "identity")

hist(sg_sim_distro, xlim=c(50,160), breaks=20)

################################ Plot 4 lengths seperately #####################################
distro_4 <- rbind(sp_sim1_distro,sp_sim2_distro,sp_sim4_distro)

p <- ggplot(distro_4, aes(x=Tree_distance)) +
geom_histogram(bins=80,alpha=0.3, color = "black", size = 0.05,fill="blue",position = "identity") +
  geom_vline(aes(xintercept=sp2cp),
             color="red", size=0.5)+
  facet_wrap(~fct_rev(ref),nrow=3,)+
  xlab("Information theoretic Robinson-Foulds distance")+
  # ggtitle("Distance from species tree")+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = "grey"))
# theme_classic()

#p + expand_limits(x=1)
p + expand_limits(x = c(0, 1))# + coord_cartesian(ylim=c(6.5,130))

ggsave("../Supplemental/Supplemental_Figure_1.png", width=185, height = 190, units = "mm", dpi = 300)

################################ T-tests #####################################
#t test
t_test4 <- t.test(cp_sim4_distro$Tree_distance, mu = sp2cp)
t_test4

t_test2 <- t.test(cp_sim2_distro$Tree_distance, mu = sp2cp)
t_test2

t_test1 <- t.test(cp_sim1_distro$Tree_distance, mu = sp2cp)
t_test1

########################### Apply gene tree IDs to distances ############################


#make a list of IDs for the gene trees
filenames <- list.files("./genetrees_IQ/", pattern="*", full.names=TRUE)
gene_tree_ID <- substr(filenames, 16,nchar(filenames)-22)

df <- data.frame(ID=gene_tree_ID,Sp_gene=gene_trees_sp_tree_dist,Cp_gene=gene_trees_cp_tree_dist)

write.csv(df,file = "generalized.RF.distances.norm.csv")

# testing to see if gene_tree_ID is in the correct order; rename.tips.phylo is from file "rename_tips.R"
test_tr1 <- read.tree(file = "genetrees_IQ/104.fasta.treefile.newick")
test_tr1 <- rename.tips.phylo(test_tr1)
rf.dist.tip.changes(test_tr1, cp_tree)

test_tr2 <- read.tree(file = "genetrees_IQ/9.fasta.treefile.newick")
test_tr2 <- rename.tips.phylo(test_tr2)
rf.dist.tip.changes(test_tr2, cp_tree)

rf.dist.tip.changes(gene_trees[[2]], cp_tree)
rf.dist.tip.changes(gene_trees[[2]], sp_tree)

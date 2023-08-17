library(phangorn)
library(ggplot2)
setwd("~/Research/plastome paper/data")

######################################################

multi_RF <- function(tree) {
  # make a variable name from the input. e.g. if sp_tree is the input, vname <- "sp_tree_dist"
  vname <- paste(deparse(substitute(tree)),"dist", sep="_")
  # assign vname to the global environment
  assign(vname, c(), envir = .GlobalEnv)
  # calculate Robinson-Fould distances for each simulated tree vs. the input tree
  for (i in 1:1000) {
    x <- RF.dist(tree, simtrees[[i]])
    new_vector <- append(get(vname) ,x)
    assign(vname, new_vector, envir = .GlobalEnv)
  }
}


# read trees
sp_tree <- read.tree("MPEST_accession.tre")
simtrees <- read.tree("simtrees.newick")
cp_tree <- read.tree("cp_tree_pruned2.tre")

# run custom function
multi_RF(sp_tree)
multi_RF(cp_tree)

# get RF dist between sp_tree and cp_tree
sp2cp<-RF.dist(sp_tree,cp_tree)

# make a data frame from the distribution of RF distances
sp_distro <- data.frame(RF=sp_tree_dist, Tree = "Species")
cp_distro <- data.frame(RF=cp_tree_dist, Tree = "Plastome")


distro <- rbind(sp_distro,cp_distro)

# plot both distributions together
ggplot(distro, aes(x=RF,fill=Tree,color=Tree)) + 
  geom_histogram(binwidth = 2) + 
  geom_vline(aes(xintercept=sp2cp),
             color="blue", linetype="dashed", size=1)

# plot distributions separate
hist(sp_tree_dist, xlim=c(60,180), breaks=20)
hist(cp_tree_dist, xlim=c(60,180), breaks=8)



# treedist(sp_tree, cp_tree)
# RF.dist(sp_tree, cp_tree)
# 
# RF.dist(sp_tree, simtrees[[1]])
# 
# dist_sp <- c()
# 
# for (i in 1:1000) {
#   x <- RF.dist(sp_tree, simtrees[[i]])
#   dist_sp <- append(dist_sp,x)
# }
# 
# hist(dist_sp)

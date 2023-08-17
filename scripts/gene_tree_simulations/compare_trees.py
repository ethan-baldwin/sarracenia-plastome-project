# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 14:37:05 2021

@author: ethan
"""

import dendropy
from dendropy.calculate import treecompare
from dendropy.model import reconcile

sp_tree = dendropy.Tree.get(path="MPEST_accession.tre", schema="newick")
simtrees = dendropy.TreeList.get(path="simtrees.newick", schema="newick", taxon_namespace=sp_tree.taxon_namespace)
cp_tree = dendropy.Tree.get(path="chloroplast_RAxML.tre", schema="newick", taxon_namespace=sp_tree.taxon_namespace)
cp_tree.retain_taxa(sp_tree.taxon_namespace)


cp_tree_pruned = dendropy.Tree.get(path="cp_tree_pruned2.tre", schema="newick", taxon_namespace=sp_tree.taxon_namespace)

distances_sp = []

#######################################################################################
print(treecompare.symmetric_difference(sp_tree, cp_tree_pruned))
print(treecompare.symmetric_difference(sp_tree, simtrees[1]))

for tree in simtrees:
    distance = treecompare.symmetric_difference(tree, sp_tree)
    distances_sp.append(distance)
 
 ########################################################################################   
 
 

# taxon set association
genes_to_species = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=sp_tree.taxon_namespace,
        num_contained=8)

# convert to containing tree
sp_tree = reconcile.ContainingTree(sp_tree,
            contained_taxon_namespace=genes_to_species.domain_taxon_namespace,
            contained_to_containing_taxon_map=genes_to_species)

for i in range(50):
    sp_tree.embed_contained_kingman(default_pop_size=1)

# embed simulated trees
for tree in simtrees:
    sp_tree.embed_tree(tree)

deep_coals = sp_tree.deep_coalescences()
stepwise_out = open("deeep_coal.txt", "w")
for tree in deep_coals:
    stepwise_out.write("%d\n" % deep_coals[tree])
print(dendropy.model.reconcile.reconciliation_discordance(simtrees[0]), sp_tree)

# for edge in cp_tree.postorder_edge_iter():
#     edge.length = 1
    
for tree in sp_tree.contained_trees:
            for edge in sp_tree.postorder_edge_iter():
                print(len(edge.tail_contained_edges[tree]))
              #  if edge.tail_node is None and sp_tree.ignore_root_deep_coalescences:
              #      continue
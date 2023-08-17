# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 20:20:47 2021

@author: ethan
"""
#&&
import dendropy
import os
from dendropy.simulate import treesim

os.chdir('C://Users//ethan//OneDrive - University of Georgia//Documents//Research//plastome paper//data')
# read in species tree
tree1 = dendropy.Tree.get(path="MPEST_accession.tre",
                          schema="newick")

# create empty TreeList objects
simtrees1 = dendropy.TreeList()
simtrees2 = dendropy.TreeList()
simtrees4 = dendropy.TreeList()

# create taxon name space common to simulated trees and species tree
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
    containing_taxon_namespace=tree1.taxon_namespace, 
    num_contained=1)

#&&
# simlate 1000 trees and store them in simtrees
for i in range(1000):
    gene_tree = treesim.contained_coalescent_tree(
        containing_tree=tree1, 
        gene_to_containing_taxon_map=gene_to_species_map)
    simtrees1.append(gene_tree)

# multiply branch lengths by 2
for edge in tree1.postorder_edge_iter():
    if type(edge.length) == float:
        edge.length = edge.length*2

#&&
# simlate 1000 trees and store them in simtrees
for i in range(1000):
    gene_tree = treesim.contained_coalescent_tree(
        containing_tree=tree1, 
        gene_to_containing_taxon_map=gene_to_species_map)
    simtrees2.append(gene_tree)

# multiply branch lengths by 2
for edge in tree1.postorder_edge_iter():
    if type(edge.length) == float:
        edge.length = edge.length*2

#&&
# simlate 1000 trees and store them in simtrees
for i in range(1000):
    gene_tree = treesim.contained_coalescent_tree(
        containing_tree=tree1, 
        gene_to_containing_taxon_map=gene_to_species_map)
    simtrees4.append(gene_tree)
    
#&&
# write simulated trees to file
simtrees1.write(
    path="simtrees.1.newick",
    schema="newick",
    )

#&&
# write simulated trees to file
simtrees2.write(
    path="simtrees.2.newick",
    schema="newick",
    )

#&&
# write simulated trees to file
simtrees4.write(
    path="simtrees.4.newick",
    schema="newick",
    )
# for tree in simtrees:
#     for edge in tree.postorder_edge_iter():
#         if edge.length<0:
#             print(edge.length)

# print(simtrees[0].as_string(schema='newick'))
# print(simtrees[0].as_ascii_plot())

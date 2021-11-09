#!bin/bash python3

from ete3 import Tree
from Bio import Phylo

#pwdTree = '02_download_ENAEBI_data_amphibian/qiime_output/intermediate/vsearch_97_denovo/ASV50-insertion-tree/tree.nwk'
#pwdOut = '02_download_ENAEBI_data_amphibian/qiime_output/intermediate/vsearch_97_denovo/ASV50-insertion-tree/tree_UM.nwk'

pwdTree = 'qiime_output_deblur/vsearch_99/tree_99_ASV50/tree.nwk'
pwdOut = 'qiime_output_deblur/vsearch_99/tree_99_ASV50/tree_UM.nwk'

#t = Phylo.read(pwdTree, "newick")
t = Tree(pwdTree)

t.convert_to_ultrametric()

t.write(format=1, outfile=pwdOut)

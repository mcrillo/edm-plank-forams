
library(ape)
library(ggtree)

# reading the text file with the raw tree (without node dates / branching times)
tree <- read.tree(file="data/phylogeny/phylogeny_21May.txt")

# plotting tree
ggtree(tree, layout ='rectangular', ladderize=TRUE) + geom_tiplab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, color = "blue") +
  xlim(0,15)

# saving the tree.txt in other formats (.nex or .tre)
# write.nexus(tree, file = "data/phylogeny/kucera_tree.nex")
# write.tree(plankforams, file = "data/phylogeny/kucera_tree.tre")

tree <- compute.brlen(tree, 1) # sets all branch lengths equal to one

ssp_dist <- cophenetic(tree) # returns a symmetrical matrix with distances for each species pair


dist.nodes(tree)
branching.times(tree)

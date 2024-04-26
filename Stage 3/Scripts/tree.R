# Install package
BiocManager::install("ggtree")

# Load package
library(ggtree)

# Load data 
tree <- read.tree("RAxML_result.raxml.bestTree")

# Create a ggtree object
ggtree_obj <- ggtree(tree)

# Plot the tree with adjusted label size
ggtree_obj + geom_tiplab(size = 1.5) + 
  geom_tree(color = "blue", size = 0.5) 
library(ape)
library(RPANDA)

library(TreeSim)

seed <- 853297604
set.seed(seed)
tree_stdiv <- sim.bd.taxa(150, 1, lambda = 1, mu = 0.2, complete= F)[[1]]
tree_ext <- sim.bd.taxa(150, 1, lambda = 1, mu = 0.8, complete= F)[[1]]
tree_rad <- sim_ClaDS(1,0.2,new_lamb_law = "uniform",new_mu_law = "uniform",condition = "taxa", taxa_stop = 150,theta = 0.1,lamb_min = 4,lamb_max = 4,mu_min = 0.2,mu_max=0.2,nShiftMax = 1)
par(mfrow = c(1,3))
plot.phylo(tree_stdiv)
plot.phylo(tree_ext)
plot.phylo(tree_rad[[1]],edge.color = tree_rad$rates)

write.tree(tree_stdiv, "tree_l1_mu02.tre")
write.tree(tree_ext, "tree_l1_mu08.tre")
write.tree(tree_rad[[1]], "tree_l1_l4_mu02.tre")

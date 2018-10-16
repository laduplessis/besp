rm(list = ls())
library(ape)

# testApeTreeLikelihood 1
cat(sprintf("\nGeneral Bayesian Skyline compared to classical skyline plot in Ape: Small tree with homochronous sampling.\n"))
tree.test <- read.tree(text = "((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);") # load tree
sk1 <- skyline(tree.test)
cat(sprintf("%.12f", sk1$logL))


# testApeTreeLikelihood 2
cat(sprintf("\nGeneral Bayesian Skyline compared to classical skyline plot in Ape: Big tree with homochronous sampling.\n"))
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree
sk1      <- skyline(tree.hiv)
#cat(hivtree.newick)
#cat("\n")
cat(sprintf("%.12f", sk1$logL))


# testApeTreeLikelihood 3
cat(sprintf("\nGeneral Bayesian Skyline compared to generalized skyline plot in Ape: Big tree with homochronous sampling.\n"))
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree
cl2      <- collapsed.intervals(coalescent.intervals(tree.hiv),0.0119)
sk2      <- skyline(cl2)
groupsizes <- c()
for (i in 1:max(cl2$collapsed.interval)) {
    groupsizes <- c(groupsizes, sum(cl2$collapsed.interval == i))
}
names(groupsizes) <- 1:length(groupsizes)
print(groupsizes)

cat(sprintf("%.12f", sk2$logL))

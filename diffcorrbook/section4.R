source("./section3.R")

source("http://bioconductor.org/biocLite.R")
biocLite(c("pcaMethods", "multtest"))
install.packages("DiffCorr")
library(DiffCorr)

# subsection 4.1
dim(eset.GSE5632.sub)
dim(eset.GSE5630.sub)

data <- cbind(eset.GSE5632.sub, eset.GSE5630.sub)
## 66 flowers samples
hc.flowers <- cluster.molecule(data[, 1:66], method = "pearson", linkage = "average")
## 60 leaves samples
hc.leaves <- cluster.molecule(data[, 67:126], method = "pearson", linkage = "average")

g1 <- cutree(hc.flowers, h = 0.4)
g2 <- cutree(hc.leaves, h = 0.4)
table(g1[table(g1)!=1])
table(g1)
table(g2)

res1 <- get.eigen.molecule(data, groups = g1, whichgroups = c(1:10), methods = "svd", n=2)
res2 <- get.eigen.molecule(data, groups = g2, whichgroups = c(11:20), methods = "svd", n=2)

## visualizing module networks
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout=layout.fruchterman.reingold(gg1))
write.modules(g1, res1, outfile = "module1_list.txt")

gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout=layout.fruchterman.reingold(gg2))
write.modules(g2, res2, outfile = "module2_list.txt")

## network file export for Cytoscape
write.graph(gg1, "tmp1.ncol", format="ncol")
write.graph(gg2, "tmp2.ncol", format="ncol")
tmp1 <- read.table("tmp1.ncol")
tmp2 <- read.table("tmp2.ncol")
tmp1$V3 <- "pp"
tmp2$V3 <- "pp"
tmp1$V4 <- "gg1forcy"
tmp2$V4 <- "gg2forcy"
tmp1 <- tmp1[, c("V4", "V1", "V3", "V2")]
tmp2 <- tmp2[, c("V4", "V1", "V3", "V2")]
write.table(tmp1, file="gg1forcy.nnf", row.names=FALSE, col.names=FALSE)
write.table(tmp2, file="gg2forcy.nnf", row.names=FALSE, col.names=FALSE)
module1_list <- read.table("module1_list.txt", skip=1)
module2_list <- read.table("module2_list.txt", skip=1)
module1_list$V1 <- sub("^", "Module", module1_list$V1)
module2_list$V1 <- sub("^", "Module", module2_list$V1 - 10)
write.table(module1_list, file="gg1forcy.nnf", append=TRUE, row.names=FALSE, col.names=FALSE)
write.table(module2_list, file="gg2forcy.nnf", append=TRUE, row.names=FALSE, col.names=FALSE)

## plotting DiffCorr groups 21 and 24
plotDiffCorrGroup(data, g1, g2, 21, 24, 1:66, 67:126,
                  scale.center=TRUE, scale.scale=TRUE,
                  ylim=c(-5,5))

## export the results (FDR < 0.05)
comp.2.cc.fdr(output.file = "Transcript_DiffCorr_res.txt",
              data[, 1:66], data[, 67:126], threshold = 0.05)


# subsection 4.2
data(AraMetLeaves)
dim(AraMetLeaves)

colnames(AraMetLeaves)
comp.2.cc.fdr(output.file = "Met_DiffCorr_res.txt",
              log10(AraMetLeaves[, 1:17]), ## Col-0 (17 samples)
              log10(AraMetLeaves[, 18:37]), ## tt4 (20 samples)
              method = "pearson",
              threshold = 1.0)
# dependencies
source("http://bioconductor.org/biocLite.R")
biocLite(c("GEOquery", "affy", "genefilter", "GOstats", "ath1121501.db"))
install.packages(c("spatstat", "igraph"))


# subsection 3.1
library(GEOquery)
data <- getGEOSuppFiles("GSE5632")
untar("GSE5632/GSE5632_RAW.tar", exdir="GSE5632")
data <- getGEOSuppFiles("GSE5630")
untar("GSE5630/GSE5630_RAW.tar", exdir="GSE5630")


# subsection 3.2
library(affy)
tgt <- list.files("./GSE5630", pattern = ".CEL.gz", full.names = TRUE)
eset.GSE5630 <- justRMA(filenames = tgt)
dim(eset.GSE5630)

tgt2 <- list.files("./GSE5632", pattern = "*.CEL.gz", full.names = TRUE)
eset.GSE5632 <- justRMA(filenames = tgt2)
dim(eset.GSE5632)

## filtering probesets
rmv <- c(grep("AFFX", rownames(eset.GSE5632)), grep("s_at", rownames(eset.GSE5632)), grep("x_at", rownames(eset.GSE5632)))
eset.GSE5632 <- eset.GSE5632[-rmv, ]
dim(eset.GSE5632)

eset.GSE5630 <- eset.GSE5630[-rmv, ]
dim(eset.GSE5630)


# subsection 3.3
library(genefilter)

e.mat <- 2^exprs(eset.GSE5632)
## keep genes with cv between .5 and 10, and where 20% of samples had exprs. > 100
ffun <- filterfun(pOverA(0.2, 100), cv(0.5, 10))
filtered <- genefilter(e.mat, ffun)
## apply the filter, and put expression back on log scale
eset.GSE5632.sub <- log2(e.mat[filtered, ])
dim(eset.GSE5632.sub)

e.mat <- 2^exprs(eset.GSE5630)
ffun <- filterfun(pOverA(0.2, 100), cv(0.5, 10))
filtered <- genefilter(e.mat, ffun)
eset.GSE5630.sub <- log2(e.mat[filtered, ])
dim(eset.GSE5630.sub)

## common probesets between GSE5632 and GSE5630
comm <- intersect(rownames(eset.GSE5632.sub), rownames(eset.GSE5630.sub))
head(comm)
length(comm)

eset.GSE5632.sub <- eset.GSE5632.sub[comm, ] # flowers
eset.GSE5630.sub <- eset.GSE5630.sub[comm, ] # leaves
dim(eset.GSE5630.sub)
dim(eset.GSE5632.sub)

GSE5632.cor <- cor(t(eset.GSE5632.sub), method = "spearman")
GSE5630.cor <- cor(t(eset.GSE5630.sub), method = "spearman")

## heatmap visualization
library(spatstat)
par(mfrow=c(1,2))
plot(im(GSE5632.cor[nrow(GSE5632.cor):1, ]), col=cm.colors(256), main="GSE5632")
plot(im(GSE5630.cor[nrow(GSE5630.cor):1, ]), col=cm.colors(256), main="GSE5630")

library(igraph)
## co-expression networks with GSE5632, SCC >= 0.95
g1 <- graph.adjacency(GSE5632.cor, weighted = TRUE, mode = "lower")
g1 <- delete.edges(g1, E(g1)[weight < 0.95])
g1 <- igraph::simplify(g1, remove.multiple = TRUE, remove.loops = TRUE)
g1 <- delete.vertices(g1, which(igraph::degree(g1)<1))
plot(g1, vertex.size=3, edge.width=3, vertex.color=ifelse(igraph::degree(g1)>20, "Magenta", "Green"), vertex.label="", layout=layout.kamada.kawai)

## co-expression networks with GSE5630, SCC >= 0.95
g2 <- graph.adjacency(GSE5630.cor, weighted = TRUE, mode = "lower")
g2 <- delete.edges(g2, E(g2)[weight < 0.95])
g2 <- igraph::simplify(g2, remove.multiple = TRUE, remove.loops = TRUE)
g2 <- delete.vertices(g2, which(igraph::degree(g2)<1))
plot(g2, vertex.size=3, edge.width=3, vertex.color=ifelse(igraph::degree(g2)>20, "Magenta", "Green"), vertex.label="", layout=layout.kamada.kawai)

## save igraph object as GML format (Cytoscape can import GML)
write.graph(g1, "g1forcy.gml", format = "gml")
write.graph(g2, "g2forcy.gml", format = "gml")


# subsection 3.4
g1.fc <- fastgreedy.community(g1)
sizes(g1.fc)

g2.fc <- fastgreedy.community(g2)
sizes(g2.fc)

## accessing module 1 in GSE5632
mod1 <- membership(g1.fc)[membership(g1.fc)==1]
## extracting probeset names in module 1
mod1.p <- names(mod1)
## accessing module 2 in GSE5632
mod2 <- membership(g1.fc)[membership(g1.fc)==2]
mod2.p <- names(mod2)
## accessing module 3 in GSE5632
mod3 <- membership(g1.fc)[membership(g1.fc)==3]
mod3.p <- names(mod3)


# subsection 3.5
library(GOstats)
library(GO.db)
library(ath1121501.db)

ls("package:ath1121501.db")
x <- ath1121501ACCNUM
mapped.probes <- mappedkeys(x)
length(mapped.probes)
geneUniv <- AnnotationDbi::as.list(x[mapped.probes])

## target probes
mod1.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod1.p])))
mod2.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod2.p])))
mod3.p.gene <- unique(unlist(AnnotationDbi::as.list(x[mod3.p])))

## mod1
hgCutoff <- 0.0001
params <- new("GOHyperGParams",
              geneIds=mod1.p.gene,
              universeGeneIds=geneUniv,
              annotation="ath1121501",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")
hgOver <- hyperGTest(params)
df <- summary(hgOver)
names(df)
pvalues(hgOver)[1:3]
htmlReport(hgOver, file = "res_mod1.html")

## mod2
params <- new("GOHyperGParams",
              geneIds=mod2.p.gene,
              universeGeneIds=geneUniv,
              annotation="ath1121501",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")
hgOver <- hyperGTest(params)
htmlReport(hgOver, file = "res_mod2.html")

## mod3
params <- new("GOHyperGParams",
              geneIds=mod3.p.gene,
              universeGeneIds=geneUniv,
              annotation="ath1121501",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")
hgOver <- hyperGTest(params)
htmlReport(hgOver, file = "res_mod3.html")
# RNAseq analysis methods
Part 2

Derived from https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

## Open RStudio 
Install edgeR. This only needs to be done once. No need to do this on the ConGen server, there is a global verion of it! 
```
BiocManager::install("edgeR")
```

## Load the edgeR package. 
We will use edgeR to identify significantly differentially expressed genes between our two groups based on counts of gene expression for each gene. 
```
library(edgeR)
```

## Set working directory
```
setwd("instructor_materials/Joanna_Kelley/2022/")
```

## Read in gene count matrix 
```
counts = read.csv("gene_count_matrix.csv", header = T, row = 1)
```

## Look at the beginning of the file (note, there are lots of zeros)
```
head(counts)
```

## Make sure the number of genes and samples is expected
```
dim(counts)
```

## Create map of samples to treatments 
```
samples = colnames(counts)
group = rep(c("Ac", "Hi", "Hy"), 6)
sample.info = data.frame(cbind(samples,group))
head(sample.info)
```

## Remove genes with zero counts

```
length(which(rowSums(counts) == 0))
gene.counts = counts[rowSums(counts) != 0,]
```

## Determine the number of genes left after removing zeros
```
dim(gene.counts)
```

## Create a DGEList of the data
```
y <- DGEList(counts=gene.counts,group=sample.info$group)

```

## Calculate the normalization factors 
From edgeR manual: The calcNormFactors function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples [30]. We call the product of the original library size and the scaling factor the effective lib

```
y <- calcNormFactors(y)
```

## Look at y
```
y
```

## Visualize the data in y, for example using an MDS plot
```
plotMDS(y)
```

## Design matrix
```
design <- model.matrix(~group)
design
```

## Estimate dispersion 
```
y <- estimateDisp(y,design)
```

## perform quasi-likelihood F-tests:
```
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
```

## summarize the results
```
summary(de <- decideTestsDGE(qlf, p.value = 0.05))
```

## plot the results
```
detags <- rownames(y)[as.logical(de)]
plotSmear(qlf, de.tags=detags)
abline(h = c(-1, 1), col = "blue")
```

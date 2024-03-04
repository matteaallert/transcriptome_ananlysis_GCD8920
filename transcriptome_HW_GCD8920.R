if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
library(DESeq2)
library(dplyr)

onetrial <- function(replicates){
  # Picks a random sample of size $replicates for wild type and compares it to
  # another sample of size $replicates for knockout mice, returns DESeq differential
  # expression analysis results
  
  # Load the files
  wild_counts <- as.matrix(read.csv("WT_raw.tsv",sep="\t",row.names=1))
  ko_counts <- as.matrix(read.csv("Snf2_raw.tsv",sep="\t",row.names=1))

  
  # Pick random samples (WITHOUT replacement) of size $replicates from each population
  wild_samples <- wild_counts[,sample(seq(from=1, to=48), replicates)]
  ko_samples <- ko_counts[,sample(seq(from=1, to=48), replicates)]
  
  # make one big matrix with them
  all <- cbind(wild_samples, ko_samples)
  colnames(all) <- NULL
  
  # label each row with which type it is
  coldata <- matrix(c(rep('wild',each=replicates),rep('knockout',each=replicates)))
  colnames(coldata) <- c('strain')

  # DESeq analysis!
  dds <- DESeqDataSetFromMatrix(countData = all, colData = coldata, design = ~strain)
  runit <- DESeq(dds)
  res <- results(runit)
  res <- results(runit, contrast=c("wild","knockout"))
  
  # add a column to the results indicating how many replicates the gene's
  # analysis is based on
  labels = matrix(c(rep(replicates,each=7125)))
  colnames(labels) <- c('replicates')
  
  # trim up the results and send them back
  answer <- cbind(labels, res[,c(2,4,6)])
  answer$gene <- rownames(res)
  return(answer)
}

sampleseqrepeats <- function(replicates, repeats){
  # Performs repeated DESeq analyses using a subsample of the
  # population
  # Inputs:
  # - replicates: How many of each population to grab for the analysis
  # - repeats: How many times to repeat the analysis
  print("First one!")
  sampled <- onetrial(replicates)
  for(val in seq(from=1, to=repeats-1)) {
    print("Beginning new repeat...")
    sampled <- rbind(sampled, onetrial(24))
  }
  return(sampled)
}

samplemeans <- function(results){
  # takes the results of the "sampleseqrepeats" function and 
  # uses the repeated tests to generate a mean log2foldChange value
  # for each gene.
  means <- 0
  for(gene in unique(results$gene)) {
    # then grab the log-fold change values that are left:
    changes = results[results$gene == gene, ]$log2FoldChange
    f <- c(gene, as.numeric(mean(changes))) # no idea why this wouldn't be a numeric automatically
    
    if(class(means) == "numeric") {
      means <- f # if it's the first entry, then the only row is this one
    } else {
      means <- rbind(means, f)
    }
  }
  rownames(means) <- NULL
  colnames(means) <- c('gene','average')
  means <- as.data.frame(means)
  return(means)
}

# Run DESeq with 24 replicates of each type, and repeat the comparison 5 times:
results <- sampleseqrepeats(24,5)

# filter results by p value:
filtered <- results[is.na(results$padj) == FALSE, ]
filtered <- filtered[filtered$padj > 0.05, ]
resOrdered <- filtered[order(filtered$padj),]
resOrdered48 <- resOrdered[1:50,]
dim(resOrdered48)
# get the means for each gene:
logfold48 <- samplemeans(resOrdered48)
head(logfold48)

##heatmaps
# top 50 genes
#install.packages("limma")
#library(limma)
install.packages("pheatmap")
library(pheatmap)
pheatmap(logfold48,
         clustering_distance_rows=euclidean,
         clustering_distance_cols=euclidean)
distance = dist(mat_data, method = "manhattan")
write.csv(logfold48, "logfold48.csv")
logfold48 = logfold48[,2:3]
install.packages("ggplot2")
library(ggplot2)

mine.heatmap <- ggplot(data = resOrdered48, mapping = aes(x = replicates ,y = gene,
                                                       fill = log2FoldChange)) +
  geom_tile() +
  xlab(label = "Sample")

mine.heatmap
heatmap((logfold48))


dds <- DESeqDataSetFromMatrix(countData = all, 
                              colData = coldata, 
                              design = ~Strain)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Strain_WildType_vs_Mutant")
plotMA(res, ylim=c(-2,2))

install.packages("pheatmap")
library(pheatmap)
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("Strain")])
all <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
all
ncol(All36) == nrow(coldata36)
dim(All36)
dim(coldata36)
dds <- DESeqDataSetFromMatrix(countData = All36, 
                              colData = coldata36, 
                              design = ~V1)
dim(coldata)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Strain_WildType_vs_Mutant")
plotMA(res, ylim=c(-2,2))


wild_counts <- as.matrix(read.csv("WT_raw.tsv",sep="\t",row.names=1))
#wild_counts36 <- wild_counts[,-2:-13]
#wild_counts24 <- wild_counts[,-2:-25]
#wild_counts12 <- wild_counts[,-2:-37]
wild_counts3 <- wild_counts[,-2:-46]
dim(wild_counts3)
ko_counts <- as.matrix(read.csv("Snf2_raw.tsv",sep="\t",row.names=1))
#ko_counts36 <- ko_counts[,-2:-13]
#ko_counts24 <- ko_counts[,-2:-25]
#ko_counts12 <- ko_counts[,-2:-37]
ko_counts3 <- ko_counts[,-2:-46]
dim(ko_counts)
head(ko_counts)
# make one big matrix with them
#all <- cbind(wild_counts, ko_counts)
#all <- cbind(wild_counts36, ko_counts36)
all <- cbind(wild_counts24, ko_counts24)
all <- cbind(wild_counts12, ko_counts12)
all <- cbind(wild_counts3, ko_counts3)
colnames(all) <- NULL
dim(all)
# label each row with which type it is
#coldata <- matrix(c(rep('wild',each=36),rep('knockout',each=36)))
coldata <- matrix(c(rep('wild',each=3),rep('knockout',each=3)))
colnames(coldata) <- c('strain')
dim(coldata)

# DESeq analysis!
dds <- DESeqDataSetFromMatrix(countData = all, colData = coldata, design = ~strain)

all <- DESeq(dds)

all.log <- rlog(all, blind=F)

library("pheatmap")
# use the log transform on the data set
topVarianceGenes <- head(order(rowVars(assay(all.log)), decreasing=T),50) 
all.matrix <- assay(all.log)[ topVarianceGenes, ]
all.matrix <- all.matrix - rowMeans(all.matrix)

# select the 'contrast' you want
all.annotation_data <- as.data.frame(colData(all.log)[c("strain")])
all.plot <- pheatmap(all.matrix, annotation_col=all.annotation_data)

all.plot36 <- pheatmap(all.matrix, annotation_col=all.annotation_data, cluster_rows = FALSE)
all.plot24 <- pheatmap(all.matrix, annotation_col=all.annotation_data, cluster_rows = FALSE)
all.plot12 <- pheatmap(all.matrix, annotation_col=all.annotation_data, cluster_rows = FALSE)
all.plot2 <- pheatmap(all.matrix, annotation_col=all.annotation_data, cluster_rows = FALSE)




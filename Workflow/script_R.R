# Charger les librairies :
library("DESeq2")
library("ggplot2")
###### Recuperer le chemin (où s'exécute le script et où sont les fichiers dont on a besoin)
dir = getwd()
setwd(dir)

ind1 = read.table(snakemake@input[[1]])
ind2 = read.table(snakemake@input[[2]])
ind3 = read.table(snakemake@input[[3]])
ind4 = read.table(snakemake@input[[4]])
ind5 = read.table(snakemake@input[[5]])
ind6 = read.table(snakemake@input[[6]])
ind7 = read.table(snakemake@input[[7]])
ind8 = read.table(snakemake@input[[8]])

###### récupération des dernières colonnes et conversion en numeric :
counts = data.frame(as.numeric(ind1[-1,7]), as.numeric(ind2[-1,7]), as.numeric(ind3[-1,7])
                   , as.numeric(ind4[-1,7]), as.numeric(ind5[-1,7]), as.numeric(ind6[-1,7]),   
                   as.numeric(ind7[-1,7]),as.numeric(ind8[-1,7]))

counts = as.matrix(counts)

###### recuperation au format data.frame des labels :
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
type = c("M", "M", "WT", "WT", "WT", "WT", "WT","M")
labels2 = data.frame(labels, type)
rowlabels = ind1[-1,1]

###### ajout des labels des colonnes et lignes de la matrice :
rownames(counts)=rowlabels
colnames(counts) = labels

###### summary :
head(counts)

###### Affichage d'une première ACP avec tous les individus :  

dds = DESeqDataSetFromMatrix(countData = counts, colData = labels2, ~type)
dds2 = DESeq(dds)
res <- results(dds2)
resLFC <- lfcShrink(dds2, coef = 2, type="normal")
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE) #normalisation rlog
pcaData <- plotPCA(rld, intgroup=c("type"), returnData=TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type)) + #construction du graphique
  geom_point(size=3) +
  geom_text(
    label=labels, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

##### Affichage d'une seconde ACP sans l'outlier : 
rld <- rld[,-8]
type <- type[-8]
labels <- labels[-8]
pcaData <- plotPCA(rld, intgroup=c("type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type)) + #construction du graphique
  geom_point(size=3) +
  geom_text(
    label=labels,
    nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

#####
#expression différentielle :
genes <- res[is.na(res$padj) == FALSE,]
genes <- genes[genes[, "padj"] < 0.05,]
int <- dim(genes)
#genes <- genes[abs(genes[, "log2FoldChange"]) > 7,]

head(res)
plotMA(resLFC, ylim=c(-5,5))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
indexes <- order(res$padj)[1:5]
#indexes <- row.names(genes)
for (value in indexes){
	plotCounts(dds, gene= value, intgroup="type") #which.min(res$padj)
}

#save :
write.table(int, file = "dim.txt")
write.table(genes, file = snakemake@output[[1]], sep = "\t",
            row.names = TRUE,col.name = TRUE)

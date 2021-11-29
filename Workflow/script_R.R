# Charger des librairies(ggplot2 est installée avec DESeq2, puisque DESeq2 lui fait appel) :
library("DESeq2")
library("ggplot2")
# Acquisition du chemin actuel
dir = getwd()
setwd(dir)
# Récupération des séquences d'ARN en tant que colonnes de la matrice de compte
ind1 = read.table(snakemake@input[[1]])
ind2 = read.table(snakemake@input[[2]])
ind3 = read.table(snakemake@input[[3]])
ind4 = read.table(snakemake@input[[4]])
ind5 = read.table(snakemake@input[[5]])
ind6 = read.table(snakemake@input[[6]])
ind7 = read.table(snakemake@input[[7]])
ind8 = read.table(snakemake@input[[8]])

# Suppression des noms et conversion en numeric :
counts = data.frame(as.numeric(ind1[-1,7]), as.numeric(ind2[-1,7]), as.numeric(ind3[-1,7])
                   , as.numeric(ind4[-1,7]), as.numeric(ind5[-1,7]), as.numeric(ind6[-1,7]),   
                   as.numeric(ind7[-1,7]),as.numeric(ind8[-1,7]))

#construction de la matrice de compte :
counts = as.matrix(counts)

# construction des dataframe de labels :
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
type = c("M", "M", "WT", "WT", "WT", "WT", "WT","M")
labels2 = data.frame(labels, type)
rowlabels = ind1[-1,1]

# ajout des labels des colonnes et lignes de la matrice :
rownames(counts)=rowlabels
colnames(counts) = labels

# affichage des premières lignes de la matrice de counts
head(counts)

# Affichage d'une première ACP avec tous les individus :  

dds = DESeqDataSetFromMatrix(countData = counts, colData = labels2, ~type)
dds2 = DESeq(dds)
res <- results(dds2)
resLFC <- lfcShrink(dds2, coef = 2, type="normal") #réduction normale des données, pour éliminer le bruit
#vsd <- vst(dds, blind=FALSE) #autre normalisation
rld <- rlog(dds, blind=FALSE) #normalisation rlog

# Affichage d'une première ACP avec tous les individus : 
pcaData <- plotPCA(rld, intgroup=c("type"), returnData=TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=type)) + #paramétrages du graphique
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

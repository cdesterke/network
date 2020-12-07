## https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA

library(WGCNA)

library(flashClust)



matrix<-read.table("matrix.txt",h=T)
pheno<-read.table("trait.txt",h=T)


# Run this to check if there are gene outliers
gsg = goodSamplesGenes(matrix, verbose = 3)
gsg$allOK 


#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
#if (!gsg$allOK)
#   {if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(matrix)[!gsg$goodGenes], collapse= ", ")));
#       if (sum(!gsg$goodSamples)>0)
#           printFlush(paste("Removing samples:", paste(rownames(matrix)[!gsg$goodSamples], collapse=", ")))
#       matrix= matrix[gsg$goodSamples, gsg$goodGenes]
#       }






table(rownames(matrix)==rownames(pheno))

save(matrix, pheno, file="SamplesAndTraits.RData")



# Choose a soft threshold power- USE A SUPERCOMPUTER IRL ------------------------------------
  
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(matrix, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
  
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
 
#from this plot, we would choose a power of 18 becuase it's the lowest power for which the scale free topology index reaches 0.90



#build a adjacency "correlation" matrix
enableWGCNAThreads()
softPower = 7
adjacency = adjacency(matrix, power = softPower, type = "signed") #specify network type
head(adjacency)
 
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM
 





# Generate Modules --------------------------------------------------------
 # Generate Modules --------------------------------------------------------
 
 library(flashClust)
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(matrix, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")


## TOM graph
TOMplot(dissTOM,geneTree,dynamicColors)


## mds plot
cmd1=cmdscale(as.dist(dissTOM),2)

plot(cmd1,col=dynamicColors,main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")




  


#plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="clusterwithoutmodulecolors.pdf")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(matrix, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()



#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

# Correlate traits --------------------------------------------------------
 
 
#Define number of genes and samples
nGenes = ncol(matrix)
nSamples = nrow(matrix)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
 
#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
 
 
#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor,
            xLabels= names(pheno),
            yLabels= names(MEs),
            ySymbols= names(MEs),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrix,
            setStdMargins= FALSE,
            cex.text= 0.5,
            zlim= c(-1,1),
            main= paste("Module-trait relationships"))
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()


## gene names from turquoise modulle
names(matrix)[moduleColors=="turquoise"]
turquoise<-data.frame(names(matrix)[moduleColors=="turquoise"])
write.table(turquoise,file="module_turquoise.txt")

################# cytoscape export
# select modules
modules = c("turquoise")
# Select module probes
inModule=is.finite(match(dynamicColors,modules))
probes = names(matrix)
modProbes=probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
weighted = TRUE, threshold = 0.2,nodeNames=modProbes, nodeAttr = dynamicColors[inModule])

cyt$nodeData$nodeName

edge<-data.frame(cyt$edgeData)
node<-data.frame(cyt$nodeData)
write.table(edge,file="turquoise_edge.txt")
write.table(node,file="turquoise_node.txt")


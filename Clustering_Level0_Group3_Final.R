################################################
########## Defi R : Clustering Level 0 #########
################################################

# Objectifs : faire 4 clusters de regroupements de gènes
# 4 méthodologies : 
## -> Euclidienne + k-mean
## -> Euclidienne + HCL
## -> Correlation + k-mean
## -> Correlation + HCL

#-----------------------
# Ouvertures et vérifications des données
#-----------------------

## Ouverture du fichier
expData = as.matrix(read.csv("cell-cycle_SCERE_DUO.txt", sep = "\t"))
#expData = read.csv("cell-cycle_SCERE_DUO.txt", sep = "\t")

## Vérification qu'il n'y a pas de valeurs manquantes : 
sum(is.na(expData))
## Visualisation des données
#boxplot(expData)


#-----------------------
# Fonction plotGenes de Gaëlle
#-----------------------

# Gene profiles
plotGenes <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    yMax = max(expData)
    
  }
  
  # Representation of the first expression profile
  plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
       ylim = c(floor(yMin), ceiling(yMax)),
       xlab = "Time point", ylab = "Gene expression level",
       main = title)
  
  # Add expression profile for other genes
  for(i in 2:nrow(expData)){
    
    lines(1:ncol(expData), expData[i,], col = "grey")
    
    # end of for()  
  }
  
  # Average expression profile
  if(meanProfile == TRUE){
    expMean = apply(expData, 2, mean)
    lines(1:ncol(expData), expMean, col = "red", 
          lwd = 1.5, lty = "dashed")
  }
  
  # end of function plotGenes()  
}


#-----------------------
## Méthode des k-means & calcul par distance euclidienne
#-----------------------

## Matrice de distance euclidienne (/!\ -> mettre des crochets pour appeler
## le tableau avec la fonction dist())
matEuc = dist(expData[,])
sum(is.na(t(matEuc)))

## Clusterisation avec k-means
resKmeans_1 <- kmeans(expData, centers = 4)
table(resKmeans_1$cluster)

cluster1 <- expData[which(resKmeans_1$cluster == 1),]
cluster2 <- expData[which(resKmeans_1$cluster == 2),]
cluster3 <- expData[which(resKmeans_1$cluster == 3),]
cluster4 <- expData[which(resKmeans_1$cluster == 4),]

## Profil d'expression des gènes
par(mfrow=c(2,2))
plotGenes(cluster1)
plotGenes(cluster2)
plotGenes(cluster3)
plotGenes(cluster4)
mtext("Kmeans & Euclidian", side = 3, line = -2, outer = TRUE)

## Heatmap
allClusters = rbind(cluster1,cluster2, cluster3, cluster4)
colnames(allClusters) = c(seq(0,245,5))
heatmap(as.matrix(allClusters), Colv = NA, Rowv = NA, scale="row", col = heat.colors(500), 
        main="Kmeans & Euclidian", xlab="Temps (min)", ylab="Genes")

#-----------------------
## Méthode des k-means & calcul par corrélation
#-----------------------
## Calcul de la matrice de distance par méthode de corrélation :
matDist <- as.dist(1 - cor(t(expData)))
  
## Clusterisation par HCL :
resKmeans_2 <- kmeans(matDist, centers = 4)
table(resKmeans_2$cluster)

cluster1 <- expData[which(resKmeans_2$cluster == 1),]
cluster2 <- expData[which(resKmeans_2$cluster == 2),]
cluster3 <- expData[which(resKmeans_2$cluster == 3),]
cluster4 <- expData[which(resKmeans_2$cluster == 4),]

## Profil d'expression des gènes
par(mfrow=c(2,2))
plotGenes(cluster1)
plotGenes(cluster2)
plotGenes(cluster3) 
plotGenes(cluster4)
mtext("Kmeans & correlation", side = 3, line = -2, outer = TRUE)

## Heatmap
allClusters = rbind(cluster1,cluster2, cluster3, cluster4)
colnames(allClusters) = c(seq(0,245,5))
heatmap(as.matrix(allClusters), Colv = NA, Rowv = NA, scale="row", col = heat.colors(500), 
        main="Kmeans & correlation", xlab="Temps (min)", ylab="Genes")


#-----------------------
## Méthode des HCL, calcul des distances par méthode euclidienne
#-----------------------

## Clusterisation
resHCL_1 <- hclust(matEuc)

cluster1 <- expData[which(cutree(resHCL_1, k = 4) == 1),]
cluster2 <- expData[which(cutree(resHCL_1, k = 4) == 2),]
cluster3 <- expData[which(cutree(resHCL_1, k = 4) == 3),]
cluster4 <- expData[which(cutree(resHCL_1, k = 4) == 4),]

## Profil d'expression des gènes
par(mfrow=c(2,2))
plotGenes(cluster1)
plotGenes(cluster2)
plotGenes(cluster3)
plotGenes(cluster4)
mtext("HCL & Euclidian", side = 3, line = -2, outer = TRUE)

## Heatmap
allClusters = rbind(cluster1,cluster2, cluster3, cluster4)
colnames(allClusters) = c(seq(0,245,5))
heatmap(as.matrix(allClusters), Colv = NA, Rowv = NA, scale="row", col = heat.colors(500), 
        main="Kmeans & Euclidian", xlab="Temps (min)", ylab="Genes")


#-----------------------
## Méthode des HCL , calcul des distances par corrélation
#-----------------------

## Clusterisation
resHCL_2 <- hclust(matDist)

cluster1 <- expData[which(cutree(resHCL_2, k = 4) == 1),]
cluster2 <- expData[which(cutree(resHCL_2, k = 4) == 2),]
cluster3 <- expData[which(cutree(resHCL_2, k = 4) == 3),]
cluster4 <- expData[which(cutree(resHCL_2, k = 4) == 4),]

## Profil d'expression des gènes
par(mfrow=c(2,2))
plotGenes(cluster1)
plotGenes(cluster2)
plotGenes(cluster3)
plotGenes(cluster4)
mtext("HCL & correlation", side = 3, line = -2, outer = TRUE)

## Heatmap
allClusters = rbind(cluster1,cluster2, cluster3, cluster4)
colnames(allClusters) = c(seq(0,245,5))
heatmap(as.matrix(allClusters), Colv = NA, Rowv = NA, scale="row", col = heat.colors(500), 
        main="HCL & correlation", xlab="Temps (min)", ylab="Genes")

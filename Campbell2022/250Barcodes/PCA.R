master <- read.csv("MasterFrequencies_PCA.csv", row.names = 1)



masterscaled <- scale(master)

pca <- prcomp(t(masterscaled))

pcs <- data.frame(pca$x)

### Export the PCs to use in graphpad
write.csv(pcs, "pcs.csv")


plot(pcs$PC1, pcs$PC3)

a1 <- data.frame(pcs$PC1[1:5], pcs$PC3[1:5])
a2 <- data.frame(pcs$PC1[6:13], pcs$PC3[6:13])
a3 <- data.frame(pcs$PC1[14:21], pcs$PC3[14:21])
a4 <- data.frame(pcs$PC1[22:29], pcs$PC3[22:29])
a5 <- data.frame(pcs$PC1[30:37], pcs$PC3[30:37])

points(a1, col = "blue")
points(a2, col = "red")
points(a3, col = "green")
points(a4, col = "purple")
points(a5, col = "black")

summary(pca)

### Loading scores
rot <- pca$rotation




####### only feces, colon and cecum
feo_master <- read.csv("MasterFrequencies_PCA_FEO.csv", row.names = 1)

feo_masterscaled <- scale(feo_master)

feo_pca <- prcomp(t(feo_masterscaled))

feo_pcs <- data.frame(feo_pca$x)

plot(feo_pcs$PC1, feo_pcs$PC2)

feo_a1 <- data.frame(feo_pcs$PC1[1:4], feo_pcs$PC2[1:4])
feo_a2 <- data.frame(feo_pcs$PC1[5:8], feo_pcs$PC2[5:8])
feo_a3 <- data.frame(feo_pcs$PC1[9:12], feo_pcs$PC2[9:12])
feo_a4 <- data.frame(feo_pcs$PC1[13:16], feo_pcs$PC2[13:16])
feo_a5 <- data.frame(feo_pcs$PC1[17:20], feo_pcs$PC2[17:20])

points(feo_a1, col = "blue")
points(feo_a2, col = "red")
points(feo_a3, col = "green")
points(feo_a4, col = "purple")
points(feo_a5, col = "black")


### Loading scores
feo_rot <- feo_pca$rotation


### Export the PCs to use in graphpad
write.csv(feo_pcs, "FEO_pcs.csv")






####### only No liver (noL)
noL_master <- read.csv("MasterFrequencies_PCA_noL.csv", row.names = 1)

noL_masterscaled <- scale(noL_master)

noL_pca <- prcomp(t(noL_masterscaled))

noL_pcs <- data.frame(noL_pca$x)

plot(noL_pcs$PC1, noL_pcs$PC2)

noL_a1 <- data.frame(noL_pcs$PC1[1:7], noL_pcs$PC2[1:7])
noL_a2 <- data.frame(noL_pcs$PC1[8:14], noL_pcs$PC2[8:14])
noL_a3 <- data.frame(noL_pcs$PC1[15:21], noL_pcs$PC2[15:21])
noL_a4 <- data.frame(noL_pcs$PC1[22:29], noL_pcs$PC2[22:29])
noL_a5 <- data.frame(noL_pcs$PC1[30:37], noL_pcs$PC2[30:37])

points(noL_a1, col = "blue")
points(noL_a2, col = "red")
points(noL_a3, col = "green")
points(noL_a4, col = "purple")
points(noL_a5, col = "black")

summary(noL_pca)

### Loading scores
noL_rot <- noL_pca$rotation


### Export the PCs to use in graphpad
write.csv(noL_pcs, "noL_pcs.csv")



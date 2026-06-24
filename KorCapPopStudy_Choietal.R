##### README ==============================================================
##### This is R scripts for analyzing population genomics data of Korean Capsella bursa-pastoris
#### Title: Local Genetic Signatures of a Global Weed: Population Structure and Ecological Drivers of Capsella bursa-pastoris Diversity in the Korean Peninsula
#### Authors: Byungwook Choi, Eunsuk Kim, Peter Tiffin, and Sang-Tae Kim
#### Co-Corresponding Authors: Byungwook Choi and Sang-Tae Kim
#### Also written by B. Choi

##### A. Packaging libraries ==============================================
library(vcfR) # Tom Jenkins
library(promises)
library(adegenet)
library(hierfstat) ## basic.stats
library(graph4lg) # mat_gen_dist
library(dplyr) # mutate
library(tidyr) # pivot_longer
library(ggplot2) # ggplot
library(egg)
library(tidyverse) # for manipulating and plotting data
library(RColorBrewer) ## Tom Jenkins ## brewer.pal
library(reshape2) ## melt
library(poppr) ## as.genclone
library(ecodist) ## Mantel's test
library(polyRAD)
library(StAMPP)
library(mmod)
library(geosphere)
library(factoextra)
library(ggrepel) ## geom_text_repel
library(pegas)
library(LEA) ## For sNMF
library(yhat) ## CA
library(PopGenReport) ## MMRR
library(scatterpie)
library(ggmap)
library(maps)
library(mapdata)
library(rnaturalearth) ## for mapping
library(rnaturalearthdata) ## for mapping
library(tmap)
library(leaflet)
library(sf)
library(ggspatial)
library(prettymapr)
library(ggplot2)
library(ggspatial)
#####=====

##### B. Importing data ===================================================
#### B-1. Population (ID) data (187 individuals from 10 populations) ======
PopMap <- read.table("Populations_187.txt",
                     header=TRUE, sep="\t", stringsAsFactors = TRUE)
PopMap %>%
  group_by(Populations) %>%
  count()

#### B-2. Genomics data constructed from Stacks2 ==========================
#### Catalog.alleles.tsv built from Statcks2
#### alleles.tsv for each individual to match with the Catalog.alleles.tsv
#### Directly imported into PolyRAD pipeline from Stacks2
#### Alleles.tsv for all individuals was in KorCapRawData.zip
RAD_DNV <- polyRAD::readStacks("catalog.alleles.tsv",
                                 "/",
                                 min.ind.with.reads = 168,
                                 min.ind.with.minor.allele = 10,
                                 possiblePloidies = list(2),
                                 contamRate = 0.001,
                                 version = 2)

RAD_DNV 
#####=====

##### C. Quality control and parameter estimation ==========================
#### C-1. Test Overdispersion ==============================================
ovdisp_RAD_DNV <- 130

hindhe_RAD_DNV <- HindHe(RAD_DNV)
hindheByLoc_RAD_DNV <- colMeans(hindhe_RAD_DNV, na.rm = TRUE)
hist(hindheByLoc_RAD_DNV, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus (RAD_DNV)")
abline(v = 0.5, col = "blue", lwd = 2) ### The peak below 0.5 indicates well-behaved diploid loci.

ALFHWE_RAD_DNV <- AddAlleleFreqHWE(RAD_DNV)
theseloci_RAN_DNV <- GetLoci(ALFHWE_RAD_DNV)[ALFHWE_RAD_DNV$alleles2loc[ALFHWE_RAD_DNV$alleleFreq >= 0.05 & ALFHWE_RAD_DNV$alleleFreq < 0.5]]
theseloci_RAN_DNV <- unique(theseloci_RAN_DNV)
myhindheByLoc_RAD_DNV <- colMeans(hindhe_RAD_DNV[ALFHWE_RAD_DNV$taxaPloidy == 2, theseloci_RAN_DNV], na.rm = TRUE)
hist(myhindheByLoc_RAD_DNV, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus, MAF >= 0.05")
abline(v = 0.5, col = "blue", lwd = 2)

#### C-2. Filtering based on the Hind/He ===================================
hh_RAD_DNV <- HindHe(RAD_DNV, ploidy = RAD_DNV$possiblePloidies[[1]])
str(hh_RAD_DNV)

TotDepthT <- rowSums(RAD_DNV$locDepth) # depth for each sample
hhByInd <- rowMeans(hh_RAD_DNV, na.rm=TRUE) # Hind/He for each sample

plot(TotDepthT, hhByInd, log = "x",
     xlab="Depth", ylab="Hind/He", main="Sample")
abline(h=0.75, lty=2)

### C-2-1. Filtering by individuals =======================================
hhByInd[hhByInd < 0.45] ## Removal of six individuals; 76, 78, 95, 137, 156, 167

hh_RAD_DNV <- hh_RAD_DNV[hhByInd > 0.45,]
rownames(hh_RAD_DNV)

RADdata_RAD_DNV <- SubsetByTaxon(RAD_DNV, rownames(hh_RAD_DNV))
RADdata_RAD_DNV ## 4,681 loci and 181 taxa

### C-2-2. Filtering by loci ==============================================
hh2_RAD_DNV <- HindHe(RADdata_RAD_DNV, ploidy = RADdata_RAD_DNV$possiblePloidies[[1]])
hhByLoc <- colMeans(hh2_RAD_DNV, na.rm = TRUE)
hist(hhByLoc, breaks=30)

mean(hhByLoc)
ExpectedHindHe(RADdata_RAD_DNV, ploidy = 2, reps = 100,  
               errorRate = 0.001, contamRate = 0.001, overdispersion = ovdisp_RAD_DNV)

## Ploidy 2
# Mean Hind/He: 0.502
# Standard deviation: 0.0374
# 95% of observations are between 0.428 and 0.575

# UNPhased SNPs 
thresh1 <- 0.428
thresh2 <- 0.575

mean(hhByLoc < thresh1) # 0.126 bottom 12.6% of SNPs would be removed
mean(hhByLoc > thresh2) # 0.824 top 82.4% of SNPs would be removed

keeploci_di <- names(hhByLoc)[hhByLoc > thresh1 & hhByLoc < thresh2]

RADdata_RAD_DNV_filtered <- SubsetByLocus(RADdata_RAD_DNV, keeploci_di)
RADdata_RAD_DNV_filtered ## 233 SNPs

RADdata_RAD_DNV_filtered <- IterateHWE(RADdata_RAD_DNV_filtered)

### C-2-3. Converting into Genind format ==================================
RAD_DNV.ind <- Export_adegenet_genind(RADdata_RAD_DNV_filtered)
#####=====

##### D. Analysing Population statistics ===================================
#### D-1. Pre-step to analyze ==============================================
### D-1-1. Importing Population information ================================
POPIND_181 <- read.table("Populations_181.txt",
                         header=TRUE, sep="\t", stringsAsFactors = TRUE)

### D-1-2. Adding Population data into Genind ==============================
RAD_DNV.ind@pop <- POPIND_181$Populations
RAD_DNV.ind
pop(RAD_DNV.ind) 

### D-1-3. Conversion genind to hierfstat ==================================
fstat.RAD_DNV <- genind2hierfstat(RAD_DNV.ind)
head(fstat.RAD_DNV)

#### D-2. Calculating Population genetics statistics =======================
### D-2-1. FST =============================================================
RAD_DNV_FST <- hierfstat::genet.dist(fstat.RAD_DNV, method = "Fst", diploid = T)
RAD_DNV_FST

RAD_DNV.ind_FST <- mat_gen_dist(x = RAD_DNV.ind, dist = "FST", null_val = TRUE)
RAD_DNV.ind_FST

RAD_DNV.ind_LinFST <- mat_gen_dist(x = RAD_DNV.ind, dist = "FST_lin", null_val = TRUE)
RAD_DNV.ind_LinFST

### D-2-2. Diversity indices; Heterozygosity ===============================
## Calculate heterozygosity per site
H.RAD_DNV.ind = basic.stats(RAD_DNV.ind, diploid = TRUE)
H.RAD_DNV.ind

## Mean observed heterozygosity per site
Ho.RAD_DNV.ind = apply(H.RAD_DNV.ind$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
Ho.RAD_DNV.ind
summary(Ho.RAD_DNV.ind)

## Mean expected heterozygosity per site
He.RAD_DNV.ind = apply(H.RAD_DNV.ind$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
He.RAD_DNV.ind
summary(He.RAD_DNV.ind)

### D-2-3. Inbreeding coeffieicnt (FIS) ====================================
apply(H.RAD_DNV.ind$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

### D-2-4. Heatmap for Population differentiation ==========================
pop_order <- c(
  "BOP","DGR","MZS","TAB","JCH","JSN","NWS","DMY","CHJ","ICH")

fst_mat <- RAD_DNV.ind_FST

fst_mat <- fst_mat[pop_order, pop_order]

fst_mat[upper.tri(fst_mat, diag = TRUE)] <- NA

fst_long <- melt(
  fst_mat,
  varnames = c("Pop1", "Pop2"),
  value.name = "FST"
) %>%
  filter(!is.na(FST)) %>%
  mutate(
    Pop1 = factor(Pop1, levels = rev(pop_order)),
    Pop2 = factor(Pop2, levels = pop_order)
  )

fg1.htm <- ggplot(fst_long,
              aes(Pop2, Pop1, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.3f", FST)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "orange",
    high = "skyblue"
  ) +
  coord_fixed() +
  theme_classic(base_size = 7, base_family = "sans") +
  labs(x = xlab, y = ylab,
       fill = expression(F[ST]),
       title = "STACKS2-PolyRAD (233 SNPs)",
       subtitle = "a) FST")+
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1, size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_blank(),
    legend.key.height = unit(.5, "cm"),
    legend.key.width  = unit(.4, "cm"),
    legend.position = c(0.85, 0.7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=11, face="bold", colour = "Navy"),
    plot.subtitle = element_text(size=11, face="bold", colour = "black"))

fg1.htm

# ggsave("Heatmap_Pop.jpg", fg1.Htm, width = 10, height = 10, units = "cm", dpi = 600)

### D-2-5. PCA for Population differentiation ==========================
set.seed(300)
# Replace missing data with the mean allele frequencies
x = tab(RAD_DNV.ind, NA.method = "mean")
# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)
# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")
# Add a column containing individuals
ind_coords$Ind = indNames(RAD_DNV.ind)
# Add a column with the site IDs
ind_coords$Site = RAD_DNV.ind$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = brewer.pal(nPop(RAD_DNV.ind), "Set3")
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
fg1.pca <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0, color = "gray50", size = 0.7)+
  geom_vline(xintercept = 0, color = "gray50", size = 0.7)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 3, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab,
       title = "",
       subtitle = "b) PCA")+
  theme_article(base_size = 7, base_family = "sans")+
  theme(axis.title = element_text(face = "bold", vjust=0.5, size=9),
        axis.text = element_text(face = "bold", size=9),
        axis.text.x = element_text(angle = 0, size=9),
        axis.text.y = element_text(size=9),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=11, face="bold", colour = "Navy"),
        plot.subtitle = element_text(size=11, face="bold", colour = "black"))

fg1.pca

# ggsave("PCA_Pop.jpg", fg1.pca , width = 10, height = 10, units = "cm", dpi = 600)

### D-2-6. DAPC for Population differentiation ==========================
set.seed(300)
x = adegenet::tab(RAD_DNV.ind, NA.method = "mean")
crossval = xvalDapc(x, RAD_DNV.ind$pop, result = "groupMean", xval.plot = TRUE)
# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
crossval$`Number of PCs Achieving Highest Mean Success`
## [1] "60"
crossval$`Number of PCs Achieving Lowest MSE`
## [1] "60"
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)
# Run a DAPC using site IDs as priors
dapc1 = adegenet::dapc(RAD_DNV.ind, RAD_DNV.ind$pop, n.pca = numPCs, n.da = 2)
# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))
## Visualize DAPC results 1
scatter(dapc1)

## Visualize DAPC results 2
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)
# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2")
# Add a column containing individuals
ind_coords$Ind = indNames(RAD_DNV.ind)
# Add a column with the site IDs
ind_coords$Site = RAD_DNV.ind$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = brewer.pal(nPop(RAD_DNV.ind), "Set3")
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
fg1.dapc <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0, color = "gray50", size = 0.7)+
  geom_vline(xintercept = 0, color = "gray50", size = 0.7)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 3, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab,
       title = "",
       subtitle = "c) DAPC")+
  theme_article(base_size = 7, base_family = "sans")+
  theme(axis.title = element_text(face = "bold", size=9),
        axis.text = element_text(face = "bold", size=9),
        axis.text.x = element_text(angle = 0, size=9),
        axis.text.y = element_text(size=9),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=11, face="bold", colour = "Navy"),
        plot.subtitle = element_text(size=11, face="bold", colour = "black"))

fg1.dapc

# ggsave("DAPC_Pop.jpg", fg1.dapc , width = 10, height = 10, units = "cm", dpi = 600)

### D-2-7. Export Figure 1 =========================================
pdf("Figure1.pdf", width = 12, height = 4)
grid.arrange(fg1.htm, fg1.pca, fg1.dapc, nrow = 1, ncol = 3)
dev.off()

dev.off()
#####=====

##### E. AMOVA =====================================================
#### E-1. Adding information (Population / Altitudinal groups)
high_elevation <- c("BOP", "DGR", "MZS", "TAB")
mid_elevation  <- c("NWS", "JSN", "JCH")
low_elevation  <- c("CHJ", "ICH", "DMY")

POPIND_181$Elevation <- with(POPIND_181, ifelse(Populations %in% high_elevation, "High",
                                                ifelse(Populations %in% mid_elevation, "Mid",
                                                       ifelse(Populations %in% low_elevation, "Low", NA))))

POPIND_181$Elevation <- factor(POPIND_181$Elevation, levels = c("Low", "Mid", "High"))

strata(RAD_DNV.ind) <- POPIND_181

agc.RAD_DNV.ind <- as.genclone(RAD_DNV.ind)
agc.RAD_DNV.ind

#### E-1. For populations =========================================
amova.RAD_DNV <- poppr::poppr.amova(agc.RAD_DNV4ind, ~Populations,
                                     method = c("ade4"), within = FALSE, nperm = 9999,
                                     filter = TRUE, threshold = 0.05)
amova.RAD_DNV

### Randtest
amova.Test.RAD_DNV <- ade4::randtest(amova.RAD_DNV, nrepet = 9999) # Test for significance
print(amova.Test.RAD_DNV)

#### E-2. For Elevations =========================================
amova.RAD_DNV_2 <- poppr::poppr.amova(agc.RAD_DNV.ind, ~Elevation/Populations,
                                       method = c("ade4"), within = FALSE, nperm = 9999,
                                       filter = TRUE, threshold = 0.05)
amova.RAD_DNV_2

amova.Test.RAD_DNV_2 <- ade4::randtest(amova.RAD_DNV_2, nrepet = 9999) # Test for significance
print(amova.Test.RAD_DNV_2)
#####=====

##### F. Mantel's test for IBD and IBE ============================
##### https://stats.oarc.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
#### F-1. Importing dataset =======================================
Env.GEO <- read.csv("Mantel_IBD+IBE.CSV", header=T)

Env.GEO$Populations <- as.factor(as.character(Env.GEO$Populations))
Env.GEO$Altitudes <- as.factor(as.character(Env.GEO$Altitudes))

#### F-1-a. Normalization for Environmental variables
Env.GEO.scale <- as.data.frame(scale(Env.GEO[, -c(1:6, 35:36)]))

Env.GEO.scale <- cbind(Env.GEO[, 1:6], Env.GEO.scale, Env.GEO[, 35:36])
str(Env.GEO.scale)
#####=====
#### F-2. Calculation of distances ================================
### F-2-1. Geographical distance (Euclidean) ======================
#### Elevation, Latitude, Longitude, Delta.distance
### Elevation [,4]
GEO.Ele.mat <- as.matrix(Env.GEO.scale[,4])
GEO.Ele.mat.dist = ecodist::distance(GEO.Ele.mat, method="euclidean") %>%
  as.matrix()
rownames(GEO.Ele.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.Ele.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.Ele.mat.dist

### Latitudes [,5]
GEO.Lat.mat <- as.matrix(Env.GEO.scale[,5])
GEO.Lat.mat.dist = ecodist::distance(GEO.Lat.mat, method="euclidean") %>%
  as.matrix()
rownames(GEO.Lat.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.Lat.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.Lat.mat.dist

### Longitudes [,6]
GEO.Lon.mat <- as.matrix(Env.GEO.scale[,6])
GEO.Lon.mat.dist = ecodist::distance(GEO.Lon.mat, method="euclidean") %>%
  as.matrix()
rownames(GEO.Lon.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.Lon.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.Lon.mat.dist

### Distance (km)
Coords <- as.data.frame(Env.GEO.scale[,6:5])
dist_mat <- distm(Coords, fun = distHaversine)  # 또는 distGeo
dist_mat_km <- dist_mat / 1000
rownames(dist_mat_km) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(dist_mat_km) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.Dis.mat.dist_km <- dist_mat_km

### F-2-2. Environmental distance (Gower) ====-====================
## BIO1 to BIO19, MT, MHT, MLT, HT, LT, MT_Fall, MT_Spring, SMP, IsoTH
### BIO1 [,16]
GEO.BIO1.mat <- as.matrix(Env.GEO.scale[,16])
GEO.BIO1.mat.dist = ecodist::distance(GEO.BIO1.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO1.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO1.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO1.mat.dist

### BIO2 [,17]
GEO.BIO2.mat <- as.matrix(Env.GEO.scale[,17])
GEO.BIO2.mat.dist = ecodist::distance(GEO.BIO2.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO2.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO2.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO2.mat.dist

### BIO3 [,18]
GEO.BIO3.mat <- as.matrix(Env.GEO.scale[,18])
GEO.BIO3.mat.dist = ecodist::distance(GEO.BIO3.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO3.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO3.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO3.mat.dist

### BIO4 [,19]
GEO.BIO4.mat <- as.matrix(Env.GEO.scale[,19])
GEO.BIO4.mat.dist = ecodist::distance(GEO.BIO4.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO4.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO4.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO4.mat.dist

### BIO5 [,20]
GEO.BIO5.mat <- as.matrix(Env.GEO.scale[,20])
GEO.BIO5.mat.dist = ecodist::distance(GEO.BIO5.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO5.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO5.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO5.mat.dist

### BIO6 [,21]
GEO.BIO6.mat <- as.matrix(Env.GEO.scale[,21])
GEO.BIO6.mat.dist = ecodist::distance(GEO.BIO6.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO6.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO6.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO6.mat.dist

### BIO7 [,22]
GEO.BIO7.mat <- as.matrix(Env.GEO.scale[,22])
GEO.BIO7.mat.dist = ecodist::distance(GEO.BIO7.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO7.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO7.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO7.mat.dist

### BIO8 [,23]
GEO.BIO8.mat <- as.matrix(Env.GEO.scale[,23])
GEO.BIO8.mat.dist = ecodist::distance(GEO.BIO8.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO8.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO8.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO8.mat.dist

### BIO9 [,24]
GEO.BIO9.mat <- as.matrix(Env.GEO.scale[,24])
GEO.BIO9.mat.dist = ecodist::distance(GEO.BIO9.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO9.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO9.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO9.mat.dist

### BIO10 [,25]
GEO.BIO10.mat <- as.matrix(Env.GEO.scale[,25])
GEO.BIO10.mat.dist = ecodist::distance(GEO.BIO10.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO10.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO10.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO10.mat.dist

### BIO11 [,26]
GEO.BIO11.mat <- as.matrix(Env.GEO.scale[,26])
GEO.BIO11.mat.dist = ecodist::distance(GEO.BIO11.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO11.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO11.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO11.mat.dist

### BIO12 [,27]
GEO.BIO12.mat <- as.matrix(Env.GEO.scale[,27])
GEO.BIO12.mat.dist = ecodist::distance(GEO.BIO12.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO12.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO12.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO12.mat.dist

### BIO13 [,28]
GEO.BIO13.mat <- as.matrix(Env.GEO.scale[,28])
GEO.BIO13.mat.dist = ecodist::distance(GEO.BIO13.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO13.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO13.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO13.mat.dist

### BIO14 [,29]
GEO.BIO14.mat <- as.matrix(Env.GEO.scale[,29])
GEO.BIO14.mat.dist = ecodist::distance(GEO.BIO14.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO14.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO14.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO14.mat.dist

### BIO15 [,30]
GEO.BIO15.mat <- as.matrix(Env.GEO.scale[,30])
GEO.BIO15.mat.dist = ecodist::distance(GEO.BIO15.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO15.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO15.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO15.mat.dist

### BIO16 [,31]
GEO.BIO16.mat <- as.matrix(Env.GEO.scale[,31])
GEO.BIO16.mat.dist = ecodist::distance(GEO.BIO16.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO16.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO16.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO16.mat.dist

### BIO17 [,32]
GEO.BIO17.mat <- as.matrix(Env.GEO.scale[,32])
GEO.BIO17.mat.dist = ecodist::distance(GEO.BIO17.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO17.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO17.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO17.mat.dist

### BIO18 [,33]
GEO.BIO18.mat <- as.matrix(Env.GEO.scale[,33])
GEO.BIO18.mat.dist = ecodist::distance(GEO.BIO18.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO18.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO18.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO18.mat.dist

### BIO19 [,34]
GEO.BIO19.mat <- as.matrix(Env.GEO.scale[,34])
GEO.BIO19.mat.dist = ecodist::distance(GEO.BIO19.mat, method="gower") %>%
  as.matrix()
rownames(GEO.BIO19.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.BIO19.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.BIO19.mat.dist

### MT [,7]
GEO.MT.mat <- as.matrix(Env.GEO.scale[,7])
GEO.MT.mat.dist = ecodist::distance(GEO.MT.mat, method="gower") %>%
  as.matrix()
rownames(GEO.MT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.MT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.MT.mat.dist

### MHT [,8]
GEO.MHT.mat <- as.matrix(Env.GEO.scale[,8])
GEO.MHT.mat.dist = ecodist::distance(GEO.MHT.mat, method="gower") %>%
  as.matrix()
rownames(GEO.MHT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.MHT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.MHT.mat.dist

### MLT [,9]
GEO.MLT.mat <- as.matrix(Env.GEO.scale[,9])
GEO.MLT.mat.dist = ecodist::distance(GEO.MLT.mat, method="gower") %>%
  as.matrix()
rownames(GEO.MLT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.MLT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.MLT.mat.dist

### HT [,10]
GEO.HT.mat <- as.matrix(Env.GEO.scale[,10])
GEO.HT.mat.dist = ecodist::distance(GEO.HT.mat, method="gower") %>%
  as.matrix()
rownames(GEO.HT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.HT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.HT.mat.dist

### LT [,11]
GEO.LT.mat <- as.matrix(Env.GEO.scale[,11])
GEO.LT.mat.dist = ecodist::distance(GEO.LT.mat, method="gower") %>%
  as.matrix()
rownames(GEO.LT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.LT.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.LT.mat.dist

### MT_Fall [,12]
GEO.MTF.mat <- as.matrix(Env.GEO.scale[,12])
GEO.MTF.mat.dist = ecodist::distance(GEO.MTF.mat, method="gower") %>%
  as.matrix()
rownames(GEO.MTF.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.MTF.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.MTF.mat.dist

### MT_Spring [,13]
GEO.MTS.mat <- as.matrix(Env.GEO.scale[,13])
GEO.MTS.mat.dist = ecodist::distance(GEO.MTS.mat, method="gower") %>%
  as.matrix()
rownames(GEO.MTS.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.MTS.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.MTS.mat.dist

### SMP [,14]
GEO.SMP.mat <- as.matrix(Env.GEO.scale[,14])
GEO.SMP.mat.dist = ecodist::distance(GEO.SMP.mat, method="gower") %>%
  as.matrix()
rownames(GEO.SMP.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.SMP.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.SMP.mat.dist

### IsoTH [,15]
GEO.IsoTH.mat <- as.matrix(Env.GEO.scale[,15])
GEO.IsoTH.mat.dist = ecodist::distance(GEO.IsoTH.mat, method="gower") %>%
  as.matrix()
rownames(GEO.IsoTH.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(GEO.IsoTH.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.IsoTH.mat.dist

### F-2-3. build a matrix =========================================
str(Env.GEO.scale)
Env.GEO.BIO.mat <- as.matrix(Env.GEO.scale[,4:34])
Env.GEO.BIO.mat

### F-2-4. Correlations between BIO- variables ====================
cor_mat <- cor(Env.GEO.scale[,16:34], method = "pearson")
cor_mat_upper <- cor_mat
cor_mat_upper[lower.tri(cor_mat_upper, diag = TRUE)] <- NA
cor_mat_upper <- round(cor_mat_upper, 2)
print(cor_mat_upper)

#####=====
#### F-3. Mantel's test ===========================================
### F-3-1. For Geographic variables ===============================
## Elevation
Mantel_DNV_Ele = vegan::mantel(RAD_DNV.ind_FST, GEO.Ele.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_Ele

## Latitudes
Mantel_DNV_Lat = vegan::mantel(RAD_DNV.ind_FST, GEO.Lat.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_Lat

## Longitudes
Mantel_DNV_Lon = vegan::mantel(RAD_DNV.ind_FST, GEO.Lon.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_Lon

## Distances
Mantel_DNV_Dist = vegan::mantel(RAD_DNV.ind_FST, GEO.Dis.mat.dist_km, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_Dist

### F-3-2. Environmental variables =================================
## BIO1
Mantel_DNV_BIO1 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO1.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO1

## BIO2
Mantel_DNV_BIO2 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO2.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO2

## BIO3
Mantel_DNV_BIO3 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO3.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO3

## BIO4
Mantel_DNV_BIO4 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO4.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO4

## BIO5
Mantel_DNV_BIO5 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO5.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO5

## BIO6
Mantel_DNV_BIO6 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO6.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO6

## BIO7
Mantel_DNV_BIO7 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO7.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO7

## BIO8
Mantel_DNVa_BIO8 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO8.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO8

## BIO9
Mantel_DNV_BIO9 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO9.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO9

## BIO10
Mantel_DNV_BIO10 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO10.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO10

## BIO11
Mantel_DNV_BIO11 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO11.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO11

## BIO12
Mantel_DNV_BIO12 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO12.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO12

## BIO13
Mantel_DNV_BIO13 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO13.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO13

## BIO14
Mantel_DNV_BIO14 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO14.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO14

## BIO15
Mantel_DNV_BIO15 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO15.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO15

## BIO16
Mantel_DNV_BIO16 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO16.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO16

## BIO17
Mantel_DNV_BIO17 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO17.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO17

## BIO18
Mantel_DNV_BIO18 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO18.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO18

## BIO19
Mantel_DNV_BIO19 = vegan::mantel(RAD_DNV.ind_FST, GEO.BIO19.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV_BIO19
#####=====
#### F-4. Visualization of IBD and IBE ============================
### F-4-1. Tabling ================================================
## F-4-1-a. Geographic variables ==================================
## Elevations
lower.tri(GEO.Ele.mat.dist)
GEO.Ele.mat.dist[lower.tri(GEO.Ele.mat.dist)] <- NA
GEO.Ele.mat.dist.melt <- melt(GEO.Ele.mat.dist, varnames = c("row", "col"))
colnames(GEO.Ele.mat.dist.melt) <- c("Pop1","Pop2","Elevation")
GEO.Ele.mat.dist.melt <- GEO.Ele.mat.dist.melt %>%
  na.omit(Elevation)
GEO.Ele.mat.dist.melt

## Latitudes
lower.tri(GEO.Lat.mat.dist)
GEO.Lat.mat.dist[lower.tri(GEO.Lat.mat.dist)] <- NA
GEO.Lat.mat.dist.melt <- melt(GEO.Lat.mat.dist, varnames = c("row", "col"))
colnames(GEO.Lat.mat.dist.melt) <- c("Pop1","Pop2","Latitude")
GEO.Lat.mat.dist.melt <- GEO.Lat.mat.dist.melt %>%
  na.omit(Latitude)
GEO.Lat.mat.dist.melt

## Longitudes
lower.tri(GEO.Lon.mat.dist)
GEO.Lon.mat.dist[lower.tri(GEO.Lon.mat.dist)] <- NA
GEO.Lon.mat.dist.melt <- melt(GEO.Lon.mat.dist, varnames = c("row", "col"))
colnames(GEO.Lon.mat.dist.melt) <- c("Pop1","Pop2","Longitude")
GEO.Lon.mat.dist.melt <- GEO.Lon.mat.dist.melt %>%
  na.omit(Longitude)
GEO.Lon.mat.dist.melt

## Distances
lower.tri(GEO.Dis.mat.dist_km)
GEO.Dis.mat.dist_km[lower.tri(GEO.Dis.mat.dist_km)] <- NA
GEO.Dis.mat.dist_km.melt <- melt(GEO.Dis.mat.dist_km, varnames = c("row", "col"))
colnames(GEO.Dis.mat.dist_km.melt) <- c("Pop1","Pop2","Distance")
GEO.Dis.mat.dist_km.melt <- GEO.Dis.mat.dist_km.melt %>%
  na.omit(Distance)
GEO.Dis.mat.dist_km.melt


## F-4-1-b. Environmental variables =================================
## BIO1
lower.tri(GEO.BIO1.mat.dist)
GEO.BIO1.mat.dist[lower.tri(GEO.BIO1.mat.dist)] <- NA
GEO.BIO1.mat.dist.melt <- melt(GEO.BIO1.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO1.mat.dist.melt) <- c("Pop1","Pop2","BIO1")
GEO.BIO1.mat.dist.melt <- GEO.BIO1.mat.dist.melt %>%
  na.omit(BIO1)
GEO.BIO1.mat.dist.melt

## BIO2
lower.tri(GEO.BIO2.mat.dist)
GEO.BIO2.mat.dist[lower.tri(GEO.BIO2.mat.dist)] <- NA
GEO.BIO2.mat.dist.melt <- melt(GEO.BIO2.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO2.mat.dist.melt) <- c("Pop1","Pop2","BIO2")
GEO.BIO2.mat.dist.melt <- GEO.BIO2.mat.dist.melt %>%
  na.omit(BIO2)
GEO.BIO2.mat.dist.melt

## BIO3
lower.tri(GEO.BIO3.mat.dist)
GEO.BIO3.mat.dist[lower.tri(GEO.BIO3.mat.dist)] <- NA
GEO.BIO3.mat.dist.melt <- melt(GEO.BIO3.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO3.mat.dist.melt) <- c("Pop1","Pop2","BIO3")
GEO.BIO3.mat.dist.melt <- GEO.BIO3.mat.dist.melt %>%
  na.omit(BIO3)
GEO.BIO3.mat.dist.melt

## BIO4
lower.tri(GEO.BIO4.mat.dist)
GEO.BIO4.mat.dist[lower.tri(GEO.BIO4.mat.dist)] <- NA
GEO.BIO4.mat.dist.melt <- melt(GEO.BIO4.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO4.mat.dist.melt) <- c("Pop1","Pop2","BIO4")
GEO.BIO4.mat.dist.melt <- GEO.BIO4.mat.dist.melt %>%
  na.omit(BIO4)
GEO.BIO4.mat.dist.melt

## BIO5
lower.tri(GEO.BIO5.mat.dist)
GEO.BIO5.mat.dist[lower.tri(GEO.BIO5.mat.dist)] <- NA
GEO.BIO5.mat.dist.melt <- melt(GEO.BIO5.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO5.mat.dist.melt) <- c("Pop1","Pop2","BIO5")
GEO.BIO5.mat.dist.melt <- GEO.BIO5.mat.dist.melt %>%
  na.omit(BIO5)
GEO.BIO5.mat.dist.melt

## BIO6
lower.tri(GEO.BIO6.mat.dist)
GEO.BIO6.mat.dist[lower.tri(GEO.BIO6.mat.dist)] <- NA
GEO.BIO6.mat.dist.melt <- melt(GEO.BIO6.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO6.mat.dist.melt) <- c("Pop1","Pop2","BIO6")
GEO.BIO6.mat.dist.melt <- GEO.BIO6.mat.dist.melt %>%
  na.omit(BIO6)
GEO.BIO6.mat.dist.melt

## BIO7
lower.tri(GEO.BIO7.mat.dist)
GEO.BIO7.mat.dist[lower.tri(GEO.BIO7.mat.dist)] <- NA
GEO.BIO7.mat.dist.melt <- melt(GEO.BIO7.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO7.mat.dist.melt) <- c("Pop1","Pop2","BIO7")
GEO.BIO7.mat.dist.melt <- GEO.BIO7.mat.dist.melt %>%
  na.omit(BIO7)
GEO.BIO7.mat.dist.melt

## BIO8
lower.tri(GEO.BIO8.mat.dist)
GEO.BIO8.mat.dist[lower.tri(GEO.BIO8.mat.dist)] <- NA
GEO.BIO8.mat.dist.melt <- melt(GEO.BIO8.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO8.mat.dist.melt) <- c("Pop1","Pop2","BIO8")
GEO.BIO8.mat.dist.melt <- GEO.BIO8.mat.dist.melt %>%
  na.omit(BIO8)
GEO.BIO8.mat.dist.melt

## BIO9
lower.tri(GEO.BIO9.mat.dist)
GEO.BIO9.mat.dist[lower.tri(GEO.BIO9.mat.dist)] <- NA
GEO.BIO9.mat.dist.melt <- melt(GEO.BIO9.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO9.mat.dist.melt) <- c("Pop1","Pop2","BIO9")
GEO.BIO9.mat.dist.melt <- GEO.BIO9.mat.dist.melt %>%
  na.omit(BIO9)
GEO.BIO9.mat.dist.melt

## BIO10
lower.tri(GEO.BIO10.mat.dist)
GEO.BIO10.mat.dist[lower.tri(GEO.BIO10.mat.dist)] <- NA
GEO.BIO10.mat.dist.melt <- melt(GEO.BIO10.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO10.mat.dist.melt) <- c("Pop1","Pop2","BIO10")
GEO.BIO10.mat.dist.melt <- GEO.BIO10.mat.dist.melt %>%
  na.omit(BIO10)
GEO.BIO10.mat.dist.melt

## BIO11
lower.tri(GEO.BIO11.mat.dist)
GEO.BIO11.mat.dist[lower.tri(GEO.BIO11.mat.dist)] <- NA
GEO.BIO11.mat.dist.melt <- melt(GEO.BIO11.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO11.mat.dist.melt) <- c("Pop1","Pop2","BIO11")
GEO.BIO11.mat.dist.melt <- GEO.BIO11.mat.dist.melt %>%
  na.omit(BIO11)
GEO.BIO11.mat.dist.melt

## BIO12
lower.tri(GEO.BIO12.mat.dist)
GEO.BIO12.mat.dist[lower.tri(GEO.BIO12.mat.dist)] <- NA
GEO.BIO12.mat.dist.melt <- melt(GEO.BIO12.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO12.mat.dist.melt) <- c("Pop1","Pop2","BIO12")
GEO.BIO12.mat.dist.melt <- GEO.BIO12.mat.dist.melt %>%
  na.omit(BIO12)
GEO.BIO12.mat.dist.melt

## BIO13
lower.tri(GEO.BIO13.mat.dist)
GEO.BIO13.mat.dist[lower.tri(GEO.BIO13.mat.dist)] <- NA
GEO.BIO13.mat.dist.melt <- melt(GEO.BIO13.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO13.mat.dist.melt) <- c("Pop1","Pop2","BIO13")
GEO.BIO13.mat.dist.melt <- GEO.BIO13.mat.dist.melt %>%
  na.omit(BIO13)
GEO.BIO13.mat.dist.melt

## BIO14
lower.tri(GEO.BIO14.mat.dist)
GEO.BIO14.mat.dist[lower.tri(GEO.BIO14.mat.dist)] <- NA
GEO.BIO14.mat.dist.melt <- melt(GEO.BIO14.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO14.mat.dist.melt) <- c("Pop1","Pop2","BIO14")
GEO.BIO14.mat.dist.melt <- GEO.BIO14.mat.dist.melt %>%
  na.omit(BIO14)
GEO.BIO14.mat.dist.melt

## BIO15
lower.tri(GEO.BIO15.mat.dist)
GEO.BIO15.mat.dist[lower.tri(GEO.BIO15.mat.dist)] <- NA
GEO.BIO15.mat.dist.melt <- melt(GEO.BIO15.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO15.mat.dist.melt) <- c("Pop1","Pop2","BIO15")
GEO.BIO15.mat.dist.melt <- GEO.BIO15.mat.dist.melt %>%
  na.omit(BIO15)
GEO.BIO15.mat.dist.melt

## BIO16
lower.tri(GEO.BIO16.mat.dist)
GEO.BIO16.mat.dist[lower.tri(GEO.BIO16.mat.dist)] <- NA
GEO.BIO16.mat.dist.melt <- melt(GEO.BIO16.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO16.mat.dist.melt) <- c("Pop1","Pop2","BIO16")
GEO.BIO16.mat.dist.melt <- GEO.BIO16.mat.dist.melt %>%
  na.omit(BIO16)
GEO.BIO16.mat.dist.melt

## BIO17
lower.tri(GEO.BIO17.mat.dist)
GEO.BIO17.mat.dist[lower.tri(GEO.BIO17.mat.dist)] <- NA
GEO.BIO17.mat.dist.melt <- melt(GEO.BIO17.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO17.mat.dist.melt) <- c("Pop1","Pop2","BIO17")
GEO.BIO17.mat.dist.melt <- GEO.BIO17.mat.dist.melt %>%
  na.omit(BIO17)
GEO.BIO17.mat.dist.melt

## BIO18
lower.tri(GEO.BIO18.mat.dist)
GEO.BIO18.mat.dist[lower.tri(GEO.BIO18.mat.dist)] <- NA
GEO.BIO18.mat.dist.melt <- melt(GEO.BIO18.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO18.mat.dist.melt) <- c("Pop1","Pop2","BIO18")
GEO.BIO18.mat.dist.melt <- GEO.BIO18.mat.dist.melt %>%
  na.omit(BIO18)
GEO.BIO18.mat.dist.melt

## BIO19
lower.tri(GEO.BIO19.mat.dist)
GEO.BIO19.mat.dist[lower.tri(GEO.BIO19.mat.dist)] <- NA
GEO.BIO19.mat.dist.melt <- melt(GEO.BIO19.mat.dist, varnames = c("row", "col"))
colnames(GEO.BIO19.mat.dist.melt) <- c("Pop1","Pop2","BIO19")
GEO.BIO19.mat.dist.melt <- GEO.BIO19.mat.dist.melt %>%
  na.omit(BIO19)
GEO.BIO19.mat.dist.melt


## F-4-1-c. Genetic distances ====================================
# FST
RAD_DNV.ind_FST
lower.tri(RAD_DNV.ind_FST)
RAD_DNV.ind_FST[lower.tri(RAD_DNV.ind_FST)] <- NA
RAD_DNV.ind_FST.melt <- melt(RAD_DNV.ind_FST, varnames = c("row", "col"))
colnames(RAD_DNV.ind_FST.melt) <- c("Pop1","Pop2","FST")
RAD_DNV.ind_FST.melt <- RAD_DNV.ind_FST.melt %>%
  na.omit(FST)
RAD_DNV.ind_FST.melt

# Linearized FST
RAD_DNV.ind_LinFST
lower.tri(RAD_DNV.ind_LinFST)
RAD_DNV.ind_LinFST[lower.tri(RAD_DNV.ind_LinFST)] <- NA
RAD_DNV.ind_LinFST.melt <- melt(RAD_DNV.ind_LinFST, varnames = c("row", "col"))
colnames(RAD_DNV.ind_LinFST.melt) <- c("Pop1","Pop2","LinFST")
RAD_DNV.ind_LinFST.melt <- RAD_DNV.ind_LinFST.melt %>%
  na.omit(LinFST)
RAD_DNV.ind_LinFST.melt

## F-4-1-d. Combining melt-data ===================================
a <- cbind(GEO.Ele.mat.dist.melt,
           GEO.Lat.mat.dist.melt[,3],
           GEO.Lon.mat.dist.melt[,3],
           GEO.Dis.mat.dist_km.melt[,3],
           GEO.BIO1.mat.dist.melt[,3],
           GEO.BIO2.mat.dist.melt[,3],
           GEO.BIO3.mat.dist.melt[,3],
           GEO.BIO4.mat.dist.melt[,3],
           GEO.BIO5.mat.dist.melt[,3],
           GEO.BIO6.mat.dist.melt[,3],
           GEO.BIO7.mat.dist.melt[,3],
           GEO.BIO8.mat.dist.melt[,3],
           GEO.BIO9.mat.dist.melt[,3],
           GEO.BIO10.mat.dist.melt[,3],
           GEO.BIO11.mat.dist.melt[,3],
           GEO.BIO12.mat.dist.melt[,3],
           GEO.BIO13.mat.dist.melt[,3],
           GEO.BIO14.mat.dist.melt[,3],
           GEO.BIO15.mat.dist.melt[,3],
           GEO.BIO16.mat.dist.melt[,3],
           GEO.BIO17.mat.dist.melt[,3],
           GEO.BIO18.mat.dist.melt[,3],
           GEO.BIO19.mat.dist.melt[,3],
           RAD_DNV4.ind_FST.melt[,3],
           RAD_DNV4.ind_LinFST.melt[,3])

colnames(a) <- c("Pop1","Pop2",
                 "Elevation",
                 "Latitude",
                 "Longitude",
                 "Distance",
                 "BIO1",
                 "BIO2",
                 "BIO3",
                 "BIO4",
                 "BIO5",
                 "BIO6",
                 "BIO7",
                 "BIO8",
                 "BIO9",
                 "BIO10",
                 "BIO11",
                 "BIO12",
                 "BIO13",
                 "BIO14",
                 "BIO15",
                 "BIO16",
                 "BIO17",
                 "BIO18",
                 "BIO19",
                 "FST",
                 "LinFST")

a$Pairwise <- paste(a$Pop1, a$Pop2, sep = "-")

a <- a %>%
  filter(Pop1 != Pop2)

### THIS a is the following input dataset named "Pairwise_IBD+IBE.csv"
#####=====
### F-4-2. Plotting ================================================
## F-4-2-a. Importing data =========================================
IBDIBE <- read.csv("Pairwise_IBD+IBE.csv")

### F-4-3. Drawing figures =========================================
## F-4-3-a. IBD_Distance ============================================
IBD_Dist_plot <- ggplot(IBDIBE, aes(Distance, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=4, alpha=.8) +
  geom_text_repel(size=3, segment.size=0.1)+
  labs(
    x = "Distances [km]", 
    y = "Linearized FST", 
    fill = "Populations",
    title = "a)",
    subtitle = bquote("Mantel r = 0.487 (" * italic(P) * " = 0.013)")) +
  theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.title.x = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.text.x = element_text(hjust=0.5, vjust=1.0, size=10),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=10, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=10, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=9, face="bold.italic", colour = "black", vjust = 0))

IBD_Dist_plot

## F-4-3-b. IBD_Latitude ============================================
IBD_Lat_plot <- ggplot(IBDIBE, aes(Latitude, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=4, alpha=.8) +
  geom_text_repel(size=3, segment.size=0.1)+
  labs(
    x = "Latitude [°N]", 
    y = "Linearized FST", 
    fill = "Populations",
    title = "b)",
    subtitle = bquote("Mantel r = 0.524 (" * italic(P) * " = 0.005)")) +
  theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.title.x = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.text.x = element_text(hjust=0.5, vjust=1.0, size=10),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=10, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=10, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=9, face="bold.italic", colour = "black", vjust = 0))

IBD_Lat_plot

## F-4-3-c. IBE_BIO3 =================================================
IBE_BIO3_plot <- ggplot(IBDIBE, aes(BIO3, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=4, alpha=.8) +
  geom_text_repel(size=3, segment.size=0.1)+
  labs(
    x = "BIO3 (Isothermality)", 
    y = "Linearized FST", 
    fill = "Populations",
    title = "c)",
    subtitle = bquote("Mantel r = 0.471 (" * italic(P) * " = 0.018)")) +
  theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.title.x = element_text(hjust=0.5, vjust=0.5, size=10),
    axis.text.x = element_text(hjust=0.5, vjust=1.0, size=10),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=10, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=10, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=9, face="bold.italic", colour = "black", vjust = 0)) 

IBE_BIO3_plot

### F-4-4. Export Figure 3 ==========================================
pdf("Figure3.pdf", width = 4, height = 12)
grid.arrange(IBD_Dist_plot, IBD_Lat_plot, IBE_BIO3_plot, nrow = 3, ncol = 1)
dev.off()

dev.off()
#####=====

##### G. MMRR with CA ==============================================
#### G-1. Building models ==========================================
### G-1-1. Scaling distance and latitudes ==========================
z.GEO.Dis.mat.dist_km <- matrix(scale(as.vector(GEO.Dis.mat.dist_km)), nrow = nrow(GEO.Dis.mat.dist_km))
rownames(z.GEO.Dis.mat.dist_km) <- rownames(GEO.Dis.mat.dist_km)
colnames(z.GEO.Dis.mat.dist_km) <- colnames(GEO.Dis.mat.dist_km)
z.GEO.Dis.mat.dist_km

z.GEO.Lat.mat.dist <- matrix(scale(as.vector(GEO.Lat.mat.dist)), nrow = nrow(GEO.Lat.mat.dist))
rownames(z.GEO.Lat.mat.dist) <- rownames(GEO.Lat.mat.dist)
colnames(z.GEO.Lat.mat.dist) <- colnames(GEO.Lat.mat.dist)
z.GEO.Lat.mat.dist

Xmats.BIO.sel <- list(BIO1 = GEO.BIO1.mat.dist,
                      BIO3 = GEO.BIO3.mat.dist,
                      BIO7 = GEO.BIO7.mat.dist,
                      BIO12 = GEO.BIO12.mat.dist,
                      BIO14 = GEO.BIO14.mat.dist,
                      BIO18 = GEO.BIO18.mat.dist)

Xmats.geo.BIO3 <- list(distance=z.GEO.Dis.mat.dist_km,
                       latitude=z.GEO.Lat.mat.dist,
                       environment=GEO.BIO3.mat.dist)

Xmats.GEO ## GEO + ELE + LAT + LON
Xmats.BIO ## BIO + BIO1 to 19 
Xmats.BIO.sel ## BIO1, BIO3, BIO7, BIO12, BIO14, BIO18
Xmats.geo.BIO3 ## GEO + BIO3

### G-1-2. Genetic distances =======================================
RAD_DNV.ind_FST

### G-2. MMRR ======================================================
PopGenReport::lgrMMRR(RAD_DNV.ind_FST, Xmats.geo.BIO3, nperm=9999)

### G-3. CA ========================================================
# Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
## G-3-1. GEO + BIO3 ===============================================
DNV_Dist_BIO3 <- yhat::commonalityCoefficients(IBDIBE, "FST", list(c("Distance","Latitude","BIO3")), imat=FALSE)
DNV_Dist_BIO3
#####=====

##### H. STRUCTURE with LEA package ================================
#### H-1. converting genind into vcf ===============================
RAD_DNV.ind
RAD_DNV.ind@pop <- NULL
loci_obj <- as.loci(RAD_DNV.ind)
write.vcf(loci_obj, file = "RAD_DNV.vcf")

#### H-2. Calling vcf2geno =========================================
LEA::vcf2geno("RAD_DNV.vcf",
              "RAD_DNV.geno")

#### H-3. Running SNMF =============================================
project_DNV = snmf("RAD_DNV.geno",
                   K = 1:20,
                   entropy = TRUE,
                   repetitions = 200,
                   project = "new", CPU = 20,
                   ploidy = 2)

jpeg("SNMF_RAD_DNV.jpeg", width = 8, height = 8, units = "cm", res = 480)
par(mar = c(4, 4, 1, 1))
par(cex.axis = 0.6, cex.lab = 0.6)
plot(project_DNV, col="grey30", pch = 19, cex = 0.6)
dev.off()

### H-3-1. K selection (K=2) ========================================
ce2 <-  cross.entropy(project_DNV, K = 2)
ce2

best_run <- which.min(ce2)
best_run


q_mat2 <- LEA::Q(project_DNV, K = 2, run = best_run) 
colnames(q_mat2) <- paste0("P", 1:2)
head(q_mat2)

q_df2 <- q_mat2 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Individual = POPIND_181$Individuals,
         Region = POPIND_181$Populations)

q_df2

write.table(q_df2, "q_df_RAD_DNV_K2.txt")

### H-3-2. K selection (K=3) ========================================
ce3 <-  cross.entropy(project_DNV, K = 3)
ce3

best_run <- which.min(ce3)
best_run


q_mat3 <- LEA::Q(project_DNV, K = 3, run = best_run) 
colnames(q_mat3) <- paste0("P", 1:3)
head(q_mat3)

q_df3 <- q_mat3 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Individual = POPIND_181$Individuals,
         Region = POPIND_181$Populations)

q_df3

write.table(q_df3, "q_df_RAD_DNV_K3.txt")

### H-3-3. K selection (K=4) ========================================
ce4 <-  cross.entropy(project_DNV, K = 4)
ce4

best_run <- which.min(ce4)
best_run


q_mat4 <- LEA::Q(project_DNV, K = 4, run = best_run) 
colnames(q_mat4) <- paste0("P", 1:4)
head(q_mat4)

q_df4 <- q_mat4 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Individual = POPIND_181$Individuals,
         Region = POPIND_181$Populations)

q_df4

write.table(q_df4, "q_df_RAD_DNV_K4.txt")

### H-3-4. K selection (K=5) ========================================
ce5 <-  cross.entropy(project_DNV, K = 5)
ce5

best_run <- which.min(ce5)
best_run


q_mat5 <- LEA::Q(project_DNV, K = 5, run = best_run) 
colnames(q_mat5) <- paste0("P", 1:5)
head(q_mat5)

q_df5 <- q_mat5 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Individual = POPIND_181$Individuals,
         Region = POPIND_181$Populations)

q_df5

write.table(q_df5, "q_df_RAD_DNV_K5.txt")


#### H-4. Plotting population structure =============================
### H-4-1. Plot (K=2) ================================================
q_df2 <- read.table("q_df_RAD_DNV_K2.txt")

q_df_long <- q_df2 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "Population", values_to = "q") 

q_df_long

region_order <- c("BOP", "DGR", "MZS", "TAB", "JCH", "JSN", "NWS", "DMY", "CHJ", "ICH")

q_df_prates <- q_df_long %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region, Individual) %>%
  mutate(Individual = forcats::fct_inorder(Individual))

q_df_prates

q_palette <- c("red","darkblue")

# region_order <- c("BOP", "CHJ", "DGR", "DMY", "ICH", "JCH", "JSN", "MZS", "NWS", "TAB")

q_df_ordered <- q_df_long %>% 
  group_by(Individual) %>%
  mutate(likely_assignment = Population[which.max(q)],
         assignment_prob = max(q)) %>%
  ungroup() %>%
  mutate(individual_num = as.numeric(sub("Cbp_", "", Individual)),
         Region = factor(Region, levels = region_order)) %>%
  arrange(Region, individual_num) %>% 
  dplyr::select(-individual_num) %>% 
  mutate(Individual = forcats::fct_inorder(factor(Individual)))

# View the ordered tibble
print(q_df_ordered)

# Remove the prefix
q_df_ordered2 <- q_df_ordered %>%
  mutate(individual_short = sub("Cbp_0*", "", Individual)) %>%
  mutate(individual_short = fct_reorder(individual_short, as.numeric(sub("Cbp_", "", Individual))))  

print(q_df_ordered2, n=200)

###
region_positions <- q_df_ordered2 %>%
  group_by(Region) %>%
  dplyr::summarize(start = min(as.numeric(sub("Cbp_", "", Individual))),
                   end = max(as.numeric(sub("Cbp_", "", Individual))),
                   max_q = max(q)) %>%
  mutate(start = start - 0.5, end = end + 0.5,
         mid = (start + end) / 2,
         label_y = 1.01) 

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "MZS", 141, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "MZS", 151.5, end))

region_positions <- region_positions %>%
  mutate(start = ifelse(Region == "NWS", 151.5, start))

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "TAB", 176.5, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "TAB", 181.5, end))

region_positions

# Create the plot
RAD_DNV_structure <- q_df_ordered2 %>% 
  ggplot(aes(x = individual_short, y = q, fill = Population)) +
  geom_col(width = 1.00) +
  geom_rect(data = region_positions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
            fill = NA, color = "grey90", linewidth = 0.1, inherit.aes = FALSE) +  # Use linewidth instead of size
  geom_text(data = region_positions, aes(x = mid, y = label_y, label = Region), 
            color = "black", size = 2.5, inherit.aes = FALSE, vjust = 0) +  # Add region labels
  scale_fill_manual(values = q_palette, labels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")) +
  labs(fill = "Region", x = "Individuals", y = "Ancestral proportions",
       title = " ", subtitle = "K = 2") + 
  theme_article() +
  coord_cartesian(ylim=c(0,1.06), expand = FALSE ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(face = "bold", size = 7, vjust = -3),
    plot.subtitle = element_text(size = 6, vjust = -2),
    axis.line = element_blank(),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 5),
    strip.background = element_rect(fill = "transparent", color = "white"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.key.height = unit(0.3, "lines"), 
    legend.key.width = unit(0.3, "lines"),  
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",              
    legend.margin = margin(0, 0, 0, -6)   
  )

RAD_DNV_structure

### H-4-2. Plot (K=3) ================================================
q_df3 <- read.table("q_df_RAD_DNV_K3.txt")

q_df_long <- q_df3 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "Population", values_to = "q") 

q_df_long

region_order <- c("BOP", "DGR", "MZS", "TAB", "JCH", "JSN", "NWS", "DMY", "CHJ", "ICH")

q_df_prates <- q_df_long %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region, Individual) %>%
  mutate(Individual = forcats::fct_inorder(Individual))

q_df_prates

q_palette <- c("red","darkblue","orange")

# region_order <- c("BOP", "CHJ", "DGR", "DMY", "ICH", "JCH", "JSN", "MZS", "NWS", "TAB")

q_df_ordered <- q_df_long %>% 
  group_by(Individual) %>%
  mutate(likely_assignment = Population[which.max(q)],
         assignment_prob = max(q)) %>%
  ungroup() %>%
  mutate(individual_num = as.numeric(sub("Cbp_", "", Individual)),
         Region = factor(Region, levels = region_order)) %>%
  arrange(Region, individual_num) %>% 
  dplyr::select(-individual_num) %>%  
  mutate(Individual = forcats::fct_inorder(factor(Individual)))

# View the ordered tibble
print(q_df_ordered)

# Remove the prefix
q_df_ordered2 <- q_df_ordered %>%
  mutate(individual_short = sub("Cbp_0*", "", Individual)) %>%
  mutate(individual_short = fct_reorder(individual_short, as.numeric(sub("Cbp_", "", Individual))))  # Ensure correct ordering

print(q_df_ordered2, n=200)

###
region_positions <- q_df_ordered2 %>%
  group_by(Region) %>%
  dplyr::summarize(start = min(as.numeric(sub("Cbp_", "", Individual))),
                   end = max(as.numeric(sub("Cbp_", "", Individual))),
                   max_q = max(q)) %>%
  mutate(start = start - 0.5, end = end + 0.5,
         mid = (start + end) / 2,
         label_y = 1.01) 

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "MZS", 141, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "MZS", 151.5, end))

region_positions <- region_positions %>%
  mutate(start = ifelse(Region == "NWS", 151.5, start))

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "TAB", 176.5, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "TAB", 181.5, end))

region_positions


# Create the plot
RAD_DNV_structure <- q_df_ordered2 %>% 
  ggplot(aes(x = individual_short, y = q, fill = Population)) +
  geom_col(width = 1.00) +
  geom_rect(data = region_positions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
            fill = NA, color = "grey90", linewidth = 0.1, inherit.aes = FALSE) +  # Use linewidth instead of size
  geom_text(data = region_positions, aes(x = mid, y = label_y, label = Region), 
            color = "black", size = 2.5, inherit.aes = FALSE, vjust = 0) +  # Add region labels
  scale_fill_manual(values = q_palette, labels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")) +
  labs(fill = "Region", x = "Individuals", y = "Ancestral proportions",
       title = " ", subtitle = "K = 3") + 
  theme_article() +
  coord_cartesian(ylim=c(0,1.06), expand = FALSE ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(face = "bold", size = 7, vjust = -3),
    plot.subtitle = element_text(size = 6, vjust = -2),
    axis.line = element_blank(),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 5),
    strip.background = element_rect(fill = "transparent", color = "white"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"), 
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.key.height = unit(0.3, "lines"),  
    legend.key.width = unit(0.3, "lines"),
    legend.text = element_text(size = 4, margin = margin(b = 0.9)),
    legend.title = element_text(size = 4, margin = margin(b = 1.5)),
    legend.position = "none",           
    legend.margin = margin(0, 0, 0, -6)       
  )

RAD_DNV_structure

### H-4-3. Plot (K=4) ================================================
q_df4 <- read.table("q_df_RAD_DNV_K4.txt")

q_df_long <- q_df4 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "Population", values_to = "q") 

q_df_long

region_order <- c("BOP", "DGR", "MZS", "TAB", "JCH", "JSN", "NWS", "DMY", "CHJ", "ICH")

q_df_prates <- q_df_long %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region, Individual) %>%
  mutate(Individual = forcats::fct_inorder(Individual))

q_df_prates

q_palette <- c("darkblue","red","orange","darkgreen")

# region_order <- c("BOP", "CHJ", "DGR", "DMY", "ICH", "JCH", "JSN", "MZS", "NWS", "TAB")

q_df_ordered <- q_df_long %>% 
  group_by(Individual) %>%
  mutate(likely_assignment = Population[which.max(q)],
         assignment_prob = max(q)) %>%
  ungroup() %>%
  mutate(individual_num = as.numeric(sub("Cbp_", "", Individual)),
         Region = factor(Region, levels = region_order)) %>%
  arrange(Region, individual_num) %>% 
  dplyr::select(-individual_num) %>%  # Optionally remove the helper column
  mutate(Individual = forcats::fct_inorder(factor(Individual)))

# View the ordered tibble
print(q_df_ordered)

# Remove the prefix
q_df_ordered2 <- q_df_ordered %>%
  mutate(individual_short = sub("Cbp_0*", "", Individual)) %>%
  mutate(individual_short = fct_reorder(individual_short, as.numeric(sub("Cbp_", "", Individual))))  # Ensure correct ordering

print(q_df_ordered2, n=200)

###
region_positions <- q_df_ordered2 %>%
  group_by(Region) %>%
  dplyr::summarize(start = min(as.numeric(sub("Cbp_", "", Individual))),
                   end = max(as.numeric(sub("Cbp_", "", Individual))),
                   max_q = max(q)) %>%
  mutate(start = start - 0.5, end = end + 0.5,
         mid = (start + end) / 2,
         label_y = 1.01) 

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "MZS", 141, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "MZS", 151.5, end))

region_positions <- region_positions %>%
  mutate(start = ifelse(Region == "NWS", 151.5, start))

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "TAB", 176.5, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "TAB", 181.5, end))

region_positions


# Create the plot
RAD_DNV_structure <- q_df_ordered2 %>% 
  ggplot(aes(x = individual_short, y = q, fill = Population)) +
  geom_col(width = 1.00) +
  geom_rect(data = region_positions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
            fill = NA, color = "grey90", linewidth = 0.1, inherit.aes = FALSE) +  # Use linewidth instead of size
  geom_text(data = region_positions, aes(x = mid, y = label_y, label = Region), 
            color = "black", size = 2.5, inherit.aes = FALSE, vjust = 0) +  # Add region labels
  scale_fill_manual(values = q_palette, labels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")) +
  labs(fill = "Region", x = "Individuals", y = "Ancestral proportions",
       title = " ", subtitle = "K = 4") + 
  theme_article() +
  coord_cartesian(ylim=c(0,1.06), expand = FALSE ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(face = "bold", size = 7, vjust = -3),
    plot.subtitle = element_text(size = 6, vjust = -2),
    axis.line = element_blank(),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 5),
    strip.background = element_rect(fill = "transparent", color = "white"),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),  
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.key.height = unit(0.3, "lines"),  
    legend.key.width = unit(0.3, "lines"),  
    legend.text = element_text(size = 4, margin = margin(b = 0.9)),
    legend.title = element_text(size = 4, margin = margin(b = 1.5)),
    legend.position = "none",                
    legend.margin = margin(0, 0, 0, -6)    
  )

RAD_DNV_structure

### H-4-4. Plot (K=5) ================================================
q_df5 <- read.table("q_df_RAD_DNV_K5.txt")

q_df_long <- q_df5 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "Population", values_to = "q") 

q_df_long

region_order <- c("BOP", "DGR", "MZS", "TAB", "JCH", "JSN", "NWS", "DMY", "CHJ", "ICH")

q_df_prates <- q_df_long %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region, Individual) %>%
  mutate(Individual = forcats::fct_inorder(Individual))

q_df_prates

q_palette <- c("red", "darkgreen","pink","darkblue","orange")

# region_order <- c("BOP", "CHJ", "DGR", "DMY", "ICH", "JCH", "JSN", "MZS", "NWS", "TAB")

q_df_ordered <- q_df_long %>% 
  group_by(Individual) %>%
  mutate(likely_assignment = Population[which.max(q)],
         assignment_prob = max(q)) %>%
  ungroup() %>%
  mutate(individual_num = as.numeric(sub("Cbp_", "", Individual)),
         Region = factor(Region, levels = region_order)) %>%
  arrange(Region, individual_num) %>% 
  dplyr::select(-individual_num) %>%  # Optionally remove the helper column
  mutate(Individual = forcats::fct_inorder(factor(Individual)))

# View the ordered tibble
print(q_df_ordered)

# Remove the prefix
q_df_ordered2 <- q_df_ordered %>%
  mutate(individual_short = sub("Cbp_0*", "", Individual)) %>%
  mutate(individual_short = fct_reorder(individual_short, as.numeric(sub("Cbp_", "", Individual))))  # Ensure correct ordering

print(q_df_ordered2, n=200)

###
region_positions <- q_df_ordered2 %>%
  group_by(Region) %>%
  dplyr::summarize(start = min(as.numeric(sub("Cbp_", "", Individual))),
                   end = max(as.numeric(sub("Cbp_", "", Individual))),
                   max_q = max(q)) %>%
  mutate(start = start - 0.5, end = end + 0.5,
         mid = (start + end) / 2,
         label_y = 1.01) 

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "MZS", 141, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "MZS", 151.5, end))

region_positions <- region_positions %>%
  mutate(start = ifelse(Region == "NWS", 151.5, start))

region_positions <- region_positions %>%
  mutate(mid = ifelse(Region == "TAB", 176.5, mid))

region_positions <- region_positions %>%
  mutate(end = ifelse(Region == "TAB", 181.5, end))

region_positions


# Create the plot
RAD_DNV_structure <- q_df_ordered2 %>% 
  ggplot(aes(x = individual_short, y = q, fill = Population)) +
  geom_col(width = 1.00) +
  geom_rect(data = region_positions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
            fill = NA, color = "grey90", linewidth = 0.1, inherit.aes = FALSE) +  # Use linewidth instead of size
  geom_text(data = region_positions, aes(x = mid, y = label_y, label = Region), 
            color = "black", size = 2.5, inherit.aes = FALSE, vjust = 0) +  # Add region labels
  scale_fill_manual(values = q_palette, labels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")) +
  labs(fill = "Region", x = "Individuals", y = "Ancestral proportions",
       title = " ", subtitle = "K = 5") +
  theme_article() +
  coord_cartesian(ylim=c(0,1.06), expand = FALSE ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(face = "bold", size = 7, vjust = -3),
    plot.subtitle = element_text(size = 6, vjust = -2),
    axis.line = element_blank(),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust = 0.5, size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 5),
    strip.background = element_rect(fill = "transparent", color = "white"),
    panel.background = element_rect(fill = "white"),  
    plot.background = element_rect(fill = "white"), 
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.key.height = unit(0.3, "lines"),
    legend.key.width = unit(0.3, "lines"), 
    legend.text = element_text(size = 4, margin = margin(b = 0.9)),
    legend.title = element_text(size = 4, margin = margin(b = 1.5)),
    legend.position = "none",               
    legend.margin = margin(0, 0, 0, -6)     
  )

RAD_DNV_structure

#####=====
#### H-5. Mapping ===================================================
### H-5-1. importing geographic information =========================
Env.GEO <- read.csv("Mantel_IBD+IBE.CSV", header = T)

Env.GEO2 <- Env.GEO %>%
  select(Populations, PopNum, Altitudes, Elevation, Latitude, Longitude)

### H-5-2. Plotting structure on map ================================
## Calculating mean q (proportion) of ancestral clusters
region_cluster_mean <- q_df_ordered2 %>%
  group_by(Region, Population) %>%
  summarise(mean_q = mean(q), .groups = "drop") %>%
  group_by(Region) %>%
  mutate(
    total_q = sum(mean_q),
    q_scaled = mean_q / total_q
  ) %>%
  ungroup()

region_pie <- region_cluster_mean %>%
  select(Region, Population, q_scaled) %>%
  tidyr::pivot_wider(names_from = Population, values_from = q_scaled) %>%
  left_join(Env.GEO2, by = c("Region" = "Populations")) 

## Plot world map using the maps package
world_map <- map_data("world")

## Plot the world map with your data overlaid
par(mar = c(1, 1, 1, 1))

Map_DNV <- ggplot() +
  # Base map of the world
  geom_map(data = world_map, map = world_map,
           aes(x = long, y = lat, map_id = region), 
           fill = "lightgray", color = "white", size = 0.1) +
  # Add data points with Latitude and Longitude
  geom_scatterpie(data = region_pie,
                  aes(x = Longitude, y = Latitude),
                  cols = c("P1", "P2", "P3", "P4", "P5"),
                  color = NA, alpha = 0.8, pie_scale = 4) +
  # Add labels for Populations near each point
  geom_text(data = region_pie, aes(x = Longitude, y = Latitude, label = Region),
            size = 2, vjust = -2.5, hjust = 0.5) +
  # Add map customization
  theme_bw() +
  labs(title = " ",
       subtitle = "Colored by Ancestral Clusters",
       color = "Clusters",
       x = "Longitude", y = "Latitude") +
  scale_fill_manual(values = q_palette,
                    name = "Ancestral clusters") +
  theme(legend.position = "right",
        plot.background  = element_rect(fill = "white", color = NA),
        axis.title.y = element_text(hjust = 0.5, vjust = 1, size = 6),
        axis.title.x = element_text(hjust = 0.5, vjust = 0, size = 6),
        axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 5),
        axis.text.y = element_text(hjust = 0, vjust = 0, size = 5),
        plot.title = element_text(size = 6, vjust = -4),   
        plot.subtitle = element_text(size = 6, vjust = -2), 
        plot.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "mm"),
        # Reduce the size of the legend
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5), 
        legend.margin = margin(0, 0, 0, -3)) + 
  coord_fixed(ratio = 1.1, xlim = c(126, 130), ylim = c(35, 38))  

Map_DNV
#####=====

##### I. Drawing map ================================================
#### I-1. Importing coordinates =====================================
cities.pop <- read.csv("Geographic_information.CSV", header=T)
head(cities.pop)
str(cities.pop)

#### I-2. Drawing Korea =============================================
register_stadiamaps(key = "01bcd91b-8243-40ce-b50b-###########")
Korea <- get_stadiamap(bbox = c(left = 125.5, bottom = 34.5, 
                                right = 130, top = 38.5), 
                       zoom = 8, maptype = c("stamen_terrain_background"), messaging = FALSE)
ggmap(Korea)

#### I-3. Adding populations ========================================
cities.pop <- cities.pop %>%
  mutate(Alt_group = case_when(
    Altitude < 100 ~ "LOW (< 100m)",
    Altitude < 400 ~ "MID (300 - 400m)",
    TRUE ~ "HIGH (> 700m)"
  )) %>%
  mutate(Alt_group = factor(
    Alt_group,
    levels = c("HIGH (> 700m)", "MID (300 - 400m)", "LOW (< 100m)")
  ))

#### I-4. Drawing a full map ========================================
a <- ggmap(Korea) +
  geom_point(
    data = cities.pop,
    aes(x = as.numeric(DD_long), y = as.numeric(DD_lat)),
    size = 1.2,
    show.legend = FALSE
  ) +
  geom_label_repel(
    data = cities.pop,
    aes(
      x = as.numeric(DD_long),
      y = as.numeric(DD_lat),
      label = Populations,
      fill = Alt_group          # ← Altitude 구간에 따라 색 구분
    ),
    nudge_y = 0.2,
    segment.size  = 0.25,
    segment.color = "grey20",
    direction     = "both",
    size = 2.0,
    label.padding = 0.25,
    point.padding = 1e-04,
    label.r = 0.25,
    show.legend = TRUE,
    min.segment.length= 0
  ) +
  scale_fill_manual(
    values = c("LOW (< 100m)" = "#3B9AB2", 
               "MID (300 - 400m)" = "#78AB46",
               "HIGH (> 700m)" = "#E5C300"),  
    name = "Elevational groups"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        label = "",   
        label.size = 0
      ),
      keyheight = unit(5, "pt"),
      keywidth  = unit(8, "pt"),
      byrow     = TRUE
    )
  ) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme_article(base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 7, face = "plain", vjust = 0.5),
    axis.title.y = element_text(size = 7, face = "plain", vjust = 1),
    axis.text.x = element_text(size = 7, face = "plain", vjust = 0.5),
    axis.text.y = element_text(size = 7, face = "plain", vjust = 1),
    legend.position = "right", 
    legend.title      = element_text(size = 5),
    legend.text       = element_text(size = 4),
    legend.key.height = unit(6, "pt"),
    legend.key.width  = unit(7, "pt"),
    legend.spacing.y  = unit(0, "pt"),
    legend.spacing.x  = unit(0, "pt"),
    legend.margin     = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,-8),
    plot.margin = margin(2, 2, 2, 2)
  )
a

#### I-5. Adding genetic clusters on the maps =======================
### I-5-1. Calculating mean q (proportion) of ancestral clusters ====
region_cluster_mean <- q_df_ordered2 %>%
  group_by(Region, Population) %>%
  summarise(mean_q = mean(q), .groups = "drop") %>%
  group_by(Region) %>%
  mutate(
    total_q = sum(mean_q),
    q_scaled = mean_q / total_q
  ) %>%
  ungroup()

### I-5-2. Pie per region ===========================================
region_pie <- region_cluster_mean %>%
  select(Region, Population, q_scaled) %>%
  tidyr::pivot_wider(names_from = Population, values_from = q_scaled) %>%
  left_join(Env.GEO2, by = c("Region" = "Populations")) 

### I-5-3. Visualization ============================================
fig2 <- ggmap(Korea) +
  geom_scatterpie(data = region_pie,
                  aes(x = Longitude, y = Latitude),
                  cols = c("P1", "P2", "P3", "P4", "P5"),
                  color = NA, alpha = 1, pie_scale = 4) +
  geom_text(data = region_pie, aes(x = Longitude, y = Latitude, label = Region),
            size = 2.5, vjust=-3, hjust=1.2) +
  labs(color = "Clusters",
       x = "Longitude", y = "Latitude") +
  scale_fill_manual(values = q_palette,
                    name = "Ancestral clusters") +
  theme_article(base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 10, face = "plain", vjust = 0),
    axis.title.y = element_text(size = 10, face = "plain", vjust = 0),
    axis.text.x = element_text(size = 10, face = "plain", vjust = 0),
    axis.text.y = element_text(size = 10, face = "plain", vjust = .5),
    legend.position = "right", 
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(10, "pt"),
    legend.key.width  = unit(10, "pt"),
    legend.spacing.y  = unit(0, "pt"),
    legend.spacing.x  = unit(0, "pt"),
    legend.margin     = margin(0,0,0,-3),
    legend.box.margin = margin(0,0,0,-3),
    plot.margin = margin(-1, 1, -2, 1)
  )

fig2
### I-5-4. Export Figure 2 ==========================================
pdf("Figure2.pdf", width = 5, height = 5)
grid.arrange(fig2, nrow = 1, ncol = 1)
dev.off()

dev.off()
#### END =============
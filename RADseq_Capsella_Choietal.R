##### PolyRAD -----------------------------------------------------------
##### Library packages --------------------------------------------------
BiocManager::install("promises")
install.packages("prettymapr")

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

#### Pre-Steps . Population (ID) data ======================================================
PopMap <- read.table("Populations.txt",
                     header=TRUE, sep="\t", stringsAsFactors = TRUE) ## 181 individuals
PopMap %>%
  group_by(Populations) %>%
  count()

### =====

##### Statistical analyses with De novo results ==============
##### Step 1. Importing dataset (De novo without filtering, directly imported with polyRAD) =======
RAD_DNV_4 <- polyRAD::readStacks("catalog.alleles.tsv",
                                 "/",
                                 min.ind.with.reads = 168,
                                 min.ind.with.minor.allele = 10,
                                 possiblePloidies = list(2),
                                 contamRate = 0.001,
                                 version = 2)

RAD_DNV_4 ## 4,681 variants

#### Step 1-1. Quality control and parameter estimation =====================================================
# overdispersionP_RAD_DNV_4 <- TestOverdispersion(RAD_DNV_4, to_test = 1:130)
# 
# sapply(overdispersionP_RAD_DNV_4[names(overdispersionP_RAD_DNV_4) != "optimal"],
#        quantile, probs = c(0.01, 0.25, 0.5, 0.75, 0.99))
# 
# ovdisp_RAD_DNV_4 <- overdispersionP_RAD_DNV_4$optimal ## 130 ## Consider higher values;; 
ovdisp_RAD_DNV_4 <- 130

hindhe_RAD_DNV_4 <- HindHe(RAD_DNV_4)
hindheByLoc_RAD_DNV_4 <- colMeans(hindhe_RAD_DNV_4, na.rm = TRUE)
hist(hindheByLoc_RAD_DNV_4, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus (RAD_DNV_4)")
abline(v = 0.5, col = "blue", lwd = 2) ### The peak below 0.5 indicates well-behaved diploid loci.

ALFHWE_RAD_DNV_4 <- AddAlleleFreqHWE(RAD_DNV_4)
theseloci_RAN_DNV_4 <- GetLoci(ALFHWE_RAD_DNV_4)[ALFHWE_RAD_DNV_4$alleles2loc[ALFHWE_RAD_DNV_4$alleleFreq >= 0.05 & ALFHWE_RAD_DNV_4$alleleFreq < 0.5]]
theseloci_RAN_DNV_4 <- unique(theseloci_RAN_DNV_4)
myhindheByLoc_RAD_DNV_4 <- colMeans(hindhe_RAD_DNV_4[ALFHWE_RAD_DNV_4$taxaPloidy == 2, theseloci_RAN_DNV_4], na.rm = TRUE)
hist(myhindheByLoc_RAD_DNV_4, col = "lightgrey",
     xlab = "Hind/He", main = "Histogram of Hind/He by locus, MAF >= 0.05")
abline(v = 0.5, col = "blue", lwd = 2)

#### Step 1-2. Filtering based on the Hind/He ===================================
hh_RAD_DNV_4 <- HindHe(RAD_DNV_4, ploidy = RAD_DNV_4$possiblePloidies[[1]])
str(hh_RAD_DNV_4)

TotDepthT <- rowSums(RAD_DNV_4$locDepth) # depth for each sample

hhByInd <- rowMeans(hh_RAD_DNV_4, na.rm=TRUE) # Hind/He for each sample

plot(TotDepthT, hhByInd, log = "x",
     xlab="Depth", ylab="Hind/He", main="Sample")
abline(h=0.75, lty=2)

hhByInd[hhByInd < 0.45] ## 76, 78, 95, 137, 156, 167 지우기

hh_RAD_DNV_4 <- hh_RAD_DNV_4[hhByInd > 0.45,]
rownames(hh_RAD_DNV_4)

## Filtering by individuals
RADdata_RAD_DNV_4 <- SubsetByTaxon(RAD_DNV_4, rownames(hh_RAD_DNV_4))
RADdata_RAD_DNV_4 ## 4,681 loci and 181 taxa

hh2_RAD_DNV_4 <- HindHe(RADdata_RAD_DNV_4, ploidy = RADdata_RAD_DNV_4$possiblePloidies[[1]]) ## ploidy, n.gen.selfing 값 바꿔도 hhByLoc 값 안 바뀜.

hhByLoc <- colMeans(hh2_RAD_DNV_4, na.rm = TRUE)
hist(hhByLoc, breaks=30)

mean(hhByLoc)

ExpectedHindHe(RADdata_RAD_DNV_4, ploidy = 2, reps = 100,  
               errorRate = 0.001, contamRate = 0.001, overdispersion = ovdisp_RAD_DNV_4)

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
#### Step 1-2-1. filter by locus ===================================================
RADdata_RAD_DNV_4_filtered <- SubsetByLocus(RADdata_RAD_DNV_4, keeploci_di)
RADdata_RAD_DNV_4_filtered ## 233 loci

RADdata_RAD_DNV_4_filtered <- IterateHWE(RADdata_RAD_DNV_4_filtered)

RAD_DNV4.ind <- Export_adegenet_genind(RADdata_RAD_DNV_4_filtered)

#### Step 1-2-2. Adding Population data into genind =============
POPIND_181 <- read.table("PopInd.txt",
                         header=TRUE, sep="\t", stringsAsFactors = TRUE)

RAD_DNV4.ind@pop <- POPIND_181$Populations
RAD_DNV4.ind
pop(RAD_DNV4.ind) 

#### Step 1-2-3. Conversion genind to hierfstat =================
fstat.RAD_DNV4 <- genind2hierfstat(RAD_DNV4.ind)
head(fstat.RAD_DNV4)

#### Step 1-2-4. Calculation Fst ================================
RAD_DNV4_FST <- hierfstat::genet.dist(fstat.RAD_DNV4, method = "Fst", diploid = T)
RAD_DNV4_FST

#            BOP        CHJ        DGR        DMY        ICH        JCH        JSN        MZS        NWS
# CHJ 0.06154473                                                                                        
# DGR 0.05657315 0.06517032                                                                             
# DMY 0.07649530 0.06955525 0.07950466                                                                  
# ICH 0.05818817 0.05328270 0.06029613 0.05221685                                                       
# JCH 0.05750654 0.05125774 0.04426715 0.06296459 0.04057695                                            
# JSN 0.04659523 0.06108026 0.04273158 0.06075160 0.04579819 0.03687662                                 
# MZS 0.06713213 0.06026828 0.06117197 0.07386671 0.05706859 0.04964539 0.04668644                      
# NWS 0.11765733 0.11942923 0.11648677 0.08686024 0.09669586 0.11423203 0.11012444 0.11633364           
# TAB 0.05874449 0.07340050 0.05296610 0.08645307 0.06379443 0.05050898 0.04695557 0.05929538 0.12560707

#### Step 1-2-5. Calculation of Heterozygosity indices ==========
## Calculate heterozygosity per site
H.RAD_DNV4.ind = basic.stats(RAD_DNV4.ind, diploid = TRUE)
H.RAD_DNV4.ind

## Mean observed heterozygosity per site
Ho.RAD_DNV4.ind = apply(H.RAD_DNV4.ind$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
Ho.RAD_DNV4.ind
summary(Ho.RAD_DNV4.ind)

## Mean expected heterozygosity per site
He.RAD_DNV4.ind = apply(H.RAD_DNV4.ind$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
He.RAD_DNV4.ind
summary(He.RAD_DNV4.ind)

## Inbreeding coefficient (FIS) / Calculate mean FIS per site.
apply(H.RAD_DNV4.ind$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

### =====

#### Step 2. Population differentiation (Fst) ===============
### Step 2-1. RAD_DNV4.ind ================================
hierfstat::pairwise.WCfst(RAD_DNV4.ind, diploid = TRUE)

# Diploid = TRUE #
#            BOP        CHJ        DGR        DMY        ICH        JCH        JSN        MZS        NWS        TAB
# BOP         NA 0.03758774 0.04269075 0.05657215 0.03150443 0.03652055 0.02217829 0.05070774 0.09988610 0.03752242
# CHJ 0.03758774         NA 0.04603157 0.04489618 0.02129207 0.02496566 0.03182426 0.03882567 0.09807642 0.04723466
# DGR 0.04269075 0.04603157         NA 0.06339357 0.03853079 0.02743990 0.02284751 0.04836573 0.10211348 0.03642650
# DMY 0.05657215 0.04489618 0.06339357         NA 0.02463838 0.04084572 0.03567372 0.05583361 0.06771487 0.06475219
# ICH 0.03150443 0.02129207 0.03853079 0.02463838         NA 0.01140060 0.01344375 0.03319400 0.07285903 0.03477631
# JCH 0.03652055 0.02496566 0.02743990 0.04084572 0.01140060         NA 0.01005679 0.03064086 0.09497383 0.02684624
# JSN 0.02217829 0.03182426 0.02284751 0.03567372 0.01344375 0.01005679         NA 0.02473747 0.08835560 0.01998113
# MZS 0.05070774 0.03882567 0.04836573 0.05583361 0.03319400 0.03064086 0.02473747         NA 0.10041265 0.04043006
# NWS 0.09988610 0.09807642 0.10211348 0.06771487 0.07285903 0.09497383 0.08835560 0.10041265         NA 0.10654692
# TAB 0.03752242 0.04723466 0.03642650 0.06475219 0.03477631 0.02684624 0.01998113 0.04043006 0.10654692         NA

RAD_DNV4.ind_FST <- mat_gen_dist(x = RAD_DNV4.ind, dist = "FST", null_val = TRUE) ## FALSE default
RAD_DNV4.ind_FST

RAD_DNV4.ind_LinFST <- mat_gen_dist(x = RAD_DNV4.ind, dist = "FST_lin", null_val = TRUE) ## FALSE default
RAD_DNV4.ind_LinFST

### =====

#### Step 3. DAPC & PCA ==============================================
### Step 3-1. DAPC RAD_DNV4.ind ==================================================================================
set.seed(300)
x = adegenet::tab(RAD_DNV4.ind, NA.method = "mean")
crossval = xvalDapc(x, RAD_DNV4.ind$pop, result = "groupMean", xval.plot = TRUE)
# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
crossval$`Number of PCs Achieving Highest Mean Success`
## [1] "60"
crossval$`Number of PCs Achieving Lowest MSE`
## [1] "60"
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)
# Run a DAPC using site IDs as priors
dapc1 = adegenet::dapc(RAD_DNV4.ind, RAD_DNV4.ind$pop, n.pca = numPCs, n.da = 2)
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
ind_coords$Ind = indNames(RAD_DNV4.ind)
# Add a column with the site IDs
ind_coords$Site = RAD_DNV4.ind$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = brewer.pal(nPop(RAD_DNV4.ind), "Set3")
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
l <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0, color = "gray50", size = 0.3)+
  geom_vline(xintercept = 0, color = "gray50", size = 0.3)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 1.5, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 1.5, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab,
       title = "DAPC",
       subtitle = "STACKS2-PolyRAD (233 snps)")+
  theme_article(base_size = 7, base_family = "sans")+
  theme(axis.title = element_text(face = "bold", vjust=0.5, size=4),
        axis.text = element_text(face = "bold", hjust = -3, vjust=0.5, size=4),
        axis.text.x = element_text(angle = 0, hjust=0, vjust=0, size=4),
        axis.text.y = element_text(hjust=0, vjust=0.5, size=4),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=5, face="bold", colour = "Navy", vjust = -2),
        plot.subtitle = element_text(size=4, face="bold", colour = "black", vjust = -2))
l

### Step 3-2. PCA RAD_DNV4.ind ===================================
set.seed(300)
# Replace missing data with the mean allele frequencies
x = tab(RAD_DNV4.ind, NA.method = "mean")
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
ind_coords$Ind = indNames(RAD_DNV4.ind)
# Add a column with the site IDs
ind_coords$Site = RAD_DNV4.ind$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = brewer.pal(nPop(RAD_DNV4.ind), "Set3")
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
m <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0, color = "gray50", size = 0.3)+
  geom_vline(xintercept = 0, color = "gray50", size = 0.3)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 1.5, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 1.5, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab,
       title = "PCA",
       subtitle = "STACKS2-PolyRAD (233 snps)")+
  theme_article(base_size = 7, base_family = "sans")+
  theme(axis.title = element_text(face = "bold", vjust=0.5, size=4),
        axis.text = element_text(face = "bold", hjust = -3, vjust=0.5, size=4),
        axis.text.x = element_text(angle = 0, hjust=0, vjust=0, size=4),
        axis.text.y = element_text(hjust=0, vjust=0.5, size=4),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=5, face="bold", colour = "Navy", vjust = -2),
        plot.subtitle = element_text(size=4, face="bold", colour = "black", vjust = -2))
m
### =====

#### Step 4. AMOVA ================================================
#### Step 4-1. RAD_DNV4.ind =======================================
# adding Elevations
high_elevation <- c("BOP", "DGR", "MZS", "TAB")
mid_elevation  <- c("NWS", "JSN", "JCH")
low_elevation  <- c("CHJ", "ICH", "DMY")

POPIND_181$Elevation <- with(POPIND_181, ifelse(Populations %in% high_elevation, "High",
                                                ifelse(Populations %in% mid_elevation, "Mid",
                                                       ifelse(Populations %in% low_elevation, "Low", NA))))

POPIND_181$Elevation <- factor(POPIND_181$Elevation, levels = c("Low", "Mid", "High"))

strata(RAD_DNV4.ind) <- POPIND_181
agc.RAD_DNV4.ind <- as.genclone(RAD_DNV4.ind)
agc.RAD_DNV4.ind

### Step 4-1-a. Populations =======================================
amova.RAD_DNV4 <- poppr::poppr.amova(agc.RAD_DNV4.ind, ~Populations,
                                          method = c("ade4"), within = FALSE, nperm = 9999,
                                          filter = TRUE, threshold = 0.05)
amova.RAD_DNV4

# $results
# Df    Sum Sq  Mean Sq
# Between samples   9  636.6877 70.74308
# Within samples  171 3478.9200 20.34456
# Total           180 4115.6077 22.86449
# 
# $componentsofcovariance
# Sigma         %
# Variations  Between samples  2.82925  12.20883
# Variations  Within samples  20.34456  87.79117
# Total variations            23.17381 100.00000
# 
# $statphi
# Phi
# Phi-samples-total 0.1220883

amova.Test.RAD_DNV4 <- ade4::randtest(amova.RAD_DNV4, nrepet = 9999) # Test for significance
print(amova.Test.RAD_DNV4)

# Observation: 2.82925 
# 
# Based on 9999 replicates
# Simulated p-value: 1e-04 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 42.0936749727 -0.0006871735  0.0045198086

plot(amova.Test.RAD_DNV4)

### Step 4-1-b. Elevations =======================================
amova.RAD_DNV4_2 <- poppr::poppr.amova(agc.RAD_DNV4.ind, ~Elevation/Populations,
                                     method = c("ade4"), within = FALSE, nperm = 9999,
                                     filter = TRUE, threshold = 0.05)
amova.RAD_DNV4_2

# $results
# Df    Sum Sq  Mean Sq
# Between Elevation                  2  157.8992 78.94960
# Between samples Within Elevation   7  478.7885 68.39836
# Within samples                   171 3478.9200 20.34456
# Total                            180 4115.6077 22.86449
# 
# $componentsofcovariance
# Sigma          %
# Variations  Between Elevation                 0.2083499   0.896605
# Variations  Between samples Within Elevation  2.6847310  11.553371
# Variations  Within samples                   20.3445614  87.550024
# Total variations                             23.2376422 100.000000
# 
# $statphi
# Phi
# Phi-samples-total     0.12449976
# Phi-samples-Elevation 0.11657896
# Phi-Elevation-total   0.00896605

amova.Test.RAD_DNV4_2 <- ade4::randtest(amova.RAD_DNV4_2, nrepet = 9999) # Test for significance
print(amova.Test.RAD_DNV4_2)

# Adjustment method for multiple comparisons:   none 
# Permutation number:   9999 
#                           Test        Obs     Std.Obs   Alter Pvalue
# 1    Variations within samples 20.3445614 -41.6989809    less 0.0001
# 2   Variations between samples  2.6847310  29.7293922 greater 0.0001
# 3 Variations between Elevation  0.2083499   0.8003741 greater 0.1920

### =====

#### Step 5. Mantel's tests =======================================
#### https://stats.oarc.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/

Env.GEO <- read.csv("Mantel_IBD+IBE.CSV", header=T)

Env.GEO$Populations <- as.factor(as.character(Env.GEO$Populations))
Env.GEO$Altitudes <- as.factor(as.character(Env.GEO$Altitudes))

### Normalization for Environmental variables
Env.GEO.scale <- as.data.frame(scale(Env.GEO[, -c(1:6, 35:36)]))

### Adding
Env.GEO.scale <- cbind(Env.GEO[, 1:6], Env.GEO.scale, Env.GEO[, 35:36])
str(Env.GEO.scale)

#### Step 5-1. Calculation of distances ===========================
#### Step 5-1-a. Geographical distance (Euclidean) ================
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

### distance [,35:36]
# GEO.Dis.mat <- as.matrix(Env.GEO.scale[,35:36])
# GEO.Dis.mat.dist = ecodist::distance(GEO.Dis.mat, method="euclidean") %>%
#   as.matrix()
# rownames(GEO.Dis.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
# colnames(GEO.Dis.mat.dist) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")

### Real distance (km)
Coords <- as.data.frame(Env.GEO.scale[,6:5])
dist_mat <- distm(Coords, fun = distHaversine)  # 또는 distGeo
dist_mat_km <- dist_mat / 1000
rownames(dist_mat_km) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
colnames(dist_mat_km) <- c("BOP","CHJ","DGR","DMY","ICH","JCH","JSN","MZS","NWS","TAB")
GEO.Dis.mat.dist_km <- dist_mat_km

#### Step 5-1-b. Environmental distance (Gower) ====-================
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

### =====

#### Step 5-2. build a matrix =====================================
str(Env.GEO.scale)
Env.GEO.BIO.mat <- as.matrix(Env.GEO.scale[,4:34])
Env.GEO.BIO.mat

### Correlations between BIO- variables
cor_mat <- cor(Env.GEO.scale[,16:34], method = "pearson")
cor_mat_upper <- cor_mat
cor_mat_upper[lower.tri(cor_mat_upper, diag = TRUE)] <- NA
cor_mat_upper <- round(cor_mat_upper, 2)
print(cor_mat_upper)


#### Step 5-2. Mantel's test ==================================
### Step 5-2-1. Geographic variables ==========================
## Elevation
Mantel_DNV4_Ele = vegan::mantel(RAD_DNV4.ind_FST, GEO.Ele.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_Ele

## Latitudes
Mantel_DNV4_Lat = vegan::mantel(RAD_DNV4.ind_FST, GEO.Lat.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_Lat

## Longitudes
Mantel_DNV4_Lon = vegan::mantel(RAD_DNV4.ind_FST, GEO.Lon.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_Lon

## Distances
Mantel_DNV4_Dist = vegan::mantel(RAD_DNV4.ind_FST, GEO.Dis.mat.dist_km, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_Dist

### Step 5-2-2. Environmental variables ==========================
## BIO1
Mantel_DNV4_BIO1 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO1.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO1

## BIO2
Mantel_DNV4_BIO2 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO2.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO2

## BIO3
Mantel_DNV4_BIO3 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO3.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO3

## BIO4
Mantel_DNV4_BIO4 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO4.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO4

## BIO5
Mantel_DNV4_BIO5 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO5.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO5

## BIO6
Mantel_DNV4_BIO6 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO6.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO6

## BIO7
Mantel_DNV4_BIO7 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO7.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO7

## BIO8
Mantel_DNV4_BIO8 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO8.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO8

## BIO9
Mantel_DNV4_BIO9 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO9.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO9

## BIO10
Mantel_DNV4_BIO10 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO10.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO10

## BIO11
Mantel_DNV4_BIO11 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO11.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO11

## BIO12
Mantel_DNV4_BIO12 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO12.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO12

## BIO13
Mantel_DNV4_BIO13 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO13.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO13

## BIO14
Mantel_DNV4_BIO14 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO14.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO14

## BIO15
Mantel_DNV4_BIO15 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO15.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO15

## BIO16
Mantel_DNV4_BIO16 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO16.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO16

## BIO17
Mantel_DNV4_BIO17 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO17.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO17

## BIO18
Mantel_DNV4_BIO18 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO18.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO18

## BIO19
Mantel_DNV4_BIO19 = vegan::mantel(RAD_DNV4.ind_FST, GEO.BIO19.mat.dist, method = "pearson", permutations = 9999, na.rm = TRUE)
Mantel_DNV4_BIO19

### =====

#### Step 5-3. Visualization of IBD and IBE =======================
### Step 5-3-1. Tabling ===========================================
## Step 5-3-1-a. Geographic variables =============================
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

### =====
## Step 5-3-1-b. Environmental variables =============================
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

#### =====
## Step 5-3-1-c. Genetic distances ================================
# FST
RAD_DNV4.ind_FST
lower.tri(RAD_DNV4.ind_FST)
RAD_DNV4.ind_FST[lower.tri(RAD_DNV4.ind_FST)] <- NA
RAD_DNV4.ind_FST.melt <- melt(RAD_DNV4.ind_FST, varnames = c("row", "col"))
colnames(RAD_DNV4.ind_FST.melt) <- c("Pop1","Pop2","FST")
RAD_DNV4.ind_FST.melt <- RAD_DNV4.ind_FST.melt %>%
  na.omit(FST)
RAD_DNV4.ind_FST.melt

# Linearized FST
RAD_DNV4.ind_LinFST
lower.tri(RAD_DNV4.ind_LinFST)
RAD_DNV4.ind_LinFST[lower.tri(RAD_DNV4.ind_LinFST)] <- NA
RAD_DNV4.ind_LinFST.melt <- melt(RAD_DNV4.ind_LinFST, varnames = c("row", "col"))
colnames(RAD_DNV4.ind_LinFST.melt) <- c("Pop1","Pop2","LinFST")
RAD_DNV4.ind_LinFST.melt <- RAD_DNV4.ind_LinFST.melt %>%
  na.omit(LinFST)
RAD_DNV4.ind_LinFST.melt

#### =====
## Step 5-3-1-d. Combining melt-data ==============================
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

head(a)

#### =====

#### Step 5-3-2. Plotting =========================================
### Step 5-3-2-a. Importing data ==================================
IBDIBE <- read.csv("IBDIBE_Jun_bc.csv")

### Step 5-3-2-b. Linear regression ===============================
### Step 5-3-2-b-1. IBD ===========================================
IBD_Dist <- lm(LinFST ~ Distance, data=IBDIBE)
IBD_Dist$coefficients
summary(IBD_Dist)

IBD_Lat <- lm(LinFST ~ Latitude, data=IBDIBE)
IBD_Lat$coefficients
summary(IBD_Lat)

### Step 5-3-2-b-2. IBE ===========================================
IBE_BIO3 <- lm(LinFST ~ BIO3, data=IBDIBE)
IBE_BIO3$coefficients
summary(IBE_BIO3)

### Step 5-3-2-c. Drawing figures =================================
IBD_Dist_plot <- ggplot(IBDIBE, aes(Distance, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=1.5, alpha=1) +
  geom_text_repel(size=2, segment.size=0.1)+
  labs(
    x = "Distances [km]", 
    y = "Linearized FST", 
    fill = "Populations",
    subtitle = bquote("Mantel r = 0.487 (" * italic(p) * " = 0.013)")) +
  theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.title.x = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.text.x = element_text(face = "bold", hjust=0.5, vjust=1.0, size=6),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=6, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=7, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=6, face="bold.italic", colour = "black", vjust = 0))

IBD_Dist_plot


IBD_Lat_plot <- ggplot(IBDIBE, aes(Latitude, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=1.5, alpha=1) +
  geom_text_repel(size=2, segment.size=0.1)+
  labs(
    x = "Latitude [°N]", 
    y = "Linearized FST", 
    fill = "Populations",
    subtitle = bquote("Mantel r = 0.524 (" * italic(p) * " = 0.005)")) +
    theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.title.x = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.text.x = element_text(face = "bold", hjust=0.5, vjust=1.0, size=6),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=6, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=7, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=6, face="bold.italic", colour = "black", vjust = 0))

IBD_Lat_plot

IBE_BIO3_plot <- ggplot(IBDIBE, aes(BIO3, LinFST, group=Pop1, label=Pairwise)) + 
  geom_point(aes(fill=Pop1), shape=20, size=1.5, alpha=1) +
  geom_text_repel(size=2, segment.size=0.1)+
  labs(
    x = "BIO3 (Isothermality)", 
    y = "Linearized FST", 
    fill = "Populations",
    subtitle = bquote("Mantel r = 0.471 (" * italic(p) * " = 0.018)")) +
  theme_bw(base_family = "sans")+
  theme(
    axis.title.y = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.title.x = element_text(face = "bold", hjust=0.5, vjust=0.5, size=6),
    axis.text.x = element_text(face = "bold", hjust=0.5, vjust=1.0, size=6),
    axis.text.y = element_text(hjust= 0.5, vjust=0.5, size=6, angle=90),
    legend.position = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(size=7, face="bold", colour = "Navy", vjust = 0),
    plot.subtitle = element_text(size=6, face="bold.italic", colour = "black", vjust = 0)) 

IBE_BIO3_plot

#### =====

#### Step 6. MMRR with CA =========================================
### Step 6-1. Building MMRR matrix ================================

MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  x <- x[!is.na(x)]  # NA 제거
  x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
  return(x)
}

### ====
### Step 6-2. Building models =====================================
## Scaling distance and latitudes
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



### Step 6-3. Genetic distances =======================================
RAD_DNV4.ind_FST

### Step 6-4. Run MMRR ================================================
PopGenReport::lgrMMRR(RAD_DNV4.ind_FST, Xmats.geo.BIO3, nperm=9999)

### Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
## GEO + BIO3
DNV_Dist_BIO3 <- yhat::commonalityCoefficients(IBDIBE, "FST", list(c("Distance","Latitude","BIO3")), imat=FALSE)
DNV_Dist_BIO3

#
Ele_Lat <- lm(Distance ~ BIO3, data=IBDIBE)
Ele_Lat$coefficients
summary(Ele_Lat)

### =============

#### Step 7. STRUCTURE with LEA ===============
### Step 7-1. converting genind into vcf ======
RAD_DNV4.ind
RAD_DNV4.ind@pop <- NULL
loci_obj <- as.loci(RAD_DNV4.ind)
write.vcf(loci_obj, file = "RAD_DNV4.vcf")

### Step 7-2. Calling vcf2geno ================
LEA::vcf2geno("RAD_DNV4.vcf",
              "RAD_DNV4.geno")

### Step 7-3. Running SNMF ===================== 
### DNV
project_DNV = snmf("RAD_DNV4.geno",
                   K = 1:20,
                   entropy = TRUE,
                   repetitions = 200,
                   project = "new", CPU = 20,
                   ploidy = 2)

jpeg("SNMF_RAD_DNV4.jpeg", width = 8, height = 8, units = "cm", res = 480)
par(mar = c(4, 4, 1, 1))
par(cex.axis = 0.6, cex.lab = 0.6)
plot(project_DNV, col="grey30", pch = 19, cex = 0.6)
dev.off()

### K selection (K=2)
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

### K selection (K=3)
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

### K selection (K=4)
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

### K selection (K=5)
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

### Step 7-4. Plotting population structure =============
## K=2
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


## K=3
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


## K=4
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

## K=5
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

### Step 7-5. Mapping =========================
## importing geograpic information
Env.GEO <- read.csv("Mantel_IBD+IBE.CSV", header = T)

Env.GEO2 <- Env.GEO %>%
  select(Populations, PopNum, Altitudes, Elevation, Latitude, Longitude)

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

### =====

#### Step 8. Drawing map =========================
### Step 8-1. Importing coordinates ==============
cities.pop <- read.csv("Geographic_information.CSV", header=T)
head(cities.pop)
str(cities.pop)

### Step 8-2. Drawing Korea ======================
register_stadiamaps(key = "01bcd91b-8243-40ce-b50b-123123123")
Korea <- get_stadiamap(bbox = c(left = 125.5, bottom = 34.5, 
                                right = 130, top = 38.5), 
                       zoom = 8, maptype = c("stamen_terrain_background"), messaging = FALSE)
ggmap(Korea)

### Step 8-3. Adding populations =================
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

### Step 8-4. Drawing a full map =================
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

### Step 8-6. Adding genetic clusters on the maps =================
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

b <- ggmap(Korea) +
  geom_scatterpie(data = region_pie,
                  aes(x = Longitude, y = Latitude),
                  cols = c("P1", "P2", "P3", "P4", "P5"),
                  color = NA, alpha = 0.9, pie_scale = 4) +
  geom_text(data = region_pie, aes(x = Longitude, y = Latitude, label = Region),
            size = 1.7, vjust = -2.5, hjust = 0.8) +
  labs(color = "Clusters",
       x = "Longitude", y = "Latitude") +
  scale_fill_manual(values = q_palette,
                    name = "Ancestral clusters") +
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
    legend.box.margin = margin(0,0,0,-9),
    plot.margin = margin(-1, 1, -2, 1)
  )

b

#### END =============
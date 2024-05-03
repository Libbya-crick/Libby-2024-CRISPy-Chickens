install.packages("factoextra")
install.packages("tidyverse")
install.packages("scatterplot3d")


library(Matrix)
library(pheatmap)
library(tidyr)
library(scatterplot3d)
library (ggplot2)
library(gridExtra)
library(plotly)
library(viridisLite)
## hiearchical clustering of all variables
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

##load csv quantification file, and clean so that coordinates and each channel has its own column in a data frame

file.HCR2 = "HCR-SOX2-TBXT-F2RL1-Bcat_3D_nuclei_features_table_fullstack.csv"
data<- read.csv(file = file.HCR5)
data<- as.data.frame(data)
data

# Separate the centers into three columns of x y z
df_split <- separate(data, Center.Mass, into = c("Center.Mass.Z", "Center.Mass.X", "Center.Mass.Y"), sep = ", ", remove = TRUE)
df_split[,4:6] <- lapply(df_split[,4:6], function(x) gsub("\\(|\\)", "", x))
# Convert the values to numeric
df_split$Center.Mass.X <- as.numeric(gsub("\\)", "", df_split$Center.Mass.X))
df_split$Center.Mass.Z <- as.numeric(gsub("\\(", "", df_split$Center.Mass.Z))
df_split$Center.Mass.Y <- as.numeric(gsub("\\)", "", df_split$Center.Mass.Y))

df_split

### Separate the nuclear averages into columns based on channels

df_split <- separate(df_split, nuclear_avgs, into = c("nuclear_avgs.DAPI", 
                                                      "nuclear_avgs.2", 
                                                      "nuclear_avgs.3",
                                                      "nuclear_avgs.4", 
                                                      "nuclear_avgs.5"), 
                     sep = ", ", remove = TRUE)

df_split[,7:11] <- lapply(df_split[,7:11], function(x) gsub("\\[|\\]", "", x))


# Convert the values to numeric
df_split$nuclear_avgs.DAPI <- as.numeric(gsub("\\(", "", df_split$nuclear_avgs.DAPI))
df_split$nuclear_avgs.2 <- as.numeric(gsub("\\)", "", df_split$nuclear_avgs.2))
df_split$nuclear_avgs.3 <- as.numeric(gsub("\\)", "", df_split$nuclear_avgs.3))
df_split$nuclear_avgs.4 <- as.numeric(gsub("\\)", "", df_split$nuclear_avgs.4))
df_split$nuclear_avgs.5 <- as.numeric(gsub("\\)", "", df_split$nuclear_avgs.5))
df_split

### Separate the hood averages into columns based on channels 
## for HCR2 1. DAPI 2. SOX2 3. TBXT 4. F2RL1 5. B-cat (discard this one)

df_split <- separate(df_split, hood_avgs, into = c("hood_avgs.DAPI", 
                                                   "hood_avgs.2", 
                                                   "hood_avgs.3",
                                                   "hood_avgs.4", 
                                                    "hood_avgs.5"), 
                     sep = ", ", remove = TRUE)

df_split[,12:16] <- lapply(df_split[,12:16], function(x) gsub("\\[|\\]", "", x))

# Convert the values to numeric
df_split$hood_avgs.DAPI <- as.numeric(gsub("\\(", "", df_split$hood_avgs.DAPI))
df_split$hood_avgs.2 <- as.numeric(gsub("\\)", "", df_split$hood_avgs.2))
df_split$hood_avgs.3 <- as.numeric(gsub("\\)", "", df_split$hood_avgs.3))
df_split$hood_avgs.4 <- as.numeric(gsub("\\)", "", df_split$hood_avgs.4))
df_split$hood_avgs.5 <- as.numeric(gsub("\\)", "", df_split$hood_avgs.5))
df_split

### Separate the cyto averages into columns based on channels 
## for HCR2 1. DAPI 2. SOX2 3. TBXT 4. F2RL1 5. B-cat (discard this one)

df_split <- separate(df_split, cyto_avgs, into = c("cyto_avgs.DAPI", 
                                                   "cyto_avgs.2", 
                                                   "cyto_avgs.3",
                                                   "cyto_avgs.4", 
                                                    "cyto_avgs.5"), 
                     sep = ", ", remove = TRUE)

df_split[,17:21] <- lapply(df_split[,17:21], function(x) gsub("\\[|\\]", "", x))

# Convert the values to numeric
df_split$cyto_avgs.DAPI <- as.numeric(gsub("\\(", "", df_split$cyto_avgs.DAPI))
df_split$cyto_avgs.2 <- as.numeric(gsub("\\)", "", df_split$cyto_avgs.2))
df_split$cyto_avgs.3 <- as.numeric(gsub("\\)", "", df_split$cyto_avgs.3))
df_split$cyto_avgs.4 <- as.numeric(gsub("\\)", "", df_split$cyto_avgs.4))
df_split$cyto_avgs.5 <- as.numeric(gsub("\\)", "", df_split$cyto_avgs.5))
df_split

## now filter on Volumes incase the segmentation messed up
hist(df_split$Volume,500)
df_split_filtered<- df_split[which(df_split$Volume >500 & df_split$Volume <3000 ),]
df_split_filtered

## generate corrected value scores by dividing by the DAPI channel to account for light diffusion in3D structure

df_split_filtered$cor.nuclear_avgs.2 <- df_split_filtered$nuclear_avgs.2 / df_split_filtered$nuclear_avgs.DAPI
df_split_filtered$cor.nuclear_avgs.3 <- df_split_filtered$nuclear_avgs.3 / df_split_filtered$nuclear_avgs.DAPI
df_split_filtered$cor.nuclear_avgs.4 <- df_split_filtered$nuclear_avgs.4 / df_split_filtered$nuclear_avgs.DAPI
df_split_filtered$cor.nuclear_avgs.5 <- df_split_filtered$nuclear_avgs.5 / df_split_filtered$nuclear_avgs.DAPI

df_split_filtered$cor.hood_avgs.2 <- df_split_filtered$hood_avgs.2 / df_split_filtered$hood_avgs.DAPI
df_split_filtered$cor.hood_avgs.3 <- df_split_filtered$hood_avgs.3 / df_split_filtered$hood_avgs.DAPI
df_split_filtered$cor.hood_avgs.4 <- df_split_filtered$hood_avgs.4 / df_split_filtered$hood_avgs.DAPI
df_split_filtered$cor.hood_avgs.5 <- df_split_filtered$hood_avgs.5 / df_split_filtered$hood_avgs.DAPI


df_split_filtered$cor.cyto_avgs.2 <- df_split_filtered$cyto_avgs.2 / df_split_filtered$cyto_avgs.DAPI
df_split_filtered$cor.cyto_avgs.3 <- df_split_filtered$cyto_avgs.3 / df_split_filtered$cyto_avgs.DAPI
df_split_filtered$cor.cyto_avgs.4 <- df_split_filtered$cyto_avgs.4 / df_split_filtered$cyto_avgs.DAPI
df_split_filtered$cor.cyto_avgs.5 <- df_split_filtered$cyto_avgs.5 / df_split_filtered$cyto_avgs.DAPI

write.csv(df_split_filtered, file = "HCR_split_df.csv")

## find distribution of cor.hood avgs to set thresholds
mean2<- mean(df_split_filtered$cor.hood_avgs.2)
sd2<- sd(df_split_filtered$cor.hood_avgs.2)

p2<- ggplot(df_split_filtered, aes(x = Image, y = cor.hood_avgs.2 ))+
  geom_violin() + 
  geom_point(size = .5, color = "blue")+
  geom_hline(yintercept=mean2, col="red")+
  geom_hline(yintercept=mean2+sd2, col="grey")+
  geom_hline(yintercept=mean2-sd2, col="grey")+
  geom_hline(yintercept=mean2+(sd2*5), col="grey")+
  geom_hline(yintercept=mean2-(sd2*5), col="grey")

mean3<- mean(df_split_filtered$cor.hood_avgs.3)
sd3<- sd(df_split_filtered$cor.hood_avgs.3)

p3<- ggplot(df_split_filtered, aes(x = Image, y = cor.hood_avgs.3 ))+
  geom_violin() + 
  geom_point(size = .5, color = "blue")+
  geom_hline(yintercept=mean3, col="red")+
  geom_hline(yintercept=mean3+sd3, col="grey")+
  geom_hline(yintercept=mean3-sd3, col="grey")+
  geom_hline(yintercept=mean3+(sd3*5), col="grey")+
  geom_hline(yintercept=mean3-(sd3*5), col="grey")

mean4<- mean(df_split_filtered$cor.hood_avgs.4)
sd4<- sd(df_split_filtered$cor.hood_avgs.4)

p4<- ggplot(df_split_filtered, aes(x = Image, y = cor.hood_avgs.4 ))+
  geom_violin() + 
  geom_point(size = .5, color = "blue")+
  geom_hline(yintercept=mean4, col="red")+
  geom_hline(yintercept=mean4+sd4, col="grey")+
  geom_hline(yintercept=mean4-sd4, col="grey")+
  geom_hline(yintercept=mean4+(sd4*5), col="grey")+
  geom_hline(yintercept=mean4-(sd4*5), col="grey")


grid<- grid.arrange(p2, p3, p4, ncol = 3)
grid

##plot corrected values in various axis or in 3d

x<- df_split_filtered$Center.Mass.X
x_range <- range(x)
y<- df_split_filtered$Center.Mass.Y
y_range<-range(y)

n2<- ggplot(df_split_filtered[which(df_split_filtered$cor.hood_avgs.3 > (mean3)&
                                      df_split_filtered$cor.hood_avgs.3 < (mean3+(sd3*5))),], aes(x = Center.Mass.X, 
                                                                                                  y = Center.Mass.Y,
                                                                                                  col = cor.hood_avgs.3))+
  geom_point()+
  scale_color_viridis_c(option = "magma") +
  xlim(x_range) +  
  ylim(y_range)
n2

#in 3d
colors <- magma(length(high.values))[factor(high.values)]
x1<- df_split_filtered$Center.Mass.X[which(df_split_filtered$cor.hood_avgs.2 > (mean2+sd2))]
y1<- df_split_filtered$Center.Mass.Y[which(df_split_filtered$cor.hood_avgs.2 > (mean2+sd2))]
z1<- df_split_filtered$Center.Mass.Z[which(df_split_filtered$cor.hood_avgs.2 > (mean2+sd2))]
z1<- (z1)

scatterplot3d(y1, x1, z1, pch = 19,main="3D Scatter Plot", xlab="y", ylab="x", zlab="z", color=colors)

## to compare expression levels between channel of interest with TBXT channel

df<- df_split_filtered[which(
  (df_split_filtered$cor.hood_avgs.3 > (mean3)&
     df_split_filtered$cor.hood_avgs.3 < (mean3+(sd3*5)))|
    (df_split_filtered$cor.hood_avgs.4 > (mean4)&
       df_split_filtered$cor.hood_avgs.4 < (mean4+(sd4*5)))),
  c('cor.hood_avgs.2','cor.hood_avgs.3','cor.hood_avgs.4')]

##normalize from 0 to 1
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#normalise data using custom function
normalised_df <- as.data.frame(lapply(df, minMax))
head(normalised_df)

# Dissimilarity matrix
d <- dist(normalised_df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.5, hang = -1)

# Cut tree into groups
sub_grp <- cutree(hc1, h = 0.6)
rect.hclust(hc1, h = 0.6, border = 2:5)

# add group to dataframe so it can be plotted

df$sub_grp <- as.character(sub_grp)
df2<- df_split_filtered[which( (df_split_filtered$cor.hood_avgs.3 > (mean3)&
     df_split_filtered$cor.hood_avgs.3 < (mean3+(sd3*5)))|
    (df_split_filtered$cor.hood_avgs.4 > (mean4)&
       df_split_filtered$cor.hood_avgs.4 < (mean4+(sd4*5)))),
  c('Center.Mass.X','Center.Mass.Y','Center.Mass.Z')]
df<- cbind(df,df2)
df

## plot as scatter plot

n7<- ggplot(df, 
            aes(x = cor.hood_avgs.4, 
                y = cor.hood_avgs.3,
                color = sub_grp))+
  geom_point()+
  scale_fill_identity()
n7

## plot group of interest back onto tbxt expression plot
x<- df_split_filtered$Center.Mass.X
z<- df_split_filtered$Center.Mass.Z
x_range <- range(x)
y<- df_split_filtered$Center.Mass.Y
y_range<-range(y)

n3<- ggplot(df, aes(x = Center.Mass.X, 
                    y = Center.Mass.Y,
                    col = cor.hood_avgs.3))+
  geom_point()+
  scale_color_viridis_c(option = "magma") +
  xlim(x_range) +  
  ylim(y_range)
#n3

overlap<- geom_point(data = df[which( df$sub_grp == 7),], aes(x = Center.Mass.X, y = Center.Mass.Y), color = "blue", size = 3)

overlap1<- n3+overlap
overlap1



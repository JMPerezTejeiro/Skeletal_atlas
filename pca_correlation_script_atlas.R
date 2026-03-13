################################################################################
######### SCRIPT TO MAKE PCA and CORRELATION MATRIX ############################
################################################################################

#Libraries required:
library(ggplot2)
library(RColorBrewer)
library(DGEobj)
library(DGEobj.utils)
library(limma)
library(edgeR)
library(devtools)
library(ggfortify)
library(gdata)
library(stats)
library(dplyr)
library(pheatmap)

#Set working directory:
setwd("")

#Load the data, counts matrix obtained after executing do_counts script in the samples obtained by featureCounts.

counts <- read.table(file = "",
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)

#PCA:

#Now, execute PCA.
factors <- colnames(counts)
factors <- gsub("_[0-9]+(\\.[0-9]+)?$", "", factors) # Eliminate _1, _2, _n tags from samples, leaving unique name and grouping them into factors.


#3. Define groups.
group <- factor(factors)

#4. Normalize.
d <- DGEList(counts = counts, group = group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

ncounts <-  cpm (d, normalized.lib.sizes = TRUE, log = TRUE) # con el log sale mejor el PCA

#5. Represent the PCA.
pca <- prcomp(t(ncounts))

#Percentage of each component.

pca.stat <- stats:::summary.prcomp(pca)$importance

#6. Vector of colours.
# Create vector making each factor unique.
uniq_factors <- unique(factors)

# Color palette.
n<- length(uniq_factors)
col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <-  unlist(mapply(brewer.pal, col_pals$maxcolors, rownames(col_pals)))
col_vector <- unique(col_vector)

# Assign one color to one factor.
colors_factors <- as.vector(factors)

#Generate colours for plots.

ggcolors <- c()

# Assign colours
for (i in 1:length(uniq_factors))
{
  colors_factors <- replace(colors_factors,
                            colors_factors == uniq_factors[i],
                            col_vector[i])
}

ggcolors <- unique(colors_factors)
names(ggcolors) <- uniq_factors



#8. Plotting with ggplot2.

axis_pca <- pca$sdev^2 / sum(pca$sdev^2)

pca_x <- as.data.frame(pca$x)

pca_x <- pca_x[, c(1,2)] # We only need the first two dimmensions.

pca_x$muestra <- factors

gg_pca <- ggplot(data = pca_x,
                 aes(x = PC1, y = PC2)) +
  geom_point(aes(color = muestra),
             size = 3) +
  #geom_text(
    #label=c(group), 
    #stat="identity",
    #position="identity",
    #check_overlap = TRUE,
    #size = 2) +
  
  scale_color_manual(values = ggcolors) + ## ggcolors is the vector that defines the colours
  labs(x = paste0("PC1 (",round(axis_pca[1]*100), "%)"),
       y = paste0("PC2 (",round(axis_pca[2]*100), "%)"),
       color = "") + # name of the legend
  guides(color = guide_legend(nrow = 2, # number of rows in the legend
                              bycol = TRUE, # ordered by columns
                              override.aes = list(size = 4))) + # size of the points in the legend
  
  ggtitle("PCA mouse basic atlas") +
  
  theme_bw() +
  theme(legend.position="bottom",
        legend.title = element_text(size = 0), # Size of the tittle of the legend
        legend.text = element_text(size = 8), # Font size of the text of the legend
        axis.text = element_text(size = 7), # Font size of the axis labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
  )

gg_pca # Para imprimirlo en pantalla


################################################################################
######################## CORRELATION MATRIX ####################################
################################################################################

# Create a mapping of samples by extracting the sample name before the underscore
sample_names <- sub("_[0-9]+(\\.[0-9]+)?$", "", colnames(counts))  

# Aggregate counts by sample (sum or mean)
merged_counts <- as.data.frame(t(sapply(unique(sample_names), function(sample) {
  rowMeans(counts[, sample_names == sample])
})))

# Normalize the counts (optional)
log_counts <- log2(merged_counts + 1)  # log1p avoids log(0) issues

# Compute Spearman correlation matrix
spearman_cor <- cor(t(log_counts), method = "spearman") 


# Save the correlation matrix
write.table(spearman_cor, "spearman_correlation_matrix.txt", sep = "\t", quote = FALSE)

# Visualize as a heatmap (optional)
spearman_hmap <- pheatmap(
  spearman_cor,
  fontsize_row = 8,  # Adjust row label size
  fontsize_col = 8,  # Adjust column label size
  angle_col = 45     # Rotate column labels to avoid cutting
)


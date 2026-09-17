# set a working directory - this is on my desktop PC
rm(list=ls()) #clean the space

###############################################################################
#
# Option 1. I will use DEP package which includes impute() function for dealing
# with missing values
###############################################################################
# load package(s)
library("DEP")
library(dplyr)
library(SummarizedExperiment)
library(preprocessCore)
library(RColorBrewer)
library(matrixStats)

# Load my data using function read csv
LTEE <- read.csv(file="Data_for_R.csv", header = TRUE) # do not use row.names=1 because the code below wont work well
LTEE

# Dimensions of data
dim(LTEE)
# Colnames of data
colnames(LTEE)
# Are there any duplicated gene names?
LTEE$Genes %>% duplicated() %>% any()

# We do not have duplicated gene names, however we can use make_unique function
# as in tutorial to add columns name and ID which are required later
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
LTEE_unique <- make_unique(LTEE, "Genes", "Protein.Names", delim = ";")

# Are there any duplicated names?
LTEE$name %>% duplicated() %>% any()

# Generate summarized experiment object
# Generate a SummarizedExperiment object using an experimental design
intensity_columns <- 4:48 # get column numbers
intensity_columns
LTEE_design <- read.csv(file="Experiment_metadata.csv", header = TRUE) 
LTEE_se <- make_se(LTEE_unique, intensity_columns, LTEE_design)
LTEE_se
mat_LTEE_se  <- assay(LTEE_se)

# Filter on missing values
# Plot a barplot of the protein identification overlap between samples
plot_frequency(LTEE_se)

# Filter for proteins that are identified in all replicates of at least one condition
LTEE_filt <- filter_missval(LTEE_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
LTEE_filt2 <- filter_missval(LTEE_se, thr = 1)

# Even less stringent
# Filter for proteins that are identified in at least one replicate of at least one condition
LTEE_filt3 <- filter_missval(LTEE_se, thr = 2)
dim(LTEE_filt3)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(LTEE_filt3)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(LTEE_filt3)

# NORMALIZE THE DATA *THINK ABOUT USING TMM DR GRIBSKOV SUGGESTION
# Normalization using VSN from DEP package
LTEE_norm <- normalize_vsn(LTEE_filt3)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(LTEE_filt3, LTEE_norm)

# Plot a heatmap of proteins with missing values
plot_missval(LTEE_filt3)
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(LTEE_norm)

# FIRST APPROACH, TRY BASIC IMPUTATION METHODS AND CHECK DISTRIBUTION
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(LTEE_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
LTEE_imp <- impute(LTEE_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
LTEE_imp_man <- impute(LTEE_norm, fun = "man", shift = 2.0, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
LTEE_imp_knn <- impute(LTEE_norm, fun = "knn", rowmax = 0.9)

# Impute missing data using the QRILC approach
LTEE_imp_QRILC <- impute(LTEE_norm, fun = "QRILC")

# Plot intensity distributions before and after imputation
plot_imputation(LTEE_norm, LTEE_imp)

###############################################################################
# Check if imputation really adjusted the values as I wanted

library(ggplot2)

# Extract the assay matrices before and after imputation
mat_norm <- assay(LTEE_norm)        # before imputation
mat_imp  <- assay(LTEE_imp_man)     # after imputation with man
mat_imp_minprob <- assay(LTEE_imp)  # after imputation with MinProb
mat_quantile <- assay(LTEE_imp_QRILC)
# Build a tidy data frame
df_obs  <- data.frame(value = as.vector(mat_norm[!is.na(mat_norm)]),  type = "Observed")
df_imp  <- data.frame(value = as.vector(mat_quantile[is.na(mat_norm)]),    type = "Imputed")
df_all  <- rbind(df_obs, df_imp)

ggplot(df_all, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values  = c("Observed" = "#457B9D", "Imputed" = "#E63946")) +
  scale_color_manual(values = c("Observed" = "#457B9D", "Imputed" = "#E63946")) +
  labs(title = "Observed vs. imputed intensity distributions",
       x = "Protein intensity", y = "Density") +
  theme_bw(base_size = 13)

##############################################################################

# SECOND APPROACH, USE A MIXED IMPUTATION METHOD AND CHECK DISTRIBUTION
# First, we have to define a logical vector defining the rows that are to be imputed with the
# missing at random (MAR) method. 
# we consider a protein to have missing values not at random (MNAR) 
# if it has missing values in all replicates of at least one condition.

# Extract protein names with missing values 
# in all replicates of at least one condition
LTEE_proteins_MNAR <- get_df_long(LTEE_filt3) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
LTEE_MNAR <- names(LTEE_filt3) %in% LTEE_proteins_MNAR

# Perform a mixed imputation
LTEE_mixed_imputation <- impute(
  LTEE_filt3, 
  fun = "mixed",
  randna = !LTEE_MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "MinProb") # imputation function for MNAR

# Plot intensity distributions before and after imputation
plot_imputation(LTEE_filt3, LTEE_mixed_imputation)
# Boxplots after mixed imputation
plot_normalization(LTEE_filt3, LTEE_mixed_imputation)
# DECIDE IF IMPUTATION IS WORTH IT OR NOT, CHECK HOW TO DECIDE. MAYBE TRY THAT IMPUTATION BY COLUMNS?
# Make more boxplots to explore what the quantile normalization did.
# TRY ONLY QUANTILE NORMALIZATION FOR NORMALIZING THE DATA, NOT FOR IMPUTATION, THIS IS DIFFERENT.
# AND READ WHAT IT IS DOING BEHING SCENES. SEE IF THTA IMPROVES THE DISTIRBUTIONS AND BOXPLOTS.

###############################################################################
# Check if mixed imputation really adjusted the values as I wanted

# Extract the assay matrices before and after imputation
mat_norm <- assay(LTEE_norm)        # before imputation
mat_imp_mixed  <- assay(LTEE_mixed_imputation)     # after imputation

# Build a tidy data frame
df_obs  <- data.frame(value = as.vector(mat_norm[!is.na(mat_norm)]),  type = "Observed")
df_imp  <- data.frame(value = as.vector(mat_imp_mixed[is.na(mat_norm)]),    type = "Imputed")
df_all  <- rbind(df_obs, df_imp)

ggplot(df_all, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values  = c("Observed" = "#457B9D", "Imputed" = "#E63946")) +
  scale_color_manual(values = c("Observed" = "#457B9D", "Imputed" = "#E63946")) +
  labs(title = "Observed vs. imputed intensity distributions",
       x = "Protein intensity", y = "Density") +
  theme_bw(base_size = 13)
#################################################################################

# DIFFERENTIAL ANALYSIS FOR THE THREE RECOMMENDED METHOS FOR MNAR VALUES (MAN, MINPROB, QRILC)
# Differential enrichment analysis  based on linear models and emphIrical Bayes statistics
# Option 1. Test every sample versus control
LTEE_data_diff <- test_diff(LTEE_imp, type = "control", control = "REL606") # Use ancentral strain REL606 as baseline
LTEE_data_diff_man <- test_diff(LTEE_imp_man, type = "control", control = "REL606") 
LTEE_data_diff_QRILC <- test_diff(LTEE_imp_QRILC, type = "control", control = "REL606") 
LTEE_data_diff_mixed <- test_diff(LTEE_mixed_imputation, type = "control", control = "REL606")

# Option 2. Make multiple contrasts of the "transition" LTEE clones to see what is changing during refinement
# Test all possible comparisons of samples
LTEE_diff_all_contrasts <- test_diff(LTEE_imp_man, type = "all")
LTEE_diff_manual_contrasts <- test_diff(LTEE_imp_man, type = "manual",
                                        test = c("CZB152_vs_ZDB564", "ZDB172_vs_ZDB564",
                                                 "ZDB143_vs_ZDB564", "CZB154_vs_ZDB564",
                                                 "ZDB96_vs_ZDB564", "ZDB107_vs_ZDB564", "REL10979_vs_ZDB564"))
                                     #   test = c("CZB152_vs_ZDB564", "ZDB564_vs_ZDB30", "ZDB172_vs_ZDB564",
                                      #           "ZDB143_vs_ZDB564", "CZB154_vs_ZDB564"))


# Denote significant proteins based on user defined cutoffs
LTEE_dep <- add_rejections(LTEE_data_diff, alpha = 0.05, lfc = log2(1.5))
LTEE_dep_man <- add_rejections(LTEE_data_diff_man, alpha = 0.05, lfc = log2(1.5))
LTEE_dep_QRILC <- add_rejections(LTEE_data_diff_QRILC, alpha = 0.05, lfc = log2(1.5))
LTEE_dep_mixed <- add_rejections(LTEE_data_diff_mixed, alpha = 0.05, lfc = log2(1.5))
  
LTEE_dep_all_man <- add_rejections(LTEE_diff_all_contrasts, alpha = 0.05, lfc = log2(2.0))
LTEE_dep_manual_man <- add_rejections(LTEE_diff_manual_contrasts, alpha = 0.05, lfc = log2(1.5))
  
# ADDITIONAL, NOT IN TUTORIAL. CHECK TOTAL NUMBER OF SIGNIFICANT PROTEINS PER CONTRAST
# Breakdown per contrast
rowData(LTEE_dep) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()

rowData(LTEE_dep_man) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()

rowData(LTEE_dep_QRILC) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()

rowData(LTEE_dep_mixed) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()

rowData(LTEE_dep_all_man) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()

rowData(LTEE_dep_manual_man) %>% 
  as.data.frame() %>% 
  dplyr::select(ends_with("_significant")) %>% 
  colSums()


# Total number of significant proteins (significant in at least one contrast)
sum(rowData(LTEE_dep)$significant)
sum(rowData(LTEE_dep_man)$significant)
sum(rowData(LTEE_dep_QRILC)$significant)
sum(rowData(LTEE_dep_mixed)$significant)

sum(rowData(LTEE_dep_all_man)$significant)
sum(rowData(LTEE_dep_manual_man)$significant)

############## END OF ADDITIONAL NOT IN TUTORIAL##########################

# VISUALIZATION OF RESULTS TO ASSIST IN DETERMINING OPTIMAT CUTOFFS
# Plot the first and second principal components
plot_pca(LTEE_dep, x = 1, y = 2, n = 500, point_size = 4)
plot_pca(LTEE_dep_man, x = 1, y = 2, n = 1000, point_size = 4)
plot_pca(LTEE_dep_QRILC, x = 1, y = 2, n = 500, point_size = 4)
plot_pca(LTEE_dep_mixed, x = 1, y = 2, n = 500, point_size = 4)

plot_pca(LTEE_dep_all_man, x = 1, y = 2, n = 500, point_size = 4)
plot_pca(LTEE_dep_manual_man, x = 1, y = 2, n = 500, point_size = 4)

# Correlation matrix
# Correlation matrix can be plotted as a heatmap to visualize Pearson correlations between samples
# Plot the Pearson correlation matrix
plot_cor(LTEE_dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(LTEE_dep_man, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(LTEE_dep_QRILC, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

plot_cor(LTEE_dep_all_man, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(LTEE_dep_manual_man, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
# This one didn't work because of internal issue with number of colors in the 
# palette of heatmap function
plot_heatmap(LTEE_dep_manual_man, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 12, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

################################################################################
# Help from Claude for this. BUT THIS HEATMAP IS TOO LARGE, VERY CONFUSING
# TRY TO REPLICATE WHAT plot_heatmap DOES SINCE I COULDN'T BYPASS THE ERROR, IT IS SOMETHING
# THAT DEP PACKAGE DOES INTERNALLY
library(pheatmap)
# 1. Extract and scale the matrix (centered = subtract row mean) This is using my significant proteins only
mat_dep_man  <- as.matrix(assay(LTEE_dep_man))
mat_dep_man  <- mat_dep_man - rowMeans(mat_dep_man, na.rm = TRUE)   # centered scaling

# 2. Run kmeans to cluster proteins into groups
set.seed(123)
km   <- kmeans(mat_dep_man, centers = 6, nstart = 25)

# 3. Sort rows by kmeans cluster so clusters appear as blocks
mat_sorted     <- mat_dep_man[order(km$cluster), ]
cluster_labels <- sort(km$cluster)

# 4. Build annotation for columns (conditions) and rows (kmeans clusters)
col_data_man  <- as.data.frame(colData(LTEE_dep_man))
condition <- col_data_man$condition   # adjust if your column name differs

col_annot <- data.frame(
  Condition = col_data_man$condition,
  Replicate = col_data_man$replicate,
  row.names = colnames(mat_dep_man)
)

row_annot <- data.frame(
  Cluster = factor(cluster_labels),
  row.names = rownames(mat_sorted)
)


# 5. Define colors
condition_levels <- unique(condition)
sidebar_colors   <- setNames(
  c("#E63946", "#457B9D", "yellow4", "chartreuse3", "darkgreen",
    "sienna3", "orangered1", "orange2", "hotpink2", "mediumvioletred", 
    "lightpink2", "royalblue3","deepskyblue3", "turquoise3", "darkblue"),   # adjust to match your my_colors
  condition_levels
)

annot_colors <- list(
  Condition = sidebar_colors,           # your existing my_colors
  Cluster   = setNames(
    RColorBrewer::brewer.pal(6, "Set2"),
    as.character(1:6)
  )
)


# 6. Plot
# Define colors for the heatmap body
hmcol <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

svg(file.path(getwd(), "heatmap_kmeans.svg"), width = 12, height = 16)

pheatmap(mat_sorted,
         color             = hmcol,
         cluster_rows      = FALSE,     # already sorted by kmeans
         cluster_cols      = TRUE,      # cluster samples
         treeheight_row    = 0,
         treeheight_col    = 30,
         annotation_col    = col_annot,
         annotation_row    = row_annot,
         annotation_colors = annot_colors,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         fontsize_col      = 10,
         main              = "Centered Protein Intensities — k=6 clusters"
)

dev.off()

# END of help from Claude
# End of trying to create a heatmap for all significant proteins
############################################################################
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(LTEE_dep_man, type = "contrast", kmeans = TRUE, 
             k = 10, col_limit = 10, show_row_names = FALSE)

plot_heatmap(LTEE_dep_manual_man, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
#############################################################################
# SOMETHING ADDITONAL. NOT IN TUTORIAL
# Get the assay data for contrasts
contrast_data <- assay(LTEE_dep_man)

contrast_manual_data <- assay(LTEE_dep_manual_man)

# Calculate variance across contrasts for each protein
protein_var <- apply(contrast_data, 1, var, na.rm = TRUE)

protein_manual_var <- apply(contrast_manual_data, 1, var, na.rm = TRUE)

# Keep top 100 most variable proteins
top100 <- names(sort(protein_var, decreasing = TRUE)[1:100])

top100_manual <- names(sort(protein_manual_var, decreasing = TRUE)[1:70])

# Subset the SummarizedExperiment object
LTEE_dep_man_top100 <- LTEE_dep_man[top100, ]

LTEE_dep_manual_man_top100 <- LTEE_dep_manual_man[top100_manual, ]

# Plot. THis is mix from tutorial and my modifications to plot only some proteins
plot_heatmap(LTEE_dep_man_top100, type = "contrast", kmeans = TRUE,
             k = 10, col_limit = 10, show_row_names = FALSE)

plot_heatmap(LTEE_dep_manual_man_top100, type = "contrast", kmeans = TRUE,
             k = 10, col_limit = 10, show_row_names = TRUE)
# TRy to modify font size
plot_heatmap(
  LTEE_dep_manual_man_top100,
  type = "contrast",
  kmeans = TRUE,
  k = 10,
  col_limit = 10,
  show_row_names = TRUE,
  row_font_size = 16,
  col_font_size = 16
)
# Save this version
# Save to file
# ----------------------------
library(svglite)
svglite(
  "./Costs_tests_poster/heatmap_properon.svg",
  width = 10,
  height = 8   # increase from 6 to 14
)

plot_heatmap(
  LTEE_dep_manual_man_top100,
  type = "contrast",
  kmeans = TRUE,
  k = 10,
  col_limit = 10,
  show_row_names = TRUE,
  row_font_size = 16,
  col_font_size = 16
)

dev.off()


#Show name of genes instead
plot_heatmap(LTEE_dep_man_top100, type = "contrast", kmeans = FALSE,
             col_limit = 10, show_row_names = TRUE,
             show_row_dend = FALSE)  # optional: hide dendrogram for cleaner look

plot_heatmap(LTEE_dep_manual_man_top100, type = "contrast", kmeans = FALSE,
             col_limit = 10, show_row_names = TRUE,
             show_row_dend = FALSE)  # optional: hide dendrogram for cleaner look

# Filter to significant proteins first, then top 100 variable
sig_proteins <- rowData(LTEE_dep_man)$significant  # logical vector
LTEE_dep_man_sig <- LTEE_dep_man[which(sig_proteins), ]

protein_var <- apply(assay(LTEE_dep_man_sig), 1, var, na.rm = TRUE)
top100 <- names(sort(protein_var, decreasing = TRUE)[1:100])
LTEE_dep_man_top100 <- LTEE_dep_man_sig[top100, ]

sig_top100_LTEE_man <- plot_heatmap(LTEE_dep_man_top100, type = "contrast", kmeans = FALSE,
            col_limit = 10, show_row_names = TRUE,
            show_row_dend = FALSE)

# For manual contrasts
sig_manual_prots <- rowData(LTEE_dep_manual_man)$significant
LTEE_dep_manualcontrasts_sig <- LTEE_dep_manual_man[which(sig_manual_prots), ]

protein_var_manual <- apply(assay(LTEE_dep_manualcontrasts_sig), 1, var, na.rm = TRUE)
top100_manual_sig <- names(sort(protein_var_manual, decreasing = TRUE)[1:100])
LTEE_dep_manual_sig_top100 <- LTEE_dep_manualcontrasts_sig[top100_manual_sig, ]

sig_top100_LTEE_manualcontrasts <- plot_heatmap(LTEE_dep_manual_sig_top100, type = "contrast", kmeans = FALSE,
                                    col_limit = 10, show_row_names = TRUE,
                                    show_row_dend = FALSE)

sig_top100_LTEE_manualcontrasts
## END OF SOMETHING ADDITIONAL
############################################################################

# ALTERNATIVE OPTION FOR HEATMAP. USE THE GENERATE OBJECT FROM DEP TO CREATE A
# HEATMAP USING HEATMAP BASIC PACKAGE
library(RColorBrewer)
library(matrixStats)
library(SummarizedExperiment)

# 1. Extract the assay matrix from the imputed SE object
LTEE_mat_man <- as.matrix(assay(LTEE_dep_man))

# FOr manual contrasts
LTEE_mat_manual_man <- as.matrix(assay(LTEE_dep_manual_man))

# 2. Subset to top 150 most variable proteins
rv     <- rowVars(LTEE_mat_man, na.rm = TRUE)
select <- order(rv, decreasing = TRUE)[seq_len(min(600, length(rv)))]
LTEE_mat_man_subset <- LTEE_mat_man[select, ]

# 3. Extract condition info for column annotation
col_data_man  <- as.data.frame(colData(LTEE_dep_man))
condition <- col_data_man$condition   # adjust if your column name differs

# 4. Define colors for the heatmap body
hmcol <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

# 5. Define colors for the condition sidebar
condition_levels <- unique(condition)
sidebar_colors   <- setNames(
  c("#E63946", "#457B9D", "yellow4", "chartreuse3", "darkgreen",
    "sienna3", "orangered1", "orange2", "hotpink2", "mediumvioletred", 
    "lightpink2", "royalblue3","deepskyblue3", "turquoise3", "darkblue"),   # adjust to match your my_colors
  condition_levels
)
col_sidebar <- sidebar_colors[condition]

# 6. Plot
# If saving to a file
svg("heatmap_top100.svg", width = 40, height = 30)

# Pre-compute the row clustering
row_clust <- hclust(dist(LTEE_mat_man_subset))


heatmap(LTEE_mat_man_subset,
        scale   = "row",         # z-score per protein across samples
        col     = hmcol,
        ColSideColors = col_sidebar,   # condition color bar on top
        Colv    = as.dendrogram(hclust(dist(t(LTEE_mat_man_subset)))),  # cluster samples
        Rowv    = as.dendrogram(row_clust),  # cluster proteins
        labCol  = colnames(LTEE_mat_man_subset),
        labRow  = rownames(LTEE_mat_man_subset),            # change to NA to hide protein names (too many)
        cexRow        = 1.0,    # shrink row labels — reduce further if still clipping
        cexCol        = 1.0,    # shrink column labels if needed
        treeheight    = 0,     # hide dendogram
        margins = c(10, 10),       # bottom and right margins
        main    = "Top 100 Variable Proteins",
        xlab    = "Sample",
        ylab    = "Protein"
)

# 7. Add a legend for the condition sidebar
# Force legend to be drawable anywhere on the device
par(xpd = NA)
legend(x = 0.85, y = 0.3, #"bottomright",  # if you want it topright and not with specific coordinates
       #inset = 0.05,
       legend = condition_levels,
       fill   = sidebar_colors,
       border = NA,
       bty    = "n",
       cex    = 0.8
)
dev.off()
dev.list()
dev.off(3)

# 6.1. Plot the whole dataset of proteins to check overall distribution and clustering
# Heatmap using the full matrix
# I think this is not a good idea because no the whole proteome changes, it is too messy
# We don't really expect the whole proteome to change, just a portion. That's what I think
# So now I'll try with the top 600 because overall there are around 550 significant proteins
svg(file.path(getwd(), "heatmap_top600_proteins.svg"), width = 12, height = 16)

par(cex = 1.2)
heatmap(LTEE_mat_man_subset,
        scale         = "row",
        col           = hmcol,
        ColSideColors = col_sidebar,
        Colv          = as.dendrogram(hclust(dist(t(LTEE_mat_man_subset)))),
        Rowv          = TRUE,              # no row clustering, too many proteins
        labCol        = colnames(LTEE_mat_man_subset),
        labRow        = NA,              # hide protein names — too many to show
        cexCol        = 0.8,
        margins       = c(8, 4),
        main          = "All Proteins",
        xlab          = "Sample",
        ylab          = "Protein"
)

par(xpd = NA)
legend(x      = 0.8,
       y      = 0.2,
       legend = condition_levels,
       fill   = sidebar_colors,
       border = NA,
       bty    = "n",
       cex    = 0.8
)

dev.off()

###############################################################################
library(reshape2)
library(forcats)
# BACK TO FOLLOWING TUTORIAL STEPS BELOW
# Volcano plots
# Plot a volcano plot for the contrast "REL606 vs ZDB564""
plot_volcano(LTEE_dep_man, contrast = "CZB152_vs_REL606", label_size = 2, add_names = TRUE)

plot_volcano(LTEE_dep_all_man, contrast = "ZDB172_vs_ZDB564", label_size = 4, add_names = TRUE)

par(mfrow=c(3,2))
plot_volcano(LTEE_dep_manual_man, contrast = "ZDB143_vs_ZDB564", label_size = 4, add_names = TRUE)
plot_volcano(LTEE_dep_manual_man, contrast = "ZDB172_vs_ZDB564", label_size = 4, add_names = TRUE)
plot_volcano(LTEE_dep_manual_man, contrast = "CZB154_vs_ZDB564", label_size = 4, add_names = TRUE)
plot_volcano(LTEE_dep_manual_man, contrast = "ZDB96_vs_ZDB564", label_size = 4, add_names = TRUE)
plot_volcano(LTEE_dep_manual_man, contrast = "ZDB107_vs_ZDB564", label_size = 4, add_names = TRUE)
plot_volcano(LTEE_dep_manual_man, contrast = "CZB152_vs_ZDB564", label_size = 4, add_names = TRUE)

# BARPLOTS of proteins of interest
plot_single(LTEE_dep_all_man, proteins = c("sucA", "sucB", "sucC"))
plot_single(LTEE_dep_manual_man, proteins = c("gltA", "icd", "sdhA"), type = "centered")

plot_single(LTEE_dep_manual_man, proteins = c("gltA", "icd", "sdhA"))
plot_single(LTEE_dep_manual_man, proteins = "gltA")


# TRy to do my own customized heatmap
# Extract all differential expression results
results <- as.data.frame(rowData(LTEE_dep_manual_man))
results

class(results)
str(results)
names(results)
results$name[1:5]
# See the available columns
colnames(results)

library(dplyr)
library(tidyr)
library(patchwork)
library(forcats)

# Check what was the error, it was not using function from dplyr but from a bioconductir package inside DEP
tmp <- results %>%
  select(name, ends_with("_diff"))

head(tmp)

tmp2 <- tmp %>%
  pivot_longer(
    cols = ends_with("_diff"),
    names_to = "comparison",
    values_to = "log2FC"
  )

head(tmp2)

tmp3 <- tmp2 %>%
  dplyr::rename(gene = name)
head(tmp3)

# TO check if two packages are conflicting
find("rename")

# Okay now full code after we check what was the error
log2fc_data <- results %>%
  dplyr::select(gene = name, dplyr::ends_with("_diff")) %>%
  tidyr::pivot_longer(
    cols = dplyr::ends_with("_diff"),
    names_to = "comparison",
    values_to = "log2FC"
  ) %>%
  dplyr::mutate(
    comparison = sub("_vs_ZDB564_diff$", "", comparison)
  )

log2fc_data

sort(unique(log2fc_data$gene))

# Order comparisons by generation, genes in operon order
generation_order <- c("ZDB143", "ZDB172", "CZB152", "CZB154", "ZDB96", "ZDB107", "REL10979")
#prp_data <- prp_data %>%
#  mutate(
#    comparison = factor(comparison, levels = generation_order),
#    gene = factor(gene, levels = c("prpB", "prpC", "prpD", "prpE"))
#  )

prp_data <- log2fc_data %>%
  dplyr::filter(gene %in% c("prpB", "prpC", "prpD", "prpE")) %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = generation_order),
    gene = factor(gene, levels = c("prpB", "prpC", "prpD", "prpE"))
  )

# Main heatmap
main_heatmap <- ggplot(prp_data, aes(x = comparison, y = fct_rev(gene), fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", log2FC)), size = 3.3,
            color = ifelse(abs(prp_data$log2FC) > 3, "white", "black")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = expression(log[2]*FC)) +
  labs(x = NULL, y = NULL, title = "prp operon: consistently upregulated except in ZDB107 and REL10979") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),   # hide here — genotype strip below will carry the labels
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )
main_heatmap
# Genotype annotation strip
genotype_data <- tibble(
  comparison = factor(generation_order, levels = generation_order),
  genotype = c("gltA1 only", "gltA1 only", "gltA1 only", "gltA1 only", "gltA2 (V152A)", "gltA2 (A124T)*", "gltA2 (A124T)*")
)

genotype_strip <- ggplot(genotype_data, aes(x = comparison, y = 1, fill = genotype)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comparison), vjust = -1.2, size = 3.2, color = "black") +
  scale_fill_manual(values = c(
    "gltA1 only"      = "grey85",
    "gltA2 (V152A)"   = "#F4A582",
    "gltA2 (A124T)*"  = "#B2182B"
  ), name = "gltA genotype") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom",
        plot.margin = margin(t = 10, b = 0))

# Stack them — genotype strip is the "annotation track" sitting right below the heatmap x-axis
main_heatmap / genotype_strip + plot_layout(heights = c(4, 1))


# another thing to try. TCA cycle genes
# Define genes in actual TCA-cycle/glyoxylate-shunt/methylcitrate-cycle order
pathway_order <- c(
  "gltA",  # citrate synthase (entry point)
  "acnB", "acnA",  # aconitase
  "icd",   # isocitrate dehydrogenase
  "sucA", "sucB",  # 2-oxoglutarate dehydrogenase
  "sucC", "sucD",  # succinyl-CoA synthetase
  "sdhA", "sdhB", "sdhC", "sdhD",  # succinate dehydrogenase
  "fumA", "fumB", "fumC",  # fumarase
  "mdh",   # malate dehydrogenase
  "aceA", "aceB",  # glyoxylate shunt (isocitrate lyase, malate synthase)
  "prpB", "prpC", "prpD", "prpE"  # methylcitrate cycle (propionate overflow branch)
)

# Example structure — replace with your actual TCA/glyoxylate/prp log2FC values
tca_data <- log2fc_data %>%   # columns: comparison, gene, log2FC
  dplyr::filter(gene %in% pathway_order) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(pathway_order)),  # rev() so gltA sits at top
    comparison = factor(comparison, levels = generation_order)
  )

# Add a pathway-module label to visually separate TCA / glyoxylate / methylcitrate blocks
tca_data <- tca_data %>%
  mutate(module = case_when(
    gene %in% c("gltA","acnB","acnA","icd","sucA","sucB","sucC","sucD",
                "sdhA","sdhB","sdhC","sdhD","fumA","fumB","fumC","mdh") ~ "TCA cycle",
    gene %in% c("aceA","aceB") ~ "Glyoxylate shunt",
    gene %in% c("prpB","prpC","prpD","prpE") ~ "Methylcitrate cycle"
  ))

ggplot(tca_data, aes(x = comparison, y = gene, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", log2FC)), size = 3,
            color = ifelse(abs(tca_data$log2FC) > 3, "white", "black")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = expression(log[2]*FC)) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Clone (vs. ZDB564)", y = NULL,
       title = "Central carbon metabolism proteomic response across Cit++ clones") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )



# Now using the normalized values instead
library(tibble)
expr <- assay(LTEE_imp_man)
sig_proteins <- results %>%
  dplyr::filter(significant) %>%
  dplyr::pull(name)

expr_sig <- expr[rownames(expr) %in% sig_proteins, ]

# Show values for all replicates
expr_long <- expr_sig %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    -gene,
    names_to = "comparison",
    values_to = "expression"
  )

expr_long

# Or collapse replicates
expr_long <- expr_long %>%
  mutate(
    comparison = sub("_\\d+$", "", comparison)
  ) %>%
  group_by(gene, comparison) %>%
  summarize(
    expression = mean(expression),
    .groups = "drop"
  )

expr_long

# Order the data
# Clones to show, ZDB564 first as the reference/baseline
clones_to_show <- c("ZDB564", "ZDB143", "ZDB172", "CZB152", "CZB154", "ZDB96", "ZDB107")

prp_normalized <- expr_long %>%                 # the tibble you described: gene, comparison, expression
  dplyr::filter(gene %in% c("prpB", "prpC", "prpD", "prpE"),
         comparison %in% clones_to_show) %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = clones_to_show),
    gene = factor(gene, levels = c("prpB", "prpC", "prpD", "prpE"))
  )

main_heatmap <- ggplot(prp_normalized, aes(x = comparison, y = fct_rev(gene), fill = expression)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", expression),
                color = expression > (min(expression) + 0.6 * diff(range(expression)))),
            size = 3.3, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Normalized\nexpression") +
  labs(x = NULL, y = NULL, title = "prp operon expression: elevated except in ZDB107") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )
main_heatmap

genotype_data_normalized <- tibble(
  comparison = factor(clones_to_show, levels = clones_to_show),
  genotype = c("ancestral (ZDB564)", "gltA1 only", "gltA1 only", "gltA1 only",
               "gltA1 only", "gltA2 (V152A)", "gltA2 (A124T)*")
)

genotype_strip_normalized <- ggplot(genotype_data_normalized, aes(x = comparison, y = 1, fill = genotype)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comparison), vjust = -1.2, size = 3.2, color = "black") +
  scale_fill_manual(values = c(
    "ancestral (ZDB564)" = "grey95",
    "gltA1 only"          = "grey75",
    "gltA2 (V152A)"       = "#F4A582",
    "gltA2 (A124T)*"      = "#B2182B"
  ), name = "gltA genotype") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(t = 10, b = 0))

main_heatmap / genotype_strip_normalized + plot_layout(heights = c(4, 1))

# For TCA cycle genes
tca_data_normalized <- expr_long %>%
  dplyr::filter(gene %in% pathway_order,
         comparison %in% clones_to_show) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(pathway_order)),
    comparison = factor(comparison, levels = clones_to_show),
    module = case_when(
      gene %in% c("gltA","acnB","acnA","icd","sucA","sucB","sucC","sucD",
                  "sdhA","sdhB","sdhC","sdhD","fumA","fumB","fumC","mdh") ~ "TCA cycle",
      gene %in% c("aceA","aceB") ~ "Glyoxylate shunt",
      gene %in% c("prpB","prpC","prpD","prpE") ~ "Methylcitrate cycle"
    )
  )

ggplot(tca_data_normalized, aes(x = comparison, y = gene, fill = expression)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", expression),
                color = expression > (min(expression) + 0.6 * diff(range(expression)))),
            size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Normalized\nexpression") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Clone", y = NULL,
       title = "Central carbon metabolism: normalized protein expression across Cit++ clones") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

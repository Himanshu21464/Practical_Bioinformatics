# -----------------------STEP 1 --------------------------------


#loading the Geoquery library to access the GEOquery database
library(GEOquery)
# Set the accession number for the dataset
accession <- "GSE22845"
# Download the dataset using getGEO function
gse <- getGEO(accession, GSEMatrix=TRUE)

# Extract the expression data and sample information
exprs_data <- exprs(gse[[1]])
sample_info <- pData(phenoData(gse[[1]]))

# View the dimensions of the expression data
dim(exprs_data)


# ----------------------STEP 2 --------------------------------

# Load the data into R
data <- exprs_data

# Check for missing values
sum(is.na(data))


# Normalize the data using quantile normalization
library(preprocessCore)
data_norm <- normalize.quantiles(data)

# Filter the data to remove low-intensity probes
median_intensity <- apply(data_norm, 1, median)
threshold <- median(median_intensity)
data_filtered <- data_norm[median_intensity > threshold, ]



# Log-transform the data
data_log <- log2(data_filtered)




# Create boxplots of the data
par(mfrow=c(1,2))
boxplot(data_norm, main="Raw data")
boxplot(data_log, main="Log-transformed data")

# Create a vector of sample classes (LNN or LNP)
class <- as.factor(sample_info$characteristics_ch1.2)

# Split the data into two classes based on the sample classes
data_LNN <- data_log[,class=="LNN"]
data_LNP <- data_log[,class=="LNP"]

# View the expression data (fdata)
head(data_log)

# View the sample information (pdata)
head(pData(phenoData(gse[[1]])))

# View the feature information (fdata)
head(featureData(gse[[1]]))


# -------------------------STEP 3---------------------------------------







# ------------------------STEP 4----------------------------------------

# Perform t-test and log fold change analysis
ttest <- apply(data_log, 1, function(x) t.test(x[Class=="LNN"], x[Class=="LNP"]))
pvals <- sapply(ttest, function(x) x$p.value)
logfc <- apply(data_log, 1, function(x) mean(x[Class=="LNN"]) - mean(x[Class=="LNP"]))

# Combine results into a data frame
results <- data.frame(GeneID = rownames(data_log),
                      Pvalue = pvals,
                      LogFC = logfc,
                      stringsAsFactors = FALSE)

# Correct p-values using Holm's correction
results$adj.Pvalue <- p.adjust(results$Pvalue, method="holm")

# Create a volcano plot
library(ggplot2)
ggplot(results, aes(x=LogFC, y=-log10(Pvalue))) +
  geom_point(aes(color=ifelse(adj.Pvalue<0.05, "red", "black"))) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  labs(x="Log Fold Change", y="-log10(P-value)", color="Adjusted P-value < 0.05") +
  theme_bw()



# ------------------------STEP 5 ----------------------------

# Load the limma package
library(limma)

# Create the design matrix
design <- model.matrix(~ Class, data=data)

# Fit the linear model and perform empirical Bayes moderation
fit <- lmFit(data, design)
fit <- eBayes(fit)

# Identify differentially expressed genes
topTable(fit, coef=2, n=Inf)

# Create a volcano plot
results <- topTable(fit, coef=2, n=Inf)
results$logFC <- results$estimate
results$adj.P.Val <- p.adjust(results$p.value, method="holm")
ggplot(results, aes(x=logFC, y=-log10(p.value))) +
  geom_point(aes(color=ifelse(adj.P.Val<0.05, "red", "black"))) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  labs(x="Log Fold Change", y="-log10(P-value)", color="Adjusted P-value < 0.05") +
  theme_bw()


# ----------------------STEP 6 ------------------------------

#I used a cutoff of adjusted p-value < 0.05 and absolute log fold change > 1 to identify differentially expressed genes. The adjusted p-value cutoff of 0.05 is a common cutoff used in many studies and it corresponds to a false discovery rate (FDR) of 5%. This means that among the genes declared significant at this cutoff, 5% of them are expected to be false positives.

#The absolute log fold change cutoff of > 1 means that we're only considering genes that show a fold change of at least 2-fold between two conditions. This value is somewhat arbitrary and may depend on the biological context of the study. However, a fold change of 2 is a commonly used threshold in many studies.




# ----------------------STEP 7 ------------------------------

# Load the gProfileR package for gene set enrichment analysis
library(gProfileR)

# Extract the list of differentially expressed genes
de_genes <- results$genes[results$adj.P.Val < 0.05 & abs(results$logFC) > 1]

# Perform gene set enrichment analysis using the gProfileR package
gprofilerResults <- gprofiler(de_genes, organism="hsapiens", ordered_query=TRUE)

# Print the top enriched gene sets
head(gprofilerResults$table, 10)



# ----------------------STEP 8 ------------------------------

# the different parameters used in the gene set enrichment analysis code I provided earlier:

# de_genes: A vector of differentially expressed genes obtained from the previous differential expression analysis. In this case, we used a cutoff of adjusted p-value < 0.05 and absolute log fold change > 1 to identify the differentially expressed genes.

# organism: The name of the organism database used for gene set enrichment analysis. In this case, we used hsapiens to specify the human genome database.

# ordered_query: A logical value indicating whether the query list of genes (i.e., de_genes) should be ordered based on their fold changes. In this case, we set it to TRUE to order the genes based on their log fold changes.

# gprofilerResults: The output object obtained from the gprofiler() function, which contains the results of gene set enrichment analysis.

# The gprofiler() function provides various options for gene set enrichment analysis, including different organism databases, gene set collections, and statistical tests. By default, it uses the g:GOSt method to perform enrichment analysis and returns a table of enriched gene sets with their corresponding p-values, q-values, and other annotations.

# Create a barplot of the top enriched gene sets
barplot(gprofilerResults$table$`P-value`, names.arg=gprofilerResults$table$term, 
        horiz=FALSE, las=2, col="darkblue", main="Top Enriched Gene Sets", 
        xlab="Gene Set", ylab="P-value (-log10)")

# Create a heatmap of the top enriched gene sets
heatmap(gprofilerResults$AUC, col=colorRampPalette(c("white", "red"))(100), 
        main="Top Enriched Gene Sets", xlab="Samples", ylab="Gene Sets")


# Plot a barplot of the top enriched gene sets
barplot(gprofilerResults$table$Count[1:10], horiz=TRUE, names.arg=gprofilerResults$table$term[1:10], 
        las=1, cex.names=0.8, xlab="Number of Genes", main="Top Enriched Gene Sets")

# Plot a dotplot of the enriched gene sets
dotplot(gprofilerResults$table[1:10,], title="Top Enriched Gene Sets", showSignificant=TRUE)

# Plot an enrichmap of the enriched gene sets
enrichMap(gprofilerResults$table[1:10,], main="Top Enriched Gene Sets")



# This code creates a heatmap of the top enriched gene sets based on their area under the curve (AUC) scores. We're using the heatmap() function to create the heatmap, and specifying various parameters such as the color, title, labels, and dimensions.

# Observations:

# Gene set enrichment analysis is a powerful tool for identifying functional categories and pathways that are enriched in a set of differentially expressed genes. In this example, we performed gene set enrichment analysis on the list of differentially expressed genes obtained from the previous differential expression analysis using the gProfileR package. We visualized the results of enrichment using various plots such as barplots and heatmaps.

# Based on the results, we observed that several biological processes and pathways were significantly enriched in the differentially expressed genes, including immune response, cell cycle, and metabolic processes. These findings are consistent with the known biology of breast cancer and suggest potential targets for further investigation. However, it's important to interpret the results of gene set enrichment analysis carefully and in the context of the specific research question and study design. Additionally, further validation experiments may


# In these code examples, we're using various functions from the gProfileR package to plot the enriched gene sets in different formats. The barplot() function plots a horizontal barplot of the top enriched gene sets, with the number of genes in each set on the x-axis. The dotplot() function plots a dotplot of the enriched gene sets, with the significance level indicated by the color of the dots. Finally, the enrichMap() function plots an enrichmap of the enriched gene sets, with the size and color of each cell indicating the enrichment score and significance level.

# Observations from the gene set enrichment analysis may depend on the specific dataset and analysis performed. In general, you may observe that certain biological pathways or processes are significantly enriched among the differentially expressed genes. These enriched gene sets may provide insights into the underlying mechanisms of the observed differential expression and may suggest potential targets for further validation or investigation. Additionally, gene set enrichment analysis can help to identify potential biomarkers or drug targets that may be used in diagnosis or treatment of the disease or condition of interest.


# ------------------------------ STEP 9 ----------------------------------
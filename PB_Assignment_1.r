
######################################### STEP 1 ###########################################################


# loading the required library
library(GEOquery)

# downloading the microarray dataset from GEO  [Note that downloading the dataset requires good internet connection. Poor internet connection may result in "Time out reached failure"]
Dataset <- getGEO("GSE37250", GSEMatrix=TRUE)

# Extract the expression data and sample information
Expressed_Data <- exprs(Dataset[[1]])

write.csv(Dataset,"PB_Assignment_1.csv")
Sample_Information <- pData(phenoData(Dataset[[1]]))


# attributes 
attr= attributes(Dataset[[1]])

#pData
p_data=pData(Dataset[[1]])

#fData
f_data=fData(Dataset[[1]])




# show all the attributes in the dataset including pdata and fdata
View(attr)
View(p_data)
View(f_data)
######################################### STEP 2 ###########################################################

# Load the data into R
data <- Expressed_Data
data <- na.omit(data) # remove rows with missing values


# Normalization
library(preprocessCore)
data_norm <- normalize.quantiles(data) # quantile normalization


# Filter the data to remove low-intensity probes
median_intensity <- apply(data_norm, 1, median)
threshold <- median(median_intensity)
data_filtered <- data_norm[median_intensity > threshold, ]

#log transformation
data_log <- log2(data_filtered-min(data_filtered))

# Visualize the data
hist(data_log[,1]) # histogram of expression values 
boxplot(data_log) # boxplot of expression values for all genes  (Log transformed data)
boxplot(data_norm) # boxplot of expression values for all genes (Normalized data)



######################################### STEP 3 ###########################################################


# Log transformation of data helps us in normalizing and stabalizing the data as the raw data contain a wide range of values and it is easier to compare log transformed dataset. linearize the relationships between gene expression and biological effects,
# and improve the interpretability of changes in gene expression levels.


######################################### STEP 4 ###########################################################


#------------------------------------------T-test--------------------------------------------------




# creating two subset based on the HIV status of the patients
HIV_Positive<- subset(p_data, characteristics_ch1.1 == "hiv status: HIV positive")
HIV_Positive<- subset(HIV_Positive, select = -c(characteristics_ch1))
HIV_Positive<- na.omit(HIV_Positive)

HIV_Negative<- subset(p_data, characteristics_ch1.1 == "hiv status: HIV negative")
HIV_Negative<- subset(HIV_Negative, select = -c(characteristics_ch1))
HIV_Negative<- na.omit(HIV_Negative)


# creating two subset based on the  geographical location of the patients
HP1<- HIV_Positive$characteristics_ch1.2=="geographical region: Malawi"
HN1<- HIV_Negative$characteristics_ch1.2=="geographical region: South Africa"



# performing t-test
t_test_perform<- apply(data_log,1, function(x) t.test(x[HP1], x[HN1], na.rm = TRUE) )
P_Values<- sapply(t_test_perform, function(x) x$p.value)


#------------------------------------------Holm Correction--------------------------------------------------
adjusted_pvalues <- p.adjust(as.vector(P_Values), method="holm")

print(t_test_perform)
print(adjusted_pvalues)
#------------------------------------------LOG fold change--------------------------------------------------

LFC <- apply(data_log, 1, function(x) mean(x[HP1]) - mean(x[HN1]))
any(is.na(LFC))
print(LFC)


####################################### VOlcano Plot ###############################################
plot(LFC, -log10(P_Values), xlab = "Log2(Fold Change)", ylab = "-log10(P-Value)", main = "Volcano Plot")

# Add a horizontal line to indicate the significance threshold
abline(h = -log10(0.05), col = "red", lty = 2)

#-----------------------------------------------------------------------------------------------------------

# Load the ggplot2 package
library(ggplot2)

data_log_df <- as.data.frame(data_log)

# Define cutoffs for significant changes
fc_cutoff <- 1
pval_cutoff <- 0.05

# Create a volcano plot
ggplot(data_log_df, aes(x = LFC, y = -log10(P_Values))) +
  geom_point(aes(color = ifelse(abs(LFC) > fc_cutoff & P_Values < pval_cutoff, "red", "black"))) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey") +
  xlab("Log2 Fold Change") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot") +
  theme_classic()


######################################### STEP 5 ###########################################################


hiv_pos <- sample_info$characteristics_ch1.1 == "hiv status: HIV positive"
hiv_neg <- sample_info$characteristics_ch1.1 == "hiv status: HIV negative"

# Define the experimental design matrix
design <- model.matrix(~ 0 + factor(sample_info$characteristics_ch1.1))

# Assign meaningful column names to the design matrix
colnames(design) <- c("HIV_pos", "HIV_neg")

# Fit a linear model to the data using the design matrix
fit <- lmFit(data, design)

# Define the contrasts of interest for differential expression analysis
cont.matrix <- makeContrasts(
  HIV_pos_vs_neg = HIV_pos - HIV_neg,
  levels = design
)

# Empirical Bayes shrinkage of the standard errors
fit <- eBayes(contrasts.fit(fit, cont.matrix))

top.table <- topTable(fit, sort.by = "P", n = Inf)
head(top.table, 20)


# Save the results to a file
write.csv(top.table, "de_results.csv", row.names = FALSE)


# Create a volcano plot of the results
with(top.table, plot(logFC, -log10(P.Value), pch=20, main="Volcano Plot", xlim=c(-2,2)))
abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1,1), col="blue", lty=2)


######################################### STEP 6 ###########################################################


#I used a cutoff of adjusted p-value < 0.05 and absolute log fold change > 1 to identify differentially expressed genes. The adjusted p-value cutoff of 0.05 is a common cutoff used in many studies and it corresponds to a false discovery rate (FDR) of 5%. This means that among the genes declared significant at this cutoff, 5% of them are expected to be false positives.

#The absolute log fold change cutoff of > 1 means that we're only considering genes that show a fold change of at least 2-fold between two conditions. This value is somewhat arbitrary and may depend on the biological context of the study. However, a fold change of 2 is a commonly used threshold in many studies.



######################################### STEP 7 ###########################################################


library(clusterProfiler)

# Load the list of differentially expressed genes
de_genes <- read.csv("de_results.csv")$gene_symbol

# Set the background gene set (e.g., all genes in the genome)
background_genes <- read.csv("all_genes.csv")$gene_symbol

# Perform enrichment analysis using the gene ontology database
enrich_result <- enrichGO(de_genes, universe = background_genes, keyType = "Symbol", ont = "ALL")

# View the top enriched terms
head(enrich_result)


######################################### STEP 8 ###########################################################


# 1. gene_set : Set of genes that is being tested for Enrichment Analysis. (by using KEGG pathway database)
# 2. Universe : Set of all genes that were tested in original DEA.
# 3. background: Set of all genes.
# 4. Pvaluecutoff and Qvaluecutoff: cutoff for determining statistically significant enriched dataset.
#    p-value cutoff = 0.05
#    q-value cutoff = 0.1
#-----------------------------------------------------------------------------------------------------------

#  Bar plot to show the top 10 enriched pathways bssed on -log10(p-value)
barplot(enrichment_results[1:10,]$pvalue, main="Top 10 Enriched Pathways",xlab="-log10(p-value)",ylab="Pathway",horiz=T)

# loading the ggplot2 library to draws the volcano plot
library(ggplot2)

# volcano plot  

# dataframe creation for volcano plot
volcano_df <- data.frame(Pathway = enrichment_results$pathway,logFC = enrichment_results$log2FoldChange,logP = -log10(enrichment_results$pvalue),stringsAsFactors = F)

# creating the volcano plot
ggplot(volcano_df, aes(x=logFC, y=logP)) +geom_point(alpha=0.8, size=2, aes(color=Pathway)) +scale_color_manual(values=rainbow(length(unique(enrichment_results$pathway)))+geom_vline(xintercept=c(-1,1), linetype=2, color="gray40") +geom_hline(yintercept=-log10(0.05), linetype=2, color="gray40") +theme_bw() +labs(x="log2 Fold Change", y="-log10(p-value)", color="Pathway")

# library for creating heatmap
library(pheatmap)

# creating the dataframes for pheatmap 
heatmap_df <- subset(diff_exp_data, select=de_genes)
heatmap_df <- heatmap_df[row.names(heatmap_df) %in% enrichment_results$gene]

# creating the heatmap
pheatmap(heatmap_df, cluster_rows=F, cluster_cols=T, color = colorRampPalette(c("blue", "white", "red"))(50))


######################################### STEP 9 ###########################################################







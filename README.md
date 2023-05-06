PB Assignment-1

Himanshu- [2021464]



STEP-1
In this assignment's first step, I selected the “GSE37250” microarray dataset and downloaded it using the GEOquery package directly in R studio. After that, I stored the data in a Variable “Dataset.”
STEP-2
After storing the data in the Dataset variable, I separated the Gene Expression data, attributes, pdata, and fdata into separate variables. Then for the EDA and preprocessing part, the expressed data is normalized, log-transformed, and filtered out the data to remove the low-intensity probes along with the removal of rows with missing values.
Also, plotting some boxplots and histograms to visualize the preprocessed data.


STEP-3
Log transformation of data helps normalize and stabilize the data as the raw data.
It contains a wide range of values, and comparing log-transformed data is more effortless. It also linearizes the relationship between gene expression and biological effects, improving the interpretability of gene expression level changes.


STEP-4
In this step, we need to perform the Differential Expression Analysis using a t-test, and for that, first, we have to create two sample groups based on the HIV status of the patients(Samples). After that, we’ll perform the two-sample t-test.
After that, we’ll find the P-value from the result and adjust the p-value by the Holm correction method and calculate the LogFC.
At last we’ll draw a volcano plot.


STEP-5
In this step, we’ll do the same by using the Limma package. And storing the top 10 differentially expressed genes based on the LogFC and P-value cutoff.
STEP-6
I used a cutoff of adjusted p-value <0.05 and absolute log fold change >1 to identify the differentially expressed gene. The adjusted p-value cutoff of 0.05 is a common cutoff value used in many studies, and it corresponds to a false discovery rate of 5%. This means that 5% of the differentially expressed genes are expected to be false positive.
The absolute logFC cutoff of >1 means that we’re only considering genes that show a fold change of at least 2-fold between two conditions.


STEP-7
For Performing Enrichment analysis using the set of genes that you have obtained using the Gene set enrichment analysis method, first, we’ll load the required library and then load the list of differentially expressed genes and background genes(all human genes).
After that, we’ll perform the enrichment analysis using the gene ontology databases.


STEP-8
Parameters:
1. Gene_set: Set of genes that are being tested for Enrichment Analysis.  [By using the KEGG pathway database]
2. Universe: Set of all genes that were tested in Original DEA.
3. Background: a set of all genes in the Homo Sapiens genome.
4. Pvaluecutoff: cutoff for determining statistically significant enriched dataset. Pvaluecutoff = 0.05
5. Qvaluecutoff: cutoff for determining statistically significant enriched dataset. Qvaluecutoff = 0.1


At  last we’ll show various plots i.e Barplot, volcano plot and pheatmap


STEP-9

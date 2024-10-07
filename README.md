# ChickpeaOmicsR
![ChickpeaOmicsR](chickpeaOmicsR.png)

**ChickpeaOmicsR** is an R package designed to facilitate multi-omics data analysis for chickpea research. The package offers tools to analyze gene expression, protein interaction networks, and genome-wide association studies (GWAS), helping chickpea breeders and researchers explore key biological processes such as flowering, salinity response, and stress adaptation.

## Installation

To install **ChickpeaOmicsR**, you can use:

```r
# Install ChickpeaOmicsR
devtools::install_github("AlsammanAlsamman/ChickpeaOmicsR")
```

## Datasets

The package comes with several preloaded datasets to assist with chickpea genome and protein studies:

1. `CaExpressionCounts`: Gene expression counts for chickpea. Generated from RNA-seq data of Gene expression counts for 28891 gene and 197 samples.
2. `CaExpressionMetadata`: Metadata associated with gene expression studies.
3. `CaGenomeAnnotation`: Genome annotations of chickpea genes.
4. `CaProteinEnrichment`: Protein enrichment analysis of all Chickpea genes.
5. `CaProteinSequenceInfo`: Protein sequence information.
6. `CaProteinInteractionNetwork`: Protein-protein interaction data.
7. `CaExpressionSignificance`: Data on significant gene expression changes.
8. `CaGWAS`: Genome-wide association study data for chickpea.

These datasets can be used to perform various analyses and visualizations within the package.
```r
# Load the datasets
data(CaExpressionCounts)
head(CaExpressionCounts)

data(CaExpressionMetadata)
head(CaExpressionMetadata)

data(CaGenomeAnnotation)
head(CaGenomeAnnotation)

data(CaProteinEnrichment)
head(CaProteinEnrichment)

data(CaProteinSequenceInfo)
head(CaProteinSequenceInfo)

data(CaProteinInteractionNetwork)
head(CaProteinInteractionNetwork)
head(CaProteinInteractionNetwork$edges)

data(CaExpressionSignificance)
head(CaExpressionSignificance)

data(CaGWAS)
head(CaGWAS)
```


## Key Functions

### Gene Expression Analysis

Retrieve gene expression data by Gene Ontology (GO) terms and plot heatmaps to visualize expression patterns.

```r
# Retrieve genes associated with specific GO terms or keywords
genes <- get_geneids_bygotermkeywords(c("Flowering"))

# Extract log fold-change values for gene expression
geneExp <- get_gene_logfc_bycategory(unlist(genes), pvalue = 0.01)

# Plot a heatmap of the gene expression
plot_gene_expression_heatmap(geneExp, plot_height = 10, plot_width = 10, file = "gene_expression_heatmap.pdf")
```

This feature helps chickpea breeders identify key genes related to traits like flowering or stress responses.

### Protein Interaction Network Analysis

Visualize protein interaction networks and overlay gene expression data for more detailed insights into molecular interactions.

```r
# Retrieve the protein interaction network for the genes
ppiTable <- get_protein_interaction_table(unlist(genes))

# Get log fold-change values for gene expression
geneLOGFC <- get_gene_logfc(unlist(genes))

# Plot the protein interaction network
plot_protein_interaction_network_with_gene_expression(ppiTable, geneLOGFC, plot_height = 20, plot_width = 20, file = "network.pdf")
```

This feature is particularly useful for chickpea breeders who want to understand how protein interactions influence important traits.

### GO Term and GWAS Integration

Combine gene expression data with GWAS results for a more comprehensive view of genotype-phenotype associations.

```r
# Generate a gene-phenotype matrix based on GO term and GWAS data
gene_pheno_matrix <- gene_pheno_GWASmatrix_byGOKeyword("Flowering", pvalue = 0.05)

# Plot the gene-phenotype matrix
plot_gene_pheno_GWASmatrix_byGOKeyword(gene_pheno_matrix, outfolder = "heatmaps", plotwidth = 5, plotheight = 8, pvalue = 0.05)
```

This allows breeders to associate specific gene functions with phenotypic traits, speeding up marker discovery and breeding decisions.

## Application for Chickpea Breeders

**ChickpeaOmicsR** offers a comprehensive toolkit for chickpea breeders, enabling them to:

- Identify key genes involved in important traits like flowering, salinity tolerance, and disease resistance.
- Visualize gene expression patterns and protein interaction networks to understand complex molecular mechanisms.
- Integrate multi-omics data, including GWAS, to discover candidate genes for breeding programs.
- Conduct gene enrichment analysis to explore biological processes and pathways relevant to chickpea improvement.
- Access preloaded datasets for quick analysis and visualization of chickpea omics data.
- Make informed breeding decisions based on genotype-phenotype associations and gene expression profiles.
- Enhance collaboration and knowledge sharing among chickpea researchers through standardized data analysis tools.

By leveraging **ChickpeaOmicsR**, breeders can make more informed decisions in chickpea breeding programs, ultimately improving crop yield and resilience.
Each function supports common workflows in chickpea research, helping breeders and biologists analyze omics data efficiently and gain insights into gene expression, protein interactions, and genetic associations.

## Examples
Here are a few scenarios illustrating how the functions you provided can be useful:

1. **Get Data for a Specific Stress Study:**
   Suppose you're a chickpea researcher interested in how different stress factors (like salinity) affect gene expression. You can use the `get_study_data("Salinity")` function to extract the gene expression counts and metadata specifically for salinity studies, allowing you to focus on this area of interest.

   ```R
   # Example usage
   salinity_data <- get_study_data("Salinity")
   head(salinity_data$counts)
   ```

2. **Retrieve Data for Specific Samples:**
   If you're conducting a follow-up analysis on specific samples collected from your experiments, you can use the `get_samples_data()` function to subset gene expression data and metadata for only those samples. For example, extracting data from SRR10331638 to SRR10331645:

   ```R
   # Example usage
   samples_data <- get_samples_data(c("SRR10331638", "SRR10331639", "SRR10331640"))
   head(samples_data$counts)
   ```

3. **Obtain Metadata for a Group of Interest:**
   To get metadata based on specific study groups, experiments, or samples, you can use `get_metadata()`. For instance, you may want metadata about an experiment titled "PRJNA232700":

   ```R
   # Example usage
   experiment_metadata <- get_metadata(experiment = "PRJNA232700")
   head(experiment_metadata)
   ```

4. **Analyze GWAS Results for Specific Genes:**
   If you're interested in a GWAS analysis to identify associations between gene variants and traits like drought tolerance, the `get_gwas_data()` function allows you to filter the data by gene and p-value. This could help you find significant associations for your genes of interest.

   ```R
   # Example usage
   gwas_data <- get_gwas_data(c("GRX3", "HIS2A"), pvalue = 0.01)
   head(gwas_data)
   ```

5. **Fetch Gene Annotations:**
   To explore gene annotations, `get_gene_annotation_bygeneid()` retrieves functional information about the genes you're interested in. This function helps integrate gene functional context into your study:

   ```R
   # Example usage
   gene_annotation <- get_gene_annotation_bygeneid(c("GRX3", "HIS2A"))
   head(gene_annotation)
   ```

6. **Gene Enrichment Analysis for GO Terms:**
   For researchers investigating gene ontology (GO) terms related to specific biological processes, the `get_geneids_bygotermkeywords()` function can extract genes associated with keywords like "Flowering" or "Salinity." It allows you to focus on genes involved in relevant pathways:

   ```R
   # Example usage
   genes_in_flowering <- get_geneids_bygotermkeywords(c("Flowering"))
   head(genes_in_flowering)
   ```

7. **Analyze Protein-Protein Interaction Networks:**
   If you want to visualize interactions among proteins of interest, you can use `get_protein_interaction_table()` to extract network data and then `plot_protein_interaction_network()` to create a network graph. This helps to identify important interactions among proteins in your study:

   ```R
   # Example usage
   ppi_table <- get_protein_interaction_table(c("LOC101501102", "LOC101507179"))
   plot_protein_interaction_network(ppi_table, file = "network.pdf")
   ```
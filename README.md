# RNA-seq Count Data Analysis for Different Types of Cancer

This repository contains an R script for analyzing RNA-seq count data to investigate different types of cancer. The main goal of this project is to use various computational tools to detect differentially expressed genes among various cancer types, identify immune-related genes, and create a network for gene-gene interaction.

## Table of Contents

- [Introduction](#introduction)
- [Preparation and Dependencies](#preparation-and-dependencies)
- [Data Analysis Workflow](#data-analysis-workflow)
- [Running the Analysis](#running-the-analysis)
- [Conclusion](#conclusion)

## Introduction

RNA-seq analysis has become an essential tool in understanding gene expression levels in various diseases, especially cancer. This project aims to comprehensively analyze the RNA-seq count data for different types of cancer, identify differentially expressed genes, and infer the underlying gene-gene interaction network.

## Preparation and Dependencies

This project utilizes several R packages. You need to ensure these packages are installed before running the script. The primary packages used include `DESeq2`, `tidyverse`, `WGCNA`, and a few others for visualization and data manipulation.

## Data Analysis Workflow

1. **Data Loading and Preparation**: The script begins by loading RNA-seq count data, metadata, and gene-related information. We filter the data to focus only on immune-related genes.

2. **Exploratory Data Analysis**: After data preparation, an exploratory data analysis is conducted. Here, we normalize the data and perform principal component analysis (PCA) to reduce dimensionality and visualize our data.

3. **Differential Gene Expression Analysis**: We perform differential gene expression analysis using DESeq2. Genes that show significant differences in expression between different types of cancer are identified.

4. **Weighted Gene Co-expression Network Analysis (WGCNA)**: This is a crucial step where we create a network of gene-gene interactions. WGCNA is a method for finding clusters (modules) of highly correlated genes. For each module, we calculate a module eigengene or 'representative' gene expression profile.

5. **Module-Trait Associations**: We associate the identified modules with clinical traits. This step can provide insights into the genes that might be critical for different cancer types.

6. **Intramodular Analysis**: Finally, we perform an intramodular analysis to identify potential driver genes within the identified modules.

## Running the Analysis

Before running the script, make sure you have the required data files and the R packages installed. Adjust the script to match your specific data file names and directories. Run the script in R, RStudio, or any R-compatible interface.

## Conclusion

The end goal of this project is to utilize the power of computational biology to contribute to our understanding of cancer. By identifying key genes and gene networks associated with different types of cancer, we can provide valuable insights into the underlying biological mechanisms. This can potentially aid in developing novel therapeutic strategies or improving existing ones.

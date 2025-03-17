# RNA-seq Count Data Analysis in Pediatric Brain Tumors


This repository contains an R script designed for the comprehensive analysis of RNA-seq count data related to different types of pediatric brain tumors including ATRT, Ependymoma, Glioblastoma, Glioma, Medulloblastoma, and Craniopharyngioma. The goal of this project is to identify immune-related gene expressions and discover their disparities among these tumor types.

## Table of Contents

- [Introduction](#introduction)
- [Preparation and Dependencies](#preparation-and-dependencies)
- [Data Analysis Workflow](#data-analysis-workflow)
- [Running the Analysis](#running-the-analysis)
- [Conclusion](#conclusion)

## Introduction

RNA-seq analysis is a powerful tool for understanding gene expression variations across different diseases, in this case, pediatric brain tumors. This project conducts a differential expression analysis to investigate the disparities in immune-related gene expressions among different types of pediatric brain tumors.

## Preparation and Dependencies

This analysis utilizes various R packages, including `DESeq2`, `tidyverse`, and `WGCNA`, among others. Ensure these packages are installed before running the script.

## Data Analysis Workflow

1. **Data Loading and Preparation**: We load the RNA-seq count data, metadata, and gene-related information, focusing on immune-related genes for our study.

2. **Exploratory Data Analysis**: The data is normalized and principal component analysis (PCA) is performed for dimensionality reduction and data visualization.

3. **Differential Gene Expression Analysis**: Using DESeq2, we conduct a differential gene expression analysis to identify the genes that show significant expression differences across the tumor types.

4. **Weighted Gene Co-expression Network Analysis (WGCNA)**: We identify clusters (modules) of highly correlated genes and compute a module eigengene for each module, representing a 'typical' gene expression profile within the module.

5. **Module-Trait Associations**: The modules identified are associated with clinical traits to identify potential gene sets that might be critical in the studied tumor types.

6. **Intramodular Analysis**: This final step identifies potential driver genes within the modules identified, providing directions for future research.

## Running the Analysis

Ensure you have the necessary data files and R packages installed before running the script. Adjust the file names and directories as per your dataset. You can run the script in R, RStudio, or any R-compatible interface.

## Conclusion

By identifying key genes and gene networks associated with different types of pediatric brain tumors, this project contributes to understanding the underlying biological mechanisms of these malignancies. This understanding could potentially aid in the development of new therapeutic strategies and improvement of existing ones.

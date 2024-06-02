# A Comprehensive Bioinformatics Analysis of Diffuse Type Gastric Cancer: Insights from WGCNA,Network and Pathway Enrichment Analysis

## Introduction

This repository contains the code for a comprehensive bioinformatics analysis of diffuse type gastric cancer (DGC), utilizing Weighted Gene Co-expression Network Analysis (WGCNA). This study aims to elucidate the molecular intricacies of DGC by identifying key gene modules and potential biomarkers.

## Abstract

This study employs an integrative approach, combining systems biology methodologies and bioinformatics analyses, to unravel the molecular intricacies of DGC. Utilizing WGCNA on gene expression data from the Gene Expression Omnibus (GEO) database (accession number GSE113255), we identified four distinct modules with unique gene expression patterns associated with DGC. Pathway enrichment analysis using Gene Ontology (GO) and KEGG pathways revealed significant involvement of specific modules in gastric cancer pathogenesis. Hub gene analysis identified key players within these modules, providing comprehensive insights into the molecular landscape of DGC and pinpointing potential therapeutic targets and biomarkers.

## Installation

To run the analysis, ensure you have R and the necessary packages installed. Follow the steps below:

1. Install R from [CRAN](https://cran.r-project.org/).
2. Install the required R packages:

```R
install.packages(c("WGCNA", "DESeq2", "flashClust", "ggplot2", "shiny", "Cytoscape", "BiocManager"))
BiocManager::install(c("GO.db", "KEGG.db", "org.Hs.eg.db"))


```markdown
# Comprehensive Bioinformatics Analysis of Diffuse Type Gastric Cancer Using WGCNA

## Introduction

This repository contains the code for a comprehensive bioinformatics analysis of diffuse type gastric cancer (DGC), utilizing Weighted Gene Co-expression Network Analysis (WGCNA). This study aims to elucidate the molecular intricacies of DGC by identifying key gene modules and potential biomarkers.

## Abstract

This study employs an integrative approach, combining systems biology methodologies and bioinformatics analyses, to unravel the molecular intricacies of DGC. Utilizing WGCNA on gene expression data from the Gene Expression Omnibus (GEO) database (accession number GSE113255), we identified four distinct modules with unique gene expression patterns associated with DGC. Pathway enrichment analysis using Gene Ontology (GO) and KEGG pathways revealed significant involvement of specific modules in gastric cancer pathogenesis. Hub gene analysis identified key players within these modules, providing comprehensive insights into the molecular landscape of DGC and pinpointing potential therapeutic targets and biomarkers.

## Installation

To run the analysis, ensure you have R and the necessary packages installed. Follow the steps below:

1. Install R from [CRAN](https://cran.r-project.org/).
2. Install the required R packages:

```R
install.packages(c("WGCNA", "DESeq2", "flashClust", "ggplot2", "shiny", "Cytoscape", "BiocManager"))
BiocManager::install(c("GO.db", "KEGG.db", "org.Hs.eg.db"))
```

## Usage

1. Data Collection and Preprocessing:
   - Gene expression data was obtained from the GEO database (accession number GSE113255).
   - Load and preprocess the data using the provided R scripts.

2. Weighted Gene Co-expression Network Construction:
   - Use the `pickSoftThreshold()` function to select an appropriate soft threshold.
   - Construct weighted adjacency matrices and topological overlap matrices (TOM).

3. Identification of Modules:
   - Apply dynamic tree-cutting methods to identify gene modules.
   - Calculate module eigengenes (ME), gene significance (GS), and module membership (MM).

4. Visualization:
   - Generate heatmaps, dendrograms, and scatter plots to visualize module-trait relationships and gene clustering.

5. Pathway Enrichment Analysis:
   - Perform GO and KEGG pathway analyses using R Shiny.
   - Visualize biomarkers and hub genes using Cytoscape.

6. Results and Discussion:
   - Analyze the identified modules and hub genes.
   - Discuss the potential implications for DGC research and therapeutic strategies.

## Features

- Comprehensive Analysis: Integrates multiple bioinformatics approaches for a thorough analysis of DGC.
- Modular Structure: Code is organized into modular scripts for easy understanding and modification.
- Visualization Tools: Includes scripts for generating various visualizations to interpret the results effectively.

Key Findings

Through our analysis, we have identified the following top genes associated with gastric cancer:

Blue Violet Module
SULF1: Down-regulated in gastric tumors, potentially linked to cancer development.
BGN: Expression correlates positively with immune cells like T-cells and macrophages.
MMP3: Associated with an increased risk of gastric cancer.
INHBA: Overexpression linked to poorer prognosis in gastric cancer patients.

Cyan Module
COL1A1: Increased expression might indicate early gastric cancer.
MMP1: Implicated in different stages of gastric cancer progression.
COL1A2: Higher levels correlate with larger tumor size and deeper invasion.
SPARC: Plays a role in various stages of cancer progression.
These genes are significantly related to gastric cancer pathogenesis and provide potential targets for further investigation and therapeutic strategies.

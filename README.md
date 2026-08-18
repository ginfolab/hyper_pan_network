# A Hypergraph-Based Pan-Network Framework for Exploring the Organizational Principles of the *Escherichia coli* Transcriptome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Language: Bash](https://img.shields.io/badge/Language-Bash-green.svg)](https://www.gnu.org/software/bash/)

This repository contains the computational framework, analysis pipelines, and reproduction scripts for our pan-transcriptomic hypergraph study in *Escherichia coli*. 

By integrating 106 diverse transcriptomic datasets, we reconstruct condition-invariant core modules and adaptive plastic network communities using frequent itemset mining (Apriori algorithm) and hypergraph modeling.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Prerequisites & Dependencies](#-prerequisites--dependencies)
- [Data Availability](#-data-availability)
- [Citation](#-citation)
- [License](#-license)

---

## 🔬 Overview

Traditional pairwise co-expression networks often struggle with noise and fail to capture multi-gene coordinated actions across diverse physiological perturbations. 

This project provides:
1. **Hypergraph Construction Pipeline**: Uses Frequent Itemset Mining (FIM / Apriori) to discover recurring multi-gene co-expression hyperedges across 106 GEO transcriptomic datasets.
2. **Bipartite Organizational Principle of the Core Network**: By systematically profiling the condition-specific activation states (modularity) of core communities across diverse datasets, the framework is capable of decoupling the network into a bipartite physiological architecture. It clearly distinguishes between an **invariant backbone** (highly independent and universally activated core units) and **adaptive plasticity** (stress-responsive communities that undergo dynamic topological rewiring and coordinated co-activation under environmental transitions).
3. **Downstream Characterization**: Functional enrichment, topological modularity analysis, and cross-mapping with RegulonDB and KEGG pathways.

---

## 📁 Repository Structure

The project is organized into five main directories to ensure full reproducibility and modularity. Each folder contains its own detailed `README.md` for specific execution guidelines:

* **[`01_pipeline/`](./01_pipeline)**: The main analytical engine containing R and Bash scripts for GEO metadata acquisition, TPM integration, hypergraph construction, power-law decoupling, and modularity evaluation.
* **[`02_figures/`](./02_figures)**: Visualization scripts utilized to generate the manuscript's main figures and supplementary charts (e.g., density distributions, auto-clustered heatmaps, and Sankey diagrams).
* **[`03_data/`](./03_data)**: Pre-processed intermediate data, expression matrices, and object files (`.RDS`) to facilitate direct inspection and reproduction.
* **[`04_Cyto/`](./04_Cyto)**: Node and edge tables explicitly formatted for network topology visualization in Gephi and Cytoscape.
* **[`05_Annotation/`](./05_Annotation)**: Functional annotation references, including mapped structures from RegulonDB (operons/regulons) and pan-genome metadata from MBGD.

---

## 📦 Prerequisites & Dependencies

The scripts are developed primarily in **R** and **Bash**. Ensure you have R (>= 4.0.0) installed along with the required packages:

```r
# Data mining & Network analysis
install.packages(c("arules", "igraph", "dplyr", "tidyr", "data.table"))
install.packages(c("future", "furrr")) # Required for parallel network scanning

# Visualization
install.packages(c("ggplot2", "patchwork", "ggraph", "tidygraph", "RColorBrewer"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "clusterProfiler"))

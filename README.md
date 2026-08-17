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
- [Quick Start / Workflow](#-quick-start--workflow)
  - [Track A: Full Pipeline (From Raw Data)](#track-a-full-pipeline-from-raw-data)
  - [Track B: Fast Reproduction (From Processed Matrices)](#track-b-fast-reproduction-from-processed-matrices)
- [Detailed Script Index](#-detailed-script-index)
  - [Part 1: Pan-Network Construction & Analysis](#part-1-pan-network-construction--analysis)
  - [Part 2: Manuscript Figures & Tables](#part-2-manuscript-figures--tables)
- [Data Availability](#-data-availability)
- [Citation](#-citation)
- [License](#-license)

---

## 🔬 Overview

Traditional pairwise co-expression networks often struggle with noise and fail to capture multi-gene coordinated actions across diverse physiological perturbations. 

This project provides:
1. **Hypergraph Construction Pipeline**: Uses Frequent Itemset Mining (FIM / Apriori) to discover recurring multi-gene co-expression hyperedges across 106 GEO transcriptomic datasets.
2. **Core vs. Adaptive Network Partitioning**: Identifies stable transcriptional backbones (high universality $U$) versus conditionally plastic responses.
3. **Downstream Characterization**: Functional enrichment, topological modularity analysis, and cross-mapping with RegulonDB and KEGG pathways.

---

## 📦 Prerequisites & Dependencies

The scripts are developed primarily in **R** and **Bash**. Ensure you have R (>= 4.0.0) installed along with the following packages:

```r
# Data mining & Network analysis
install.packages(c("arules", "igraph", "dplyr", "tidyr", "data.table"))

# Visualization
install.packages(c("ggplot2", "pheatmap", "ggalluvial", "RColorBrewer"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "clusterProfiler"))

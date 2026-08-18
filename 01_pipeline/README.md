# Pan-Network Analytical Pipeline: Chronological Tutorial

This directory contains the core analytical engine of the *E. coli* hypergraph pan-network project. The scripts provided here form a complete, end-to-end pipeline—from raw data acquisition to hypergraph construction and biological benchmarking.

## ⚠️ Global Precautions & Setup
1. **Hardcoded Paths**: The scripts currently contain environment-specific absolute paths (e.g., `/Users/jiangzhenbo/...` or `/lustre/home/...`). **You must modify `setwd()` and index paths** in the scripts to match your local or server directory structure before execution.
2. **Conda Environments**: The Shell scripts assume the presence of a Conda environment named `salmon`. Ensure you have Conda installed and the environment properly configured (`source ~/anaconda3/bin/activate salmon`).
3. **Hardware Requirements**: The pipeline utilizes heavy parallel computing (`mclapply`, `furrr`). A machine with **>= 32GB RAM** and a multi-core CPU is strongly recommended.

---

## Phase 1: Data Acquisition & Preprocessing

### Step 1: Pan-genome Reference & Indexing
Extracts the *E. coli* pan-genome sequences from the MBGD database, formats them into a standard FASTA file, and builds the mapping index using Salmon.
> 💡 **Precaution:** Ensure MySQL is accessible for querying the MBGD database.

```bash
bash Ecoli_pangenome_new.sh
```

### Step 2: GEO Metadata Retrieval
Batch-downloads and parses experimental metadata (strain, genotype, growth conditions) for all 106 GSE datasets using the GEOquery package, building a robust metadata registry.

```bash
Rscript "Backup and download the description information of the GEO data.R"
```

### Step 3: Automated Transcript Quantification (Single-end & Paired-end)
Executes the core quantification process (prefetching SRA, converting to FASTQ, and Salmon quantification).
💡 Crucial Operational Guide: > 1. For each target GSE directory (e.g., GSE12345/), you must ensure a corresponding SRR_Acc_List.txt file (containing all target SRR accession IDs) is pre-placed inside.
2. You must distribute the respective core script (single_Ecoli.sh or paired_Ecoli.sh) into each GSE subfolder before running the batch wrappers.
3. The wrapper scripts feature a fail-safe size-check mechanism (re-running if TPM.txt is missing or < 1MB). To conserve disk space, intermediate .sra and .fastq files are automatically deleted after successful quantification.

```bash
# Run for single-end datasets
bash New_Run_single_TPM_size.sh

# Run for paired-end datasets
bash New_Run_paired_TPM_size.sh
```

### Step 4: TPM Integration & Quality Control
Aggregates the individual `TPM.txt` outputs into a unified expression matrix. This script applies critical quality controls: filtering out low-expressed genes, capping overly large datasets to 30 samples to prevent size bias, and discarding datasets with fewer than 15 valid samples.

```bash
Rscript Combine_TPM.R
```

---

## Phase 2: Hypergraph Construction Engine

### Step 5: Network Inference & Hyperedge Mining
*This is the most computationally intensive step.* It calculates FDR-corrected Pearson correlations, partitions the networks using Affinity Propagation (AP) clustering, and employs frequent itemset mining (Eclat/Apriori) to discover multi-gene hyperedges. It utilizes `data.table` and multi-core `mclapply` for extreme memory optimization.

```bash
Rscript "From integrating the TPM files to Apriori (histogram-based statistical network features).R"
```

---

## Phase 3: Benchmarking & Biological Optimization

### Step 6: Internal Topological Benchmarking
Calculates network topological properties (Modularity, Clustering Coefficient, Density) across varying Universality ($U$) cutoffs. Identifies the empirical stabilization point ($U=15$) where the core network architecture forms.

```bash
Rscript "Internal index compares differences in network types (such as modularity, etc.).R"
```

### Step 7: Regulatory Ground Truth Preparation
Parses standard transcriptional regulation definitions from RegulonDB (v12.0). It constructs hierarchical reference dictionaries for Operons, `RC_only`, `RC_muti`, and Regulons to serve as biological ground truth for downstream validation.

```bash
Rscript "Create the Regulon and Regulon_combination files.R"
```

### Step 8: External Biological Benchmarking
Benchmarks the derived core modules against the RegulonDB Operon reference. It computes Precision, Recall, and Mapped Counts across dynamic $U$ cutoffs, demonstrating the superior biological fidelity of the hyperedge-based method over traditional pairwise thresholds.

```bash
Rscript "External index module mapping of Operon's statistical graph.R"
```

### Step 9: Core Module Mapping & Threshold Optimization
Systematically maps the finalized core communities against the full hierarchical RegulonDB tiers. It dynamically scans F-score thresholds (0.1 to 1.0) to optimize boundaries and outputs comprehensive statistical visuals (ridge plots, stacked bar charts) and regulatory tree graphs via `ggraph`.

```bash
Rscript "Query_pan_network(Creating the Core Network Module).R"
```

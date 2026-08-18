# Core Dataset Repository (`03_data`)

This directory stores the essential intermediate and background reference datasets required for the pan-transcriptomic and hypergraph network analyses of *E. coli*. 

> **Note on File Size Management:** To comply with GitHub's file size limits and maintain repository efficiency, large raw expression files (such as bulky TPM matrices) have been omitted from this version control path. Only core structural and annotation summaries are hosted here.

---

## 📂 Included Files Description

* **[AP_info_all_0.5.txt](./AP_info_all_0.5.txt)**
  * **Description:** Contains the comprehensive Affinity Propagation (AP) clustering results derived from dataset-specific Pearson correlation matrices (at quality/parameter cutoff $q = 0.5$). 
  * **Purpose:** Serves as the fundamental input mapping table linking individual gene identifiers to their respective batch-specific AP clusters, forming the basis for subsequent itemset and hypergraph mining.

* **[The latest annotation of RegulonDB genes.RDS](./The%20latest%20annotation%20of%20RegulonDB%20genes.RDS)**
  * **Description:** The latest standardized background annotation database parsed from RegulonDB (v12.0) and custom combinatorial dictionaries (including Operons, Regulons, `RC_only`, and `RC_muti` tiers).
  * **Purpose:** Acts as the primary biological ground-truth reference for functional module matching, Precision/Recall benchmarking, and regulatory network mapping.

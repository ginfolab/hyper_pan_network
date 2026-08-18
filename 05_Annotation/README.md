# Functional Annotation and Reference Database (`05_Annotation`)

This directory houses comprehensive functional annotation results, gene ontology mappings, pathway references, and specialized sequence sets for downstream biological interpretation and enrichment analysis of *E. coli* network modules.

---

## 📂 Included Files Description

* **[ALL_module_(Count >= 2 & qvalue < 0.05)_annotation.RDS](./ALL_module_%28Count%20%3E%3D%202%20%26%20qvalue%20%3C%200.05%29_annotation.RDS)**
  * **Description:** An R serialized data file containing the aggregated functional enrichment and annotation results across all network modules (stringently filtered by overlap count $\ge 2$ and adjusted $q$-value $< 0.05$).
  * **Purpose:** Serves as the primary data source for querying module-specific biological functions and functional enrichment summaries.

* **[Ec_gene2go](./Ec_gene2go)**
  * **Description:** The *E. coli*-specific mapping table linking individual gene identifiers to their corresponding Gene Ontology (GO) terms.
  * **Purpose:** Provides the backbone association dictionary required for functional enrichment and gene-set over-representation tests.

* **[Ec_go2term](./Ec_go2term)**
  * **Description:** The ontological dictionary defining GO terms alongside their descriptive functional category definitions and ontologies (Biological Process, Molecular Function, Cellular Component).
  * **Purpose:** Acts as the reference lookup dictionary to translate GO IDs into readable functional descriptions during visualization and reporting.

* **[Ecoli_genename.fa](./Ecoli_genename.fa)**
  * **Description:** A FASTA-formatted sequence database containing standard *E. coli* gene names and corresponding nucleotide/protein sequences.
  * **Purpose:** Used for sequence-level validation, identifier translations, and homology matching across the pipeline.

* **[Kmodule4Ec05.txt](./Kmodule4Ec05.txt)**
  * **Description:** A curated text file recording KEGG module mappings and functional pathway associations tailored for the project's analytical framework.
  * **Purpose:** Facilitates metabolic pathway reconstruction and functional module profiling.

* **[mge_genes.txt](./mge_genes.txt)**
  * **Description:** A specialized list identifying mobile genetic element (MGE) related genes (such as transposons, phages, and plasmids) within the *E. coli* genome scope.
  * **Purpose:** Enables the identification and filtering of horizontally acquired or mobile genomic components during core network evaluation.

* **[target_gene_list_final.RDS](./target_gene_list_final.RDS)**
  * **Description:** An R serialized list containing the finalized, high-confidence target gene subsets used across core analytical workflows.
  * **Purpose:** Acts as the definitive gene universe or filtered target pool for subsequent network topology and functional mapping scripts.

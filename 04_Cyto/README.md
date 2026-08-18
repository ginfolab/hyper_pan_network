# Cytoscape Network Visualization Directory (`04_Cyto`)

This directory contains the essential node and edge files used for network visualization, topological layout, and module-group interaction analysis in **Cytoscape**.

---

## 📂 Included Files Description

* **[Module(Cyto)_relation_node.txt](./Module(Cyto)_relation_node.txt)**
  * **Description:** The comprehensive network node attribute table. It records detailed metadata for each network module (e.g., `Node`, `Node_name`, `Size`, `Group`, core/non-core gene distribution metrics), functional enrichment summaries (`GO_term`), custom color mappings (`Module_color_mapping`, `Function_color`, `Region_color`), and hierarchical regional classifications (`Region`, `SubRegion`).
  * **Purpose:** Serves as the primary node attribute source in Cytoscape for mapping structural sizes, functional annotations, and multi-tier regional color schemes.

* **[Module(Cyto)_relation_edge.txt](./Module(Cyto)_relation_edge.txt)**
  * **Description:** The comprehensive network edge attribute table. It captures relationships between modules (e.g., `Module1`, `Module2`), interaction categorizations (`Edge_type` such as *Between_group* or *Within_group*), evolutionary dynamics (`Evolutionary_Pattern`, `Dominant_Clade`), specificity strengths, clade-specific percentages (`Pct_Edge_Clade_*`), itemset counts (`Itemset_num`), and shared gene intersections (`common_Itemset_genes`).
  * **Purpose:** Serves as the edge attribute source in Cytoscape for defining network connectivity, evolutionary trajectories, and inter-module hyperedge sharing characteristics.

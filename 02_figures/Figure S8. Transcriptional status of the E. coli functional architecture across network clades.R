# ==============================================================================
# Figure S8: 核心转录组网络拓扑可塑性与功能分区全景图 (CV-Based Plasticity)
# ==============================================================================

# ==============================================================================
# 0. 加载必要的 R 包
# ==============================================================================
suppressMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(magick)
  library(RColorBrewer)
})

# ==============================================================================
# 1. 设置路径与加载外部注释与拓扑数据
# ==============================================================================
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"

core_gene_path    <- file.path(base_path, "PCC/EC_coregene.txt")
network_gene_path <- file.path(base_path, "Gephi/clustid2Module2Step.RDS")
excel_file_path   <- file.path(base_path, "PCC/New_pangenome_Escherichia_coli_TPM5.xlsx")
Mobile_and_Defense_Genetic_Elements_path <- file.path(base_path, "MGE_table.txt")

New_Gene_Anno_path     <- file.path(base_path, "Module_Gene_Region_Annotation.RDS")
New_Dataset_Clade_path <- file.path(base_path, "Dataset_Clade_Mapping.txt")

# --- 加载基础数据 ---
Core_gene_raw <- read_delim(core_gene_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)$X1
Core_genes    <- paste0("Ec_", Core_gene_raw)

Core_network_data <- readRDS(network_gene_path)
Core_network_genes <- Core_network_data$id 

Gene_anno_df <- readRDS(New_Gene_Anno_path)
subregion_lookup <- setNames(Gene_anno_df$SubRegion, Gene_anno_df$Ec_clustid)
region_lookup    <- setNames(Gene_anno_df$Region, Gene_anno_df$Ec_clustid)

# ==============================================================================
# 2. 读取 TPM 矩阵并计算【转录可塑性 (CV)】- 融合 Presence 逻辑
# ==============================================================================
sheet_names <- excel_sheets(excel_file_path)
sheet_cv_list <- list()

message(">>> 正在读取 Excel：存在即计算 CV，缺失则留空 (深灰)...")
pb <- txtProgressBar(min = 0, max = length(sheet_names), style = 3)

for (i in seq_along(sheet_names)) {
  sheet_id <- sheet_names[i]
  raw_df <- suppressMessages(read_xlsx(excel_file_path, sheet = sheet_id))
  
  # 提取原始基因名并清洗 (使用你的原始逻辑)
  gene_names_raw <- raw_df[[1]]
  clean_genes <- str_replace(gene_names_raw, "^[^&]+&[^&]+&", "Ec_")
  
  # 提取数值矩阵并计算逐行 CV
  expr_data <- as.matrix(raw_df[, -1])
  class(expr_data) <- "numeric"
  
  row_cv <- apply(expr_data, 1, function(x) {
    x <- x[!is.na(x)]
    if(length(x) < 2) return(0)  # 如果只有1个样本，无法波动，记为 0 (最稳定)
    m <- mean(x)
    if(m == 0) return(0)         # 均值为0，无波动，记为 0
    return(sd(x) / m)            # 计算变异系数
  })
  
  # 汇总：合并清洗后的基因名与 CV (去重取最大变异)
  temp_df <- data.frame(Gene = clean_genes, CV = row_cv) %>%
    group_by(Gene) %>%
    summarise(CV = max(CV, na.rm = TRUE), .groups = "drop")
  
  # 将结果存入字典 (只有在这个表里出现的基因才会有记录)
  sheet_cv_list[[sheet_id]] <- setNames(temp_df$CV, temp_df$Gene)
  setTxtProgressBar(pb, i)
}
close(pb)

all_excel_genes <- sort(unique(unlist(lapply(sheet_cv_list, names))))
target_genes <- intersect(all_excel_genes, Core_network_genes)
message(paste0("\n>>> 最终绘图基因数 (交集): ", length(target_genes)))

# 【核心修改】：构建全 NA 矩阵。只要基因没在这个 sheet 出现，它就永远是 NA
cv_matrix <- matrix(NA, nrow = length(sheet_names), ncol = length(target_genes))
rownames(cv_matrix) <- sheet_names  
colnames(cv_matrix) <- target_genes

# 填装 CV 值：只有存在的基因才会被赋值
for (sheet in sheet_names) {
  current_cvs <- sheet_cv_list[[sheet]]
  intersect_genes <- intersect(names(current_cvs), target_genes)
  if(length(intersect_genes) > 0) {
    cv_matrix[sheet, intersect_genes] <- current_cvs[intersect_genes]
  }
}
filtered_cv_df <- as.data.frame(cv_matrix)


# ==============================================================================
# 3. 统一行排序 (对接 Figure 5 最新 Dataset / Clade 顺序)
# ==============================================================================
message(">>> 正在同步 Figure 5 演化分支结构...")
dataset_clade_df <- read_tsv(New_Dataset_Clade_path, show_col_types = FALSE)

ordered_datasets <- dataset_clade_df$Dataset
valid_datasets <- intersect(ordered_datasets, rownames(filtered_cv_df))
sorted_final_df <- filtered_cv_df[valid_datasets, ]

plot_matrix <- as.matrix(sorted_final_df)


# ==============================================================================
# 4. 精准锁定列分组信息 (SubRegion 级别 - 带有前缀免疫)
# ==============================================================================
all_plot_genes <- colnames(plot_matrix)

clean_plot_genes  <- str_remove(all_plot_genes, "^Ec_")
clean_lookup_keys <- str_remove(names(subregion_lookup), "^Ec_")

safe_subregion_lookup <- setNames(unname(subregion_lookup), clean_lookup_keys)
safe_region_lookup    <- setNames(unname(region_lookup), clean_lookup_keys)

current_subregions <- safe_subregion_lookup[clean_plot_genes]
current_subregions[is.na(current_subregions)] <- "Other"

current_parent_regions <- safe_region_lookup[clean_plot_genes]
current_parent_regions[is.na(current_parent_regions)] <- "Other"

unique_valid_subregions <- sort(unique(current_subregions[current_subregions != "Other"]))
subregion_levels <- c(unique_valid_subregions, "Other")
column_split_factor <- factor(current_subregions, levels = subregion_levels)


# ==============================================================================
# 5. 组装全局多维注释 (Top & Left) - 【绝对纯净命名版】
# ==============================================================================
message(">>> 正在生成高维联合注释系统...")

# --- 5.1 SubRegion 颜色映射 ---
region_base_colors <- c("Region_1"="#7876B1", "Region_2"="#E18727", "Region_3"="#0072B5", "Region_4"="#35A595", "Other"="grey80")

subregion_colors <- sapply(subregion_levels, function(sr) {
  if (sr == "Other") return("grey80")
  parent_reg <- str_extract(sr, "Region_\\d+") 
  return(as.character(region_base_colors[parent_reg])) 
})

subregion_counts <- table(column_split_factor)
subregion_with_n <- paste0(names(subregion_counts), " (n=", as.numeric(subregion_counts), ")")
names(subregion_with_n) <- names(subregion_counts)

anno_vec_subregion <- subregion_with_n[as.character(column_split_factor)]
anno_col_subregion <- setNames(as.character(subregion_colors), subregion_with_n[names(subregion_colors)])
anno_col_subregion <- anno_col_subregion[!is.na(names(anno_col_subregion))]

# --- 5.2 MGE & Defense 标记 ---
mge_df <- read.table(Mobile_and_Defense_Genetic_Elements_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mge_lookup <- setNames(mge_df$type, mge_df$Ec_clustid)

anno_mge_type <- mge_lookup[all_plot_genes]
anno_mge_type[is.na(anno_mge_type) | anno_mge_type == ""] <- "None"

raw_mge_colors <- c("GI"="#33A02C", "phage"="#E31A1C", "None"="#F0F0F0")
present_types <- setdiff(unique(anno_mge_type), "None")
extra_types <- setdiff(present_types, names(raw_mge_colors))
if(length(extra_types) > 0) {
  extra_cols <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(extra_types)), extra_types)
  raw_mge_colors <- c(raw_mge_colors, extra_cols)
}

used_mge_types <- unique(anno_mge_type)
final_mge_colors <- raw_mge_colors[used_mge_types]
mge_colors_vec <- setNames(as.character(final_mge_colors), names(final_mge_colors))

# --- 5.3 Core Genome 标记 ---
anno_vec_core <- ifelse(all_plot_genes %in% Core_genes, "Core", "Non-Core")
core_colors_vec <- setNames(c("#D95F02", "#F0F0F0"), c("Core", "Non-Core"))

# --- 5.4 组装 Top Annotation ---
anno_colors <- list(
  Core_genome          = core_colors_vec,
  Clustering_SubRegion = anno_col_subregion, 
  MGE_Defense          = mge_colors_vec               
)

top_anno <- HeatmapAnnotation(
  Core_genome          = anno_vec_core, 
  Clustering_SubRegion = anno_vec_subregion,
  MGE_Defense          = anno_mge_type,
  col                  = anno_colors, 
  simple_anno_size     = unit(0.3, "cm"), 
  show_annotation_name = TRUE, 
  annotation_name_side = "left",
  annotation_label     = c("Core genome", "Clustering SubRegion", "MGE & Defense Type"),
  annotation_name_gp   = gpar(fontsize = 6, fontface = "bold")
)

# --- 5.5 组装 Left Annotation (Clade 物理染色切分) ---
tree_colors_row <- c("#D9A441", "#734C7A", "#3D726D")
clade_names <- sort(unique(dataset_clade_df$Clade)) 
clade_colors_vec <- setNames(tree_colors_row, clade_names)

row_clades <- dataset_clade_df$Clade[match(rownames(plot_matrix), dataset_clade_df$Dataset)]

left_anno <- rowAnnotation(
  Clade = row_clades,
  col = list(Clade = clade_colors_vec),
  simple_anno_size = unit(0.3, "cm"),
  show_annotation_name = FALSE
)


# ==============================================================================
# 6. 准备 Right Annotation (按 Region 汇总的平均可塑性条形图)
# ==============================================================================
valid_regions <- sort(unique(current_parent_regions[current_parent_regions != "Other"]))

anno_bar_list <- list()
for(reg in valid_regions) {
  reg_genes <- all_plot_genes[current_parent_regions == reg]
  if(length(reg_genes) > 0) {
    # na.rm = TRUE 忽略那些深灰色 (不表达) 的情况
    row_means <- as.numeric(rowMeans(plot_matrix[, reg_genes, drop = FALSE], na.rm = TRUE))
    row_means[is.nan(row_means)] <- 0 
    
    anno_bar_list[[reg]] <- anno_barplot(
      row_means,
      which = "row", 
      gp = gpar(fill = as.character(region_base_colors[reg]), col = NA),
      width = unit(1.5, "cm"),
      axis_param = list(labels_rot = 0, gp = gpar(fontsize = 5)), # 条形图坐标轴数字也同步调小
      border = FALSE
    )
  }
}

if(length(anno_bar_list) > 0) {
  # 【关键修改】：在这里通过 annotation_name_gp 修改底部 Region 标签字体
  right_anno <- do.call(rowAnnotation, c(anno_bar_list, list(
    annotation_name_side = "bottom",
    annotation_name_rot = 45,           # 建议旋转45度，防止字体挤在一起
    annotation_name_gp = gpar(fontsize = 10, fontface = "plain") # 控制底部文字大小
  )))
} else {
  right_anno <- NULL
}


# ==============================================================================
# 7. 渲染终极拓扑可塑性热图 (保持 2.5 上限与 RdYlBu 配色)
# ==============================================================================
message(">>> 正在渲染 S8 终极可塑性热图...")

# 11 级 RdYlBu 色带
cv_colors <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
               "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

# 根据你之前的分位数结果，设为 2.5 依然是最合理的
cv_breaks <- seq(0, 2.5, length.out = 11)
col_fun <- colorRamp2(cv_breaks, cv_colors)

ht_final <- Heatmap(
  plot_matrix,
  name = "Plasticity (CV)",
  col = col_fun,
  
  # 不表达的基因设为高级深灰色
  na_col = "grey40", 
  
  # --- 注释挂载 ---
  top_annotation = top_anno,
  left_annotation = left_anno,
  right_annotation = right_anno, 
  
  # --- 行设置 (Clade 分组) ---
  cluster_rows = FALSE, 
  row_split = row_clades,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8, fontface = "plain"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  
  # --- 列设置 (SubRegion 聚类) ---
  column_split = column_split_factor, 
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,              
  clustering_distance_columns = "euclidean", 
  clustering_method_columns = "ward.D2", 
  
  column_gap = unit(1.5, "mm"),          
  show_column_dend = FALSE,
  show_column_names = FALSE,
  column_title_rot = 45,
  column_title_gp = gpar(fontsize = 8, fontface = "plain"),
  
  border = TRUE,
  use_raster = TRUE
)

# 打印查看效果
ht_final

# ==============================================================================
# 8. 导出高清 PDF
# ==============================================================================
output_pdf <- file.path(base_path, "Figure_S8_Plasticity_Landscape.pdf")
pdf(output_pdf, width = 20, height = 12)
draw(ht_final, merge_legend = TRUE)
dev.off()

message(paste0(">>> 任务圆满完成！PDF 已保存至：\n    ", output_pdf))
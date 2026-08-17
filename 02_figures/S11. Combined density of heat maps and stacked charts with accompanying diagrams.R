# ==============================================================================
# Figure S11 终极完整版: 模块效能全景热图 + 分布统计图 
# (基于 Log-MinMax标准化密度 + 🌟RdYlGn翻转高级绿红发散色带 + RStudio实时预览)
# ==============================================================================

# ==============================================================================
# 0. 环境加载
# ==============================================================================
suppressMessages({
  library(readxl); library(dplyr); library(stringr); library(tibble)
  library(readr); library(ComplexHeatmap); library(circlize)
  library(RColorBrewer); library(purrr); library(tidyr)
  library(ggplot2); library(patchwork); library(grid)
})

# ==============================================================================
# 1. 路径设置与参数
# ==============================================================================
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
PCC_power <- 1 

s10_rds_input <- file.path(base_path, paste0("PCC/Raw_Weighted_7State_Top5PCT_withModularity_PCC", PCC_power, ".RDS"))
Dataset_Mapping_path <- file.path(base_path, "Dataset_Clade_Mapping.txt")
Module_Mapping_path  <- file.path(base_path, "Module_Region_Mapping.txt")

# ==============================================================================
# 2. 提取最原始的物理骨架 (绝不改变任何物理位置和 Order)
# ==============================================================================
message(">>> 1. 正在加载最原始物理骨架...")

dataset_df <- read_tsv(Dataset_Mapping_path, show_col_types = FALSE)
exact_row_order <- dataset_df$Dataset
original_clade_levels <- unique(dataset_df$Clade) 

# 建立原始明文映射字典
clade_dict_raw <- setNames(dataset_df$Clade, dataset_df$Dataset)

# 建立一个只换名字、不换位置的标签映射字典
label_rename_dict <- c("Clade_1" = "Clade_3", "Clade_2" = "Clade_2", "Clade_3" = "Clade_1")

module_df <- read_tsv(Module_Mapping_path, show_col_types = FALSE) %>%
  mutate(Module = str_replace(Module, "^S1_M", "Module_"))
exact_col_order <- module_df$Module
col_subregions <- factor(module_df$SubRegion, levels = unique(module_df$SubRegion))

# ==============================================================================
# 4. 从 RDS 加载效能矩阵 (执行 Log-MinMax 联合标准化)
# ==============================================================================
message(">>> 2. 正在执行 Log 对数压缩与全局 0-1 重新标准化...")
raw_data <- readRDS(s10_rds_input)

df_wide <- raw_data %>%
  select(Dataset, Community, Density_Weighted) %>%
  mutate(Community = str_replace(Community, "^S1_M", "Module_")) %>%
  pivot_wider(names_from = Community, values_from = Density_Weighted, values_fill = 0)

mat_raw <- as.matrix(df_wide[, -1])
rownames(mat_raw) <- df_wide$Dataset

# Log 标准化算法：先 ln(x + 1) 压缩极端高值，再全局归一化
mat_log <- log1p(mat_raw * 15) 
min_v   <- min(mat_log, na.rm = TRUE)
max_v   <- max(mat_log, na.rm = TRUE)
mat_norm <- (mat_log - min_v) / (max_v - min_v)

# 严格对齐物理骨架
existing_ds  <- intersect(exact_row_order, rownames(mat_norm))
existing_mod <- intersect(exact_col_order, colnames(mat_norm))
plot_matrix  <- mat_norm[existing_ds, existing_mod]

# 将无激活的极小值置为 NA 以凸显灰黑色背景
plot_matrix[plot_matrix < 1e-5] <- NA

# ==============================================================================
# 6. 配色与高精热图渲染 
# ==============================================================================
region_base_colors <- c("Region_1"="#7876B1", "Region_2"="#E18727", "Region_3"="#0072B5", "Region_4"="#35A595")
sub_to_reg_lookup <- module_df %>% select(SubRegion, Region) %>% distinct()
subregion_colors <- setNames(region_base_colors[sub_to_reg_lookup$Region], sub_to_reg_lookup$SubRegion)

top_anno <- HeatmapAnnotation(
  SubRegion = as.character(col_subregions[match(existing_mod, exact_col_order)]),
  col = list(SubRegion = subregion_colors), show_annotation_name = FALSE, simple_anno_size = unit(0.3, "cm")
)

aligned_row_clades_raw <- factor(clade_dict_raw[existing_ds], levels = original_clade_levels)
clade_colors <- setNames(c("#D9A441", "#734C7A", "#3D726D"), original_clade_levels)

left_anno <- rowAnnotation(
  Clade = as.character(aligned_row_clades_raw),
  col = list(Clade = clade_colors), show_annotation_name = FALSE, simple_anno_size = unit(0.3, "cm")
)

# 核心调色板变换
cv_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")))(100)
valid_vals <- plot_matrix[!is.na(plot_matrix)]
cv_breaks <- seq(min(valid_vals), max(valid_vals), length.out = 100)
col_fun <- colorRamp2(cv_breaks, cv_colors)

ht_final <- Heatmap(
  plot_matrix, 
  name = "Efficiency Score\n(Log-Norm Density)", 
  col = col_fun, na_col = "grey10",
  top_annotation = top_anno, left_annotation = left_anno,
  cluster_rows = FALSE, 
  row_split = aligned_row_clades_raw, 
  row_title = label_rename_dict[levels(aligned_row_clades_raw)],
  row_title_rot = 0, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_columns = FALSE, column_split = col_subregions[match(existing_mod, exact_col_order)],
  row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
  
  # 🌟 修改点 1：将 column_names_gp 和 row_names_gp 的 fontsize 从 6 放大到 10
  show_column_names = TRUE, column_names_gp = gpar(fontsize = 10), column_names_rot = 45, 
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 10),
  
  column_title_gp = gpar(fontsize = 10, fontface = "bold"), column_title_rot = 45,
  border = TRUE, use_raster = TRUE
)

# ==============================================================================
# 7. 构建基于 Log 标准化数据的各自独立内部比例堆叠图
# ==============================================================================
message(">>> 3. 正在基于 Log 标准化数据计算各 Clade 内部独立比例并执行无偏堆叠...")

df_long <- as.data.frame(plot_matrix) %>%
  rownames_to_column(var = "Dataset") %>%
  pivot_longer(cols = -Dataset, names_to = "Module", values_to = "Score")

df_all <- df_long %>%
  mutate(
    Clade_Raw = clade_dict_raw[Dataset],
    Clade_Display = factor(label_rename_dict[Clade_Raw], levels = c("Clade_1", "Clade_2", "Clade_3"))
  )

clade_totals <- df_all %>% 
  group_by(Clade_Display) %>% 
  summarise(Total_Cells = n(), .groups = "drop")

# --- 图A：静默状态各自内部占比堆叠 ---
df_na_stacked <- df_all %>% 
  filter(is.na(Score) | Score == 0) %>%
  group_by(Clade_Display) %>%
  summarise(NA_Count = n(), .groups = "drop") %>%
  right_join(clade_totals, by = "Clade_Display") %>%
  mutate(NA_Count = replace_na(NA_Count, 0), Prop_Own_Clade = NA_Count / Total_Cells) 

p_na <- ggplot(df_na_stacked, aes(x = "Inactive", y = Prop_Own_Clade, fill = Clade_Display)) +
  geom_col(position = "stack", color = "grey20", linewidth = 0.3, width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c("Clade_1"="#3D726D", "Clade_2"="#734C7A", "Clade_3"="#D9A441")) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() + labs(x = "", y = "Sum of Clade Internal Proportions") +
  theme(axis.text.x = element_text(face = "bold"), legend.position = "none")

# --- 图B：20个区间，每个 Clade 内部占比后再堆叠 ---
breaks_seq_20 <- seq(0, 1, length.out = 21)

df_valid_binned <- df_all %>% 
  filter(!is.na(Score) & Score > 0) %>%
  mutate(Bin = cut(Score, breaks = breaks_seq_20, include.lowest = TRUE, right = FALSE)) %>%
  group_by(Bin, Clade_Display) %>%
  summarise(Bin_Count = n(), .groups = "drop") %>%
  complete(Bin, Clade_Display, fill = list(Bin_Count = 0)) %>%
  filter(!is.na(Bin)) %>%
  left_join(clade_totals, by = "Clade_Display") %>%
  mutate(
    Prop_Own_Clade = Bin_Count / Total_Cells,
    Bin_Mid = as.numeric(str_extract(Bin, "[0-9.]+")) + 0.025
  )

p_dist <- ggplot(df_valid_binned, aes(x = Bin_Mid, y = Prop_Own_Clade, fill = Clade_Display)) +
  geom_col(position = "stack", color = "grey20", linewidth = 0.2, width = 0.046, alpha = 0.9) +
  scale_fill_manual(values = c("Clade_1"="#3D726D", "Clade_2"="#734C7A", "Clade_3"="#D9A441")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), "%"), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() + 
  labs(x = "Efficiency Score (Log-Standardized Density) - 20 Bins", 
       y = "Sum of Clade Internal Proportions", fill = "Environmental Clade")

p_bottom <- p_na + p_dist + plot_layout(widths = c(1, 6)) + 
  plot_annotation(title = "Unbiased Clade Profiling: Log-Standardized Density Stacked Across 20 Intervals",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

# ==============================================================================
# 8. 双层渲染引擎 (先画到 RStudio 窗口展示，再同步导出 PDF)
# ==============================================================================

plot_figure_s11 <- function() {
  grid.newpage()
  # 🌟 修改点 2：将 heights 的比例从 c(3, 1) 修改为 c(5, 1)，大幅度压缩下方柱状图的高度
  pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(5, 1), "null"))))
  
  # 渲染上半区：热图矩阵
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_final, merge_legend = TRUE, newpage = FALSE)
  upViewport()
  
  # 渲染下半区：无偏自我占比堆积图
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p_bottom, newpage = FALSE)
  upViewport()
}

# --- 步骤 1：直接在 RStudio Plots 画布上激活显示 ---
message(">>> 4. [RStudio 优先渲染] 正在把大图打在屏幕上，请查看 RStudio 右下角 Plots 面板预览...")
plot_figure_s11()

# --- 步骤 2：建立本地 PDF 管道无损写盘 ---
message(">>> 5. [PDF 设备同步] 正在向本地路径无损写盘导出...")
output_pdf <- file.path(base_path, paste0("Figure_S11_Ultimate_Landscape_PCC", PCC_power, "_LogNorm_RdYlGn_Style.pdf"))

pdf(output_pdf, width = 20, height = 16)
plot_figure_s11()
dev.off()

message(paste0(">>> [大功告成] RStudio 实时高亮画面与本地 PDF 文件均已完美就绪！\nPDF 已保存至：\n ", output_pdf))



# ==============================================================================
# 提取 Figure S11 核心 Density 数据并整合保留分类骨架信息 (Clade, Region, SubRegion)
# ==============================================================================

message(">>> 正在启动带完整分类信息的 11 密度数据提取程序...")

# 1. 设定输入与输出路径
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
Dataset_Mapping_path <- file.path(base_path, "Dataset_Clade_Mapping.txt")
Module_Mapping_path  <- file.path(base_path, "Module_Region_Mapping.txt")

output_dir  <- file.path(base_path, "paper图片")
output_file <- file.path(output_dir, "Figure_S11_Density_with_Metadata.csv")

# 如果输出目录不存在，则自动创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. 加载行与列的物理骨架元数据
message(">>> 正在读取元数据映射表...")
dataset_meta <- read_tsv(Dataset_Mapping_path, show_col_types = FALSE) %>%
  select(Dataset, Clade)

module_meta  <- read_tsv(Module_Mapping_path, show_col_types = FALSE) %>%
  mutate(Module_Clean = str_replace(Module, "^S1_M", "Module_")) %>%
  select(Module_Clean, Region, SubRegion)

# 3. 提取原始 RDS 密度数据并进行全信息打标签流式整合
message(">>> 正在融合 Density 数值与元数据骨架...")

# a. 将原始矩阵转换为标准长表
df_long_raw <- raw_data %>%
  select(Dataset, Community, Density_Weighted) %>%
  mutate(Module_Clean = str_replace(Community, "^S1_M", "Module_")) %>%
  select(-Community)

# b. 左右双向无损连接元数据
final_density_output <- df_long_raw %>%
  # 内联样本对应的进化环境分群 (Clade)
  left_join(dataset_meta, by = "Dataset") %>%
  # 内联模块对应的空间网络分区 (Region & SubRegion)
  left_join(module_meta, by = "Module_Clean") %>%
  # 规范列名并重排顺序，让元数据走在前面，数值列走在最后
  rename(Module = Module_Clean, Density = Density_Weighted) %>%
  select(Dataset, Clade, Module, Region, SubRegion, Density)

# 4. 写入本地磁盘
write_csv(final_density_output, output_file)

# ==============================================================================
# 验证与大功告成提示
# ==============================================================================
if (file.exists(output_file)) {
  message(">>> [成功] 包含 Clade/Region/SubRegion 完整信息的数据已成功导出！")
  message(">>> 文件保存路径为：")
  message("    ", output_file)
  
  # 打印数据摘要，方便实时核对
  total_records <- nrow(final_density_output)
  unique_ds     <- length(unique(final_density_output$Dataset))
  unique_mod    <- length(unique(final_density_output$Module))
  
  message(paste0(">>> [数据摘要] 成功生成总计 ", total_records, " 行的标准长表数据。"))
  message(paste0("    (覆盖了 ", unique_ds, " 个不同数据集与 ", unique_mod, " 个核心模块)"))
  message(">>> [表格列架构预阅]: Dataset | Clade | Module | Region | SubRegion | Density")
} else {
  warning(">>> [错误] 文件写入失败，请检查路径或权限。")
}
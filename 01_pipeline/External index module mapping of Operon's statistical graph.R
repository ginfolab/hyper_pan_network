# ==============================================================================
# 1. 加载包和设置路径
# ==============================================================================
library(igraph)
library(dplyr)
library(tidyverse)
library(stringr)

# 基础路径
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"

# ==============================================================================
# 2. 准备 "标准答案" (只保留 Operon)
# ==============================================================================
# 读取数据库
Function_term_table_raw <- readRDS(file.path(base_path, "RegulonDB/最新的RegulonDB基因注释.RDS"))

# 【关键修改】：只保留 Operon
Operon_table <- Function_term_table_raw %>%
  filter(Type == "Operon")

# ==============================================================================
# 3. 预读取网络数据 (在循环外读取，极大提高速度)
# ==============================================================================
# 读取原始边列表

Single_edge_001_clean <- readRDS(file.path(base_path, "PCC/final_old001_method_gene_pair_U_clean.RDS"))
Single_edge_005_clean <- readRDS(file.path(base_path, "PCC/final_old005_method_gene_pair_U_clean.RDS"))
Hyper_edge_clean <- readRDS(file.path(base_path, "PCC/final_Hyper_edge_U_clean.RDS"))



# ==============================================================================
# 4. 定义通用的计算函数
# ==============================================================================
calculate_module_scores <- function(edge_data, network_name, u_val, reference_table) {
  
  # 1. 建图与聚类
  g <- graph_from_data_frame(edge_data, directed = FALSE)
  if (ecount(g) == 0) return(data.frame())
  
  E(g)$Weight_s <- E(g)$Universality^2
  clusters <- cluster_louvain(g, resolution = 1, weights = E(g)$Weight_s)
  
  # 2. 整理模块信息
  module_df <- data.frame(id = clusters$names, membership = clusters$memberships[1,]) %>%
    mutate(module_info = paste0("Module_", membership)) %>%
    group_by(module_info) %>%
    summarise(Genes = paste(id, collapse = ", "), .groups = "drop") %>%
    mutate(Module_size = str_count(Genes, ",") + 1)
  
  # 3. 与 Reference (Operon) 进行匹配计算
  matched_df <- expand_grid(module_df, reference_table) %>%
    mutate(
      overlap_genes = map2_chr(
        str_split(Genes, ", "),
        str_split(Term_clustid, ", "),
        ~ intersect(.x, .y) %>% paste(collapse = ", ")
      )
    ) %>%
    filter(overlap_genes != "") %>% 
    mutate(overlap_genes_count = str_count(overlap_genes, ",") + 1) %>%
    mutate(
      Precision = overlap_genes_count / Module_size,
      Recall    = overlap_genes_count / Term_size,
      F_score   = if_else(Precision + Recall == 0, 0, 2 * Precision * Recall / (Precision + Recall))
    ) %>%
    # 筛选条件：保留 F_score >= 0.5
    filter(F_score >= 0.5) %>%
    select(module_info, Function_term, Precision, Recall, F_score, Module_size, Term_size) %>%
    mutate(
      Network = network_name,
      Universality_cutoff = u_val
    )
  
  return(matched_df)
}

# ==============================================================================
# 5. 主循环 (Universality Loop)
# ==============================================================================
start_U = 5
end_U = 70
results_list <- list()

message("开始循环计算...")
pb <- txtProgressBar(min = start_U, max = end_U, style = 3)

for (U_cut in start_U:end_U) {
  
  # --- 1. Top_0.01 (原 S1) ---
  current_edge_s1 <- Single_edge_001_clean %>% filter(Universality >= U_cut)
  if(nrow(current_edge_s1) > 0) {
    # 【修改】：直接使用新名称 "Top_0.01"
    res_s1 <- calculate_module_scores(current_edge_s1, "Top_0.01", U_cut, Operon_table)
    results_list[[paste0("S1_", U_cut)]] <- res_s1
  }
  
  # --- 2. Top_0.05 (原 S5) ---
  current_edge_s5 <- Single_edge_005_clean %>% filter(Universality >= U_cut)
  if(nrow(current_edge_s5) > 0) {
    # 【修改】：直接使用新名称 "Top_0.05"
    res_s5 <- calculate_module_scores(current_edge_s5, "Top_0.05", U_cut, Operon_table)
    results_list[[paste0("S5_", U_cut)]] <- res_s5
  }
  
  # --- 3. hyeperedge-based (原 Hyper) ---
  current_edge_hyper <- Hyper_edge_clean %>% filter(Universality >= U_cut)
  if(nrow(current_edge_hyper) > 0) {
    # 【修改】：直接使用新名称 "hyeperedge-based"
    res_hyper <- calculate_module_scores(current_edge_hyper, "hyeperedge-based", U_cut, Operon_table)
    results_list[[paste0("Hyper_", U_cut)]] <- res_hyper
  }
  
  setTxtProgressBar(pb, U_cut)
}
close(pb)

# 合并结果
final_operon_statistics <- bind_rows(results_list)

# 保存中间结果 (备份)
output_file <- file.path(base_path, paste0("Operon_Detail_Scores_U", start_U, "_", end_U, ".RDS"))
saveRDS(final_operon_statistics, output_file)
message("统计计算完成，已保存至: ", output_file)

# ==============================================================================
# 6. 数据准备：分为“质量”和“数量”两组
# ==============================================================================

# 确保 Network 因子顺序正确 (用于图例排序)
plot_df <- final_operon_statistics %>%
  mutate(Network = factor(Network, levels = c("Top_0.01", "Top_0.05", "hyeperedge-based")))

# --- A. 质量数据 (Precision & Recall) ---
quality_data <- plot_df %>%
  select(Network, Universality_cutoff, Precision, Recall) %>% 
  pivot_longer(
    cols = c("Precision", "Recall"),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Precision", "Recall")))

# --- B. 数量数据 (Mapped Count) ---
quantity_data <- plot_df %>%
  group_by(Network, Universality_cutoff) %>%
  summarise(Mapped_Count = n(), .groups = 'drop')

# ==============================================================================
# 7. 绘图 (Left: P/R, Right: Count)
# ==============================================================================

# 定义颜色 (你指定的颜色)
my_colors <- c(
  "Top_0.01"         = "#1d9bf7",
  "Top_0.05"         = "navy",
  "hyeperedge-based" = "#ff7f0e"
)

# 通用主题设置
common_theme <- theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 11)
  )

# --- 绘制左图 (Precision & Recall) ---
p_left <- ggplot(quality_data, aes(x = Universality_cutoff, y = Score, color = Network, fill = Network)) +
  # 误差带 (95% CI)
  stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", alpha = 0.2, color = NA) +
  # 均值线
  stat_summary(geom = "line", fun = "mean", size = 1) +
  # U=15 深绿色虚线
  geom_vline(xintercept = 15, linetype = "dashed", color = "darkgreen", size = 0.8) +
  # 分面展示 Precision 和 Recall
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  # 颜色应用
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "A. Quality: Precision & Recall",
    y = "Score (Mean +/- 95% CI)",
    x = "Universality Cutoff (U)"
  ) +
  common_theme

# --- 绘制右图 (Mapped Count) ---
p_right <- ggplot(quantity_data, aes(x = Universality_cutoff, y = Mapped_Count, color = Network)) +
  # 折线
  geom_line(size = 1.2) +
  # U=15 深绿色虚线
  geom_vline(xintercept = 15, linetype = "dashed", color = "darkgreen", size = 0.8) +
  # 颜色应用
  scale_color_manual(values = my_colors) +
  labs(
    title = "B. Quantity: Mapped Count",
    subtitle = "(Modules with F-score > 0.5)",
    y = "Number of Mapped Operons",
    x = "Universality Cutoff (U)"
  ) +
  common_theme

# ==============================================================================
# 8. 拼图与保存
# ==============================================================================
# 左右拼接，共用图例
final_plot <- p_left + p_right +
  plot_layout(widths = c(1, 1), guides = "collect") & 
  theme(legend.position = "top")

print(final_plot)
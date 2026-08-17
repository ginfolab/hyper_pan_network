#####数据制作####
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli")

library(igraph)
library(dplyr)
library(stringr)

cor_cut <- c("001","005")

# 初始化最终表格
Final_statistics <- data.frame(
  data_set = character(),
  U = numeric(),
  Number_of_edges = numeric(),
  Number_of_nodes = numeric(),
  Average_degree = numeric(),
  Modularity = numeric(),
  Clustering_coefficient = numeric(),
  Density = numeric(),
  ScaleFree_R2 = numeric(),
  Number_of_communities = numeric(),       # 新增
  Community_size_mean = numeric(),         # 新增
  Community_size_median = numeric()        # 新增
)

# ----------------- 单边网络循环 -----------------
for (i in seq_along(cor_cut)) {
  
  data_U_cut <- readRDS(paste0("./PCC/final_old", cor_cut[i], "_method_gene_pair_U.RDS")) %>%
    mutate(Genepair = str_replace_all(Genepair, "\\b(\\w+)\\b", "Ec_\\1")) %>%
    rename(gene_pairs = Genepair)
  
  colnames(data_U_cut) <- c("gene_pairs", "Universality")
  
  U_set <- sort(unique(data_U_cut$Universality))
  
  for (u in U_set) {
    
    data_U <- data_U_cut[data_U_cut$Universality >= u, ]
    edge_list <- do.call(rbind, strsplit(as.character(data_U$gene_pairs), ","))
    edge_list <- as.data.frame(edge_list, stringsAsFactors = FALSE)
    colnames(edge_list) <- c("from", "to")
    
    graph <- graph_from_data_frame(edge_list, directed = FALSE)
    
    dataName <- paste0("single_net_CorCut", cor_cut[i])
    
    Number_of_edges <- ecount(graph)
    Number_of_nodes <- vcount(graph)
    Average_degree <- mean(degree(graph))
    
    # 社区检测
    community <- cluster_louvain(graph)
    Modularity <- modularity(community)
    Clustering_coefficient <- transitivity(graph, type = "global")
    Density <- edge_density(graph)
    
    # 计算社区信息
    community_sizes <- sizes(community)  # 各社区节点数
    Number_of_communities <- length(community_sizes)
    Community_size_mean <- mean(community_sizes)
    Community_size_median <- median(community_sizes)
    
    # ----------------- scale-free R² -----------------
    deg <- degree(graph)
    deg_table <- as.data.frame(table(deg))
    deg_table$deg <- as.numeric(as.character(deg_table$deg))
    deg_table$Freq <- as.numeric(deg_table$Freq)
    deg_table <- deg_table[deg_table$deg > 0, ]
    
    if(nrow(deg_table) > 1){
      fit <- lm(log10(Freq) ~ log10(deg), data = deg_table)
      ScaleFree_R2 <- summary(fit)$r.squared
    } else {
      ScaleFree_R2 <- NA
    }
    
    # ----------------- 保存到表格 -----------------
    current_table <- data.frame(
      data_set = dataName,
      U = u,
      Number_of_edges = Number_of_edges,
      Number_of_nodes = Number_of_nodes,
      Average_degree = Average_degree,
      Modularity = Modularity,
      Clustering_coefficient = Clustering_coefficient,
      Density = Density,
      ScaleFree_R2 = ScaleFree_R2,
      Number_of_communities = Number_of_communities,
      Community_size_mean = Community_size_mean,
      Community_size_median = Community_size_median
    )
    
    Final_statistics <- rbind(Final_statistics, current_table)
  }
}

# ----------------- 超边网络循环 -----------------
hyper_U_original <- readRDS("./Escherichia_coli_2items_AP0.5PCC_Uall.RDS")[, -3]
colnames(hyper_U_original) <- c("gene_pairs", "Universality")
U_set <- sort(unique(hyper_U_original$Universality))

for (u in U_set) {
  
  data_U <- hyper_U_original[hyper_U_original$Universality >= u, ]
  edge_list <- do.call(rbind, strsplit(as.character(data_U$gene_pairs), ","))
  edge_list <- as.data.frame(edge_list, stringsAsFactors = FALSE)
  colnames(edge_list) <- c("from", "to")
  
  graph <- graph_from_data_frame(edge_list, directed = FALSE)
  dataName <- "Hyperedge_net"
  
  Number_of_edges <- ecount(graph)
  Number_of_nodes <- vcount(graph)
  Average_degree <- mean(degree(graph))
  
  community <- cluster_louvain(graph)
  Modularity <- modularity(community)
  Clustering_coefficient <- transitivity(graph, type = "global")
  Density <- edge_density(graph)
  
  # 社区信息
  community_sizes <- sizes(community)
  Number_of_communities <- length(community_sizes)
  Community_size_mean <- mean(community_sizes)
  Community_size_median <- median(community_sizes)
  
  # scale-free R²
  deg <- degree(graph)
  deg_table <- as.data.frame(table(deg))
  deg_table$deg <- as.numeric(as.character(deg_table$deg))
  deg_table$Freq <- as.numeric(deg_table$Freq)
  deg_table <- deg_table[deg_table$deg > 0, ]
  
  if(nrow(deg_table) > 1){
    fit <- lm(log10(Freq) ~ log10(deg), data = deg_table)
    ScaleFree_R2 <- summary(fit)$r.squared
  } else {
    ScaleFree_R2 <- NA
  }
  
  current_table <- data.frame(
    data_set = dataName,
    U = u,
    Number_of_edges = Number_of_edges,
    Number_of_nodes = Number_of_nodes,
    Average_degree = Average_degree,
    Modularity = Modularity,
    Clustering_coefficient = Clustering_coefficient,
    Density = Density,
    ScaleFree_R2 = ScaleFree_R2,
    Number_of_communities = Number_of_communities,
    Community_size_mean = Community_size_mean,
    Community_size_median = Community_size_median
  )
  
  Final_statistics <- rbind(Final_statistics, current_table)
}

# ----------------- 保存最终结果 -----------------
saveRDS(Final_statistics, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Final_statistics.RDS")
##超边网络


######绘制Figure 2
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
#画图
# 加载必要的包
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(scales) # 用于生成渐变色
Final_statistics <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Final_statistics.RDS")
# 将数据转换为长格式以便于 ggplot2 处理
long_data <- Final_statistics %>%
  
  pivot_longer(
    cols = c(Number_of_nodes, Number_of_edges,  Average_degree, Modularity, Clustering_coefficient ), # 添加 Number_of_nodes
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(value > 0 & !is.na(value)) %>%
  filter(metric %in% c(
    # "Number_of_edges", "Number_of_nodes", 
    "Modularity", "Clustering_coefficient")) %>%
  mutate(metric = factor(metric, 
                         levels = c(
                           "Number_of_nodes", 
                           "Number_of_edges", 
                           "Average_degree",
                           "Density",
                           "ScaleFree_R2"
                         )))
# 自动检测 data_set 是 single_net 还是 hyperedge_net
long_data$data_group <- ifelse(grepl("single_net", long_data$data_set), "single_net", "hyperedge_net")
unique_data_sets <- unique(long_data$data_set)

# 为 single_net 和 hyperedge_net 分别生成颜色梯度
num_single_net <- length(unique_data_sets[grepl("single_net", unique_data_sets)])
num_hyperedge_net <- length(unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])

single_net_colors <- colorRampPalette(c("#1d9bf7", "navy"))(num_single_net)
hyperedge_net_colors <- colorRampPalette(c("#ff7f0e"))(num_hyperedge_net)

# 创建颜色映射表
color_mapping <- c(
  setNames(single_net_colors, unique_data_sets[grepl("single_net", unique_data_sets)]),
  setNames(hyperedge_net_colors, unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])
)
library(patchwork)

# 绘图
library(dplyr)

# 筛选并按节点数量从大到小排序（从右往左）
hyperedge_data <- Final_statistics %>%
  filter(data_set == "Hyperedge_net") %>%
  arrange(desc(Number_of_nodes))  # 从右往左

# 计算差分（近似一阶导数）
hyperedge_data <- hyperedge_data %>%
  mutate(delta_mod = c(NA, diff(Modularity)))

# 设置一个“平稳”阈值（模块度变化非常小即视为平稳）
threshold <- 0.01

# 找到第一个满足“变化幅度小于阈值”的点
plateau_index <- which(abs(hyperedge_data$delta_mod) < threshold)[1]
plateau_point <- hyperedge_data[plateau_index, ]

# 输出最优 U 值
cat("最早达到平衡的 U 为：", plateau_point$U, "\n")

# ==================== 1. 数据准备 (保持不变) ====================
long_data <- Final_statistics %>%
  pivot_longer(
    cols = c(Number_of_edges, Number_of_nodes, Average_degree, Modularity, Clustering_coefficient),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::filter(value > 0 & !is.na(value))

long_data$data_group <- ifelse(grepl("single_net", long_data$data_set), "single_net", "hyperedge_net")
unique_data_sets <- unique(long_data$data_set)
num_single_net <- length(unique_data_sets[grepl("single_net", unique_data_sets)])
num_hyperedge_net <- length(unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])
single_net_colors <- colorRampPalette(c("#1d9bf7", "navy"))(num_single_net)
hyperedge_net_colors <- colorRampPalette(c("#ff7f0e"))(num_hyperedge_net)
color_mapping <- c(
  setNames(single_net_colors, unique_data_sets[grepl("single_net", unique_data_sets)]),
  setNames(hyperedge_net_colors, unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])
)

# ==================== 2. 主题设置 (关键修改) ====================

# 修改说明：
# 1. plot.title face="plain": 为了让 expression 里的 regular 生效
# 2. axis.title face="plain": 将横纵轴标签改为 Regular
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "plain", hjust = 0), 
    plot.title.position = "plot", 
    axis.title.x = element_text(size = 14, face = "plain"), # 改为 plain
    axis.title.y = element_text(size = 14, face = "plain"), # 改为 plain
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.position = "none"
  )

# ==================== 3. 绘图 ====================

# --- 图 A: Modularity ---
p_modularity <- long_data %>%
  filter(metric == "Modularity") %>%
  ggplot(aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  labs(
    x = "Universality change",
    y = "", 
    title = expression(paste(bold("A"), "   Modularity")) 
  ) +
  my_theme

# --- 图 B: Clustering_coefficient ---
p_clustering <- long_data %>%
  filter(metric == "Clustering_coefficient") %>%
  ggplot(aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  labs(
    x = "Universality change",
    y = "",
    title = expression(paste(bold("B"), "   Clustering_coefficient"))
  ) +
  my_theme

# --- 图 C: U where Modularity First Stabilizes ---
p_bottom <- ggplot(Final_statistics, aes(x = Number_of_nodes, y = Modularity, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(data = head(Final_statistics[Final_statistics$data_set == "Hyperedge_net", ], 10), 
            aes(label = paste("U =", U)), 
            vjust = -1, size = 3, color = "grey60") +
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = plateau_point$Number_of_nodes, linetype = "dashed", color = "darkgreen") +
  geom_point(data = plateau_point, aes(x = Number_of_nodes, y = Modularity), color = "darkgreen", size = 3) +
  annotate("text", x = plateau_point$Number_of_nodes, y = plateau_point$Modularity + 0.05, 
           label = paste("Stable U =", plateau_point$U), color = "darkgreen", size = 5, fontface = "bold") +
  labs(
    x = "Number of Nodes",
    y = "Modularity",
    color = "Data Set",
    title = expression(paste(bold("C"), "   U where Modularity First Stabilizes (from Right to Left)"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "plain", hjust = 0),
    plot.title.position = "plot",
    legend.position = "bottom", 
    legend.title = element_text(size = 12, face = "bold"), # 图例标题通常保留粗体好看，如果也要regular请改为plain
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "plain"), # 改为 plain
    axis.title.y = element_text(size = 14, face = "plain"), # 改为 plain
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

# ==================== 4. 组合输出 ====================

layout <- (p_modularity | p_clustering) / p_bottom
layout + plot_layout(heights = c(1, 1.2))








##################Figure S3####
library(tidyverse)
library(patchwork)

# ==================== 1. 数据准备 ====================
# 筛选需要的指标
long_data_aux <- Final_statistics %>%
  pivot_longer(
    cols = c(Density, ScaleFree_R2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(value > 0 & !is.na(value)) %>%
  mutate(metric = factor(metric, levels = c("Density", "ScaleFree_R2")))

# 设置颜色映射 (保持一致)
long_data_aux$data_group <- ifelse(grepl("single_net", long_data_aux$data_set), "single_net", "hyperedge_net")
unique_data_sets <- unique(long_data_aux$data_set)
num_single_net <- sum(grepl("single_net", unique_data_sets))
num_hyperedge_net <- sum(grepl("Hyperedge_net", unique_data_sets))

single_net_colors <- colorRampPalette(c("#1d9bf7", "navy"))(num_single_net)
hyperedge_net_colors <- colorRampPalette(c("#ff7f0e"))(num_hyperedge_net)

color_mapping <- c(
  setNames(single_net_colors, unique_data_sets[grepl("single_net", unique_data_sets)]),
  setNames(hyperedge_net_colors, unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])
)

# ==================== 2. 通用主题设置 ====================
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "plain", hjust = 0), # 标题左对齐，常规字体
    plot.title.position = "plot", 
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.position = "none" # 默认不显示图例，只在最后一张图显示
  )

# ==================== 3. 绘图 ====================

# --- 图 A: Density ---
p_density <- long_data_aux %>%
  filter(metric == "Density") %>%
  ggplot(aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() +
  scale_color_manual(values = color_mapping) +
  # 添加辅助线
  geom_vline(xintercept = 15, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = 15, y = Inf, label = "U = 15", color = "darkgreen", 
           vjust = 2, hjust = -0.1, size = 4) +
  labs(
    x = "Universality change",
    y = "",
    title = expression(paste(bold("A"), "   Density"))
  ) +
  my_theme

# --- 图 B: ScaleFree_R2 ---
p_scalefree <- long_data_aux %>%
  filter(metric == "ScaleFree_R2") %>%
  ggplot(aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() +
  scale_color_manual(values = color_mapping) +
  # 添加辅助线
  geom_vline(xintercept = 15, linetype = "dashed", color = "darkgreen") +
  # ScaleFree 这边如果不需要文字标签可以注释掉下面这句，或者保留以保持对称
  # annotate("text", x = 15, y = Inf, label = "U = 15", color = "darkgreen", vjust = 2, hjust = -0.1, size = 4) +
  labs(
    x = "Universality change",
    y = "",
    title = expression(paste(bold("B"), "   ScaleFree_R2"))
  ) +
  my_theme

# --- 图 C: Nodes vs Communities ---
p_communities <- ggplot(Final_statistics, aes(x = Number_of_nodes, y = Number_of_communities, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  # 标记 U 值
  geom_text(data = head(Final_statistics[Final_statistics$data_set == "Hyperedge_net", ], 10), 
            aes(label = paste("U =", U)), 
            vjust = -1, size = 2.5, color = "grey60") +
  scale_color_manual(values = color_mapping) +
  # 辅助线 (来自 plateau_point)
  geom_vline(xintercept = plateau_point$Number_of_nodes, linetype = "dashed", color = "darkgreen") +
  
  labs(
    x = "Number of Nodes",
    y = "Number of Communities",
    color = "Data Set",
    title = expression(paste(bold("C"), "   Nodes vs. Communities under changing Universality"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "plain", hjust = 0),
    plot.title.position = "plot",
    legend.position = "bottom", # 图例放底部
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

# ==================== 4. 组合输出 ====================

layout_aux <- (p_density | p_scalefree) / p_communities

# 设置高度比例 (上1 下1.5)
layout_aux + plot_layout(heights = c(1, 1.5))

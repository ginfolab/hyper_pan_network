######################################绘制网络特性比较折线图########################################
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


# 保存两个 ggplot 对象
p1 <- ggplot(long_data, aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, ncol = 2) +     
  scale_y_log10() +
  scale_color_manual(values = color_mapping) +
  labs(
    x = "Universality change",
    y = "",
    color = "Data Set",
    title = "Network characteristics statistics"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "None",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p2 <- ggplot(Final_statistics, aes(x = Number_of_nodes, y = Modularity, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(data = head(Final_statistics[Final_statistics$data_set == "Hyperedge_net", ], 10), 
            aes(label = paste("U =", U)), 
            vjust = -1, size = 2.5, color = "grey60") +
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = plateau_point$Number_of_nodes, linetype = "dashed", color = "darkgreen") +
  geom_point(data = plateau_point, aes(x = Number_of_nodes, y = Modularity), color = "darkgreen", size = 3) +
  annotate("text", x = plateau_point$Number_of_nodes, y = plateau_point$Modularity + 0.03, 
           label = paste("Stable U =", plateau_point$U), color = "darkgreen", size = 4) +
  labs(
    x = "Number of Nodes",
    y = "Modularity",
    color = "Data Set",
    title = "U where Modularity First Stabilizes (from Right to Left)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# 三行两列组合
p1 / p2 + plot_layout(heights = c(1, 2))





######绘制辅助图  固定辅助线
library(tidyverse)

long_data <- Final_statistics %>%
  pivot_longer(
    cols = c(
      Average_degree,
      Modularity,
      Clustering_coefficient,
      Density,
      ScaleFree_R2
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(value > 0 & !is.na(value)) %>%
  filter(metric %in% c("Density", "ScaleFree_R2")) %>%
  mutate(metric = factor(metric, levels = c("Density", "ScaleFree_R2")))

# 为分组标签
long_data$data_group <- ifelse(grepl("single_net", long_data$data_set),
                               "single_net", "hyperedge_net")

unique_data_sets <- unique(long_data$data_set)

# 颜色
num_single_net <- sum(grepl("single_net", unique_data_sets))
num_hyperedge_net <- sum(grepl("Hyperedge_net", unique_data_sets))

single_net_colors <- colorRampPalette(c("#1d9bf7", "navy"))(num_single_net)
hyperedge_net_colors <- colorRampPalette(c("#ff7f0e"))(num_hyperedge_net)

color_mapping <- c(
  setNames(single_net_colors, unique_data_sets[grepl("single_net", unique_data_sets)]),
  setNames(hyperedge_net_colors, unique_data_sets[grepl("Hyperedge_net", unique_data_sets)])
)


# 绘图
p3 <- ggplot(long_data, aes(x = U, y = value, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, ncol = 2) +
  scale_y_log10() +
  scale_color_manual(values = color_mapping) +
  
  # ---- 这里是你要的虚线和标签 ----
geom_vline(xintercept = 15, linetype = "dashed", color = "darkgreen") +
  annotate(
    "text",
    x = 15, y = Inf,
    label = "U = 15",
    color = "darkgreen",
    vjust = 6.5,
    hjust = -0.1,
    size = 4
  ) +
  
  
  labs(
    x = "Universality change",
    y = "",
    title = "Network characteristics statistics"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

p4 <- ggplot(Final_statistics, aes(x = Number_of_nodes, y = Number_of_communities, color = data_set, group = data_set)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(data = head(Final_statistics[Final_statistics$data_set == "Hyperedge_net", ], 10), 
            aes(label = paste("U =", U)), 
            vjust = -1, size = 2.5, color = "grey60")+
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = plateau_point$Number_of_nodes, linetype = "dashed", color = "darkgreen") +
  
  labs(
    x = "Number of Nodes",
    y = "Number of Communities",
    color = "Data Set",
    title = "Nodes vs. Communities under changing Universality "
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# 三行两列组合
p3 / p4 + plot_layout(heights = c(1, 2))


######################################绘制网络特性比较折线图########################################
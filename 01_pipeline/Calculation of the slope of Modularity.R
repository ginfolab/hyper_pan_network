FDR_THRESHOLD <- 0.05    
TOP_PROP      <- 0.05   
suffix <- paste0("FDR", gsub("\\.", "", as.character(FDR_THRESHOLD)), 
                 "_Top", gsub("\\.", "", as.character(TOP_PROP)))
rds_dir <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/"
rds_path <- paste0(rds_dir, "Community_Performance_", suffix, ".RDS")


Module_size <- read.table('/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新修改label和颜色).txt', 
                          header = TRUE, 
                          sep = "\t", 
                          quote = "", 
                          comment.char = "", 
                          stringsAsFactors = FALSE) %>%
  select(Community = Node_name, Size)


all_network_results <- readRDS(rds_path) %>%
  mutate(Community = str_replace(Community, "^S1_M", "Module_")) %>%
  filter(!str_detect(Community, "Group")) %>% # 仅保留有效模块
  left_join(Module_size, by = "Community")



# 加载所需的包
library(ggplot2)
library(dplyr)
library(scales)
library(ggpubr)
library(patchwork)

# 确保数据已过滤掉 NA
plot_data <- all_network_results %>%
  filter(!is.na(Size) & !is.na(Q_c)) 

# ---------------------------------------------------------
# P1: 原始 Size vs 原始 Q_c (不取 Log，不平方)
# ---------------------------------------------------------
p1 <- ggplot(plot_data, aes(x = Size, y = Q_c)) +
  geom_point(alpha = 0.5, color = "#1f77b4", size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
  geom_smooth(method = "loess", color = "#2ca02c", se = TRUE) +
  theme_bw() +
  labs(title = "Size vs Q_c (Raw)", 
       x = "Module Size", y = "Q_c") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5)

# ---------------------------------------------------------
# P2: 原始 Size vs 对数 Q_c (取 Pseudo-log，不平方)
# ---------------------------------------------------------
p2 <- ggplot(plot_data, aes(x = log(Size), y = log(Q_c))) +
  geom_point(alpha = 0.5, color = "#1f77b4", size = 1.5) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(-1, 0, 1, 10, 100, 1000, 5000)) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
  geom_smooth(method = "loess", color = "#2ca02c", se = TRUE) +
  theme_bw() +
  labs(title = "Module Size(log) vs Q_c (log)", 
       x = "Module Size (log)", y = "Q_c (log)") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5)

# ---------------------------------------------------------
# P3: 平方 Size^2 vs 原始 Q_c (不取 Log，平方)
# ---------------------------------------------------------
p3 <- ggplot(plot_data, aes(x = Size^2, y = Q_c)) +
  geom_point(alpha = 0.5, color = "#ff7f0e", size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
  geom_smooth(method = "loess", color = "#2ca02c", se = TRUE) +
  theme_bw() +
  labs(title = "Size^2 vs Q_c (Raw)", 
       x = bquote(Size^2), y = "Q_c") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5)

# ---------------------------------------------------------
# P4: 平方 Size^2 vs 对数 Q_c (取 Pseudo-log，平方)
# ---------------------------------------------------------
p4 <- ggplot(plot_data, aes(x = log(Size^2), y = log(Q_c))) +
  geom_point(alpha = 0.5, color = "#ff7f0e", size = 1.5) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(-1, 0, 1, 10, 100, 1000, 5000)) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
  geom_smooth(method = "loess", color = "#2ca02c", se = TRUE) +
  theme_bw() +
  labs(title = "Size^2(log) vs Q_c (log)", 
       x = paste0(bquote(Size^2), "(log)"), y = "Q_c (log)") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5)

# ---------------------------------------------------------
# 组合 4 张图 (2x2 布局)
# ---------------------------------------------------------
combined_plot <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(tag_levels = 'A', # 自动加上 A, B, C, D 标签
                  title = "Comprehensive Comparison: Size vs Modularity (Q_c)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

# 打印最终结果
print(combined_plot)



# 加载必要的包
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 确保计算了两个新的比值指标
plot_data <- plot_data %>%
  mutate(
    Ratio_Size = Q_c / Size,
    Ratio_Size2 = Q_c / (Size^2)
  )

# ---------------------------------------------------------
# 2. 绘制原始 Q_c 的直方图 (新增加在最上方的图)
# ---------------------------------------------------------
h0 <- ggplot(plot_data, aes(x = Q_c)) +
  geom_histogram(bins = 50, fill = "#2ca02c", color = "black", alpha = 0.7) +
  stat_bin(bins = 50, geom = "text", 
           aes(label = ifelse(after_stat(count) > 0, after_stat(count), "")), 
           vjust = -0.5, size = 3, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(title = "Distribution of Q_c",
       subtitle = "Raw modularity score",
       x = "Q_c (Modularity)", 
       y = "Frequency (Count)")

# ---------------------------------------------------------
# 3. 绘制 Q_c / Size 的直方图
# ---------------------------------------------------------
h1 <- ggplot(plot_data, aes(x = Ratio_Size)) +
  geom_histogram(bins = 50, fill = "#1f77b4", color = "black", alpha = 0.7) +
  stat_bin(bins = 50, geom = "text", 
           aes(label = ifelse(after_stat(count) > 0, after_stat(count), "")), 
           vjust = -0.5, size = 3, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(title = "Distribution of Q_c / Size",
       subtitle = "Average modularity contribution per node",
       x = "Q_c / Size", 
       y = "Frequency (Count)")

# ---------------------------------------------------------
# 4. 绘制 Q_c / Size^2 的直方图
# ---------------------------------------------------------
h2 <- ggplot(plot_data, aes(x = Ratio_Size2)) +
  geom_histogram(bins = 50, fill = "#ff7f0e", color = "black", alpha = 0.7) +
  stat_bin(bins = 50, geom = "text", 
           aes(label = ifelse(after_stat(count) > 0, after_stat(count), "")), 
           vjust = -0.5, size = 3, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(title = "Distribution of Q_c / Size^1.4439",
       subtitle = "Normalized modularity density",
       x = bquote("Q_c /"~Size^1.4439), 
       y = "Frequency (Count)")

# ---------------------------------------------------------
# 5. 合并并显示图形 (上、中、下排列)
# ---------------------------------------------------------
# 直接使用 / 符号将三个图从上到下排列
combined_hist <- h0 / h1 / h2 + 
  plot_annotation(tag_levels = 'A',
                  title = "Distributions of Raw and Normalized Modularity (with Counts)")

print(combined_hist)





# ==============================================================================
# 完整代码：模块大小与模块度 Q_c 的对数幂律拟合与可视化
# ==============================================================================

# 加载必要的包
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

# 1. 路径与参数定义
FDR_THRESHOLD <- 0.05    
TOP_PROP      <- 0.05   
suffix <- paste0("FDR", gsub("\\.", "", as.character(FDR_THRESHOLD)), 
                 "_Top", gsub("\\.", "", as.character(TOP_PROP)))
rds_dir <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/"
rds_path <- paste0(rds_dir, "Community_Performance_", suffix, ".RDS")

# 2. 读取模块大小数据
Module_size <- read.table('/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新修改label和颜色).txt', 
                          header = TRUE, 
                          sep = "\t", 
                          quote = "", 
                          comment.char = "", 
                          stringsAsFactors = FALSE) %>%
  select(Community = Node_name, Size)

# 3. 读取网络结果数据并合并
all_network_results <- readRDS(rds_path) %>%
  mutate(Community = str_replace(Community, "^S1_M", "Module_")) %>%
  filter(!str_detect(Community, "Group")) %>% 
  left_join(Module_size, by = "Community")

# 4. 过滤有效数据（排除 NA 以及非正数，确保对数变换合法）
plot_data <- all_network_results %>%
  filter(!is.na(Size) & !is.na(Q_c) & Q_c > 0 & Size > 0) 

# 5. 拟合 log-log 线性模型并提取参数
fit_model <- lm(log(Q_c) ~ log(Size), data = plot_data)
fit_summary <- summary(fit_model)
slope_k <- coef(fit_model)[2]
intercept_b <- coef(fit_model)[1]
r_squared <- fit_summary$r.squared
p_val <- fit_summary$coefficients[2, 4]

# 打印拟合结果到控制台
message(sprintf(">>> Empirical log-log slope (k): %.4f (R² = %.3f)", slope_k, r_squared))

# 6. 绘制发表级的高颜值 log(Size) vs log(Q_c) 拟合图
p_log_scaling <- ggplot(plot_data, aes(x = log(Size), y = log(Q_c))) +
  # 散点
  geom_point(alpha = 0.6, color = "#2b5c8f", size = 2) +
  # 线性回归线与置信区间
  geom_smooth(method = "lm", color = "#d95f02", fill = "#fdb462", linetype = "solid", se = TRUE, linewidth = 1) +
  # 主题美化
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold")
  ) +
  # 坐标轴与标题
  labs(
    title = "Empirical Power-Law Decoupling of Module Modularity",
    subtitle = sprintf("Fitted Scaling Factor: k = %.3f (R² = %.3f, p < 2.2e-16)", slope_k, r_squared),
    x = expression(bold(log(Module~Size))),
    y = expression(bold(log(Raw~Modularity~Q[c])))
  ) +
  # 添加规范的数学公式标注
  annotate(
    "text", 
    x = min(log(plot_data$Size)) + 0.5, 
    y = max(log(plot_data$Q_c)) - 0.5,
    label = paste0("log(Q[c]) == ", round(slope_k, 4), " * log(Size) ", 
                   ifelse(intercept_b >= 0, "+ ", "- "), 
                   abs(round(intercept_b, 3))),
    parse = TRUE, size = 4.5, color = "#d95f02", fontface = "bold", hjust = 0
  )

# 显示图形
print(p_log_scaling)
# 显示图形
print(p_log_scaling)


#######开始画图############## Figure S1
# 自定义顺序
library(dplyr)
library(ggplot2)
library(tidyr)
sample_count <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/PCC_statistics.RDS") %>%
  select(gene_number,edge_number,sample_count,AP_cluster_count)
# 转换为长格式并设定 factor 顺序
metric_order <- c("sample_count","gene_number", "edge_number",  "AP_cluster_count")

metric_labels <- c("Sample Count","Gene Number", "Edge Number",  "AP cluster number")

data_long <- sample_count %>%
  pivot_longer(
    cols = all_of(metric_order),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = metric_order, labels = metric_labels))


# 绘图
ggplot(data_long, aes(x = value)) +
  geom_histogram(bins = 15, fill = "#4C9F70", color = "white") +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  labs(x = "Value Range", y = "Number of Networks") +
  theme_minimal(base_size = 14)
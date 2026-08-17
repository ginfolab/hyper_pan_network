################################################################################
#####################################################整合所有的TPM 并进行初步预处理
################################################################################
#############这部分要修改
cutoff  <- 5
Species <- "Escherichia_coli"
path <- paste0("/lustre/home/users/exr/ncbi/public/sra/", Species)
setwd(path)

folder_paths <- list.dirs(path, recursive = FALSE)
folder_paths <- folder_paths[grep("^GSE", basename(folder_paths))]
folder_names <- basename(folder_paths)

library(readxl)
library(openxlsx)
wb <- createWorkbook()

# 创建数据框
final_table <- data.frame(
  dataname = character(),
  gene_count = integer(),
  sample_count = integer(),
  shared_columns = integer()
)

# 存储所有数据框的列名
column_names_list <- list()

for (variable in 1:length(folder_names)) {
  current_path <- paste(folder_paths[variable], "/TPM.txt", sep = "")
  
  if (!file.exists(current_path)) {
    next  # 文件不存在，跳到下一个循环
  }
  
  TPM <- read.table(current_path, header = TRUE, sep = "\t", row.names = 1)
  
  filtered_TPM <- TPM[rowSums(TPM >= cutoff) >= 1, ]
  filtered_TPM <- filtered_TPM[, colSums(filtered_TPM > 0) > 1000]
  
  if (ncol(filtered_TPM) < 15) {
    next
  }
  
  # 记录当前 sheet 的列名
  column_names_list[[folder_names[variable]]] <- colnames(filtered_TPM)
  
  # 添加数据到 Excel
  addWorksheet(wb, sheetName = folder_names[variable])
  writeData(wb, sheet = folder_names[variable], filtered_TPM, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
  
  # 统计
  dataname <- folder_names[variable]
  gene_count <- nrow(filtered_TPM)
  sample_count <- ncol(filtered_TPM)
  
  # 创建数据框
  summary_table <- data.frame(
    dataname = dataname,
    gene_count = gene_count,
    sample_count = sample_count,
    shared_columns = 0  # 占位符
  )
  
  final_table <- rbind(final_table, summary_table)
}

# 计算列名重复数
all_column_names <- unlist(column_names_list)
column_freq <- table(all_column_names)  # 统计每个列名的出现次数

# 计算每个数据框的 shared_columns
for (i in seq_along(final_table$dataname)) {
  dataset_name <- final_table$dataname[i]
  if (dataset_name %in% names(column_names_list)) {
    dataset_columns <- column_names_list[[dataset_name]]
    final_table$shared_columns[i] <- sum(column_freq[dataset_columns] > 1)  # 统计当前数据框的列名在其他数据框中出现的次数
  }
}
# 保存statistics数据
saveRDS(final_table, paste0("New_pangenome_", Species, "_dataSizeStatistics.RDS"))
# 写入 Excel
existing_excel_path <- paste0("New_pangenome_", Species, "_TPM", cutoff, ".xlsx")
file.remove(existing_excel_path)
saveWorkbook(wb, existing_excel_path, overwrite = TRUE)

#####删除sheet多余的累赘名字
# 读取 Excel 文件
file_path <- "New_pangenome_Escherichia_coli_TPM5.xlsx"
wb <- loadWorkbook(file_path)
# 获取所有 sheet 名称
sheets <- names(wb)
# 仅修改包含 `_` 的 sheet
modified_sheets <- sheets[grepl("_", sheets)]
base_names <- sub("_.*", "", modified_sheets)
# 确保新名称唯一
existing_names <- setdiff(sheets, modified_sheets)  # 现有不变的 sheet 名称
all_names <- c(existing_names, base_names)  # 先合并
unique_names <- make.unique(all_names, sep = "_")  # 生成唯一名称
unique_names <- unique_names[(length(existing_names) + 1):length(unique_names)]  # 取修改部分
# 确保 unique_names 生成正确
print(data.frame(modified_sheets, unique_names))
# 修改 sheet 名称
for (i in seq_along(modified_sheets)) {
  renameWorksheet(wb, modified_sheets[i], unique_names[i])
}
# 保存 Excel 文件
saveWorkbook(wb, file_path, overwrite = TRUE)
######   ######   ############   ######   ######





######################Apriori 构建  创建 ######
###################先进行AP聚类。
q_value = 0.5
library(cluster)
library(apcluster)
library(tidyr)
library(readxl)
library(openxlsx)
library(readr)
library(tibble)
library(dplyr)
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/")
filenames <- "New_pangenome_Escherichia_coli_TPM5.xlsx"
sheet_names <- excel_sheets(filenames)

# Create a function for parallel processing
for (number in 1:length(sheet_names)) {
  mat <- read_excel(filenames, sheet = sheet_names[number]) %>%
    mutate(GENE_names = ...1 %>% sub("^[^&]*&", "", .) %>% sub("^[^&]*&", "", .)) %>%
    column_to_rownames("GENE_names") %>%
    select(-...1) %>%
    as.matrix()
  
  # 计算 Pearson 相关性
  cor_res <- Hmisc::rcorr(t(mat), type = "pearson")
  cor_mat <- as.matrix(cor_res$r)
  p_raw <- cor_res$P
  
  # FDR 矫正 P 值
  p_vec <- p_raw[upper.tri(p_raw)]
  p_adj_vec <- p.adjust(p_vec, method = "BH")
  
  # 重构对称的矫正 P 值矩阵
  p_adj_mat <- matrix(NA, nrow = nrow(p_raw), ncol = ncol(p_raw),
                      dimnames = dimnames(p_raw))
  p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  diag(p_adj_mat) <- 0
  
  # 过滤：保留 FDR < 0.05 的相关系数，其余设为 0
  cor_mat[p_adj_mat >= 0.05] <- 0
  diag(cor_mat) <- 1  # 对角线设为1（自相关）
  
  # Affinity Propagation 聚类
  apres_p100 <- apcluster(cor_mat, q = q_value)
  apres_p100@sim <- cor_mat  # 记录相似矩阵
  rawcluster <- labels(apres_p100, type = "exemplars")
  rcTransform <- labels(apres_p100, type = "enum")
  
  # 整理聚类结果
  clusterdata <- as.data.frame(rcTransform)
  colnames(clusterdata) <- "cluster"
  clusterdata$batch <- sheet_names[number]
  clusterdata <- unite(clusterdata, "batch_cluster", batch, cluster)
  clusterdata$ID <- rownames(cor_mat)
  rownames(clusterdata) <- clusterdata$ID
  
  # 保存结果
  RESULT <- paste0("./APcluster_info/", sheet_names[number], "_q=", q_value, ".txt")
  write.table(clusterdata, RESULT, quote = FALSE, sep = "\t", row.names = FALSE)
}
#整合到一起
setwd("./APcluster_info/")
path <- paste("_q=",q_value,".txt",sep = "")
file_names <- list.files(pattern = path)
AP_info <- read.table(file_names[1],header = T,sep = "\t")
for (variable in 2:length(file_names)) {
  candidate_data <- read.table(file_names[variable],header = T,sep = "\t")
  AP_info <- rbind(AP_info,candidate_data)
} 
length(unique(AP_info$batch_cluster))
path1 <- paste("./AP_info_all_",q_value,".txt",sep = "")
AP_info$ID <-  paste0("Ec_",AP_info$ID)
write.table(AP_info,path1,sep = "\t",row.names = F,quote = F)
#################



#####################Apriori 进行分析
######### 准备工作
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/APcluster_info")
library(arules)
library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
Species = c("Escherichia_coli")
cutoff = 0.5
inputfile = "PCC"
corenet_cutoff = 15
path1 <- paste("./AP_info_all_0.5.txt")
groceries <- read.transactions(path1, format="single", sep="\t",header = T,cols=c("batch_cluster", "ID"))
groceries_filter_Items <- groceries[, itemFrequency(groceries, type = "absolute") >= corenet_cutoff]
txn_size <- size(groceries_filter_Items)
groceries_final <-  groceries_filter_Items[txn_size > 1]
######### 准备工作

#########只创建pairwise itemsets#########
# 定义 eclat 参数# 运行 eclat 算法
support_threshold <- 1 / dim(groceries@itemsetInfo)[1]
grocery_rules <- eclat(data = groceries,
                       parameter = list(support = support_threshold,
                                        minlen = 2,
                                        maxlen = 2))
# 将结果转换为数据框

# 使用 data.table 读取并优化处理 # 保存结果
resultapriori <- as.data.table(as(grocery_rules, "data.frame"))
resultapriori[, gene_count := str_count(items, ",") + 1]
resultapriori[, items := gsub("\\{|\\}", "", items)]
resultapriori <- resultapriori[, .(items, count, gene_count)]
names(resultapriori) <- c("Itemsets", "Universality", "Gene_count")
resultapriori  <- as.data.frame(resultapriori)

resultapriori_final <- resultapriori %>%
  mutate(Itemsets = str_remove_all(Itemsets, "Ec_"))%>%    # 删除“EC_”
  separate(Itemsets, into = c("gene1", "gene2"), sep = ",", convert = TRUE) %>%  # 拆分成两列并转换为数值
  mutate(Itemsets = pmap_chr(list(gene1, gene2), ~ paste(sort(c(..1, ..2)), collapse = ","))) %>%  # 排序后合并 
  dplyr::select(Itemsets, Universality,Gene_count)%>%
mutate(Itemsets = str_replace_all(Itemsets, "\\b(\\w+)\\b", "Ec_\\1"))  # 添加“EC_”

path_normal <-  paste("../../",Species[1],"_2items_AP",cutoff,inputfile,"_Uall.RDS",sep = "")
saveRDS(resultapriori_final, path_normal)
#########只创建pairwise itemsets#########




#########创建closed itemsets (多核极速优化版)
###先定义闭集的阈值：
support_threshold <- corenet_cutoff / dim(groceries@itemsetInfo)[1]

message(">>> [1/4] 正在运行 Eclat 算法并追踪 tidLists (这步取决于数据集大小)...")
grocery_rules <- eclat(data = groceries_final,
                       parameter = list(support = support_threshold,
                                        minlen = 2,
                                        maxlen = 100,
                                        target = "closed frequent itemsets",
                                        tidLists = TRUE),
                       control = list(verbose = FALSE)) 

message(">>> [2/4] 正在提取并转换追踪列表...")
tids <- tidLists(grocery_rules)           
tids_list <- as(tids, "list")             

# ==============================================================================
# 🌟 核心提速方案：多核并行字符串拼接 (仅限 Mac/Linux)
# ==============================================================================
message(">>> [3/4] 正在启动多核并行引擎处理字符串拼接...")
library(parallel)
# 自动检测你的电脑有几个 CPU 核心，留出 1 个给系统，剩下的全功率运算
no_cores <- detectCores() - 1 

# 使用 mclapply 替代原有的 vapply，速度提升 N 倍 (N = 你的 CPU 核心数)
tids_strings <- unlist(mclapply(tids_list, function(x) paste(x, collapse = ";"), mc.cores = no_cores))

# ==============================================================================
# 🌟 核心提速方案：全 data.table 原地内存修改
# ==============================================================================
message(">>> [4/4] 正在使用 data.table 与 stringi 进行内存级极速清洗...")
library(stringi)

# 转换 data.table
resultapriori <- as.data.table(as(grocery_rules, "data.frame"))

# 极速赋值（不产生内存复制）
resultapriori[, Dataset_List := tids_strings]

# 使用最高效的 stri_ 系列函数进行字符替换和计数 (比 gsub 快得多)
resultapriori[, items := stri_replace_all_fixed(items, c("{", "}"), "", vectorize_all = FALSE)] 
resultapriori[, gene_count := stri_count_fixed(items, ",") + 1]

# 提取并重命名列
resultapriori <- resultapriori[, .(Itemsets = items, 
                                   Universality = count, 
                                   Gene_count = gene_count, 
                                   Dataset_List)] 

output_path <- sprintf("../../%s_closed_U%s.RDS", Species[1], corenet_cutoff)
saveRDS(resultapriori, file = output_path) 
message(paste0("[+] 计算完成！结果已极速保存至：\n    ", output_path))
#########创建closed itemsets




#########创建maximal itemsets
support_threshold <- corenet_cutoff / dim(groceries@itemsetInfo)[1]

grocery_rules <- eclat(data = groceries,
                       parameter = list(support = support_threshold,
                                        minlen = 2,
                                        maxlen = 100,
                                        target = "maximally frequent itemsets"),
control = list(verbose = FALSE)) # 减少输出信息以加速

# 将结果转换为数据框
# 使用 data.table 读取并优化处理 # 保存结果

resultapriori <- as.data.table(as(grocery_rules, "data.frame"))
resultapriori[, gene_count := str_count(items, ",") + 1]
resultapriori[, items := gsub("\\{|\\}", "", items)]
resultapriori <- resultapriori[, .(items, count, gene_count)]
names(resultapriori) <- c("Itemsets", "Universality", "Gene_count")
resultapriori  <- as.data.frame(resultapriori)
path_normal <-  paste("../../",Species[1],"_maximal_U",corenet_cutoff,".RDS",sep = "")
saveRDS(resultapriori, path_normal)




########制作pairwise edge 和node 在动态U下的统计########
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

df  <- readRDS( "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Escherichia_coli_2items_AP0.5PCC_Uall.RDS")

final_df <- data.frame(Universality = integer(),
                       n_pairs = integer(),
                       n_genes = integer())


for (i in unique(sort(df$Universality))) {
  current_data <- df %>% filter(Universality >= i)
  current_statistics_df <- data.frame(Universality = paste0("≥",i),
                                      n_pairs = dim(current_data)[1],
                                      n_genes = current_data %>% separate_rows(Itemsets, sep = ",") %>% pull(Itemsets) %>% unique() %>% length())
  
  final_df <- rbind(final_df, current_statistics_df)
}

final_df$propotion_genes <- paste0(round( final_df$n_genes / final_df$n_genes[1],4) *100,"%")
final_df$propotion_pairs <- paste0(round( final_df$n_pairs / final_df$n_pairs[1],4) *100,"%")

write.csv(final_df,"pairwise_edge_node_statistics.csv",row.names = F)

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ==== 数据预处理 ====
final_df <- final_df %>%
  mutate(
    Universality_num = as.numeric(str_remove(Universality, "≥")),
    propotion_genes_num = as.numeric(str_remove(propotion_genes, "%")),
    propotion_pairs_num = as.numeric(str_remove(propotion_pairs, "%"))
  ) %>%
  select(Universality_num, propotion_genes_num, propotion_pairs_num)

# 只保留 "Pairs" 数据
pair_df <- final_df %>%
  select(Universality_num, propotion_pairs_num) %>%
  rename(Percentage = propotion_pairs_num)

# 只保留 "Genes" 数据
gene_df <- final_df %>%
  select(Universality_num, propotion_genes_num) %>%
  rename(Percentage = propotion_genes_num)

# ==== 左图（Nodes / Genes，无断轴） ====
p_nodes <- ggplot(gene_df, aes(x = Universality_num, y = Percentage)) +
  geom_line(color = "#1f78b4", size = 1.1) +
  geom_point(color = "#1f78b4", size = 2) +
  scale_x_continuous(breaks = seq(0, max(gene_df$Universality_num), by = 5)) +
  labs(x = "Universality threshold (≥U)",
       y = "Percentage of genes (%)",
       title = "Nodes (Genes)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

# ==== 上半部分（Pairs，100–50%） ====
p_top <- ggplot(pair_df, aes(x = Universality_num, y = Percentage)) +
  geom_line(color = "#e31a1c", size = 1.1) +
  geom_point(color = "#e31a1c", size = 2) +
  coord_cartesian(ylim = c(50, 100)) +
  scale_x_continuous(breaks = seq(0, max(pair_df$Universality_num), by = 5)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "", title = "Edges (Gene pairs)") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

# ==== 下半部分（Pairs，0–1%） ====
p_bottom <- ggplot(pair_df, aes(x = Universality_num, y = Percentage)) +
  geom_line(color = "#e31a1c", size = 1.1) +
  geom_point(color = "#e31a1c", size = 2) +
  coord_cartesian(ylim = c(0, 3)) +
  scale_x_continuous(breaks = seq(0, max(pair_df$Universality_num), by = 5)) +
  theme_minimal(base_size = 14) +
  labs(x = "Universality threshold (≥U)", y = "")

# ==== 拼接右图（断轴上下拼） ====
p_pairs <- p_top / p_bottom + plot_layout(heights = c(1, 3))

# ==== 最终左右拼接（左:Nodes；右:Pairs） ====
final_plot <- p_nodes | p_pairs +
  plot_layout(widths = c(1, 1.5))  # 控制左右宽度比例

final_plot

########制作pairwise edge 和node 在动态U下的统计########

####统计每个数据出现了多少个边，点， APcluster， sample count####
library(dplyr)
library(cluster)
library(apcluster)
library(tidyr)
library(readxl)
library(openxlsx)
library(readr)
library(tibble)
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/")
filenames <- "New_pangenome_Escherichia_coli_TPM5.xlsx"
sheet_names <- excel_sheets(filenames)

# Create a function for parallel processing
data_statistics <- data.frame(
  data_name = character(),
  edge_number = integer(),
  edge_number_positive = integer(),
  edge_number_negative =  integer(),
  minimal_PCC_abs_correlation = numeric()
)
data_c_df_final <- data.frame()

for (number in 1:length(sheet_names)) {
  mat <- read_excel(filenames, sheet = sheet_names[number]) %>%
    mutate(GENE_names = ...1 %>% sub("^[^&]*&", "", .) %>% sub("^[^&]*&", "", .)) %>%
    column_to_rownames("GENE_names") %>%
    select(-...1) %>%
    as.matrix()
  
  # 计算 Pearson 相关性
  cor_res <- Hmisc::rcorr(t(mat), type = "pearson")
  cor_mat <- as.matrix(cor_res$r)
  p_raw <- cor_res$P
  
  # FDR 矫正 P 值
  p_vec <- p_raw[upper.tri(p_raw)]
  p_adj_vec <- p.adjust(p_vec, method = "BH")
  
  # 重构对称的矫正 P 值矩阵
  p_adj_mat <- matrix(NA, nrow = nrow(p_raw), ncol = ncol(p_raw),
                      dimnames = dimnames(p_raw))
  p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  diag(p_adj_mat) <- 0
  
  # 过滤：保留 FDR < 0.05 的相关系数，其余设为 0
  cor_mat[p_adj_mat >= 0.05] <- 0
  diag(cor_mat) <- 1  # 对角线设为1（自相关）
  
  data_c_df <- as.data.frame(as.table(cor_mat)) %>%
    setNames(c("Gene1", "Gene2", "Correlation")) %>%
    filter(Correlation != 0,
           Gene1 != Gene2)
  data_c_df_save <- data_c_df %>% select(Correlation)
  data_c_df_final <- rbind(data_c_df_final,data_c_df_save)
  
  
  current_data <- data.frame(
    data_name = sheet_names[number],
    edge_number = dim(data_c_df)[1],
    gene_number =  length(unique(c(data_c_df$Gene1, data_c_df$Gene2))),
    edge_number_positive = sum(data_c_df$Correlation > 0),
    edge_number_negative = sum(data_c_df$Correlation < 0),
    minimal_PCC_abs_correlation = abs(data_c_df$Correlation[which.min(abs(data_c_df$Correlation))])
  )
  data_statistics <- rbind(data_statistics,current_data)
}
saveRDS(data_statistics,"data_statistics.RDS")
saveRDS(data_c_df_final,"data_c_df_final.RDS")

#改文件夹 设置
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/")
data_statistics <- readRDS("data_statistics.RDS")
data_c_df_final <-  readRDS("data_c_df_final.RDS")
#读取每个数据的AP_cluster 数量并整合
AP_count <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/APcluster_info/AP_info_all_0.5.txt", header = TRUE, sep = "\t") %>%
  mutate(data = str_remove(batch_cluster, "_.*")) %>%
  select(-ID) %>%
  distinct() %>%
  count(data, name = "AP_cluster_count") %>%
  rename(data_name = data)


# 读取每个数据的sample 数量并整合
sample_count <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/New_pangenome_Escherichia_coli_dataSizeStatistics.RDS") %>%
  select(dataname, sample_count) %>%
  mutate(dataname = str_remove(dataname, "_.*")) %>%
  rename(data_name = dataname) %>%
  left_join(data_statistics,by = "data_name" ) %>%
  left_join(AP_count, "data_name")

saveRDS(sample_count,"/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/PCC_statistics.RDS")




#######开始画图##############
# 自定义顺序
library(ggplot2)
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
####统计每个数据出现了多少个边和点 APcluster sample count####

####统计统计所有数据的 PCC cor####
data_c_df_final  <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/data_c_df_final.RDS")
sample_count <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/PCC_statistics.RDS")
library(tidyverse)
# 构建长格式数据
plot_data <- tibble(
  value = c(data_c_df_final$Correlation, sample_count$minimal_PCC_abs_correlation),
  metric = c(
    rep("Correlation", length(data_c_df_final$Correlation)),
    rep("Minimal PCC (abs)", length(sample_count$minimal_PCC_abs_correlation))
  )
)

# 绘图
ggplot(plot_data, aes(x = value)) +
  geom_histogram(bins = 15, fill = "#4C9F70", color = "white") +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  labs(x = "Value Range", y = "Number of Networks") +
  theme_minimal(base_size = 14)
####统计统计所有数据的 PCC cor####

####一些额外的统计#######
table(plot_data$metric)

#Correlation Minimal PCC (abs) 
#832037914               106 
 positive_count <- sum(plot_data$value > 0)
negative_count <- sum(plot_data$value < 0)

 # 计算比例
 total_count <- positive_count + negative_count
 positive_ratio <- positive_count / total_count
negative_ratio <- negative_count / total_count

 # 打印结果
cat("正数比例：", round(positive_ratio, 3), "\n")
#正数比例： 0.86 
cat("负数比例：", round(negative_ratio, 3), "\n")
#负数比例： 0.14 
####一些额外的统计#######





######################韦恩图 对比新老方法的边的关系######
library(cluster)
library(apcluster)
library(tidyr)
library(readxl)
library(openxlsx)
library(readr)
library(tibble)
library(dplyr)
setwd("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/")
filenames <- "New_pangenome_Escherichia_coli_TPM5.xlsx"
sheet_names <- excel_sheets(filenames)
final_old001_method_gene_pair <- data_frame()
final_old005_method_gene_pair <- data_frame()
final_new_method_gene_pair <- data_frame()

# 提取所有的 FDR 显著的边 保存起来
for (number in 1:length(sheet_names)) {
  mat <- read_excel(filenames, sheet = sheet_names[number]) %>%
    mutate(GENE_names = ...1 %>% sub("^[^&]*&", "", .) %>% sub("^[^&]*&", "", .)) %>%
    column_to_rownames("GENE_names") %>%
    select(-...1) %>%
    as.matrix()
  
  # 计算 Pearson 相关性
  cor_res <- Hmisc::rcorr(t(mat), type = "pearson")
  cor_mat <- as.matrix(cor_res$r)
  p_raw <- cor_res$P
  
  # FDR 矫正 P 值
  p_vec <- p_raw[upper.tri(p_raw)]
  p_adj_vec <- p.adjust(p_vec, method = "BH")
  
  # 重构矫正后的 P 值矩阵
  p_adj_mat <- matrix(NA, nrow = nrow(p_raw), ncol = ncol(p_raw))
  rownames(p_adj_mat) <- rownames(p_raw)
  colnames(p_adj_mat) <- colnames(p_raw)
  p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  diag(p_adj_mat) <- 0  # 可选：设为0 或 NA
  
  # 提取基因对、相关系数和矫正后的 P 值为三列
  gene_names <- rownames(cor_mat)
  combn_index <- which(upper.tri(cor_mat), arr.ind = TRUE)
  
  genepair_df <- data.frame(
    Genepair = paste(pmin(as.numeric(gene_names[combn_index[, 1]]), as.numeric(gene_names[combn_index[, 2]])),
                     pmax(as.numeric(gene_names[combn_index[, 1]]), as.numeric(gene_names[combn_index[, 2]])),
                     sep = ","),
    Correlation = cor_mat[upper.tri(cor_mat)],
    FDR = p_adj_mat[upper.tri(p_adj_mat)],
    dataset = sheet_names[number]
  ) %>%
    unique()
  
  saveRDS(genepair_df,paste0("./correlation/",sheet_names[number],"_all_gene_pairs.RDS"))
  
  
  #提取前百分之一 cor 的基因对 旧方法
  top1_percent_threshold <- quantile(genepair_df$Correlation, probs = 0.99, na.rm = TRUE)
  current_old001_method_gene_pair <- genepair_df %>%
    filter(Correlation >= top1_percent_threshold) %>%
    select(Genepair)
  final_old001_method_gene_pair <- final_old001_method_gene_pair %>%
    bind_rows(current_old001_method_gene_pair) 
  
  #提取前百分之五 cor 的基因对 旧方法
  top5_percent_threshold <- quantile(genepair_df$Correlation, probs = 0.95, na.rm = TRUE)
  current_old005_method_gene_pair <- genepair_df %>%
    filter(Correlation >= top5_percent_threshold) %>%
    select(Genepair)
  final_old005_method_gene_pair <- final_old005_method_gene_pair %>%
    bind_rows(current_old005_method_gene_pair) 
  
  #提取显著的 cor 的基因对 新方法
  current_new_method_gene_pair <- genepair_df %>%
    filter(FDR < 0.05)%>%
    select(Genepair)
  final_new_method_gene_pair <- final_new_method_gene_pair %>%
    bind_rows(current_new_method_gene_pair) 
}

final_old005_method_gene_pair_U <- final_old005_method_gene_pair %>%
  count(Genepair, name = "Universality")
saveRDS(final_old005_method_gene_pair_U,"final_old005_method_gene_pair_U.RDS")

final_old001_method_gene_pair_U <- final_old001_method_gene_pair %>%
  count(Genepair, name = "Universality")
saveRDS(final_old001_method_gene_pair_U,"final_old001_method_gene_pair_U.RDS")

#绘图
library(stringr)
library(dplyr)
final_old005_method_gene_pair <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/final_old005_method_gene_pair_U.RDS")
final_old001_method_gene_pair <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/final_old001_method_gene_pair_U.RDS")
final_new_method_gene_pair <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Escherichia_coli_2items_AP0.5PCC_Uall.RDS") %>%
 rename(Genepair = Itemsets) %>%
mutate(Genepair = str_remove_all(Genepair,"Ec_"))
  
new_pairs <- as.character(na.omit(final_new_method_gene_pair$Genepair))
old001_pairs <- as.character(na.omit(final_old001_method_gene_pair$Genepair))
old005_pairs <- as.character(na.omit(final_old005_method_gene_pair$Genepair))
library(eulerr)

# 准备数据
venn_data <- list(Significant_cor_edge = new_pairs,
                  Top1_percentage_cor_edge = old001_pairs ,
                  Top5_percentage_cor_edge =  old005_pairs)

# 构建欧拉图对象
fit <- euler(venn_data)
plot(fit,
     fills = list(fill = c("orange", "lightblue", "navy"), alpha = 0.6),
     quantities = TRUE)
######################韦恩图 对比新老方法的边的关系######


#########创建PCC_correlation_all(TPM5_FDR0.05).rds
rds_files <- list.files("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/correlation", 
                        pattern = "\\.RDS$", full.names = TRUE)

library(purrr)  # 确保加载 purrr
library(dplyr)
library(stringr)
combined_data <- lapply(rds_files, function(f) {
  tryCatch({
    readRDS(f)
  }, error = function(e) {
    message("无法读取文件: ", f)
    return(NULL)
  })
}) %>% purrr::compact() %>% bind_rows() %>%   mutate(Genepair = str_replace_all(Genepair, "\\b(\\w+)\\b", "Ec_\\1"))

saveRDS(combined_data,"PCC_correlation_all(TPM5_FDR0.05).rds")



####################################################################################################################################
############################################用旧方法计算基因edge 的Universality############################################
################################################################################################################################################################################
# 加载必要的包
library(dplyr)
library(purrr)

# 1. 设置文件夹路径
folder_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/correlation/"

# 2. 获取该文件夹下所有 RDS 文件的完整路径
rds_files <- list.files(path = folder_path, pattern = "\\.RDS$", full.names = TRUE)

# 3. 批量读取、筛选并合并
message("开始读取并筛选显著边，请稍候...")

sig_edges_combined <- map_dfr(rds_files, function(file) {
  # 读取单个 RDS 文件
  df <- readRDS(file)
  
  # 立即筛选并丢弃多余列，释放内存压力
  df %>%
    filter(FDR < 0.05) %>%   # 【关键】筛选 FDR 小于 0.05
    select(Genepair)         # 【关键】只保留 Genepair 列
})

# 4. 统计 Universality (计算每个 Genepair 出现的总次数)
message("开始统计 Universality...")

final_universality_df <- sig_edges_combined %>%
  count(Genepair, name = "Universality") %>%
  arrange(desc(Universality)) # 按出现次数从高到低排序，方便查看

# 查看最终结果的前几行
head(final_universality_df)

# 5. (可选) 保存最终结果
# saveRDS(final_universality_df, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/Gene_Universality_Summary.RDS")

message("全部处理完成！")
library(stringr)

# 给每个基因 ID 加上 Ec_ 前缀
final_universality_df <- final_universality_df %>%
  mutate(
    # "(\\d+)" 匹配任何连续的数字，"Ec_\\1" 表示在匹配到的数字前加上 Ec_
    Genepair = str_replace_all(Genepair, "(\\d+)", "Ec_\\1")
  )
saveRDS(final_universality_df, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/Old_method_pairs_Universality_Summary.RDS")


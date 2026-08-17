# ==============================================================================
# 0. 加载核武器级别的网络计算包
# ==============================================================================
library(dplyr)
library(tidyr)
library(stringr)
library(igraph)   
library(purrr)    
# --- 【新增并行计算包】 ---
library(future)
library(furrr)

FDR_THRESHOLD <- 0.05    # FDR 显著性阈值
TOP_PROP      <- 0.05  # 保留前百分之多少的强边 (0.05 即 5%)
MIN_NODES     <- 3       # 模块幸存节点最少数量

# ==============================================================================
# 1. 清洗与整合超级字典 (双轨制逻辑)
# ==============================================================================
message(">>> 正在从 Module_anno 构建超级社区字典...")







Module_anno <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_Gene_Region_Annotation.RDS")
module_dict <- split(Module_anno$Gene, Module_anno$Module)

clean_group_anno <- Module_anno %>%
  filter(!is.na(Group) & str_detect(Group, "^Group_"))
group_dict <- split(clean_group_anno$Gene, clean_group_anno$Group)

# 1.3 合体为原始的超级大字典 (包含所有大大小小的模块)
original_community_list <- c(module_dict, group_dict)

# --- 【逻辑修正点 1：在这里提取原始、未过滤的全局核心户口本！】 ---
global_core_genes <- unique(unlist(original_community_list))
message(paste0(">>> 全局核心基因户口本已确立！共包含 ", length(global_core_genes), " 个核心基因。"))

# 1.4 【执行铁律】：剔除初始基因数量少于 5 个的社区，这批人不作为独立队伍接受考核
eval_community_list <- original_community_list[lengths(original_community_list) >= 5]  
message(paste0(">>> 过滤完成！共计 ", length(eval_community_list), " 支精锐队伍 (Size >= 5) 即将投入 106 重宇宙进行考核。"))

# ==============================================================================
# 2. 定义核心扫描引擎 (修正补全版)
# ==============================================================================
evaluate_communities_in_network <- function(rds_file, comm_list, core_gene_pool) {
  
  net_df <- readRDS(rds_file)
  dataset_name <- net_df$dataset[1]
  
  # --- Top 5% 强边过滤 ---
  significant_edges <- net_df %>% filter(FDR < FDR_THRESHOLD)
  if(nrow(significant_edges) < 10) return(NULL) 
  
  valid_edges <- significant_edges %>%
    arrange(desc(abs(Correlation))) %>%
    slice_head(prop = TOP_PROP)
  
  message(paste0(">>> [", dataset_name, "] 显著边数: ", nrow(significant_edges), 
                 " | 保留 Top ",TOP_PROP," 边数: ", nrow(valid_edges)))
  
  # --- 构建边表 ---
  edge_list <- valid_edges %>%
    separate(Genepair, into = c("Gene1", "Gene2"), sep = ",") %>%
    mutate(
      Gene1 = paste0("Ec_", Gene1),
      Gene2 = paste0("Ec_", Gene2),
      weight = abs(Correlation)
    ) %>% 
    select(Gene1, Gene2, weight)     
  
  # 【⚠️ 核心补全：创建图对象及计算全局指标】
  g <- graph_from_data_frame(edge_list, directed = FALSE)
  total_m <- sum(E(g)$weight)              
  node_strengths_all <- strength(g)          
  
  # --- 逻辑修正点 2：使用外部传入的原始户口本 ---
  valid_all_core_in_net <- intersect(core_gene_pool, V(g)$name)
  
  # 铸造一个只有核心基因的封闭小宇宙
  g_core_only <- induced_subgraph(g, valid_all_core_in_net)
  node_strengths_core_only <- strength(g_core_only)
  # --------------------------------------------------
  
  res_df <- map_dfr(names(comm_list), function(comm_name) {
    
    comm_genes <- comm_list[[comm_name]]
    valid_genes <- intersect(comm_genes, V(g)$name)
    n_nodes <- length(valid_genes)
    
    # 至少需要 3 个幸存节点才有统计意义
    if(n_nodes < MIN_NODES) {
      return(data.frame(
        Community = comm_name, Q_c = NA, Density = NA, Mean_PCC = NA,
        Ratio_Internal = NA, Ratio_CrossCore = NA, Ratio_Accessory = NA, 
        Valid_Nodes = n_nodes
      ))
    }
    
    sub_g <- induced_subgraph(g, valid_genes)
    E_c <- sum(E(sub_g)$weight)               
    K_c <- sum(node_strengths_all[valid_genes])   
    Q_c <- (E_c / total_m) - (K_c / (2 * total_m))^2 
    
    density <- ecount(sub_g) / (n_nodes * (n_nodes - 1) / 2)
    mean_pcc <- ifelse(ecount(sub_g) > 0, mean(E(sub_g)$weight), 0)
    
    # --- 三分天下计算 ---
    strength_internal <- 2 * E_c 
    strength_all_core <- sum(node_strengths_core_only[valid_genes])
    strength_cross_core <- strength_all_core - strength_internal
    strength_accessory <- K_c - strength_all_core
    
    ratio_internal  <- ifelse(K_c > 0, strength_internal / K_c, 0)
    ratio_crosscore <- ifelse(K_c > 0, strength_cross_core / K_c, 0)
    ratio_accessory <- ifelse(K_c > 0, strength_accessory / K_c, 0)
    
    data.frame(
      Community = comm_name, 
      Q_c = Q_c, 
      Density = density, 
      Mean_PCC = mean_pcc,
      Ratio_Internal = ratio_internal,     
      Ratio_CrossCore = ratio_crosscore,   
      Ratio_Accessory = ratio_accessory,   
      Valid_Nodes = n_nodes
    )
  })
  
  res_df$Dataset <- dataset_name
  return(res_df)
}


# ==============================================================================
# 3. 启动 106 个数据集的星际扫荡 (多核并发版)
# ==============================================================================
message(">>> 引擎点火，正在检测系统可用核心数...")

# 检测你的电脑有几个 CPU 核心，并保留 1 个核心给操作系统以免电脑卡死
# 如果你的内存较小（比如只有 8GB），建议强制设定为 workers = 4
num_cores <- availableCores() - 3 
plan(multisession, workers = num_cores)

message(paste0(">>> 成功唤醒 ", num_cores, " 个并发线程！开启狂暴扫描模式..."))

rds_dir <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/correlation/"
rds_files <- list.files(rds_dir, pattern = "\\.RDS$", full.names = TRUE)

# 【核心修改】：把 map_dfr 换成了 future_map_dfr
# .options = furrr_options(seed = TRUE) 是为了保证多线程下随机数的稳定性（虽然我们没用到随机数，但这是个好习惯）
# 为了能实时看到输出，建议将 .progress 设为 TRUE
all_network_results <- future_map_dfr(
  rds_files, 
  function(x) {
    # 提取文件名作为标识
    current_file <- basename(x)
    
    # 运行核心函数
    res <- evaluate_communities_in_network(x, eval_community_list, global_core_genes)
    
    # 【关键】：运行完一个，报一条信息
    # 这里的 message 会在每个 worker 完成任务时尝试写回主控制台
    message(paste0("Done: [", current_file, "] - ", Sys.time()))
    
    return(res)
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE  # 开启 furrr 自带的进度条
)

# 整理列顺序
all_network_results <- all_network_results %>%
  select(Dataset, Community, Valid_Nodes, Q_c, Density, Mean_PCC, 
         Ratio_Internal, Ratio_CrossCore, Ratio_Accessory)

message(">>> 扫荡完毕！大一统矩阵已生成！")

# 导出终极数据表
suffix <- paste0("FDR", gsub("\\.", "", as.character(FDR_THRESHOLD)), 
                 "_Top", gsub("\\.", "", as.character(TOP_PROP)))
# --- 自动生成带参数的文件名 ---
output_base <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/Community_Performance_"
csv_path <- paste0(output_base, suffix, ".csv")
rds_path <- paste0(output_base, suffix, ".RDS")

write.csv(all_network_results, csv_path, row.names = FALSE)
saveRDS(all_network_results, rds_path)
message(paste0(">>> 结果已安全着陆于: ", output_base))

# 跑完后关闭多线程，释放内存
plan(sequential)












# ==============================================================================
# 终极自主聚类热图：Euclidean距离 + 0-1标准化 + 自适应颜色分布
# ==============================================================================
# --- 参数设置 ---
FDR_THRESHOLD <- 0.05    
TOP_PROP      <- 0.05   
suffix <- paste0("FDR", gsub("\\.", "", as.character(FDR_THRESHOLD)), 
                 "_Top", gsub("\\.", "", as.character(TOP_PROP)))

output_base <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/PCC/Community_Performance_"
rds_path <- paste0(output_base, suffix, ".RDS")
#Top module的颜色
ref_color_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_color_mapping.RDS") %>% 
     dplyr::select(Module, color_mapping) %>%
     filter(str_detect(Module, "^S1_M")) %>% 
     mutate(Module = str_replace(Module, "^S1_M", "Module_"))
# ==============================================================================
# 0. 加载必要的包
# ==============================================================================
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
})

message(">>> 1. 正在加载数据并重构矩阵...")

# 读取数据，过滤并转为宽矩阵
all_network_results <- readRDS(rds_path) %>%
  mutate(Community = str_replace(Community, "^S1_M", "Module_")) %>%
  filter(!str_detect(Community, "Group")) # 确保提取有效的模块

df_wide <- all_network_results %>%
  select(Dataset, Community, Q_c) %>%
  pivot_wider(names_from = Community, values_from = Q_c)

mat_qc <- as.matrix(df_wide[, -1])

# 强力清洗行列名称（去除空格和括号）
rownames(mat_qc) <- str_trim(str_remove(df_wide$Dataset, "\\s*\\(.*?\\)"))
colnames(mat_qc) <- str_trim(str_remove(colnames(mat_qc), "\\s*\\(.*?\\)"))

# ==============================================================================
# 终极自主聚类热图：Spearman相关性 + 负数剔除 + 维度对齐
# ==============================================================================

# --- 1. 数据加载与重构 ---
message(">>> 1. 正在加载数据并重构矩阵...")
all_network_results <- readRDS(rds_path) %>%
  mutate(Community = str_replace(Community, "^S1_M", "Module_")) %>%
  filter(!str_detect(Community, "Group")) 

df_wide <- all_network_results %>%
  select(Dataset, Community, Q_c) %>%
  pivot_wider(names_from = Community, values_from = Q_c)

mat_qc <- as.matrix(df_wide[, -1])
rownames(mat_qc) <- str_trim(str_remove(df_wide$Dataset, "\\s*\\(.*?\\)"))
colnames(mat_qc) <- str_trim(str_remove(colnames(mat_qc), "\\s*\\(.*?\\)"))

# --- 2. 矩阵预处理：剔除负数与标准化 ---
message(">>> 2. 正在执行负数剔除与 0-1 标准化...")

# 【关键点】：先将负数设为 NA
#mat_qc[mat_qc < 0] <- NA
is_na_idx <- is.na(mat_qc)

# 备份用于计算的矩阵 (NA 转 0)
mat_for_calc <- mat_qc
mat_for_calc[is_na_idx] <- 0

# 全局 Min-Max 标准化
min_val <- min(mat_for_calc, na.rm = TRUE)
max_val <- max(mat_for_calc, na.rm = TRUE)
mat_scaled <- (mat_for_calc - min_val) / (max_val - min_val)

# 备份用于画图的矩阵 (还原 NA)
mat_plot_final <- mat_scaled
mat_plot_final[is_na_idx] <- NA

# --- 3. 颜色映射与 Region 准备 ---
message(">>> 3. 准备颜色与 Region 注释...")
colors_11 <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
               "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
valid_vals <- mat_plot_final[!is.na(mat_plot_final)]
dynamic_breaks <- unique(quantile(valid_vals, probs = seq(0, 1, length.out = 11)))
col_fun <- colorRamp2(dynamic_breaks, colorRampPalette(colors_11)(length(dynamic_breaks)))

Module_anno <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_Gene_Region_Annotation.RDS")
community_region_map <- Module_anno %>%
  select(Community = Module, Region) %>%
  mutate(Community = str_replace(str_trim(Community), "^S1_M", "Module_"), Region = str_trim(Region)) %>%
  filter(!is.na(Region)) %>% distinct() 

col_regions <- community_region_map$Region[match(colnames(mat_plot_final), community_region_map$Community)]
col_regions[is.na(col_regions)] <- "Other"
region_colors <- c("Region_1" = "#E41A1C", "Region_2" = "#377EB8", "Region_3" = "#4DAF4A", "Region_4" = "#984EA3", "Other" = "grey")

# --- 4. 维度对齐与聚类计算 ---
message(">>> 4. 正在对齐维度并计算 Spearman 聚类树...")
# 剔除标准差为 0 的行列
row_sds <- apply(mat_scaled, 1, sd)
col_sds <- apply(mat_scaled, 2, sd)
valid_rows <- row_sds > 0
valid_cols <- col_sds > 0

mat_calc_aligned <- mat_scaled[valid_rows, valid_cols]
mat_plot_aligned <- mat_plot_final[valid_rows, valid_cols]

# 计算树
hc_rows <- hclust(as.dist(1 - cor(t(mat_calc_aligned), method = "spearman")), method = "ward.D2")
hc_cols <- hclust(as.dist(1 - cor(mat_calc_aligned, method = "spearman")), method = "ward.D2")

# ==============================================================================
# 5. 渲染热图 (修正变量名)
# ==============================================================================
ht_qc_auto <- Heatmap(
  mat_plot_aligned,        # <--- 这里由 mat_render 改为 mat_plot_aligned
  name = "Norm Q_c",
  col = col_fun,
  na_col = "#505050",      # NA 渲染为深灰色
  
  # 顶部注释也要对应过滤
  top_annotation = HeatmapAnnotation(
    Region = col_regions[valid_cols], 
    col = list(Region = region_colors),
    annotation_name_side = "left"
  ),
  
  cluster_rows = hc_rows,
  cluster_columns = hc_cols,
  
  # 自动切割
  row_split = 3,           
  column_split = 4,        
  
  row_title = "Cluster %s", 
  column_title = "Group %s",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  row_gap = unit(3, "mm"),
  column_gap = unit(3, "mm"),
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5), 
  
  heatmap_legend_param = list(
    title = "Norm Q_c",
    title_gp = gpar(fontface = "bold"),
    color_bar = "continuous",
    legend_height = unit(5, "cm")
  )
)

# 打印热图
draw(ht_qc_auto)





# 导出为 PDF
output_pdf <- paste0("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Heatmap_AutoCluster_Qc_FDR", FDR_THRESHOLD, "_TOP", TOP_PROP, ".pdf")
pdf(output_pdf, width = 20, height = 12)

dev.off()

message(paste0(">>> 大功告成！全自动聚类热图已导出至：\n", output_pdf))



# ==============================================================================
# 提取热图聚类与分组信息
# ==============================================================================
# 1. 强制渲染并获取对象
ht_obj <- draw(ht_qc_auto)

# 2. 提取行分组 (Dataset)
r_order_list <- row_order(ht_obj)

df_rows <- data.frame() # 先初始化
for(i in seq_along(r_order_list)) {
  indices <- r_order_list[[i]]
  tmp_df <- data.frame(
    Dataset = rownames(mat_render)[indices],
    Row_Cluster = paste0("Cluster_", i),
    stringsAsFactors = FALSE
  )
  df_rows <- rbind(df_rows, tmp_df)
}

# 3. 提取列分组 (Community)
c_order_list <- column_order(ht_obj)

df_columns <- data.frame() # 先初始化
for(j in seq_along(c_order_list)) {
  indices <- c_order_list[[j]]
  tmp_df <- data.frame(
    Community = colnames(mat_render)[indices],
    Column_Group = paste0("Group_", j),
    stringsAsFactors = FALSE
  )
  df_columns <- rbind(df_columns, tmp_df)
}

Module2region_info <-  df_columns  %>%
  rename(Module = Community,
         Region = Column_Group) %>%
  mutate(Region = str_replace(Region, "Group_", "Region_")) %>%
  
# ==============================================================================
# 4. 立即检查
# ==============================================================================
if(nrow(df_rows) > 0){
  message(">>> 成功抓取！行数：", nrow(df_rows))
  print(head(df_rows))
} else {
  # 如果还是空，唯一的可能就是 mat_render 根本没有 rownames
  message(">>> 严重错误：mat_render 的行名或列名为 NULL！")
  print(head(rownames(mat_render)))
}

Dataset2cluster_info <- df_rows %>%
  rename(Dataset = Dataset, Clade = Row_Cluster) %>%
  mutate(Clade = str_replace(Clade, "Cluster_", "Clade_"))








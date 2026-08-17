suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(readr)
  library(stringr)
  library(purrr)    
  library(future)
  library(furrr)
  library(progressr) 
})

handlers(handler_txtprogressbar(char = "=", style = 3))

message(">>> 🚀 正在启动 E. coli 泛网络 [Top 5% 强骨架 + 7相态全计数(全量原始 Acc) + 模块密度] 终极流水线...")

# ==============================================================================
# 0. 全局路径与配置
# ==============================================================================
base_path         <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
rds_dir           <- file.path(base_path, "PCC/correlation/")
table2_clade_file <- file.path(base_path, "paper图片/Table2 GSE_info_sorted.xlsx")
output_target_dir <- file.path(base_path, "PCC") 

# 🌟 设置参数
FDR_THRESHOLD <- 0.05
TOP_PROP      <- 0.05 
PCC_power     <- 1    
MIN_NODES     <- 3      

suffix <- paste0("FDR", gsub("\\.", "", as.character(FDR_THRESHOLD)), "_Top", gsub("\\.", "", as.character(TOP_PROP)))

# 定位 Modularity 结果文件
modularity_rds_path <- file.path(base_path, "PCC", paste0("Community_Performance_", suffix, ".RDS"))
if (!file.exists(modularity_rds_path)) {
  stop("[-] 找不到 Modularity RDS 文件，请检查路径。")
}
mod_data <- readRDS(modularity_rds_path)

# 修复列名 Q_c
if ("Q_c" %in% names(mod_data)) {
  mod_data <- mod_data %>% rename(Modularity = Q_c)
}

# 消除命名冲突
if ("Community" %in% names(mod_data)) {
  mod_data <- mod_data %>% 
    mutate(Community = str_replace(Community, "^S1_M", "Module_"))
}

# ------------------------------------------------------------------------------
# 准备字典 1：核心基因与原始 15-universality 网络 (用于 Signed 边)
# ------------------------------------------------------------------------------
# 读取全量文件 (这里不加 filter，为了后续共用)
all_pairwise_edges <- readRDS(file.path(base_path, "Escherichia_coli_2items_AP0.5PCC_Uall.RDS"))

pairwise_Edge <- all_pairwise_edges %>% filter(Universality >= 15)

Module_anno <- readRDS(file.path(base_path, "Module_Gene_Region_Annotation.RDS")) %>%
  mutate(Raw_Gene_ID = str_remove(Ec_clustid, "^Ec_"))

anno_clean <- Module_anno %>% select(Raw_Gene_ID, Module, Group) %>% distinct(Raw_Gene_ID, .keep_all = TRUE) 
global_core_genes <- unique(Module_anno$Ec_clustid[!is.na(Module_anno$Module)])

pairwise_Edge_classified <- pairwise_Edge %>%
  separate(Itemsets, into = c("Gene1", "Gene2"), sep = ",", remove = FALSE) %>%
  mutate(Gene1_clean = str_remove(Gene1, "^Ec_"), Gene2_clean = str_remove(Gene2, "^Ec_")) %>%
  left_join(anno_clean, by = c("Gene1_clean" = "Raw_Gene_ID")) %>% rename(Module1 = Module, Group1 = Group) %>%
  left_join(anno_clean, by = c("Gene2_clean" = "Raw_Gene_ID")) %>% rename(Module2 = Module, Group2 = Group) %>%
  mutate(Edge_type = case_when(
    is.na(Module1) | is.na(Module2) ~ NA_character_,
    Module1 == Module2 ~ "IM", Group1 == Group2 ~ "IG", TRUE ~ "CG"
  )) %>% filter(!is.na(Edge_type))

# 输出 Table 6
pairwise_Edge_final <- pairwise_Edge_classified %>% select(Edge = Itemsets, Gene1,Gene2, Universality, Module_Gene1 = Module1 , Group_Gene1 = Group1,Module_Gene2 = Module2 , Group_Gene2 = Group2,Edge_type)
write_csv(pairwise_Edge_final, file.path(base_path, "paper图片/Table6 Core network gene pairs.csv"))

edge_class_hash <- pairwise_Edge_classified %>%
  separate(Itemsets, into = c("G1", "G2"), sep = ",") %>%
  mutate(
    G1_pure = str_remove(G1, "^Ec_"), G2_pure = str_remove(G2, "^Ec_"),
    Ordered_Pair = if_else(G1_pure < G2_pure, paste0(G1_pure, "_", G2_pure), paste0(G2_pure, "_", G1_pure))
  ) %>% select(Ordered_Pair, Dict_Edge = Edge_type) %>% distinct(Ordered_Pair, .keep_all = TRUE)

clade_df <- read_excel(table2_clade_file) %>% select(Dataset, Clade) %>% distinct()

# ==============================================================================
# 2. 核心大循环：全景 7相态提取引擎 (已移除 Acc 过滤限制)
# ==============================================================================
Module_anno_with_prefix <- Module_anno %>% mutate(Gene = paste0("Ec_", Raw_Gene_ID))
clean_mod_anno <- Module_anno_with_prefix %>% filter(!is.na(Region))
mod_dict_all <- split(clean_mod_anno$Gene, clean_mod_anno$Module)
valid_mods <- names(mod_dict_all)[lengths(mod_dict_all) >= 5]
mod_dict <- mod_dict_all[valid_mods]

evaluate_network_by_edge_type <- function(rds_file, p) {
  on.exit(p())
  
  net_df <- readRDS(rds_file)
  dataset_name <- net_df$dataset[1]
  sig_edges <- net_df %>% filter(FDR < FDR_THRESHOLD)
  if(nrow(sig_edges) < 10) return(NULL) 
  
  valid_edges <- sig_edges %>% 
    arrange(desc(abs(Correlation))) %>% slice_head(prop = TOP_PROP) %>%
    separate(Genepair, into = c("Gene1", "Gene2"), sep = ",") %>%
    mutate(
      G1_pure = str_remove(Gene1, "^Ec_"), G2_pure = str_remove(Gene2, "^Ec_"),
      Ordered_Pair = if_else(G1_pure < G2_pure, paste0(G1_pure, "_", G2_pure), paste0(G2_pure, "_", G1_pure)),
      Gene1_full = paste0("Ec_", G1_pure), Gene2_full = paste0("Ec_", G2_pure),
      edge_weight = (abs(Correlation))^PCC_power 
    )
  
  valid_edges <- valid_edges %>%
    left_join(anno_clean, by = c("G1_pure" = "Raw_Gene_ID")) %>% rename(Module1 = Module, Group1 = Group) %>%
    left_join(anno_clean, by = c("G2_pure" = "Raw_Gene_ID")) %>% rename(Module2 = Module, Group2 = Group) %>%
    left_join(edge_class_hash, by = "Ordered_Pair") %>%
    mutate(
      # 🌟 核心修复：纯粹的 Acc 判定，无条件接受所有游离边
      Final_Edge_Type = case_when(
        !is.na(Dict_Edge) ~ paste0("Signed_", Dict_Edge),
        is.na(Module1) | is.na(Module2) ~ "Acc",
        Module1 == Module2 ~ "Unsigned_IM",
        Group1 == Group2 ~ "Unsigned_IG",
        TRUE ~ "Unsigned_CG"
      )
    )
  
  all_nodes_in_g <- unique(c(valid_edges$Gene1_full, valid_edges$Gene2_full))
  
  res <- map_dfr(names(mod_dict), function(mod_name) {
    v_genes <- intersect(mod_dict[[mod_name]], all_nodes_in_g)
    n <- length(v_genes)
    
    if(n < MIN_NODES) {
      return(data.frame(
        Community = mod_name, Valid_Nodes = n,
        Possible_IM_Edges = 0,
        Count_Signed_IM = 0, Count_Unsigned_IM = 0, 
        Count_Signed_IG = 0, Count_Unsigned_IG = 0,
        Count_Signed_CG = 0, Count_Unsigned_CG = 0,
        Count_Acc = 0,
        d_S_IM = 0, d_U_IM = 0, d_S_IG = 0, d_U_IG = 0, d_S_CG = 0, d_U_CG = 0, d_Acc = 0
      ))
    }
    
    mod_edges <- valid_edges %>% filter(Gene1_full %in% v_genes | Gene2_full %in% v_genes)
    
    # 边数量拓扑统计 
    Count_Signed_IM   <- sum(mod_edges$Final_Edge_Type == "Signed_IM", na.rm=TRUE)
    Count_Unsigned_IM <- sum(mod_edges$Final_Edge_Type == "Unsigned_IM", na.rm=TRUE)
    Count_Signed_IG   <- sum(mod_edges$Final_Edge_Type == "Signed_IG", na.rm=TRUE) * 0.5
    Count_Unsigned_IG <- sum(mod_edges$Final_Edge_Type == "Unsigned_IG", na.rm=TRUE) * 0.5
    Count_Signed_CG   <- sum(mod_edges$Final_Edge_Type == "Signed_CG", na.rm=TRUE) * 0.5
    Count_Unsigned_CG <- sum(mod_edges$Final_Edge_Type == "Unsigned_CG", na.rm=TRUE) * 0.5
    # 🌟 Acc 不参与模块间双重遍历
    Count_Acc         <- sum(mod_edges$Final_Edge_Type == "Acc", na.rm=TRUE) 
    
    # 能量(权重)统计
    d_S_IM <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Signed_IM"], na.rm=T)
    d_U_IM <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Unsigned_IM"], na.rm=T)
    d_S_IG <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Signed_IG"], na.rm=T) * 0.5
    d_U_IG <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Unsigned_IG"], na.rm=T) * 0.5
    d_S_CG <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Signed_CG"], na.rm=T) * 0.5
    d_U_CG <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Unsigned_CG"], na.rm=T) * 0.5
    d_Acc  <- sum(mod_edges$edge_weight[mod_edges$Final_Edge_Type == "Acc"], na.rm=T)  
    
    data.frame(
      Community = mod_name, Valid_Nodes = n,
      Possible_IM_Edges = ifelse(n > 1, n * (n - 1) / 2, 0),
      Count_Signed_IM = Count_Signed_IM, Count_Unsigned_IM = Count_Unsigned_IM,
      Count_Signed_IG = Count_Signed_IG, Count_Unsigned_IG = Count_Unsigned_IG,
      Count_Signed_CG = Count_Signed_CG, Count_Unsigned_CG = Count_Unsigned_CG,
      Count_Acc = Count_Acc,
      d_S_IM = d_S_IM, d_U_IM = d_U_IM, 
      d_S_IG = d_S_IG, d_U_IG = d_U_IG, 
      d_S_CG = d_S_CG, d_U_CG = d_U_CG, 
      d_Acc = d_Acc
    ) 
  })
  res$Dataset <- dataset_name
  return(res)
}

num_cores <- availableCores() - 2 
plan(multisession, workers = num_cores)
rds_files <- list.files(rds_dir, pattern = "\\.RDS$", full.names = TRUE)

message(">>> [Phase 1/3] 正在分布式计算 Top 5% 网络全景 (提取全量原始 Acc)...")
with_progress({
  p <- progressor(steps = length(rds_files)) 
  raw_results <- future_map_dfr(rds_files, evaluate_network_by_edge_type, p = p, .options = furrr_options(seed = TRUE))
})
plan(sequential)

# ==============================================================================
# 3. 结果合并与密度计算
# ==============================================================================
message(">>> [Phase 2/3] 正在缝合数据并计算网络拓扑密度...")

full_region_map <- Module_anno_with_prefix %>% 
  select(Community = Module, Region, SubRegion) %>% 
  distinct()

df_raw_weighted <- raw_results %>%
  left_join(full_region_map, by = "Community") %>%
  left_join(clade_df, by = "Dataset") %>%
  left_join(mod_data %>% select(Dataset, Community, Modularity), by = c("Dataset", "Community")) %>%
  filter(!is.na(Region), !is.na(Clade)) %>%
  
  mutate(
    Ratio_Silenced = if_else(Valid_Nodes < MIN_NODES | Modularity < 0, 1, 0),
    Total = if_else(Ratio_Silenced == 1, 1, d_S_IM + d_U_IM + d_S_IG + d_U_IG + d_S_CG + d_U_CG + d_Acc + 1e-15),
    
    # 拓扑密度
    Density_Topological = if_else(Ratio_Silenced == 1 | Possible_IM_Edges == 0, 0, 
                                  (Count_Signed_IM + Count_Unsigned_IM) / Possible_IM_Edges),
    
    # 加权密度
    Density_Weighted = if_else(Ratio_Silenced == 1 | Possible_IM_Edges == 0, 0, 
                               (d_S_IM + d_U_IM) / Possible_IM_Edges),
    
    Ratio_Signed_IM   = if_else(Ratio_Silenced == 1, 0, d_S_IM / Total),
    Ratio_Unsigned_IM = if_else(Ratio_Silenced == 1, 0, d_U_IM / Total),
    Ratio_Signed_IG   = if_else(Ratio_Silenced == 1, 0, d_S_IG / Total),
    Ratio_Unsigned_IG = if_else(Ratio_Silenced == 1, 0, d_U_IG / Total),
    Ratio_Signed_CG   = if_else(Ratio_Silenced == 1, 0, d_S_CG / Total),
    Ratio_Unsigned_CG = if_else(Ratio_Silenced == 1, 0, d_U_CG / Total),
    Ratio_Acc         = if_else(Ratio_Silenced == 1, 0, d_Acc / Total)
  ) %>% 
  select(-Total) %>%
  select(Dataset, Community, Clade, Region, SubRegion, Valid_Nodes, Modularity, Ratio_Silenced, 
         Possible_IM_Edges, starts_with("Density_"), starts_with("Count_"), starts_with("Ratio_"), starts_with("d_"))

# ==============================================================================
# 4. 终极双格式导出
# ==============================================================================
final_csv_name <- paste0("Raw_Weighted_7State_Top5PCT_withModularity_PCC", PCC_power, ".csv")
final_rds_name <- paste0("Raw_Weighted_7State_Top5PCT_withModularity_PCC", PCC_power, ".RDS")

write_csv(df_raw_weighted, file.path(output_target_dir, final_csv_name))
saveRDS(df_raw_weighted, file.path(output_target_dir, final_rds_name))

message(paste0("===========================================================\n",
               "🏆 全网最强底层数据生成完毕 (提取全量原始 Acc！)：\n",
               "CSV 保存至: ", final_csv_name, "\n",
               "RDS 保存至: ", final_rds_name, "\n",
               "==========================================================="))
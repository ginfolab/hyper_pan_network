library(dplyr)
library(stringr)
library(clusterProfiler)
library(readr)
library(writexl)

# ==============================================================================
# 1. 环境设置与路径定义
# ==============================================================================
base_path  <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
gephi_path <- file.path(base_path, "Gephi")
annot_path <- file.path(base_path, "Annotation")

gene_list <- readRDS(file.path(gephi_path, "Escherichia_coli_single_network_node_U15Gephi.RDS"))
names(gene_list) <- c("id","Module","Weight")

cut_type        <- "qvalue"
cutoff      <- 0.05
ego_module_size <- 5   
min_mod_size    <- 1   
U               <- 15

# ==============================================================================
# 2. 准备背景文件 (Annotation)
# ==============================================================================
term2gene_original <- read.csv(file.path(annot_path, 'Ec_gene2go'), header = F, sep = "\t", quote = "", fill = TRUE)
term2gene <- term2gene_original[term2gene_original$V2 %in% gene_list$id, ]    
term2name <- read.csv(file.path(annot_path, 'Ec_go2term'), header=F, sep="\t", quote="", fill=T)
colnames(term2gene) <- c("term", "gene")
colnames(term2name) <- c("term", "name")

# ==============================================================================
# 3. 构建 "Master List" (源头打标签，彻底避免串台 Bug)
# ==============================================================================
all_genes_df <- data.frame()

# --- A. 处理 Gephi 模块 ---
clusters <- readRDS(file.path(gephi_path, "Gephi_cluster_info.RDS"))
Final_step <- dim(clusters$memberships)[1]
module_stats <- readRDS(file.path(gephi_path, "Ec_module_modularities.RDS")) 

for (step in 1:Final_step) {
  if (step == 2) {
    cat(paste("Skipping Step", step, "as requested.\n"))
    next
  }
  
  g_file <- file.path(gephi_path, paste0("Escherichia_coli_single_network_node_U", U, "Gephi.RDS"))
  if(!file.exists(g_file)) next
  
  curr_genes <- readRDS(g_file)[, c("id", paste0("Step_", step))]
  colnames(curr_genes) <- c("Gene", "ID") 
  
  valid_ids <- module_stats %>% 
    filter(Step == paste0("S", step), Module_sizes >= min_mod_size) %>%
    pull(tag_final)
  
  curr_genes <- curr_genes %>% filter(ID %in% valid_ids)
  
  if(nrow(curr_genes) > 0) {
    # 【修复重点】：在此处直接根据 Step 赋予明确的 Type 和清理 ID
    if(step == 1) {
      curr_genes$Type <- "Module"
      curr_genes$ID <- str_replace(curr_genes$ID, "^S1_M", "Module_")
    } else if(step == 3) {
      curr_genes$Type <- "Group"
      curr_genes$ID <- str_replace(curr_genes$ID, "^S3_M", "Group_")
    } else {
      # 如果存在 S4, S5 等高阶结构，保留独立 Type，防止混入 Region
      curr_genes$Type <- paste0("Gephi_Step_", step)
    }
    all_genes_df <- rbind(all_genes_df, curr_genes)
  }
}

# --- B. 处理 Region ---
region_data <- readRDS(file.path(base_path, "Module_Gene_Region_Annotation.RDS")) %>%
  filter(!is.na(Region)) %>%
  dplyr::select(Gene = Ec_clustid, ID = SubRegion) %>% # 注：如果你想用大区，这里换成 ID = Region
  distinct() %>%
  mutate(Type = "Region") # 【修复重点】：直接打上 Region 标签

# --- C. 合并最终 Master List ---
master_input <- rbind(all_genes_df, region_data) %>% distinct()

# ==============================================================================
# 4. 执行富集分析
# ==============================================================================
enrichment_results <- data.frame()
unique_ids <- unique(master_input$ID)

cat("开始富集分析，共需处理 ID 数量:", length(unique_ids), "\n")

for (uid in unique_ids) {
  target_genes <- master_input %>% filter(ID == uid) %>% pull(Gene)
  target_genes <- target_genes[target_genes %in% term2gene$gene] 
  
  if(length(target_genes) < 1) next 
  
  ego <- tryCatch({
    enricher(target_genes, 
             TERM2GENE = term2gene, 
             TERM2NAME = term2name, 
             pAdjustMethod = "BH", 
             pvalueCutoff = 1,      
             qvalueCutoff = 1,
             minGSSize = 1,         
             maxGSSize = 10000)
  }, error = function(e) NULL)
  
  if (!is.null(ego) && nrow(ego@result) > 0) {
    res <- ego@result %>% 
      filter(Count >= ego_module_size, .data[[cut_type]] < cutoff) %>%
      filter(!grepl("Unknown|Others|unknown|Other", Description))
    
    if(nrow(res) > 0) {
      # 【语义优化】：原代码 res$Module，现在改为 Node_ID，更通用
      res$Node_ID <- uid 
      enrichment_results <- rbind(enrichment_results, res)
    }
  }
}

# ==============================================================================
# 5. 合并元数据与最终清洗
# ==============================================================================
if (nrow(enrichment_results) > 0) {
  
  region_map <- readRDS(file.path(base_path, "Module_region_map.RDS")) %>%
    rename(Module_info = Module)
  region_map$Module_info <- str_replace(region_map$Module_info, "^S1_M", "Module_")
  
  mod2grp_map <- readRDS(file.path(gephi_path, "Module_relation_node1.RDS")) %>% 
    dplyr::select(Module = Id, Group) %>% 
    filter(Group != "Regulon") %>%
    mutate(
      Module = str_replace(Module, "^S1_M", "Module_"),
      Group = str_replace(Group, "^S3_M", "Group_"),
      Group = str_replace(Group, "^S1_M", "Module_")
    )
  
  final_output <- enrichment_results %>%
    left_join(unique(master_input[, c("ID", "Type")]), by = c("Node_ID" = "ID")) %>%
    # 注意：这里的映射仅对 Type == "Module" 有意义
    left_join(region_map, by = c("Node_ID" = "Module_info")) %>%
    left_join(mod2grp_map, by = c("Node_ID" = "Module")) %>%
    mutate(
      Group  = ifelse(Type == "Module", Group, NA),
      Region = ifelse(Type == "Module", Region, NA)
    ) %>%
    mutate(Term_type = case_when(
      # 【建议】：如果 term2name 的 ID 列有数据库前缀，最好改用 grepl(..., ID)
      grepl("tigr", Description, ignore.case = TRUE)    ~ "TIGR",
      grepl("Regulon", Description, ignore.case = TRUE) ~ "Regulon",
      grepl("cog", Description, ignore.case = TRUE)     ~ "COG",
      grepl("kegg", Description, ignore.case = TRUE)    ~ "KEGG",
      grepl("mbgd", Description, ignore.case = TRUE)    ~ "MBGD",
      grepl("Operon", Description, ignore.case = TRUE)  ~ "Operon",
      grepl("Kmodule", Description, ignore.case = TRUE) ~ "Kmodule",
      TRUE ~ "GO"
    )) %>%
    # 调整列顺序
    dplyr::select(Node_ID, Type, Description, Term_type, qvalue, Count, Group, Region, geneID, everything()) %>%
    arrange(Type, Node_ID, qvalue)
  
  output_rds <- file.path(base_path, paste0("ALL_module_(Count_ge_", ego_module_size, "_",cut_type,"_lt_", cutoff, ")_S2_Removed.RDS"))
  saveRDS(final_output, output_rds)
  
  output_xlsx <- file.path(base_path, paste0("ALL_module_(Count_ge_", ego_module_size, "_",cut_type,"_lt_", cutoff, ")_S2_Removed.xlsx"))
  write_xlsx(final_output, output_xlsx)
  
  cat("分析完成！\n")
  print(table(final_output$Type)) 
  
} else {
  cat("未找到符合条件的富集结果。\n")
}




################################################################################
################################################################################
############################# Table 3 summarize ################################
################################################################################
################################################################################
library(dplyr)
library(stringr)
library(writexl)
library(readr)

# ==============================================================================
# 0. 载入外部物理骨架骨架与顺序文件
# ==============================================================================
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
module_order_file <- file.path(base_path, "Module_Region_Mapping.txt")

# 转换为标准明文格式并读取
module_order_df <- read.table(module_order_file, header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(Module = str_replace(Module, "^S1_M", "Module_"))

# 提取严格的物理模块排序向量
exact_module_sequence <- module_order_df$Module

# ==============================================================================
# 1. 准备工作：确保 final_output 包含 Term_type 列，并格式化 Node_ID 为 Module_x
# ==============================================================================
table3_data <- final_output %>%
  mutate(
    # 确保 Node_ID 的格式与外部骨架完全一致
    Node_ID = str_replace(Node_ID, "^S1_M", "Module_"),
    Term_type = case_when(
      grepl("tigr", Description, ignore.case = TRUE)    ~ "TIGR",
      grepl("Regulon", Description, ignore.case = TRUE) ~ "Regulon",
      grepl("cog", Description, ignore.case = TRUE)     ~ "COG",
      grepl("kegg", Description, ignore.case = TRUE)    ~ "KEGG",
      grepl("mbgd", Description, ignore.case = TRUE)    ~ "MBGD",
      grepl("Operon", Description, ignore.case = TRUE)  ~ "Operon",
      grepl("Kmodule", Description, ignore.case = TRUE) ~ "Kmodule",
      TRUE ~ "Other" # 主要是 GO
    )
  )

# ==============================================================================
# 2. 核心逻辑：排序、汇总与 Representation 计算
# ==============================================================================

# 定义优先级列表
priority_1 <- c("KEGG", "Kmodule", "TIGR", "COG") # 按你注释中要求的准确代谢到大分类排序
priority_2 <- c("Regulon", "Operon", "MBGD")

summary_table_raw <- table3_data %>%
  # --- 第一步：组内 qvalue 排序 ---
  arrange(Node_ID, qvalue) %>%
  
  # --- 第二步：分组汇总 ---
  group_by(Node_ID, Type, Group) %>% # 暂时不放 Region，后面用外部文件严格对齐
  summarise(
    # 1. 保留你原始的各个分类汇总 (用分号连接)
    TIGR_Info    = paste(Description[Term_type == "TIGR"],    collapse = "; "),
    COG_Info     = paste(Description[Term_type == "COG"],     collapse = "; "),
    KEGG_Info    = paste(Description[Term_type == "KEGG"],    collapse = "; "),
    MBGD_Info    = paste(Description[Term_type == "MBGD"],    collapse = "; "),
    Regulon_Info = paste(Description[Term_type == "Regulon"], collapse = "; "),
    Operon_Info  = paste(Description[Term_type == "Operon"],  collapse = "; "),
    Kmodule_Info = paste(Description[Term_type == "Kmodule"], collapse = "; "),
    
    # 2. 保留你原始的 Representation (核心难点计算)
    Representation = {
      curr_type <- Term_type
      curr_desc <- Description
      idx_p1 <- which(curr_type %in% priority_1)
      
      if (length(idx_p1) > 0) {
        curr_desc[idx_p1[1]]
      } else {
        idx_p2 <- which(curr_type %in% priority_2)
        if (length(idx_p2) > 0) {
          curr_desc[idx_p2[1]]
        } else {
          NA_character_
        }
      }
    },
    
    Best_qvalue = min(qvalue),
    .groups = "drop" 
  ) %>%
  
  # --- 第三步：清理空字符串 ---
  mutate(across(ends_with("_Info"), ~ifelse(. == "", NA, .))) %>% 
  
  # ==============================================================================
# 【附加步骤】：智能合成 Final_Summary (面向读者的最终短句)
# ==============================================================================
mutate(
  # 提取纯粹的核心功能 (去除数字前缀)
  Core_Function = case_when(
    !is.na(KEGG_Info)    ~ str_extract(KEGG_Info, "(?<=:)[^;]+"),
    !is.na(Kmodule_Info) ~ str_extract(Kmodule_Info, "(?<=:)[^;]+") %>% str_remove(", prokaryotes and chloroplasts"),
    !is.na(TIGR_Info)    ~ str_extract(TIGR_Info, "(?<=:)[^;]+"),
    !is.na(COG_Info)     ~ str_extract(COG_Info, "(?<=:)[^;]+"),
    TRUE ~ "Uncharacterized/Unknown Function"
  ),
  Core_Function = str_trim(Core_Function),
  
  # 提取转录调控因子
  Main_Regulator = if_else(
    !is.na(Regulon_Info), 
    str_extract(Regulon_Info, "(?<=Regulon:)[^;]+"), 
    NA_character_
  ),
  
  # 合成短句
  Final_Summary = case_when(
    !is.na(Main_Regulator) & Core_Function != "Uncharacterized/Unknown Function" ~ 
      paste0(Core_Function, " (Regulated by ", Main_Regulator, ")"),
    
    !is.na(Main_Regulator) & Core_Function == "Uncharacterized/Unknown Function" ~ 
      paste0("Unknown Function (Regulated by ", Main_Regulator, ")"),
    
    TRUE ~ Core_Function
  )
)

# ==============================================================================
# 3. 严格执行外部骨架对齐与排序 (🌟核心修改点)
# ==============================================================================
summary_table <- module_order_df %>%
  # 外部文件的 Module 与汇总表的 Node_ID 严格合并，顺带补齐 SubRegion 和 Region
  left_join(summary_table_raw, by = c("Module" = "Node_ID")) %>%
  rename(Node_ID = Module) %>%
  
  # 按照物理文件中的严格行顺序进行重排
  arrange(match(Node_ID, exact_module_sequence)) %>%
  
  # 如果 enrichment 结果中不存在某个 module，在这里会被保留，其余列填 NA。如果需要剔除未注到的，取消注释下面这行：
  # filter(!is.na(Type)) %>%
  
  # ==============================================================================
# --- 第五步：最后的列排序 (完美保留了你的心血，并插入 SubRegion) ---
# ==============================================================================
dplyr::select(
  Node_ID, Type, Group, Region, SubRegion,   # <-- 🌟 顺序第一梯队，加入了 SubRegion
  Representation,                           # <-- 这是你原本提取的代表性功能 (保留)
  Final_Summary,                            # <-- 这是新增的一句话短句总结
  Best_qvalue, 
  TIGR_Info, COG_Info, KEGG_Info, MBGD_Info,                 # <-- 详细分类 (保留)
  Regulon_Info, Operon_Info, Kmodule_Info                  # <-- 详细分类 (保留)
)

# ==============================================================================
# 4. 保存结果
# ==============================================================================

# 保存为 RDS
saveRDS(summary_table, file.path(base_path, "Table3_Module_Functional_Summary.RDS"))

# 保存为 Excel
write_xlsx(summary_table, file.path(base_path, "Table3_Module_Functional_Summary.xlsx"))

cat("Table 3 生成完毕！已严格按照外部 Mapping 序列排序并补齐 SubRegion！\n")

# 打印查看最重要的五个对比列验证顺序与增列
print(head(summary_table %>% dplyr::select(Node_ID, Region, SubRegion, Representation, Final_Summary)))







#####给cytoscape 的网络上色node 文件添加 functional annotation 信息 #####
####################################################################################
####################################################################################给cytoscape添加 functional annotation 信息
####################################################################################

# ==============================================================================
# . 生成最终表
# ==============================================================================

MBGD_anno <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/MBGD_anno.txt", 
                        header = T, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE) %>%
  mutate(anno_short = ifelse(nchar(Descri) > 10, 
                             paste0(substr(Descri, 1, 10), "..."), 
                             Descri)) %>%
  filter(!is.na(color_name) & color_name != "") %>%
  # 3. 新增步骤：删除 color_name 后面的数字部分
  mutate(color_name = gsub("[0-9]+", "", color_name))%>% 
  rename(Annotation = category) %>%
  select(-Descri) %>% 
  distinct(Annotation, .keep_all = TRUE) %>%
  ######## 只选择mbgd anno
  filter(grepl("mbgd", Annotation, ignore.case = TRUE))
  
  
Ec_anno <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Ec_gene_anno.txt", header = T, sep = "\t", quote = "", fill = TRUE)
Gene_Module <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Escherichia_coli_single_network_node_U15Gephi.RDS") %>%
  dplyr::select(Gene = id, Module = Step_1)
Module_Group <- read.table('/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新).txt', header = TRUE, sep = "\t", check.names = FALSE, quote = "") %>%  dplyr::select(Module = Node, Group )
Module_region <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_region.rds")  %>% dplyr::select(Module = Module_Name , Region) %>% distinct()


Gene_Module_anno <- Gene_Module %>% left_join(Module_Group, by = c("Module")) %>% left_join(Module_region, by = c( "Module")) %>%
  left_join(Ec_anno, by = c("Gene" = "Ec_clustid")) %>%
  select(Gene, Module ,Group,Region, KEGG = kegg_anno, COG = cog_anno, TIGR = tigr_anno, MBGD = mbgd_anno)
saveRDS(Gene_Module_anno, "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_Gene_Region_Annotation.RDS")

#########选取module最大的 注释信息 ##########
library(dplyr)
library(tidyr)
library(stringr)
# ==============================================================================
# 生成 Module Summarize
# ==============================================================================

# 1. 数据重塑：将四列注释合并为一列，方便统计
long_anno <- Gene_Module_anno %>%
  pivot_longer(
    cols = c(KEGG, COG, TIGR, MBGD), 
    names_to = "Source", 
    values_to = "Annotation"
  ) %>%
  # 去除空值和空字符串
  filter(!is.na(Annotation) & Annotation != "") %>%
  left_join(MBGD_anno, by = "Annotation") %>%
  filter(!is.na(anno_short))  # 只保留在 MBGD_anno 中有对应注释的行)

# 2. 统计每个 Module 内部，各 Annotation 的出现次数
module_anno_counts <- long_anno %>%
  count(Module, Group, Region, Annotation,anno_short,color_name, name = "Freq") %>%
  arrange(Module, desc(Freq)) %>%
  filter(!is.na(color_name))
# 3. 找出每个 Module 的 Top Candidates (出现次数最多的注释)
# 注意：这一步可能会让一个 Module 保留多行（如果有平局）
top_candidates <- module_anno_counts %>%
  group_by(Module) %>%
  filter(Freq == max(Freq)) %>%
  ungroup()

# 4. 计算 Group 级别的“流行度” (用于 Tie-breaking)
# 逻辑：统计某个 Annotation 在该 Group 下多少个 Module 中都是 Top Candidate
group_popularity <- top_candidates %>%
  # 确保每个 Module 对每个 Annotation 只投一票 (虽然 filter 后通常唯一，但保险起见)
  distinct(Group, Module, Annotation) %>% 
  count(Group, Annotation, name = "Group_Score")

# 5. 最终筛选：合并 Group 分数并决出唯一的 Summarize_anno
# 5. 最终筛选：合并 Group 分数并决出唯一的 Summarize_anno
Module_summarize <- top_candidates %>%
  # 关联 Group 分数
  left_join(group_popularity, by = c("Group", "Annotation")) %>%
  group_by(Module) %>%
  # 添加一个随机数列，用于最后的随机二选一
  mutate(Random_Tie_Breaker = runif(n())) %>%
  # 排序优先级：Module(Cyto)_relation_node(最新修改label和颜色)
  # 1. Freq (Module内出现次数，核心标准)
  # 2. Group_Score (Group内流行度，您的主要决胜标准)
  # 3. Random (随机兜底)
  arrange(desc(Freq), desc(Group_Score), Random_Tie_Breaker) %>%
  # 取第一名
  slice(1) %>%
  ungroup() %>%
  # 选择最终列 (此时 Summarize_anno 是完全原始的、带小数点的)
  select(Module, Region, anno_short, Function_color = color_name)

# 先读取
Cyto_node <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新).txt", 
                              header = TRUE, 
                              sep = "\t", 
                              check.names = FALSE,
                              quote = "") # 加上 quote="" 防止引号导致的读取错误

# 强制把最后一列改名为 Module_color_mapping
# (假设它是最后一列，根据你之前的 names 输出)
colnames(Cyto_node)[ncol(Cyto_node)] <- "Module_color_mapping"


##手动选择的module
manual_modules <- c("S1_M1", "S1_M307", "S1_M44", "S1_M179", "S1_M11")
# 然后再合并
Cyto_node_final <- Cyto_node %>%
  left_join(Module_summarize, by = c("Node" = "Module")) %>%
 # rename(Function_color = color_name) %>%
  mutate(
    Region_color = case_when(
      Region == "Region_1" ~ "hotpink" ,  # 紫色
      Region == "Region_2" ~ "darkolivegreen",  # 粉色/洋红
      Region == "Region_3" ~ "mediumpurple" ,  # 橄榄绿
      Region == "Region_4" ~ "turquoise",  # 青色/蓝绿
      TRUE ~ "white"                      # 如果有其他Region，默认标灰
    ),
    # 新增逻辑：如果 Module_color_mapping 是 grey，则改为 white，否则保持原样
    Module_color_mapping = if_else(Module_color_mapping == "grey", "white", Module_color_mapping),
    Module_color_mapping = if_else(Node %in% manual_modules, "lightgrey", Module_color_mapping),
    Function_color = if_else(is.na(Function_color), "white", Function_color), # 如果没有注释，Function_color 也标灰
   # Function_color = if_else(Function_color == "gray", "white", Function_color),
    Node_Id = str_remove(Node_name, "Module_"),
    anno_short = if_else(is.na(anno_short), Node_Id, anno_short)
  )  

write.table(Cyto_node_final,'/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新修改label和颜色).txt', sep = "\t", row.names = F, quote = F)  
  


  
##### 画图显示颜色图例 #####
library(ggplot2)
library(dplyr)
MBGD_color <- MBGD_anno %>% filter(!str_detect(Annotation, "\\.")) 
# Node_name,Module_color_mapping   
# Region,Region_color   
# anno_short,Function_color

# 1. 提取唯一的 Region 和对应的颜色
legend_data <- MBGD_color %>%
  #filter(!Region_color == "white") %>%    ##修改
  #select(color_name,anno_short) %>%   ##修改  
  distinct() %>% 
  arrange(Annotation)     ##修改 

# 2. 画一个只有图例的图
p_legend <- ggplot(legend_data, aes(x = 1, y = anno_short, fill = color_name)) +        #修改
  geom_tile(width = 0.5, height = 0.5) +  # 画色块
  scale_fill_identity() +                 # 关键：告诉 ggplot 直接使用列里的颜色名称，不要自动分配
  geom_text(aes(label = anno_short), x = 1.6, hjust = 0, size = 5) + # 添加文字标签                     #修改
  theme_void() +                          # 去掉背景坐标轴
  xlim(0.5, 3) +                          # 调整画布范围
  labs(title = "MBGD Annotation Color Legend" )

# 3. 显示
print(p_legend)


  
  
  
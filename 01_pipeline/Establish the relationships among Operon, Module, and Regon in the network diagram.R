# ==============================================================================
# 0. 准备数据 (示例数据，实际运行请使用你的变量)
# ==============================================================================
# 设定阈值
Strict_Simp_Cut <- 0.8 
Perfect_F_Cut   <- 0.8 
#Overlap_F_Cut   <- 0.2 
minModuleSize <- 5  # 最小模块尺寸阈值
minOverlapSize <- 2  # 最小重叠基因数阈值

library(dplyr)
library(stringr) # 如果需要 str_detect
library(tidyverse)
library(stringr)
Core_gene_module <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/clustid2Module2Step.RDS") %>%
  select(Ec_id = id, Module = Step_1) %>%
  mutate(Module = str_replace_all(Module, "S1_M", "Module_")) %>%
  group_by(Module) %>%
  summarize(Ec_id = paste(Ec_id, collapse = ",")) 

Core_gene_other_info <- read.table('/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新).txt',header = T,sep = "\t") %>%
  select(Module = Node, Size, Group) %>%
  mutate(Module = str_replace_all(Module, "S1_M", "Module_")) 
Core_network_node_table <- Core_gene_module %>% left_join(Core_gene_other_info, by = "Module")


Function_term_table <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Regulon&Operon_Node.txt",header = T,sep = "\t") %>%
  filter(Type %in% c("Regulon","Operon")) %>%
  rename(Term = Node) 


Hyperedge <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Escherichia_coli_closed_U15.RDS")



# ==============================================================================
# 1. 展开数据 & 准备全局尺寸
# ==============================================================================
message("Step 1: 数据展开与预处理...")

# 展开 Module 表
module_genes <- Core_network_node_table %>%
  filter(Size >= minModuleSize) %>%  # 【修改点】在这里加入过滤，只保留尺寸 >= 5 的模块
  select(Module, Module_Size = Size, Ec_id) %>%
  separate_rows(Ec_id, sep = ",\\s*") %>% 
  mutate(Ec_id = str_trim(Ec_id))

# 展开 Term 表
term_genes <- Function_term_table %>%
  select(Term, Type, Term_Size = Size, members) %>%
  separate_rows(members, sep = ",\\s*") %>%
  rename(Ec_id = members) %>%
  mutate(Ec_id = str_trim(Ec_id))

# 准备全局尺寸查找表
Global_Size_Map <- term_genes %>%
  distinct(Term, Type, Term_Size)

# 展开超边数据 (为后面做准备)
hyper_long <- Hyperedge %>%
  mutate(HID = row_number()) %>% 
  select(HID, Universality, Itemsets) %>%
  separate_rows(Itemsets, sep = ",\\s*") %>% 
  rename(Ec_id = Itemsets) %>%
  mutate(Ec_id = str_trim(Ec_id))

# ==============================================================================
# 2. 计算筛选指标 (生成 Result_Mapping)
# ==============================================================================
message("Step 2: 计算 F-score 并筛选 (生成 Valid Operons)...")



# 计算初步映射
Result_Mapping_Raw <- module_genes %>%
 
  inner_join(term_genes, by = "Ec_id") %>% 
  group_by(Module, Term, Type, Module_Size, Term_Size) %>%
  summarise(Overlap_Count = n(), .groups = "drop") %>% 
  #决定了 “谁有资格成为一个独立的点（Node）”
  filter(Overlap_Count >= minOverlapSize) %>%      
  mutate(
    Precision = Overlap_Count / Module_Size,
    Recall = Overlap_Count / Term_Size,
    F_score = (2 * Precision * Recall) / (Precision + Recall)
  ) %>%
  mutate(Mapping_state = case_when(
    F_score >= Perfect_F_Cut ~ "Perfect_Match",
    Precision >= Strict_Simp_Cut ~ "Module_is_Subset",
    Recall >= Strict_Simp_Cut ~ "Module_Contains_Term",
    Overlap_Count >= minOverlapSize ~ "Partial Overlap",
    TRUE ~ "No_mapping"
  )) %>%
  filter(Mapping_state != "No_mapping") %>%
  select(Module, Term, Term_Size, Type, Overlap_Count, F_score, Precision, Recall, Mapping_state)

# ==============================================================================
# 3. 【核心修改】重构目标集合 (Valid Operons vs Other_genes)
# ==============================================================================
message("Step 3: 重构 Target 集合 (将无效 Operon 归入 Other_genes)...")

# 3.1 提取 "有效 Operon" (Valid Operons)
# 只有在 Result_Mapping_Raw 中出现的 Operon 才保留
Valid_Operons_List <- Result_Mapping_Raw %>%
  filter(Type == "Operon") %>%
  select(Module, Target_Name = Term) %>%
  distinct() %>%
  mutate(Target_Type = "Operon")

# 找回有效 Operon 包含的基因
Set_Valid_Operons <- Valid_Operons_List %>%
  inner_join(term_genes, by = c("Target_Name" = "Term"),relationship = "many-to-many") %>%
  select(Module, Target_Name, Target_Type, Ec_id)

# 3.2 定义 "Other_genes" (原 No_Operon + 被淘汰的 Operon)
# 逻辑：Module 中的基因，如果不在 Set_Valid_Operons 里，就归为 Other_genes
Genes_In_Valid <- Set_Valid_Operons %>%
  select(Module, Ec_id) %>%
  distinct() %>%
  mutate(Is_Valid = TRUE)

Set_Other_Genes <- module_genes %>%
  left_join(Genes_In_Valid, by = c("Module", "Ec_id")) %>%
  filter(is.na(Is_Valid)) %>% # 没匹配上的
  select(Module, Ec_id) %>%
  mutate(
    Target_Name = "Other_genes",  # 【改名】这里改为 Other_genes
    Target_Type = "Other_genes"
  )

# 3.3 合并生成新的 All_Targets (这将是后续画图的唯一依据)
New_All_Targets <- bind_rows(Set_Valid_Operons, Set_Other_Genes)

# 计算新的 Local Overlap (用于 Node ID 和 Size)
New_Target_Local_Stats <- New_All_Targets %>%
  group_by(Module, Target_Name, Target_Type) %>%
  summarise(Local_Overlap_Count = n(), .groups = "drop")

# ==============================================================================
# 4. 准备 Cytoscape 节点映射 (基于 New_All_Targets)
# ==============================================================================
message("Step 4: 准备 Cytoscape ID 映射...")

# 4.1 定义 ID 生成规则
Generate_Node_ID <- function(name, type, module, overlap) {
  if_else(type == "Other_genes",
          paste0(module, "_", name), # Module_1_Other_genes
          paste0(module, "_", name, "(", overlap, ")")) # Module_1_lacZYA(3)
}

# 4.2 构建节点-基因字典
# A. Regulon (来自 Result_Mapping 中出现的 Regulon)
# 注意：Result_Mapping_Raw 里包含 Regulon 的对应关系，我们需要提取出来
Valid_Regulons <- Result_Mapping_Raw %>% filter(Type == "Regulon") %>% select(Regulon_Name = Term) %>% distinct()

Dict_Regulon <- Valid_Regulons %>%
  inner_join(term_genes, by = c("Regulon_Name" = "Term")) %>%
  select(Node_ID = Regulon_Name, Ec_id) %>%
  distinct()

# B. Target (基于 New_All_Targets)
Dict_Target <- New_All_Targets %>%
  left_join(New_Target_Local_Stats, by = c("Module", "Target_Name", "Target_Type")) %>%
  mutate(Node_ID = Generate_Node_ID(Target_Name, Target_Type, Module, Local_Overlap_Count)) %>%
  select(Node_ID, Ec_id) %>%
  distinct()

# 合并字典
Dict_All_Nodes <- bind_rows(Dict_Regulon, Dict_Target) %>% distinct()

# 4.3 映射超边
Node_Hyper_Map <- Dict_All_Nodes %>%
  inner_join(hyper_long, by = "Ec_id") %>%
  select(Node_ID, HID, Universality) %>%
  distinct()

# ==============================================================================
# 5. 构建边表 (Edge Table)
# ==============================================================================
message("Step 5: 构建网络边表 (仅包含 Valid Operons 和 Other_genes)...")

# --- 5.1 构建 Regulon 调控边 ---
# 我们需要重新计算 Regulon 和 New_All_Targets 之间的 F-score
# 因为之前的 Result_Mapping 是针对所有 Term 的，现在我们要把它映射到具体的 Node ID 上

# 提取 Result_Mapping 中关于 Regulon 的部分
Regulon_Edges_Base <- Result_Mapping_Raw %>%
  filter(Type == "Regulon") %>%
  select(Module, Regulon_Name = Term, F_score, Overlap_Count)

# 关键：我们需要知道这个 Regulon 连的是哪个 Valid Operon (或 Other_genes?)
# 其实 Result_Mapping 本身并没有记录 Regulon -> Operon 的直接关系（它是 Regulon -> Module 的关系吗？不对）
# 等等，Result_Mapping 的每一行是一个 (Module, Term) 对。如果 Term 是 Regulon，这只说明 Regulon 和 Module 有重叠。
# 它并没有说明这个 Regulon 调控了 Module 里的哪个 Operon。
# 【必须重新计算 Regulon -> New_All_Targets 的具体连线】

# 重新计算 Regulon -> New_All_Targets 的重叠
Edges_Regulon_Target <- Valid_Regulons %>%
  inner_join(term_genes, by = c("Regulon_Name" = "Term")) %>% # 展开 Regulon 基因
  select(Regulon_Name, Ec_id) %>%
  inner_join(New_All_Targets, by = "Ec_id",relationship = "many-to-many") %>% # 匹配到 Module 内部的具体 Target
  group_by(Module, Regulon_Name, Target_Name, Target_Type) %>%
  summarise(Overlap_Count = n(), .groups = "drop") %>%
  # 关联尺寸计算 F-score
  left_join(Global_Size_Map %>% select(Regulon_Name = Term, R_Size = Term_Size), by = "Regulon_Name") %>%
  left_join(New_Target_Local_Stats, by = c("Module", "Target_Name", "Target_Type")) %>% # Target Size 用 Local
  left_join(Global_Size_Map %>% select(Target_Name = Term, O_Size_Global = Term_Size), by = "Target_Name") %>%
  mutate(
    # 分母：Regulon 用全局，Operon 用全局，Other_genes 用局部
    Denom_Target = if_else(Target_Type == "Operon", as.numeric(O_Size_Global), as.numeric(Local_Overlap_Count)),
    Precision = Overlap_Count / as.numeric(R_Size),
    Recall = Overlap_Count / Denom_Target,
    F_score = (2 * Precision * Recall) / (Precision + Recall)
  ) %>%
  # 再次筛选 (确保连线质量) 决定了 “谁有资格画上一条连线（Edge）”
  filter(
   # F_score >= Overlap_F_Cut, 
    Overlap_Count >= minOverlapSize) %>%           
  mutate(
    Source_Node = Regulon_Name,
    Target_Node = Generate_Node_ID(Target_Name, Target_Type, Module, Local_Overlap_Count),
    Interaction = "Regulates",
    Module_Context = Module
  ) %>%
  select(Source_Node, Target_Node, F_score, Overlap_Count, Module_Context, Interaction)

# 计算 Regulon 边的超边权重
Regulon_Weights <- Edges_Regulon_Target %>%
  select(Source_Node, Target_Node) %>% distinct() %>%
  inner_join(Node_Hyper_Map, by = c("Source_Node" = "Node_ID"), relationship = "many-to-many") %>%
  rename(HID_S = HID, U_S = Universality) %>%
  inner_join(Node_Hyper_Map, by = c("Target_Node" = "Node_ID", "HID_S" = "HID")) %>%
  group_by(Source_Node, Target_Node) %>%
  summarise(Median_Universality = median(U_S), .groups = "drop")

Cyto_Edge_Regulon_Weighted <- Edges_Regulon_Target %>%
  left_join(Regulon_Weights, by = c("Source_Node", "Target_Node")) %>%
  mutate(Median_Universality = replace_na(Median_Universality, 0))

# --- 5.2 挖掘 Module 内部横向连接 ---
Internal_Targets_Info <- New_All_Targets %>%
  left_join(New_Target_Local_Stats, by = c("Module", "Target_Name", "Target_Type")) %>%
  mutate(Node_ID = Generate_Node_ID(Target_Name, Target_Type, Module, Local_Overlap_Count)) %>%
  select(Node_ID, Parent_Module = Module) %>% distinct() %>%
  inner_join(Node_Hyper_Map, by = "Node_ID")

Internal_Edges_Calc <- Internal_Targets_Info %>%
  rename(Node_A = Node_ID, U_Val = Universality) %>%
  inner_join(Internal_Targets_Info %>% select(Node_B = Node_ID, Parent_Module, HID, U_Val = Universality), 
             by = c("Parent_Module", "HID", "U_Val"), relationship = "many-to-many") %>%
  filter(Node_A < Node_B) %>% 
  group_by(Node_A, Node_B, Parent_Module) %>%
  summarise(Median_Universality = median(U_Val), .groups = "drop")

Cyto_Edge_Internal <- Internal_Edges_Calc %>%
  mutate(
    Source_Node = Node_A, Target_Node = Node_B,
    F_score = 0, Overlap_Count = 0,
    Module_Context = Parent_Module, Interaction = "Hyperedge_Connection"
  ) %>%
  select(colnames(Cyto_Edge_Regulon_Weighted))

# 合并所有边
Cyto_Edge_Combined <- bind_rows(Cyto_Edge_Regulon_Weighted, Cyto_Edge_Internal)

# ==============================================================================
# 6. 构建节点表 (Node Table)
# ==============================================================================
message("Step 6: 构建节点表 (Target=LocalSize, Regulon=GlobalSize)...")

# A. Target 节点 (Valid Operons + Other_genes)
Nodes_Targets <- New_All_Targets %>%
  distinct(Module, Target_Name, Target_Type) %>%
  left_join(New_Target_Local_Stats, by = c("Module", "Target_Name", "Target_Type")) %>%
  mutate(
    Node_ID = Generate_Node_ID(Target_Name, Target_Type, Module, Local_Overlap_Count),
    Node_Type = Target_Type,
    Parent_Module = Module,
    Display_Label = if_else(Target_Type == "Other_genes", Target_Name, paste0(Target_Name, "(", Local_Overlap_Count, ")")),
    Size = Local_Overlap_Count # 使用局部大小
  ) %>%
  select(Node_ID, Node_Type, Parent_Module, Display_Label, Size)

# B. Regulon 节点
Nodes_Regulons <- data.frame(Node_ID = unique(Cyto_Edge_Combined$Source_Node)) %>%
  filter(Node_ID %in% Cyto_Edge_Regulon_Weighted$Source_Node) %>% 
  mutate(Node_Type = "Regulon", Parent_Module = NA, Display_Label = Node_ID) %>%
  left_join(Global_Size_Map %>% select(Node_ID = Term, R_Size = Term_Size), by = "Node_ID") %>%
  mutate(Size = coalesce(as.numeric(R_Size), 5)) %>%
  select(Node_ID, Node_Type, Parent_Module, Display_Label, Size)

Color_Map <- c(
  "Regulon"     = "#377EB8", 
  "Operon"      = "#C62828", 
  "Other_genes" = "#B0BEC5"  
)

Cyto_Node_Final <- bind_rows(Nodes_Targets, Nodes_Regulons) %>% 
  distinct(Node_ID, .keep_all = TRUE) %>%
  # 【关键修改】添加 Color 列
  mutate(
    Node_Color = case_when(
      Node_Type == "Regulon" ~ Color_Map["Regulon"],
      Node_Type == "Operon"  ~ Color_Map["Operon"],
      Node_Type == "Other_genes" ~ Color_Map["Other_genes"],
      TRUE ~ "#000000" # 默认黑色 (防错)
    )
  )
# ==============================================================================
# 7. 处理孤立点与导出
# ==============================================================================
message("Step 7: 处理孤立点并导出...")

Connected_Nodes <- unique(c(Cyto_Edge_Combined$Source_Node, Cyto_Edge_Combined$Target_Node))
Isolated_Nodes <- Cyto_Node_Final %>% filter(!Node_ID %in% Connected_Nodes)

if(nrow(Isolated_Nodes) > 0) {
  Dummy_Edges <- Isolated_Nodes %>%
    mutate(
      Source_Node = Node_ID, Target_Node = Node_ID,
      F_score = 0, Overlap_Count = 0, Median_Universality = 0,
      Module_Context = Parent_Module, Interaction = "Isolated"
    ) %>%
    select(colnames(Cyto_Edge_Combined))
  Cyto_Edge_Final_Export <- bind_rows(Cyto_Edge_Combined, Dummy_Edges)
} else {
  Cyto_Edge_Final_Export <- Cyto_Edge_Combined
}

Out_Path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/"
if(!dir.exists(Out_Path)) dir.create(Out_Path, recursive = TRUE)

write.table(Cyto_Node_Final, paste0(Out_Path, "Cyto_Node_Table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Cyto_Edge_Final_Export, paste0(Out_Path, "Cyto_Edge_Table_Complete.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

message("完成！网络已简化：仅保留符合筛选条件的 Operon，其余归为 Other_genes。")






# ==============================================================================
# 8. [新增] 基于 color_mapping_df 过滤模块  TOP11
# ==============================================================================
message("Step 8: 基于 color_mapping_df 过滤最终网络...")
color_mapping_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Top_Module_color_S1.RDS") %>%
  select(Module_N,color_mapping)
  
# 假设 color_mapping_df 已经存在，并且包含 Module_N 列
# 如果没有 Module_N，请先运行: 
# color_mapping_df <- color_mapping_df %>% mutate(Module_N = str_replace_all(Module, "S1_M", "Module_"))

# 1. 获取目标模块列表
Target_Modules <- unique(color_mapping_df$Module_N)
message(paste("目标保留模块数量:", length(Target_Modules)))

# 2. 过滤 Target 节点 (Operon 和 Other_genes)
# 逻辑：只保留 Parent_Module 在目标列表里的节点
Filtered_Nodes_Targets <- Cyto_Node_Final %>%
  filter(Node_Type != "Regulon") %>%
  filter(Parent_Module %in% Target_Modules)

# 3. 过滤 Regulon 节点
# 逻辑：Regulon 本身没有模块归属，但如果它连接了被保留的 Target，它就应该被保留
# 先找出被保留 Target 的 ID
Keep_Target_IDs <- Filtered_Nodes_Targets$Node_ID

# 从边表中找出连接这些 Target 的 Regulon
Related_Regulons <- Cyto_Edge_Final_Export %>%
  filter(Target_Node %in% Keep_Target_IDs) %>%
  filter(Source_Node %in% Cyto_Node_Final$Node_ID[Cyto_Node_Final$Node_Type == "Regulon"]) %>%
  pull(Source_Node) %>%
  unique()

Filtered_Nodes_Regulons <- Cyto_Node_Final %>%
  filter(Node_ID %in% Related_Regulons)

# 4. 合并最终节点表
Cyto_Node_Filtered <- bind_rows(Filtered_Nodes_Targets, Filtered_Nodes_Regulons)

# 5. 过滤边表
# 逻辑：Source 和 Target 必须都在过滤后的节点表中
Valid_Node_IDs <- Cyto_Node_Filtered$Node_ID

Cyto_Edge_Filtered <- Cyto_Edge_Final_Export %>%
  filter(Source_Node %in% Valid_Node_IDs & Target_Node %in% Valid_Node_IDs)





# ==============================================================================
# 9. 导出过滤后的文件
# ==============================================================================
Out_Path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/"
if(!dir.exists(Out_Path)) dir.create(Out_Path, recursive = TRUE)

# 导出节点表
write.table(Cyto_Node_Filtered, 
            paste0(Out_Path, "Cyto_Node_Table_Filtered.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 导出边表
write.table(Cyto_Edge_Filtered, 
            paste0(Out_Path, "Cyto_Edge_Table_Filtered.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

message(paste("过滤完成！"))
message(paste("原节点数:", nrow(Cyto_Node_Final), "-> 新节点数:", nrow(Cyto_Node_Filtered)))
message(paste("原边数:", nrow(Cyto_Edge_Final_Export), "-> 新边数:", nrow(Cyto_Edge_Filtered)))

"#377EB8"
"#E41A1C"
"#B0BEC5"
"#2A9D8F"
"#E41A1C"
"#BDC3C7"
"#F1C40F"


library(tidyverse)
library(tidygraph)
library(ggraph)
# ==============================================================================
# 0. 数据准备 (修改部分)
# ==============================================================================

# 给节点表加上颜色，同时新增一列【节点文字颜色】
Node_With_Color <- Cyto_Node_Filtered %>%
  left_join(color_mapping_df %>% select(Module_N, Module_Color = color_mapping), 
            by = c("Parent_Module" = "Module_N")) %>%
  mutate(
    # 1. 节点本身的颜色 (保持不变)
    Final_Color = if_else(is.na(Module_Color), "#333", Module_Color),
    Size = as.numeric(Size),
    
    # 2. 【新增】节点文字的颜色 (Operon浅灰，Regulon黑)
    # 提示：'gray60' 比 'lightgray' 在白底上更容易看清，如果觉得太深可以换回 'gray80'
    Node_Text_Color = if_else(Node_Type == "Regulon", "black", "gray60") 
  )

# 边的颜色处理 (保持不变)
Edge_With_Color <- Cyto_Edge_Filtered %>%
  left_join(Node_With_Color %>% select(Node_ID, Target_Node_Color = Final_Color), 
            by = c("Target_Node" = "Node_ID"))

# ==============================================================================
# 1. 构建图对象 & 计算布局
# ==============================================================================

tbl_g <- tbl_graph(nodes = Node_With_Color, edges = Edge_With_Color) %>%
  mutate(degree = centrality_degree()) 

set.seed(123) 
Layout_Data <- create_layout(tbl_g, layout = "fr", weights = Median_Universality + 0.1)

# 【修改】计算 Label 位置和颜色
Module_Labels <- Layout_Data %>%
  filter(Node_Type != "Regulon", !is.na(Parent_Module)) %>% 
  group_by(Parent_Module) %>%
  summarise(
    x = mean(x),             
    y = min(y) - 1,          
    # 【关键】这里一定要取 Final_Color，确保和节点颜色一致
    Label_Color = first(Final_Color) 
  )

# ==============================================================================
# 2. 绘图 (修改部分)
# ==============================================================================

p <-ggraph(Layout_Data) + 
  
  # --- 1. 画边 (保持不变) ---
  geom_edge_link(
    aes(
      width = F_score,
      color = Target_Node_Color, 
      filter = Interaction == "Regulates" 
    ),
    alpha = 0.6,
    arrow = arrow(length = unit(2, 'mm'), type = "closed"), 
    end_cap = circle(2, 'mm') 
  ) +
  
  # --- 2. 画点 (保持不变) ---
  geom_node_point(
    aes(
      size = Size,            
      color = Final_Color,    
      filter = Node_Type != "Regulon" 
    ),
    alpha = 0.9
  ) +
  
  geom_node_point(
    aes(
      size = Size,
      filter = Node_Type == "Regulon"
    ),
    color = "#333333", 
    shape = 18,        
    alpha = 1
  ) +
  
  # --- 3. 画模块名称 (修改颜色) ---
  geom_text(
    data = Module_Labels,      
    aes(
      x = x, 
      y = y, 
      label = Parent_Module,
      color = Label_Color  # 【修改】这里绑定到数据中的颜色，而不是写死 "black"
    ), 
    size = 5,                  
    fontface = "bold",         
    vjust = -0.2,                
    alpha = 0.8 # 稍微调高一点不透明度，因为彩色字有时候比黑字难看清
  ) +
  
  # --- 4. 画具体节点标签 (修改颜色) ---
  geom_node_text(
    aes(
      label = Display_Label,
      color = Node_Text_Color # 【修改】这里绑定我们刚才新建的颜色列
    ),
    size = 2,
    repel = TRUE, 
    bg.color = "white", 
    bg.r = 0.1
  ) +
  
  # --- 5. 调整外观 ---
  scale_edge_width(range = c(0.2, 1.5)) + 
  scale_size(range = c(2, 8)) +           
  
  # 【重要】因为我们在aes里映射了颜色字符串，scale_color_identity 是必须的
  # 它会让 ggplot 直接读取 "red", "#333", "gray60" 这些字符串作为颜色
  scale_color_identity() +                 
  scale_edge_color_identity() +
  
  theme_graph() +                         
  theme(legend.position = "none")


p

ggsave(filename = "前十个模块的operon信息和regulon调控信息.pdf", plot = p, width = 12, height = 10, device = "pdf")




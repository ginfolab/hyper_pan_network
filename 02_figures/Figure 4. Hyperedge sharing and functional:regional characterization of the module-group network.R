library(dplyr)
library(tidyr)
library(stringr)

word_length <- 15
Category_name <- "cog"

# ==============================================================================
# 1. 基础数据准备：颜色映射与 MBGD 注释库
# ==============================================================================
# (修复了结尾多余的 %>%，并加上了我们之前说的 distinct 去重)
color_table <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/category_color_clean.txt", 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(Category = dbname, MBGD_color = color) %>% 
  filter(grepl(Category_name, Category, ignore.case = TRUE)) %>%
  mutate(MBGD_color = paste0("#", MBGD_color)) %>% 
  distinct(Category, .keep_all = TRUE) 

MBGD_anno <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/MBGD_anno.txt", 
                        header = TRUE, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE) %>%
  mutate(anno_short = ifelse(nchar(Descri) > word_length, 
                             paste0(substr(Descri, 1, word_length), "..."), 
                             Descri)) %>%
  filter(grepl(Category_name, category, ignore.case = TRUE)) %>%
  select(Annotation = category, Descri, anno_short) %>%
  left_join(color_table, by = c("Annotation" = "Category")) %>%
  rename(color_name = MBGD_color) %>%
  filter(!is.na(color_name) & color_name != "") %>%
  filter(!grepl("Unknown|Others|unknown|Other|General|Hypothetical", Descri)) %>%
  distinct(Annotation, .keep_all = TRUE) 
# ==============================================================================
# 2. 获取 Module 与 Region 以及 Gene_count (供后续 80% 阈值使用)
# ==============================================================================
Gene_Module_anno <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Module_Gene_Region_Annotation.RDS")

# 【修改点】：提前统一命名格式，防止后续 Join 失败
Module_Region_map <- Gene_Module_anno %>%
  distinct(Module, Region) %>%
  mutate(Module = str_replace(Module, "Module_", "S1_M"))

# 计算每个模块的总基因数，用于分母
module_gene_counts <- Gene_Module_anno %>%
  select(Module,Ec_clustid ) %>%
   distinct() %>%
  count(Module, name = "Total_Genes")

# ==============================================================================
# 3. 策略 A：基于富集分析显著性 (qvalue) 提取 Annotation (高优先级)
# ==============================================================================
Module_anno_qval <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Final_All_Enrichment.RDS") %>%
  filter(grepl("cog", Description, ignore.case = TRUE),
         Type == "Module") %>%
  group_by(Module) %>%
  slice_min(order_by = qvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Category = sub(":.*", "", Description)) %>% 
  select(Module, Category) 

Module_anno_color_A <- Module_anno_qval %>%
  left_join(MBGD_anno, by = c("Category" = "Annotation")) %>%
  filter(!is.na(color_name) & color_name != "") %>%
  distinct(Module, .keep_all = TRUE) %>% 
  select(Module, anno_short, Function_color = color_name) %>% 
  mutate(Module = str_replace(Module, "Module_", "S1_M"))

# ==============================================================================
# 4. 策略 B：严格频率打捞 (频次最高 && 占比 >= 80% && 组内最流行)
# ==============================================================================
long_anno <- Gene_Module_anno %>%
  pivot_longer(cols = c(cog_anno), names_to = "Source", values_to = "Annotation") %>%
  # mutate(Annotation = str_remove(Annotation, "\\..*")) %>% 
  left_join(MBGD_anno, by = "Annotation")

# 统计频率并严格过滤 80% 阈值
top_candidates_B <- long_anno %>%
  filter(!is.na(color_name)) %>%
  count(Module, Group, Annotation, anno_short, color_name, name = "Freq") %>%
  left_join(module_gene_counts, by = "Module") %>%
  mutate(Proportion = Freq / Total_Genes) %>%
  filter(Proportion >= 0.8) %>%               # 【核心】：基因占比必须达到 80% 以上
  group_by(Module) %>%
  filter(Freq == max(Freq)) %>%               # 依然要保证是模块内出现频次最高的
  ungroup()

# 组内流行度决胜机制
group_popularity_B <- top_candidates_B %>%
  distinct(Group, Module, Annotation) %>% 
  count(Group, Annotation, name = "Group_Score")

Module_anno_color_B <- top_candidates_B %>%
  left_join(group_popularity_B, by = c("Group", "Annotation")) %>%
  group_by(Module) %>%
  mutate(Random_Tie_Breaker = runif(n())) %>%
  arrange(desc(Freq), desc(Group_Score), Random_Tie_Breaker) %>%
  slice(1) %>%
  ungroup() %>%
  select(Module, anno_short, Function_color = color_name) %>%
  mutate(Module = str_replace(Module, "Module_", "S1_M"))

# ==============================================================================
# 5. 合并策略 A 与 B 并整合 Region 数据
# ==============================================================================
# 找出已经被 策略A (q-value) 成功注释的模块
annotated_by_A <- unique(Module_anno_color_A$Module)

# 将策略B中那些没有被策略A覆盖的模块补充进来
Module_anno_combined <- bind_rows(
  Module_anno_color_A,
  Module_anno_color_B %>% filter(!Module %in% annotated_by_A)
)

# 与 Region 结合，生成最终汇总表
Module_summarize <- Module_Region_map %>%
  left_join(Module_anno_combined, by = "Module")



#创建Module-subregion
Module_subregion <- Gene_Module_anno %>%
  distinct(Module, SubRegion) %>%
  mutate(Module = str_replace(Module, "Module_", "S1_M")) 
# ==============================================================================
# 6. 读取 Cyto 节点信息并生成最终颜色网络映射
# ==============================================================================
Cyto_node <- read.table("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新).txt", 
                        header = TRUE, sep = "\t", check.names = FALSE, quote = "")

#distribution_file <- read.csv("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/paper图片/Module_Topological_IM_Dynamics_Master.csv", header = TRUE, stringsAsFactors = FALSE) %>%

#  mutate(Node = str_replace(Module, "Module_", "S1_M")) %>%
#select(Node, Topological_Pattern,Overall_Dispersion, Dispersion_Level,  Max_Prop_Diff,Highest_IM_Clade, Highest_CM_Clade )  
  
  
colnames(Cyto_node)[ncol(Cyto_node)] <- "Module_color_mapping"

manual_modules <- c("S1_M1", "S1_M307", "S1_M44", "S1_M179", "S1_M11")

Cyto_node_final <- Cyto_node %>%
  left_join(Module_summarize, by = c("Node" = "Module")) %>%
   left_join(Module_subregion, by = c("Node" = "Module")) %>%
  mutate(
    Region_color = case_when(
      Region == "Region_1" ~ "#7876B1",  
      Region == "Region_2" ~ "#E18727",  
      Region == "Region_3" ~ "#0072B5" ,  
      #Region == "Region_4" ~ "turquoise",  
      TRUE ~ "white"                      
    ),
    Module_color_mapping = if_else(Module_color_mapping == "grey", "white", Module_color_mapping),
    Module_color_mapping = if_else(Node %in% manual_modules, "lightgrey", Module_color_mapping),
    Function_color = if_else(is.na(Function_color), "white", Function_color), 
    Node_Id = str_remove(Node_name, "Module_"),
    anno_short = if_else(is.na(anno_short), Node_Id, anno_short)
  ) %>% 
  mutate(
    SubRegion_color = case_when(
      # 🌟 Region 1：紫色/荧光霓虹家族 (引入多色相跳变，反差极大)
      SubRegion == "Region_1.1" ~ "#E8D5FF",  # 1.1: 明亮清透的浅丁香紫 (白底基础色)
      SubRegion == "Region_1.2" ~ "#D946EF",  # 1.2: 高饱和荧光洋红/ Fuchsia (视觉极强冲撞点)
      SubRegion == "Region_1.3" ~ "#6B21A8",  # 1.3: 纯正浓郁的深罗兰紫
      SubRegion == "Region_1.4" ~ "#2E0854",  # 1.4: 极深墨葡萄紫 (近乎黑紫，与1.2形成断层对比)
      
      # 🌟 Region 2：太阳暖色家族 (从高亮鲜黄直接切入烈火橙、再沉淀至暗红)
      SubRegion == "Region_2.1" ~ "#FDE047",  # 2.1: 极其耀眼的明亮高光黄 (在白底上瞬间锁定)
      SubRegion == "Region_2.2" ~ "#EA580C",  # 2.2: 强烈的炽热烈火橙 (与2.1黄和2.3红界限极其分明)
      SubRegion == "Region_2.3" ~ "#7F1D1D",  # 2.3: 极其深沉的暗夜血红 (高质感暗暖色)
      
      # 🌟 Region 3：冷艳冰川与深海家族 (从明亮青天蓝跨越到科技深蓝、再到极地黑蓝)
      SubRegion == "Region_3.1" ~ "#06B6D4",  # 3.1: 鲜艳夺目的冰川青/湖蓝色 (跳出普通蓝色的高亮色)
      SubRegion == "Region_3.2" ~ "#1D4ED8",  # 3.2: 严谨权威的经典皇家蓝
      SubRegion == "Region_3.3" ~ "#030712",  # 3.3: 极深邃的深渊黑蓝 (在白底上几乎呈现黑色效果)
      
      # 兜底防御（若有未定义区域，在白底上显示为极具个性的高质感钛灰色）
      TRUE ~ "white"                      
    )
  ) #%>%
  #left_join(distribution_file, by = c("Node"))




write.table(Cyto_node_final,'/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Module(Cyto)_relation_node(最新修改label和颜色).txt', sep = "\t", row.names = F, quote = F)  



##### 画图显示颜色图例 #####
##### 画图显示颜色图例 #####
library(ggplot2)
library(dplyr)

# 1. 从合并了策略A和策略B的 Module_anno_combined 中提取唯一用到的注释和颜色
legend_data <- Module_anno_combined %>%
  select(anno_short, Function_color) %>%  # 提取注释和对应的颜色
  filter(!is.na(anno_short) & !is.na(Function_color)) %>% # 双保险：过滤掉可能存在的空值
  distinct() %>%                          # 去重，保留唯一的注释和颜色对应关系
  arrange(anno_short)                     # 按字母顺序排序，让图例看起来更有条理

# 2. 画一个只有图例的图
p_legend <- ggplot(legend_data, aes(x = 1, y = anno_short, fill = Function_color)) +        
  geom_tile(width = 0.5, height = 0.5) +  # 画色块
  scale_fill_identity() +                 # 关键：告诉 ggplot 直接使用列里的十六进制颜色名称
  geom_text(aes(label = anno_short), x = 1.4, hjust = 0, size = 5) + # 添加文字标签 (x=1.4 让文字离色块近一点)
  theme_void() +                          # 去掉背景坐标轴
  xlim(0.5, 3) +                          # 调整画布范围，防止文字被截断
  labs(title = "COG Annotation Color Legend") + # 标题改为 COG
  theme(plot.title = element_text(hjust = 0.2, face = "bold", margin = margin(b = 20))) # 稍微调整一下标题位置

# 3. 显示
print(p_legend)

# 4. 保存 PDF
ggsave(
  filename = "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Cyto/Function_legend.pdf", 
  plot = p_legend, 
  width = 8,           # 宽度
  height = 6,          # 高度
  device = "pdf",      # 明确指定保存格式
  useDingbats = FALSE  # 防止在某些 PDF 查看器中出现乱码
)
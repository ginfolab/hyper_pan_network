
# ==============================================================================
# 1. 环境设置与库加载
# ==============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)

# --- 基础路径设置 ---
base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"

# --- 参数设置 ---
Target_module <- 1        
Top_selection <- 11       
ego_qvalue    <- 0.05     
ego_module_size <- 5      
max_gravity_limit <- 7   
word_len = 20
# 1. 定义你想手动添加的模块
manual_modules <- c("S1_M1", "S1_M307", "S1_M44", "S1_M179", "S1_M11")
# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 读取模块信息
# ------------------------------------------------------------------------------
# 读取原始数据并过滤 Step (先不进行 slice_max，以免丢失数据)
raw_module_data <- readRDS(file.path(base_path, "Gephi/Ec_module_modularities.RDS")) %>%
  dplyr::filter(Step %in% paste0("S", Target_module))

# 3. 找出原本 Top 11 的模块 ID
top_module_ids <- raw_module_data %>%
  slice_max(order_by = Module_sizes, n = Top_selection) %>%
  pull(tag_final)

# 4. 合并 ID (Top 11 + 手动 5 个)，并去重
final_target_ids <- unique(c(top_module_ids, manual_modules))

# 5. 根据合并后的 ID 筛选数据，并整理格式
target_module_info_raw <- raw_module_data %>%
  dplyr::filter(tag_final %in% final_target_ids) %>% # 筛选
  dplyr::rename(Module = tag_final) %>%
  dplyr::select(Module, Module_sizes) %>%
  arrange(desc(Module_sizes)) # 可选：按大小重新排序

# ------------------------------------------------------------------------------
# 2.2 读取富集结果并提取最显著描述 (核心修改)
# ------------------------------------------------------------------------------
enrichment_file <- file.path(base_path, paste0("Gephi/ALL_module_(Count >= ", ego_module_size, " & qvalue < ", ego_qvalue, ")_annotation.RDS"))

ego_all_final_Top <- readRDS(enrichment_file) %>%
  filter(Module %in% target_module_info_raw$Module) %>%
  filter(!grepl("RC_only|RC_muti", Description)) %>%
  filter(!grepl("Unknown|Others|unknown|Other|General|Hypothetical", Description)) %>% # 二次过滤无意义描述
  mutate(Description = str_remove(Description, "^.*[:~] *")) %>%    # "^.*:" 意思是：从开头 (^) 开始匹配任意字符 (.*) 直到遇到第一个冒号 (:)
  mutate(Description = str_trunc(Description, word_len ))  # 限制描述长度，太长会挤压图形

# 【核心修改】：按 qvalue 排序取最显著的一个 + 清洗字符串
module_descriptions <- ego_all_final_Top %>%
  group_by(Module) %>%
  # 1. 按 qvalue 从小到大排序，取第一个 (最显著)
  arrange(qvalue) %>% 
  slice(1) %>% 
  ungroup() %>%
  dplyr::select(Module, Description) %>%
  mutate(
    # 2. 仅保留冒号后面的内容 
    # 正则解释: ^[^:]+:\\s* 匹配开头直到第一个冒号和后面的空格
    Clean_Desc = sub("^[^:]+:\\s*", "", Description),
    # 3. 加上括号
    Label_Sub = paste0("(", Clean_Desc, ")")
  )

# ------------------------------------------------------------------------------
# 2.3 构建新的名称映射 (Old_ID -> New_Base_Name)
# ------------------------------------------------------------------------------
# 这里的逻辑是将 S1_M1 映射为 Module_1，并关联描述
node_metadata <- target_module_info_raw %>%
  left_join(module_descriptions, by = "Module") %>%
  mutate(
    # 生成基础新名字: S1_M1 -> Module_1
    New_Name_Base = str_replace(Module, "^.*_M", "Module_"),
    # 处理描述: 如果没有富集结果(NA)，则为空字符串，否则保留括号内容
    Label_Sub = ifelse(is.na(Label_Sub), "", Label_Sub)
  )

# 创建 ID 映射表: Key=S1_M1, Value=Module_1
# 注意：网络图内部还是用 Module_1 做 ID
id_map <- setNames(node_metadata$New_Name_Base, node_metadata$Module)


# ------------------------------------------------------------------------------
# 2.4 颜色映射 (修改版：使用预定义颜色 + 手动模块深灰)
# ------------------------------------------------------------------------------

# 1. 读取预定义的颜色表
# 假设该文件包含列: Module (e.g., S1_M1), color_mapping (e.g., #FF0000 或 red)
ref_color_df <- readRDS("/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli/Gephi/Module_color_mapping.RDS") %>% 
  dplyr::select(Module, color_mapping)

# 2. 初始化最终颜色向量
# final_target_ids 包含了所有要画的模块 (Top 11 + Manual 5)
# 默认将所有模块颜色设为 "grey30" (深灰色)
final_colors <- setNames(rep("grey30", length(final_target_ids)), final_target_ids)

# 3. 匹配并覆盖颜色
# 找出哪些目标模块在参考表中存在
valid_modules <- intersect(final_target_ids, ref_color_df$Module)

# 提取这些存在的模块对应的颜色
matched_colors <- ref_color_df %>%
  filter(Module %in% valid_modules) %>%
  pull(color_mapping, name = Module)

# 用参考表中的颜色覆盖默认的灰色
# 剩下的 Manual 模块因为不在 ref_color_df 中，所以保持 "grey30"
final_colors[names(matched_colors)] <- matched_colors

# 4. 关键步骤：转换键名 (Old_ID -> New_Base_Name)
# 因为画图时节点的 name 是 "Module_1", "Module_44" 等
# 而目前的 final_colors 的 names 是 "S1_M1", "S1_M44" 等
# 我们需要利用 2.3 中生成的 id_map 进行转换

# id_map 的结构是: names = S1_M1, value = Module_1
# 我们根据 final_colors 的旧名字找到对应的新名字
new_names_for_color <- id_map[names(final_colors)]

# 更新颜色向量的名称
names(final_colors) <- new_names_for_color

# 5. 赋值给 color_mapping 供画图使用
color_mapping <- final_colors

# ==============================================================================
# 3. 构建网络数据
# ==============================================================================

# 准备绘图数据 (Top 10 connection)
plot_data <- ego_all_final_Top %>%
  group_by(Module) %>%
  arrange(qvalue) %>%
  slice_head(n = 10) %>% 
  ungroup() %>%
  mutate(Significance = -log10(qvalue)) %>%
  # 【关键】：把原始 Module 列换成新名字
  mutate(Module_New = id_map[Module])

# --- A. 构建节点 (Nodes) ---
# Module 节点
nodes_module <- node_metadata %>%
  mutate(name = New_Name_Base, type = "Module") %>% # name 是 Module_1
  rename(size = Module_sizes) %>%
  select(name, size, type, Label_Sub) # 保留 Label_Sub 供画图用

# Term 节点
nodes_term <- plot_data %>%
  group_by(Description) %>%
  summarise(size = sum(Count)) %>%
  rename(name = Description) %>%
  mutate(type = "Term", Label_Sub = NA) 

nodes <- bind_rows(nodes_module, nodes_term)

# --- B. 构建边 (Edges) ---
edges <- plot_data %>%
  select(Module_New, Description, Significance) %>% 
  rename(from = Module_New, to = Description) %>%
  mutate(visual_width = Significance,
         weight = pmin(Significance, max_gravity_limit), 
         Source_Module = from)

# ==============================================================================
# 4. 创建图对象
# ==============================================================================
igraph_obj <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
graph <- as_tbl_graph(igraph_obj)

# ==============================================================================
# 5. 绘图 (双层标签法)
# ==============================================================================
p_network <- ggraph(graph, layout = "fr", niter = 5000) + 
  
  # 1. 画边
  geom_edge_link(aes(width = visual_width, color = Source_Module), alpha = 0.5) +
  scale_edge_width(range = c(0.5, 3), name = "-log10(qvalue)") +
  scale_edge_color_manual(values = color_mapping) +
  
  # 2. 画 Module 节点
  geom_node_point(aes(size = size, filter = type == "Module", color = name)) +
  scale_color_manual(values = color_mapping) +
  
  # 3. 画 Term 节点
  geom_node_point(aes(size = size, filter = type == "Term"), color = "grey80", alpha = 0.8) +
  scale_size(range = c(2, 18), name = "Gene Count") +
  
  # 4. 加 Term 标签
  geom_node_text(aes(label = name, filter = type == "Term"), 
                 repel = TRUE, size = 4, color = "black",alpha = 0.5, bg.r = 0.1, max.overlaps = 30) +
  
  # 5. 加 Module 标签 - 第一层：模块名 (大号字体)
  geom_node_text(aes(label = name, filter = type == "Module"), 
                 fontface = "bold", 
                 color = "black", 
                 size = 8, 
                # alpha = 0.5,
                 nudge_y = 0.2, # 向上提
                 repel = FALSE) + 
  
  # 6. 加 Module 标签 - 第二层：功能描述 (小号字体)
  geom_node_text(aes(label = Label_Sub, filter = type == "Module"), 
                 fontface = "plain", 
                 color = "black", 
                 size = 5,        # 大小减半
                # alpha = 0.5,
                 nudge_y = -0.25, # 向下压，避开上面的大字
                 repel = FALSE) +
  
  theme_void() +
  theme(legend.position = "none") + 
  labs(title = paste0("Module-Function Network"))

# ==============================================================================
# 6. 保存
# ==============================================================================
print(p_network)

save_file <- file.path(base_path, paste0("paper图片/Figure S5. Functional connectivity network of 16 representative modules.pdf"))
ggsave(save_file, p_network, width = 14, height = 12)

cat("绘图完成！文件已保存至:", save_file, "\n")
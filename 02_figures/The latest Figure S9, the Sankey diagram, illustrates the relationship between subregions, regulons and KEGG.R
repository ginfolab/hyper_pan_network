# ==============================================================================
# Figure S10 终极绝美版: 渐变色映射 + 【最小交叉拓扑排序】 + 完美排版 (Top-N 显著版)
# ==============================================================================

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

suppressMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel) 
  library(scales) 
})

base_path <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
excel_file <- file.path(base_path, 'ALL_module_(Count_ge_5_qvalue_lt_0.05)_S2_Removed.RDS')

# 1. 读取并基础过滤
filtered_data <- readRDS(excel_file) %>% filter(Type == "Region" )

# ==============================================================================
# 2. 【核心修改】提取上游：每个 SubRegion 保留最显著的前 N 个 Regulon
# ==============================================================================
top_n_num <- 3  # 👉 【参数调节区】你想保留前几个？可以在这里直接改成 3 或 5

message(paste0(">>> 正在提取每个 SubRegion 最显著的前 ", top_n_num, " 个 Regulon..."))

df_regulon <- filtered_data %>% filter(Term_type == "Regulon") %>%
  mutate(Regulon = str_remove(Description, "^Regulon:"), Reg_Score = -log10(p.adjust)) %>%
  select(Node_ID, Regulon, Reg_Score) %>% distinct()

# 【新逻辑】：按 SubRegion (Node_ID) 分组，根据 Reg_Score 降序取前 N 个
df_regulon_top <- df_regulon %>%
  group_by(Node_ID) %>%
  slice_max(order_by = Reg_Score, n = top_n_num, with_ties = FALSE) %>%
  ungroup()

# ==============================================================================
# 3. 下游智能降级补位 (Smart Fallback)
# ==============================================================================
message(">>> 正在提取核心功能...")

df_function_smart <- filtered_data %>%
  filter(Term_type %in% c("KEGG", "Kmodule", "TIGR", "COG")) %>%
  group_by(Node_ID) %>%
  mutate(Best_Type = case_when(
    any(Term_type == "Kmodule")    ~ "Kmodule",
    any(Term_type == "KEGG") ~ "KEGG",
    any(Term_type == "TIGR")    ~ "TIGR",
    TRUE ~ "COG"
  )) %>%
  filter(Term_type == Best_Type) %>%
  ungroup() %>%
  mutate(
    Function = str_extract(Description, "(?<=:)[^;]+"),
    Function = str_remove(Function, ",.*"),
    Function = str_trim(Function),
    Function = str_trunc(Function, width = 30, ellipsis = "..."),
    Func_Score = -log10(p.adjust) 
  ) %>%
  arrange(Node_ID, desc(Func_Score)) %>%
  distinct(Node_ID, Function, .keep_all = TRUE) %>%
  select(Node_ID, Function, Func_Score)

# ==============================================================================
# 4. 联结网络，计算压缩动态线宽
# ==============================================================================
# 【修改】：将 df_regulon_shared 替换为 df_regulon_top
sankey_data <- df_regulon_top %>%
  inner_join(df_function_smart, by = "Node_ID", relationship = "many-to-many") %>%
  mutate(
    SubRegion = Node_ID,
    Parent_Region = str_extract(Node_ID, "Region_\\d+"),
    Combined_Weight_Raw = Reg_Score + Func_Score,
    Adjusted_Weight = scales::rescale(Combined_Weight_Raw, to = c(1, 5))
  )

# ==============================================================================
# 5. 【全新黑魔法】：根据 SubRegion 拓扑关系进行「最小交叉分层排序」
# ==============================================================================
message(">>> 正在进行最小交叉拓扑重排 (消除乱麻连线)...")

# 1. 中轴 (SubRegion) 的绝对顺序：自然数字排序
subregion_levels <- str_sort(unique(sankey_data$SubRegion), numeric = TRUE)

# 2. 右轴 (Function) 排序逻辑：寻找它的"主导 SubRegion"，并跟随其物理位置！
function_levels <- sankey_data %>%
  group_by(Function, SubRegion) %>%
  summarise(Weight = sum(Combined_Weight_Raw), .groups = "drop") %>%
  group_by(Function) %>%
  arrange(desc(Weight)) %>%
  slice(1) %>% 
  ungroup() %>%
  mutate(SubRegion_Rank = match(SubRegion, subregion_levels)) %>% 
  arrange(SubRegion_Rank, desc(Weight)) %>% 
  pull(Function)

# 3. 左轴 (Regulon) 排序逻辑：同理，跟随主导 SubRegion 排序以减少交叉
regulon_levels <- sankey_data %>%
  group_by(Regulon, SubRegion) %>%
  summarise(Weight = sum(Combined_Weight_Raw), .groups = "drop") %>%
  group_by(Regulon) %>%
  arrange(desc(Weight)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(SubRegion_Rank = match(SubRegion, subregion_levels)) %>%
  arrange(SubRegion_Rank, desc(Weight)) %>%
  pull(Regulon)

sankey_data <- sankey_data %>%
  mutate(
    Regulon = factor(Regulon, levels = regulon_levels),
    SubRegion = factor(SubRegion, levels = subregion_levels),
    Function = factor(Function, levels = function_levels)
  )

# ==============================================================================
# 6. 动态调配亚区渐变色 (Gradient Colors)
# ==============================================================================
base_region_colors <- c("Region_1" = "#7876B1", "Region_2" = "#E18727", 
                        "Region_3" = "#0072B5", "Region_4" = "#35A595")

lighten_col <- function(col, frac) {
  rgb_v <- col2rgb(col); rgb_new <- rgb_v + (255 - rgb_v) * frac
  rgb(rgb_new[1], rgb_new[2], rgb_new[3], maxColorValue=255)
}
darken_col <- function(col, frac) {
  rgb_v <- col2rgb(col); rgb_new <- rgb_v * (1 - frac)
  rgb(rgb_new[1], rgb_new[2], rgb_new[3], maxColorValue=255)
}

subregion_colors <- list()
for (pr in names(base_region_colors)) {
  subs <- str_sort(unique(sankey_data$SubRegion[sankey_data$Parent_Region == pr]), numeric = TRUE)
  n_subs <- length(subs)
  if (n_subs > 0) {
    if (n_subs == 1) {
      subregion_colors[[subs[1]]] <- base_region_colors[pr]
    } else {
      c_light <- lighten_col(base_region_colors[pr], 0.5) 
      c_dark  <- darken_col(base_region_colors[pr], 0.3)
      pal <- colorRampPalette(c(c_light, base_region_colors[pr], c_dark))(n_subs)
      names(pal) <- subs
      subregion_colors <- c(subregion_colors, pal)
    }
  }
}
subregion_colors <- unlist(subregion_colors)

# ==============================================================================
# 7. 渲染终极桑基图 (完美动态排版)
# ==============================================================================
message(">>> 正在渲染带渐变色与防重叠标签的终极大图...")

p_sankey <- ggplot(sankey_data, 
                   aes(y = Adjusted_Weight, axis1 = Regulon, axis2 = SubRegion, axis3 = Function)) +
  
  geom_flow(aes(fill = SubRegion), width = 1/5, alpha = 0.8, curve_type = "cubic") + 
  geom_stratum(width = 1/5, color = "grey40", fill = "grey95", alpha = 0.9) +
  
  geom_text_repel(
    stat = "stratum", 
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)), 
    size = 3.5, color = "black", direction = "y", 
    nudge_x = -0.25, hjust = 1,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3
  ) +
  
  geom_text_repel(
    stat = "stratum", 
    aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA)), 
    size = 3.5, color = "black", direction = "y", 
    nudge_x = 0, hjust = 0.5,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3
  ) +
  
  geom_text_repel(
    stat = "stratum", 
    aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA)), 
    size = 3.5, color = "black", direction = "y", 
    nudge_x = 0.25, hjust = 0,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3
  ) +
  
  # 【修改】：将 X 轴左侧标签从 "Shared Regulators" 改为 "Top Regulators"
  scale_x_discrete(limits = c("Top Regulators", 
                              "Structural Sub-Regions", 
                              "KEGG pathway"), 
                   expand = c(0.25, 0.35)) + 
  
  scale_fill_manual(values = subregion_colors) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", color = "black", vjust = 5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0.5, margin = margin(b = 20)),
    legend.position = "none" 
  ) +
  labs(
    y = "", fill = "",
    # 【修改】：主标题相应修改
    title = "Information Flow: From Top Regulons to Phenotypic Functions"
  )
p_sankey



output_pdf <- file.path(base_path, "Figure_S9_Ultimate_Information_Flow_TopN.pdf")

pdf(output_pdf, width = 20, height = max(10, length(unique(sankey_data$SubRegion)) * 0.9))
print(p_sankey)
dev.off()

message(paste0(">>> 绝美水平直流版生成完毕！节点已按拓扑归类，交叉降至最低！\n>>> 文件已保存至: ", output_pdf))
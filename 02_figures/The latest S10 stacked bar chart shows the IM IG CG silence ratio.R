suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork) 
  library(ggrepel)   
})

message(">>> 正在启动 Figure S13 饱满平铺版: (精准标记 Group_14，修复合并Bug)...")

# ==============================================================================
# 0. 全局路径与参数配置
# ==============================================================================
base_path     <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
PCC_power     <- 1
max_cap_value <- 1000  

# 读取底层文件
input_raw_csv <- file.path(base_path, "PCC", paste0("Raw_Weighted_7State_Top5PCT_withModularity_PCC", PCC_power, ".csv"))
if (!file.exists(input_raw_csv)) stop("[-] 找不到底层的 Raw CSV 文件。")
df_raw <- read_csv(input_raw_csv, show_col_types = FALSE)

# 读取注释字典
anno_path <- file.path(base_path, "Module_Gene_Region_Annotation.RDS")
if (!file.exists(anno_path)) stop("[-] 找不到 Module_Gene_Region_Annotation.RDS")
Module_anno <- readRDS(anno_path)

region_counts <- Module_anno %>%
  filter(!is.na(Region)) %>%
  group_by(Region) %>%
  summarise(Gene_Count = n_distinct(Ec_clustid), .groups = "drop") %>%
  mutate(Width_A = 0.2 + 0.7 * (Gene_Count / max(Gene_Count)))

# ==============================================================================
# 🌟 第1部分：【图 A 数据准备】
# ==============================================================================
df_clean <- df_raw %>% filter(!is.na(Region), !is.na(Clade))

positive_mod <- df_clean %>% filter(Ratio_Silenced == 0, Modularity > 0) %>% pull(Modularity)
quantiles <- quantile(positive_mod, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

df_global <- df_clean %>%
  mutate(
    Modularity_Group = case_when(
      Ratio_Silenced == 1 | is.na(Modularity) | Modularity <= 0 ~ "Silenced (NA)", 
      Modularity <= quantiles[2] ~ "Low (0-25%)",
      Modularity <= quantiles[3] ~ "Medium-Low (25-50%)",
      Modularity <= quantiles[4] ~ "Medium-High (50-75%)",
      TRUE ~ "High (75-100%)"
    )
  ) %>%
  mutate(Modularity_Group = factor(Modularity_Group, levels = c("Silenced (NA)", "Low (0-25%)", "Medium-Low (25-50%)", "Medium-High (50-75%)", "High (75-100%)"))) %>%
  group_by(Clade, Region, Modularity_Group) %>%
  summarise(Module_Count = n(), .groups = "drop") %>%   
  group_by(Clade, Region) %>%
  mutate(
    Total_Modules = sum(Module_Count),                  
    Proportion = Module_Count / Total_Modules,          
    Label_Text = ifelse(Proportion > 0.05, sprintf("%.1f", Proportion * 100), ""),
    Dummy_X = "" 
  ) %>%
  ungroup() %>%
  left_join(region_counts, by = "Region")

df_global$Region <- factor(df_global$Region, levels = str_sort(unique(df_global$Region), numeric = TRUE))
df_global$Clade  <- factor(df_global$Clade, levels = str_sort(unique(df_global$Clade), numeric = TRUE)) 

# ==============================================================================
# 🌟 第2部分：【图 B 数据准备】(含 Group 标签合并)
# ==============================================================================
df_processed <- df_raw %>%
  filter(Ratio_Silenced == 0, !is.na(Region), Clade %in% c("Clade_1", "Clade_2")) %>%
  mutate(
    Region_Clean = paste0("Region_", str_remove(Region, "^Region_")),
    Merged_IM  = d_S_IM + d_U_IM,
    Merged_IG  = d_S_IG + d_U_IG,
    Merged_CG  = d_S_CG + d_U_CG,
    Merged_Acc = d_Acc,
    IM_per_gene  = if_else(Valid_Nodes > 0, (2 * Merged_IM) / Valid_Nodes, 0),
    IG_per_gene  = if_else(Valid_Nodes > 0, (2 * Merged_IG) / Valid_Nodes, 0),  
    CG_per_gene  = if_else(Valid_Nodes > 0, (2 * Merged_CG) / Valid_Nodes, 0),  
    Acc_per_gene = if_else(Valid_Nodes > 0, Merged_Acc / Valid_Nodes, 0)
  ) %>%
  filter(Region_Clean %in% c("Region_1", "Region_2"))

df_stats <- df_processed %>%
  group_by(Community, Clade, Region_Clean) %>%
  summarise(
    IM_mean  = mean(IM_per_gene, na.rm = TRUE), IM_sd  = sd(IM_per_gene, na.rm = TRUE),
    IG_mean  = mean(IG_per_gene, na.rm = TRUE), IG_sd  = sd(IG_per_gene, na.rm = TRUE),
    CG_mean  = mean(CG_per_gene, na.rm = TRUE), CG_sd  = sd(CG_per_gene, na.rm = TRUE),
    Acc_mean = mean(Acc_per_gene, na.rm = TRUE), Acc_sd = sd(Acc_per_gene, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_sd"), ~ if_else(is.na(.), 0, .)))

df_axis_ready <- df_stats %>%
  pivot_wider(
    id_cols = c(Community, Region_Clean),
    names_from = Clade,
    values_from = c(IM_mean, IM_sd, IG_mean, IG_sd, CG_mean, CG_sd, Acc_mean, Acc_sd)
  ) %>%
  filter(!is.na(IM_mean_Clade_1) & !is.na(IM_mean_Clade_2))

# 🌟 终极修复：用 Module 去和 Community 匹配，而不是 Ec_clustid！
df_axis_ready <- df_axis_ready %>%
  left_join(
    Module_anno %>% select(Module, Group) %>% distinct(),
    by = c("Community" = "Module")
  )

# ==============================================================================
# 🌟 第3部分：【图 A 渲染】
# ==============================================================================
color_palette_mod <- c("Silenced (NA)" = "#BDC3C7", "Low (0-25%)" = "#AED6F1", "Medium-Low (25-50%)" = "#5DADE2", "Medium-High (50-75%)" = "#2874A6", "High (75-100%)" = "#154360")
base_theme <- theme_minimal() + theme(
  strip.text.x = element_text(size = 14, face = "bold", color = "#2C3E50", margin = margin(b = 10)),
  strip.text.y.left = element_text(size = 14, face = "bold", color = "#2C3E50", angle = 90, margin = margin(r = 10)),
  axis.title.y = element_text(size = 14, face = "bold", color = "#2C3E50", margin = margin(r = 15)),
  axis.text.y = element_text(size = 11, color = "grey30"),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "#2C3E50", fill = NA, linewidth = 1), 
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20), color = "#2C3E50"),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 12),
  plot.margin = margin(20, 10, 20, 20)
)

p_global <- ggplot(df_global, aes(x = Dummy_X, y = Proportion, fill = Modularity_Group)) +
  geom_col(aes(width = Width_A), position = "stack", color = "#2C3E50", linewidth = 0.6) +
  geom_text(aes(label = Label_Text, color = Modularity_Group), position = position_stack(vjust = 0.5), fontface = "bold", size = 4) +
  scale_color_manual(values = c("Silenced (NA)"="white", "Low (0-25%)"="black", "Medium-Low (25-50%)"="black", "Medium-High (50-75%)"="white", "High (75-100%)"="white"), guide = "none") +
  facet_grid(Clade ~ Region, switch = "y") +
  scale_fill_manual(name = "Modularity Quantiles:", values = color_palette_mod, drop = FALSE, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(labels = function(x) x * 100, expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of Modules (%)", title = "A. Composition of Modules by Modularity") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

# ==============================================================================
# 🌟 第4部分：【图 B 渲染】 
# ==============================================================================
color_palette_region <- c("Region_1" = "#7876B1", "Region_2" = "#E18727")

build_scatter_plot <- function(df, edge_type, show_labels = FALSE, target_group = "Group_14") {
  x_mean_col <- paste0(edge_type, "_mean_Clade_1")
  x_sd_col   <- paste0(edge_type, "_sd_Clade_1")
  y_mean_col <- paste0(edge_type, "_mean_Clade_2")
  y_sd_col   <- paste0(edge_type, "_sd_Clade_2")
  
  plot_df <- df %>%
    select(Community, Region_Clean, Group, 
           X_mean = !!sym(x_mean_col), X_sd = !!sym(x_sd_col),
           Y_mean = !!sym(y_mean_col), Y_sd = !!sym(y_sd_col)) %>%
    mutate(X_mean_capped = pmin(X_mean, max_cap_value), Y_mean_capped = pmin(Y_mean, max_cap_value))
  
  max_center_value <- max(c(plot_df$X_mean_capped, plot_df$Y_mean_capped, 5), na.rm = TRUE)
  axis_max         <- max_center_value * 1.15 
  
  p <- ggplot(plot_df, aes(x = X_mean_capped, y = Y_mean_capped)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, Y_mean_capped - Y_sd), ymax = pmin(axis_max, Y_mean_capped + Y_sd)), 
                  color = "grey75", width = 0, linewidth = 0.4, alpha = 0.5) +
    geom_errorbarh(aes(xmin = pmax(0, X_mean_capped - X_sd), xmax = pmin(axis_max, X_mean_capped + X_sd)), 
                   color = "grey75", height = 0, linewidth = 0.4, alpha = 0.5) +
    geom_point(aes(fill = Region_Clean), color = "#2C3E50", shape = 21, size = 3.8, stroke = 0.6, alpha = 0.85)
  
  # 🌟 给指定 Group 打标签
  if (show_labels) {
    # 不再限制只能是黄点(Region_2)，只要是这个 Group 全标出来
    label_data <- plot_df %>% filter(!is.na(Group) & Group == target_group)
    
    cat("👉 [雷达诊断] 在", edge_type, "图中，成功定位到", nrow(label_data), "个", target_group, "的模块。\n")
    
    if (nrow(label_data) > 0) {
      p <- p + geom_text_repel(
        data = label_data,
        aes(label = Community),        
        size = 4,                 
        color = "#C0392B",          
        fontface = "bold",
        max.overlaps = 100,         
        box.padding = 0.8,          
        point.padding = 0.5, 
        segment.color = "grey30",
        segment.size = 0.5,
        min.segment.length = 0      
      )
    }
  }
  
  p <- p +
    scale_fill_manual(name = "Network Regions:", values = color_palette_region) +
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    labs(title = edge_type, x = "Clade 1 (Mean Strength)", y = "Clade 2 (Mean Strength)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "#2C3E50", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "grey92", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", color = "#2C3E50", hjust = 0.5, margin = margin(b=6)),
      axis.title = element_text(size = 10, face = "bold", color = "#2C3E50"),
      axis.text = element_text(size = 9, color = "grey30"),
      plot.margin = margin(6, 6, 6, 6) 
    )
  return(p)
}

# ==============================================================================
# 调用函数绘图 (仅对 IG 开启特定 Group 的标签引擎)
# ==============================================================================
p_im  <- build_scatter_plot(df_axis_ready, "IM")
p_ig  <- build_scatter_plot(df_axis_ready, "IG", show_labels = TRUE, target_group = "Group_14") 
p_cg  <- build_scatter_plot(df_axis_ready, "CG")
p_acc <- build_scatter_plot(df_axis_ready, "Acc")

p_b_matrix <- (p_im + p_ig) / (p_cg + p_acc) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    plot.margin = margin(2, 2, 2, 2)
  )

p_b_final <- p_b_matrix + plot_annotation(
  title = "B.Region 1 vs Region 2 (Clade 1 to Clade 2)",
  subtitle = "Points represent individual modules; dashed lines indicate the y=x equilibrium boundary.",
  theme = theme(
    plot.title = element_text(size = 16, face = "bold", color = "#2C3E50", hjust = 0.5, margin = margin(b=5)),
    plot.subtitle = element_text(size = 11, face = "italic", color = "grey40", hjust = 0.5, margin = margin(b=10))
  )
)

# ==============================================================================
# 🌟 第5部分：【终极全景拼接与比例拉伸】
# ==============================================================================
p_final <- p_global + wrap_elements(p_b_final) + plot_layout(ncol = 2, widths = c(1, 3))

output_pdf <- file.path(base_path, "paper图片", paste0("FigureS10_Global_AandB_Expanded_PCC", PCC_power, ".pdf"))
ggsave(output_pdf, plot = p_final, width = 26, height = 13, device = cairo_pdf)

message(paste0("\n[+] 🏆 完美平铺！散点全面释放！请查看：\n    ", output_pdf))
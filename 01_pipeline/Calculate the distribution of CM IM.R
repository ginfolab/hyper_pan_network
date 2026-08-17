suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(patchwork) # 🌟 拼图神器
})

message(">>> 正在启动 Figure S12 (极简图注版): 25/50/75% 面积切分对比图...")

# ==============================================================================
# 0. 全局路径与参数配置
# ==============================================================================
base_path         <- "/Users/jiangzhenbo/Desktop/Rwork/test/New_Escherichia_coli"
PCC_power         <- 1 

unweighted_csv <- file.path(base_path, "PCC", paste0("Raw_Unweighted_Topology_PCC", PCC_power, ".csv"))
weighted_csv   <- file.path(base_path, "PCC", paste0("Raw_Weighted_PureEnergy_PCC", PCC_power, ".csv"))

if (!file.exists(unweighted_csv) | !file.exists(weighted_csv)) {
  stop("[-] 找不到底层的 Raw CSV 文件，请检查路径。")
}

# ==============================================================================
# 🌟 核心函数：1D 核密度图与 25/50/75 颜色映射解算流水线
# ==============================================================================
build_density_plot <- function(input_file, plot_title) {
  
  df_raw <- read_csv(input_file, show_col_types = FALSE)
  
  # 1. 提取并合并 Integration State (CM = IG + CG)
  df_density <- df_raw %>%
    filter(!is.na(Region), !is.na(Clade), !is.na(Ratio_IM)) %>%
    mutate(CM_Proportion = Ratio_IG + Ratio_CG) 
  
  df_density$Region <- factor(df_density$Region, levels = str_sort(unique(df_density$Region), numeric = TRUE))
  df_density$Clade  <- factor(df_density$Clade,  levels = rev(str_sort(unique(df_density$Clade), numeric = TRUE)))
  
  # 2. 精确计算面积分布的 25%, 50%, 75% 截断点 (🌟 移除文本标签逻辑)
  cutoffs_df <- df_density %>%
    group_by(Clade, Region) %>%
    summarise(
      Q25 = quantile(CM_Proportion, 0.25, na.rm = TRUE),
      Q50 = quantile(CM_Proportion, 0.50, na.rm = TRUE),
      Q75 = quantile(CM_Proportion, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "Quantile", values_to = "Cutoff_Value")
  
  # 3. 精确计算 1D 核密度曲线
  density_df <- df_density %>%
    group_by(Clade, Region) %>%
    summarise(
      dens = list(density(CM_Proportion, from = 0, to = 1, n = 512, adjust = 1.2)),
      .groups = "drop"
    ) %>%
    mutate(
      x = purrr::map(dens, ~.$x),
      y = purrr::map(dens, ~.$y)
    ) %>%
    select(-dens) %>%
    unnest(cols = c(x, y))
  
  # 4. 寻找密度曲线与三条垂直切分线的精准交点
  intersect_dots <- cutoffs_df %>%
    left_join(density_df, by = c("Clade", "Region"), relationship = "many-to-many") %>%
    group_by(Clade, Region, Quantile) %>%
    slice_min(abs(x - Cutoff_Value), n = 1) %>% 
    ungroup() %>%
    select(Clade, Region, Quantile, Cutoff_Value, intersect_y = y)
  
  # 5. 渲染高阶 1D 密度分布图
  p_density <- ggplot() +
    # 背景冰蓝色山峰
    geom_area(data = density_df, aes(x = x, y = y), fill = "#A9D6E5", alpha = 0.7, color = "#2C3E50", linewidth = 0.8) +
    
    # 🌟 绘制 3 条垂直虚线 (将颜色映射到 Quantile)
    geom_segment(data = intersect_dots, aes(x = Cutoff_Value, xend = Cutoff_Value, y = 0, yend = intersect_y, color = Quantile), 
                 linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
    
    # 🌟 绘制 3 个交点 (边缘颜色同样映射到 Quantile)
    geom_point(data = intersect_dots, aes(x = Cutoff_Value, y = intersect_y, color = Quantile), 
               size = 2, shape = 21, fill = "white", stroke = 1.2) +
    
    # 🌟 彻底移除了 geom_text，画面极其干净
    
    facet_grid(Clade ~ Region, switch = "y") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1), labels = scales::percent_format()) +
    
    # 🌟 添加统一的高级颜色映射图注 (Legend)
    scale_color_manual(
      name = "Cumulative Density Area Cutoffs:",
      values = c("Q25" = "#27AE60", "Q50" = "#E74C3C", "Q75" = "#8E44AD"),
      labels = c("Q25" = "25% Area (1st Quartile)", "Q50" = "50% Area (Median)", "Q75" = "75% Area (3rd Quartile)")
    ) +
    
    labs(
      x = "Integration State Proportion (CM = Ratio_CG + Ratio_IG)",
      y = "Kernel Density Estimation (Module Density)",
      title = plot_title
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20), color = "#2C3E50"),
      strip.text.x = element_text(size = 14, face = "bold", color = "#2C3E50", margin = margin(b = 10)),
      strip.text.y.left = element_text(size = 14, face = "bold", color = "#2C3E50", angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(size = 14, face = "bold", color = "#2C3E50", margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, face = "bold", color = "#2C3E50", margin = margin(r = 15)),
      axis.text = element_text(size = 11, color = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.border = element_rect(color = "#2C3E50", fill = NA, linewidth = 0.8),
      
      # 图注的字体设置
      legend.title = element_text(size = 14, face = "bold", color = "#2C3E50"),
      legend.text = element_text(size = 13, color = "#2C3E50"),
      legend.key.width = unit(2, "cm"),
      
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p_density)
}

# ==============================================================================
# 🌟 第4部分：生成两张密度图并使用 Patchwork 左右无缝拼接
# ==============================================================================
message(">>> [Phase 4] 正在计算并渲染 Unweighted (未加权) AUC 面积切分图...")
p_unweighted_density <- build_density_plot(
  input_file = unweighted_csv, 
  plot_title = "A. Unweighted Topology (Quartile Distribution)"
)

message(">>> [Phase 4] 正在计算并渲染 Weighted (加权能量) AUC 面积切分图...")
p_weighted_density <- build_density_plot(
  input_file = weighted_csv, 
  plot_title = "B. Weighted Energy (Quartile Distribution)"
)

# 🌟 使用 Patchwork 将两张图左右拼接，并智能合并底部图注 (guides = "collect")
message(">>> [Phase 5] 正在执行无缝合并渲染并导出超宽 PDF...")
p_combined_density <- (p_unweighted_density + p_weighted_density) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

output_pdf_combined <- file.path(base_path, "paper图片", paste0("FigureS12_Quartile_Legend_Comparison_PCC", PCC_power, ".pdf"))

# 画布宽度拉满，完美容纳两张图 (28x10)
ggsave(output_pdf_combined, plot = p_combined_density, width = 28, height = 10, device = cairo_pdf)

message(paste0("[+] 史诗级双核极简图注版对比图渲染完成！保存至：\n    ", output_pdf_combined))